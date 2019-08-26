#--------------------
# Stats
#--------------------

F_test <- function(stat, v1, v2) 
{
  # see rrBLUP::GWAS (local function `score.calc`)
  pbeta(v2 / (v2 + v1 * stat), v2/2, v1/2)
}

#--------------------
# Formula
#--------------------

#' @export
parseFormula <- function(f)
{
  # @ http://stackoverflow.com/questions/14049440/extract-components-from-mixed-model-lme4-formula
  
  ### args
  f <- as.formula(f)
  
  ### process
  terms <- terms(f)
  
  # vars
  vars <- all.vars(terms) # character
  
  response.ind <- attributes(terms)$response
  response <- vars[response.ind]
  
  # effects
  labels <- attributes(terms)$term.labels
  labels <- gsub(" ", "", labels)

  ind <- grepl("\\|", labels)
  ranef <- labels[ind]
  fixef <- labels[!ind]
  
  splits <- strsplit(ranef, "\\|")
  stopifnot(all(laply(splits, length) == 2))
  
  ranef1 <- laply(splits, function(x) x[1])
  ranef2 <- laply(splits, function(x) x[2])
  
  ### return
  list(vars = vars, labels = labels,
    response = response, ranef = ranef, ranef1 = ranef1, ranef2 = ranef2,
    fixef = fixef)
}

#--------------------
# Linear Alegbra
#--------------------

#' @export
relfac <- function(mat, method.relfac = "auto", tol.relfac.svd = 1e-10, tol.relfac.evd = 1e-10) 
{
  # ?Matrix::chol
  # Returned value: a matrix of class Cholesky, i.e., upper triangular: R such that R'R = x.
  # Note that another notation is equivalent x = L L', where L is a lower triangular 
  # @ http://en.wikipedia.org/wiki/Cholesky_decomposition
  #
  # If the substitution is Z* = Z L, then Z*' = L' Z' = R Z'      
  # @ http://www.journalofanimalscience.org/content/88/2/497.long%5BVazquez%20et%20al.,%202010%5D

  # inc
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  
  rmat <- switch(method.relfac, 
    "chol" = relfac.chol(mat),
    "svd" = relfac.svd(mat, tol.relfac.svd),
    "evd" = relfac.evd(mat, tol.relfac.evd),
    "auto" = {
      # filters
      mat.duplicated <- any(duplicated(llply(1:ncol(mat), function(i) mat[, i])))
      
      # auto
      if(mat.duplicated) {
        rmat <- relfac.evd(mat, tol.relfac.evd)
      } else {
        rmat <- suppressWarnings(try({
          relfac.chol(mat)
        }, silent = TRUE))
        
        if(class(rmat) == "try-error") {
          rmat <- relfac.evd(mat, tol.relfac.evd)
        }
      }
      
      rmat
    },
    stop())
    
  rownames(rmat) <- rownames(mat)
  colnames(rmat) <- colnames(mat)
  
  return(rmat)
}

relfac.chol <- function(mat)
{
  Matrix::chol(mat)
}

relfac.svd <- function(mat, tol.relfac.svd)
{
  mat <- as.matrix(mat)
  out <- La.svd(mat)
  
  d <- out$d
  
  D <- diag(out$d)
  D2 <- sqrt(D)

  Lt <- D2 %*% out$vt
  
  ind <- (abs(Lt) < tol.relfac.svd)
  Lt[ind] <- 0
  
  #Matrix::Matrix(Lt)
  as(Lt, "dgCMatrix")
}

relfac.evd <- function(mat, tol)
{
  out <- eigen(mat, symmetric = TRUE)
  
  # clean eigen values
  ind <- (abs(out$values) < tol)
  if(any(ind)) {
    out$values[ind] <- 0
  }

  # clean eigen vectors only if eigen values were previously cleaned
  if(any(ind)) {
    ind_vec <- (abs(out$vectors) < tol)
    if(any(ind_vec)) {
      out$vectors[ind_vec] <- 0
    }
  }
  
  # return
  R <- diag(sqrt(out$values)) %*% t(out$vectors)
  
  #return(Matrix::Matrix(R))
  return(as(R, "dgCMatrix"))
  
  ### old code
  if(length(ind)) {
    R <- diag(sqrt(out$values[-ind])) %*% t(out$vectors[, -ind]) 
    #L <- out$vectors[, -ind] %*% diag(sqrt(out$values[-ind])) 
  } else {
    R <- diag(sqrt(out$values)) %*% t(out$vectors)
    #L <- out$vectors %*% diag(sqrt(out$values))
  }
      
  out$values[ind] <- 0    
      
  Matrix::Matrix(R)
}

#--------------------
# LogLik functions
#--------------------

#' @export
logLikNum <- function(object) -devCritFun2(object) / 2

devCritFun2 <- function(object, REML) 
{
  if(missing(REML)) {
    REML <- isREML(object)
  }
  cmp <- object@devcomp$cmp
   
  switch(as.character(REML),
    "TRUE" = {
      ## adjust ML results to REML
	    lnum <- log(2*pi*cmp[["pwrss"]])
	    n <- object@devcomp$dims[["n"]]
	    nmp <- n - length(object@beta)
      ldW <- sum(log(weights(object, method = "prior")))
      
      - ldW + cmp[["ldL2"]] + cmp[["ldRX2"]] + nmp*(1 + lnum - log(nmp))
    },
    "FALSE" = {
      ## adjust REML results to ML
      n <- object@devcomp$dims[["n"]]
      lnum <- log(2*pi*cmp[["pwrss"]])
      ldW <- sum(log(weights(object, method = "prior")))
      
      - ldW + cmp[["ldL2"]] + n*(1 + lnum - log(n))
    },
    as.numeric(NA))
}

#--------------------
# Profile functions
#--------------------

#' Profile variance proportions
#'
#' Adapted/copied from \code{varianceProf} / \code{logProf} function in \code{lme4}.
#'
#' @export
varpropProf <- function(x, prof.scale = c("sdcor", "varcov"), 
  ranef = TRUE, na.action = na.omit) 
{
  stopifnot(requireNamespace("splines"))

  prof.scale <- match.arg(prof.scale)

  stopifnot(inherits(x, "thpr"))
  cn <- colnames(x)
  
  if(length(sigs <- grep(paste0("^\\.", if (ranef) "sig" else "sigma"), cn))) {
    sigP <- paste0("^(\\.sig", 
      if (ranef) "(ma)?"
      else "ma", ")")
   
    repP <- "\\1prop"
    colnames(x) <- cn <- sub(sigP, repP, cn)
    levels(x[[".par"]]) <- sub(sigP, repP, levels(x[[".par"]]))
    names(attr(x, "backward")) <- names(attr(x, "forward")) <- 
      sub(sigP, repP, names(attr(x, "forward")))
    
    # convert estimates to proportions    
    x2_sum <- 0
    for(nm in cn[sigs]) {
      x2_sum <- switch(prof.scale,
        "sdcor" = x2_sum + x[[nm]]^2,
        "varcov" = x2_sum + x[[nm]],
        stop("prof.scale"))
    }

    for(nm in cn[sigs]) {
      x[[nm]] <- switch(prof.scale,
        "sdcor" = x[[nm]]^2 / x2_sum,
        "varcov" = x[[nm]] / x2_sum,
        stop("prof.scale"))
    }

    # fit splines
    for(nm in cn[sigs]) {
      fr <- x[x[[".par"]] == nm, TRUE, drop = FALSE]
      form <- eval(substitute(.zeta ~ nm, list(nm = as.name(nm))))
      
      print(form)
      print(head(fr))
      
      attr(x, "forward")[[nm]] <- isp <- splines::interpSpline(form, fr, na.action = na.action)
      attr(x, "backward")[[nm]] <- splines::backSpline(isp)
    }
  }
  
  return(x)
}


#--------------------
# Print functions
#--------------------

#' @export
VarProp <- function(x) 
{
  vf <- VarCorr(x)
  vf <- as.data.frame(vf)
  #vf <- vf[, c("grp", "vcov")]
  
  vf <- within(vf, prop <- vcov / sum(vcov))
  
  return(vf)
}

#--------------------
# ranef functions
#--------------------

#' @export
relmatRanef <- function(x, whichel, ...) 
{
  if(missing(whichel)) {
    stop(paste("Specify a single name of random effect. Available names:", 
      paste(names(ranef(x)), collapse = ", "))) 
  }
  
  pred <- ranef(x, whichel = whichel)[[whichel]]
  stopifnot(ncol(pred) == 1)
  nms <- rownames(pred)
  pred <- pred[, 1]

  # derive the transformed prediction (if necessary)
  relfac <- x@optinfo$relmat$relfac
  
  if(whichel %in% names(relfac)) {
    relfac <- relfac[[whichel]]

    relfac <- relfac[nms, nms] # to make sure names are the same
  
    # see https://github.com/cran/pedigreemm/blob/master/R/pedigree.R#L367 
    pred <- crossprod(pred, relfac) %>% as.numeric
    names(pred) <- nms
  }
  
  return(pred)
}



