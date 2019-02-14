#---------------------
# Function `relmatGlmer`
#---------------------

#' Fit a generalized linear mixed model (GLMM) with relative matrices.
#'
#' Its implementation follows that of `glmer` function in lme4 R package.
#  See https://github.com/lme4/lme4/blob/master/R/lmer.R 
#'
#' @export
relmatGlmer <- function(...)
{
  mc <- match.call()
  
  env <- parent.frame(1)
  env$relmat_glmer <- relmat_glmer
  
  mc[[1]] <- quote(relmat_glmer)
  eval(mc, env)
}

relmat_glmer <- function(formula, data = NULL, family = gaussian,
  control = glmerControl(), start = NULL, 
  verbose = 0L, nAGQ = 1L,
  subset, weights, na.action, offset,
  contrasts = NULL, mustart, etastart, devFunOnly = FALSE, 
  ...,
  # relLmer-specific argument in comparison with lmer function
  method.relfac = "auto",
  debug = 0,
  modOnly = FALSE, controlOnly = FALSE,
  vcControl = list(),
  relmat = list(), check.nobs = "ignore", calc.derivs = TRUE 
)
{
  mc <- mcout <- match.call()
  missCtrl <- missing(control)
  
  if (!inherits(control, "glmerControl")) {
  	if(!is.list(control)) {
  	  stop("'control' is not a list; use glmerControl()")
  	}
  	
  	## back-compatibility kluge
  	msg <- "Use control=glmerControl(..) instead of passing a list"
  	if(length(cl <- class(control))) {
  	  msg <- paste(msg, "of class", dQuote(cl[1]))
  	}
  	warning(msg, immediate. = TRUE)
  	
  	control <- do.call(glmerControl, control)
  }

  ## relmatGlmer-specific part
  #relnms <- names(relmat)
  #stopifnot(all(relnms %in% names(data)))
  #relcls <- laply(relnms, function(x) class(data[[x]]))
  #if(!all(relcls == "factor")) {
  #  stop("all variables related to 'relmat' to be factors")
  #}
  
  if(check.nobs == "ignore") {
    control$checkControl$check.nobs.vs.rankZ <- "ignore"
    control$checkControl$check.nobs.vs.nlev <- "ignore"
    control$checkControl$check.nobs.vs.nRE <- "ignore"
  }
  control$calc.derivs <- calc.derivs
    
  ## family-checking code duplicated here and in glFormula (for now) since
  ## we really need to redirect at this point; eventually deprecate formally
  ## and clean up
  if(is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame(2))
  }
  if(is.function(family)) {
    family <- family()
  }
  if (isTRUE(all.equal(family, gaussian()))) {
    ## redirect to lmer (with warning)
    warning("calling glmer() with family=gaussian (identity link) as a shortcut to lmer() is deprecated;",
      " please call lmer() directly")
    
    mc[[1]] <- quote(lme4::lmer)
    mc["family"] <- NULL            # to avoid an infinite loop
    
    return(eval(mc, parent.frame()))
  }
  
  if(controlOnly) return(control)
  
  mc$control <- control ## update for back-compatibility kluge
  ## https://github.com/lme4/lme4/issues/50

  ## parse relmat, etc arguments
  nrelmat <- length(relmat)
  
  ## see https://github.com/lme4/lme4/issues/50
  ## parse the formula and data
  mc[[1]] <- quote(lme4::glFormula)
  ##glmod <- eval(mc, parent.frame(1L)) # Commented line 1
  ##mcout$formula <- glmod$formula # Commented line 2
  ##glmod$formula <- NULL # Commented line 3

  #----
  # Start of relmatLmer-specific part
  #----
  # remove call arguments which are specific to `relmat_glmer`
  # before calling lme4::lFormula
  args <- c("method.relfac", "modOnly", "controlOnly", "debug", "vcControl", 
    "relmat", "check.nobs", "calc.derivs")
  for(arg in args) {
    if(arg %in% names(mc)) {
      mc[[arg]] <- NULL
    }
  }
  glmod <- eval(mc, parent.frame(1L)) 

  stopifnot(is.list(relmat), length(names(relmat)) == length(relmat))
  relfac <- list()
  
  # list of random term factor
  flist <- glmod$reTrms[["flist"]]
  fnmns <- names(flist) # e.g. (1|gr) + (1|id) 
  
  # list of random term factors with custom covariances
  # - specified in `relmat` argument
  # - need special processing
  relnms <- names(relmat)

  if(any(fnmns %in% relnms)) {
    # random-effects design matrix components
    Ztlist <- glmod$reTrms[["Ztlist"]]
    
    ind <- which(fnmns %in% relnms)
    for(i in ind) {
      fn <- fnmns[i]
      
      if(debug) {
        cat(" - Zlist element", i, "/", length(fnmns), ":", fn, "\n")
      }

      zn <- glmod$fr[, fn]
      
      if(class(zn) != "factor") {
        zn <- as.factor(zn)
      }
      zn.unique <- levels(zn) 
      
      stopifnot(!is.null(rownames(relmat[[fn]])))
      rn <- rownames(relmat[[fn]])

      stopifnot(all(zn.unique %in% rn))
      
      mat <- Matrix::Matrix(relmat[[fn]][zn.unique, zn.unique], sparse = TRUE)
      relfac[[fn]] <- relfac(mat, method.relfac)
      
      # ?Matrix::chol
      # Returned value: a matrix of class Cholesky, i.e., upper triangular: R such that R'R = x.
      # Note that another notation is equivalent x = L L', where L is a lower triangular 
      # @ http://en.wikipedia.org/wiki/Cholesky_decomposition

      # If the substitution is Z* = Z L, then Z*' = L' Z' = R Z'
      # @ http://www.journalofanimalscience.org/content/88/2/497.long%5BVazquez%20et%20al.,%202010%5D
      #Ztlist[[i]] <-  relfac[[i]] %*% Ztlist[[i]]

      if(debug) {
        cat(" - relfac[[fn]]:", nrow(relfac[[fn]]), "x", ncol(relfac[[fn]]), "\n")
        cat(" - Ztlist[[i]]:", nrow(Ztlist[[i]]), "x", ncol(Ztlist[[i]]), "\n")
        print(head(rownames(Ztlist[[i]])))
      }
      
      ncol_relac <- ncol(relfac[[fn]])
      nrow_Zt <- nrow(Ztlist[[i]])
      ncol_Zt <- ncol(Ztlist[[i]])
      
      # scenario 1
      if(ncol_relac == nrow_Zt) {
        Ztlist[[i]] <- relfac[[fn]] %*% Ztlist[[i]]
      } else if(nrow_Zt > ncol_Zt) {
        # scenario 2
        stopifnot(ncol_Zt == ncol_relac)
        
        if((nrow_Zt %% ncol_Zt) == 0) {
          ngr <- length(zn.unique)
          if(debug) {
            cat(" - ngr:", ngr, "\n")
          }    
          # check: sequence `zn.unique` exatcly match `zn` inside a block
          nblocks.x <- length(zn) / length(zn.unique)
          for(j in 1:nblocks.x) {
            jnd <- ngr * (j - 1) + seq(1, ngr)
                  
            stopifnot(all(zn[jnd] %in% zn.unique))
          }
      
          nblocks.y <- nrow(Ztlist[[i]]) / ncol(Ztlist[[i]])
          nblocks <- nblocks.x * nblocks.y
          stopifnot(nblocks != 1)
  
          if(debug) {
            cat(" - nblocks.x:", nblocks.x, "\n")
            cat(" - nblocks.y:", nblocks.y, "\n")
            cat(" - nblocks:", nblocks, "\n")
          }  
          
          # divide by blocks
          Ztlisti <- list()
              
          for(j in 1:nblocks.y) {
            if(debug) {
              cat("  --  block", j, "/", nblocks.y, "\n")
            }
          
            jnd <- seq(j, by = nblocks, length = ngr)
            Ztlisti[[j]] <- list()
    
            for(k in 1:nblocks.x) {
              knd <- ngr * (k - 1) + seq(1, ngr)
      
              Ztlisti[[j]][[k]] <- Ztlist[[i]][jnd, knd]
            }
          }
        
          # update per block
          for(j in 1:nblocks.y) {
            jnd <- seq(j,  by = nblocks, length = ngr)
    
            for(k in 1:nblocks.x) {
              knd <- ngr * (k - 1) + seq(1, ngr)
       
              if(debug) {
                cat(" - relfac[[fn]]:", nrow(relfac[[fn]]), "x", ncol(relfac[[fn]]), "\n")
                cat(" - Ztlisti[[j]][[k]]:", nrow(Ztlisti[[j]][[k]]), "x", 
                  ncol(Ztlisti[[j]][[k]]), "\n")
                print(head(rownames(Ztlisti[[j]][[k]])))
              }
      
              Ztlisti[[j]][[k]] <- relfac[[fn]] %*% Ztlisti[[j]][[k]]
            }
          }
  
          for(j in 1:nblocks.y) {
            jnd <- seq(j,  by = nblocks, length = ngr)
    
            for(k in 1:nblocks.x) {
              knd <- ngr * (k - 1) + seq(1, ngr)
      
              Ztlist[[i]][jnd, knd] <- Ztlisti[[j]][[k]]
            }
          }
        } else {
          stop("nrow_Zt > ncol_Zt, but (nrow_Zt %% ncol_Zt) is not zero") 
        }
      } else {
        stop("nrow_Zt < ncol_Zt")
      }
    }
  
    glmod$reTrms[["Ztlist"]] <- Ztlist
    glmod$reTrms[["Zt"]] <- do.call(rBind, Ztlist)
  }
  #-------------------------------
  # end of relmatGlmer-specific code
  #------------------------------- 
  
  mcout$formula <- glmod$formula # Uncommented line 2
  glmod$formula <- NULL # Uncommented line 3
  
  if(modOnly) return(glmod)
  
  ## create deviance function for covariance parameters (theta)
  nAGQinit <- if(control$nAGQ0initStep) 0L else 1L
  devfun <- do.call(mkGlmerDevfun, c(glmod, 
    list(verbose = verbose, control = control, nAGQ = nAGQinit)))
    
  if (nAGQ==0 && devFunOnly) return(devfun)

  ## optimize deviance function over covariance parameters
  ## FIXME: perhaps should be in glFormula instead??
  if(is.list(start)) {
    start.bad <- setdiff(names(start),c("theta","fixef"))
    if(length(start.bad) > 0) {
      stop(sprintf("bad name(s) for start vector (%s); should be %s and/or %s",
        paste(start.bad, collapse = ", "),
        shQuote("theta"),
        shQuote("fixef")), call. = FALSE)
    }
    if (!is.null(start$fixef) && nAGQ == 0)
      stop("should not specify both start$fixef and nAGQ==0")
  }
  
  ## FIX ME: allow calc.derivs, use.last.params etc. if nAGQ=0
  if(control$nAGQ0initStep) {
    opt <- optimizeGlmer(devfun,
      optimizer = control$optimizer[[1]],
      ## DON'T try fancy edge tricks unless nAGQ=0 explicitly set
      restart_edge=if (nAGQ==0) control$restart_edge else FALSE,
      boundary.tol=if (nAGQ==0) control$boundary.tol else 0,
      control = control$optCtrl,
      start=start,
      nAGQ = 0,
      verbose=verbose,
      calc.derivs=FALSE)
  }

  if(nAGQ > 0L) {
    start <- lme4:::updateStart(start,theta = opt$par)

    ## update deviance function to include fixed effects as inputs
    devfun <- updateGlmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ)

    if (devFunOnly) return(devfun)
    ## reoptimize deviance function over covariance parameters and fixed effects
    opt <- optimizeGlmer(devfun,
      optimizer = control$optimizer[[2]],
      restart_edge=control$restart_edge,
      boundary.tol=control$boundary.tol,
      control = control$optCtrl,
      start=start,
      nAGQ=nAGQ,
      verbose = verbose,
      stage=2,
      calc.derivs=control$calc.derivs,
      use.last.params=control$use.last.params)
  }

  cc <- if (!control$calc.derivs) NULL else {
    if (verbose > 10) cat("checking convergence\n")

    lme4:::checkConv(attr(opt,"derivs"),opt$par,
      ctrl = control$checkConv,
      lbound=environment(devfun)$lower)
  }

  ## prepare output
  mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr, 
    mc = mcout, lme4conv=cc)      
}

