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

  if(!all(relnms %in% fnmns)) {
    warning("not all relmat ID variables match RE ID terms")
  }
  
  for(i in seq_along(fnmns)) {
    fn <- fnmns[i]
    
    if(fn %in% relnms) { # update only when matching; otherwise, ignore
      if(debug) {
        cat(" - Zlist element", fn, "\n")
      }
      
      # check names
      zn <- glmod$fr[, fn]
      if(class(zn) != "factor") { # convert to factor
        zn <- as.factor(zn) 
      }
      zn.unique <- levels(zn)

      stopifnot(!is.null(rownames(relmat[[fn]])))
      rn <- rownames(relmat[[fn]])

      stopifnot(all(zn.unique %in% rn))
      
      # compute a relative factor R: K = R'R
      # See lme4qtl:::relfac
      K <- Matrix::Matrix(relmat[[fn]][zn.unique, zn.unique], sparse = TRUE)
      R <- relfac(K, method.relfac)
      relfac[[fn]] <- R

      # compute a substitution Z*
      # See similar code comments in lme4qtl::relmatLmer      
      pi <- length(glmod$reTrms$cnms[[i]])
      Zi_t <- glmod$reTrms$Ztlist[[i]] 
      Zi_t <- kronecker(R, diag(1, pi)) %*% Zi_t # t(Z*)

      # put the new t(Z*) back into the appropriate slot `Ztlist`
      glmod$reTrms$Ztlist[[i]] <- Zi_t
    }
  }
  # update the full Zt matrix (the slot `Zt`) by combining all Zt matrices (the slot `Ztlist`)
  glmod$reTrms[["Zt"]] <- do.call(rbind, glmod$reTrms$Ztlist)
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

