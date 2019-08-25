#---------------------
# Function `relmatLmer`
#---------------------

#' Fit a linear mixed model (LMM) with relative matrices.
#'
#' Its implementation follows that of lmer function in lme4 R package.
#'
#' @export
relmatLmer <- function(...)
{
  mc <- match.call()
  
  env <- parent.frame(1)
  env$relmat_lmer <- relmat_lmer
  
  mc[[1]] <- quote(relmat_lmer)
  eval(mc, env)
}

relmat_lmer <- function(formula, data = NULL, REML = TRUE,
  control = lmerControl(), start = NULL,
  verbose = 0L, subset, weights, na.action, offset,
  contrasts = NULL, devFunOnly = FALSE, 
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
  
  ## see functions in modular.R for the body ...
  if(!missCtrl && !inherits(control, "lmerControl")) {
    if(!is.list(control)) {
      stop("'control' is not a list; use lmerControl()")
    }
  
    ## back-compatibility kluge
    warning("passing control as list is deprecated: please use lmerControl() instead", immediate. = TRUE)
    control <- do.call(lmerControl, control)
  }
  
  #----
  # relmatLmer-specific part
  #----
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
  
  if (!is.null(list(...)[["family"]])) {
    stop("calling lmer with 'family' is deprecated; please use glmer() instead")
  }
  
  if(controlOnly) return(control)
  
  mc$control <- control ## update for back-compatibility kluge
  ## https://github.com/lme4/lme4/issues/50

  ## parse relmat, etc arguments
  nrelmat <- length(relmat)
  
  ## parse data and formula
  mc[[1]] <- quote(lme4::lFormula)
  #lmod <- eval(mc, parent.frame(1L)) ## parse data and formula
  
  #----
  # relmatLmer-specific part
  #----
  # remove call arguments which are specific to `relmat_lmer` and to be removed before calling lme4::lFormula
  args <- c("method.relfac", "modOnly", "controlOnly", "debug", "vcControl", 
    "relmat", "check.nobs", "calc.derivs")
  for(arg in args) {
    if(arg %in% names(mc)) {
      mc[[arg]] <- NULL
    }
  }

  # @ https://github.com/lme4/lme4/commit/21dfb0114562a36c59ab6b90e4dde7311c7e4c70
  #lmod <- eval(mc, environment(formula)) ## parse data and formula
  lmod <- eval(mc, parent.frame(1L))
  
  #-------------------------------
  # start of solaris-specific code
  #-------------------------------
  stopifnot(is.list(relmat), length(names(relmat)) == length(relmat))
  relfac <- list()
  
  # list of random term factor
  flist <- lmod$reTrms[["flist"]]
  fnmns <- names(flist) # e.g. (1|gr) + (1|id) 
  
  # list of random term factors with custom covariances
  # - specified in `relmat` argument
  # - need special processing
  relnms <- names(relmat)

   for(i in seq_along(fnmns)) {
    fn = fnmns[i]
    if(fn %in% relnms) {
      if(debug) {
        cat((" - Zlist element", fn, "\n")
      }
      
      # check names
      zn <- lmod$fr[, fn]
      if(class(zn) != "factor") {
        zn <- as.factor(zn)
      }
      zn.unique <- levels(zn)

      stopifnot(!is.null(rownames(relmat[[fn]])))
      rn <- rownames(relmat[[fn]])

      stopifnot(all(zn.unique %in% rn))
      
      # compute a relative factor R: K = R'R
      # See lme4qtl:::relfac
      K <- Matrix::Matrix(relmat[[fn]][zn.unique, zn.unique], sparse = TRUE)
      R <- relfac(Ki, method.relfac)
      relfac[[fn]] <- R

      # compute a substitution Z*
      # 1. for random effects ~(1|f): Z* = R Z'
      # 2. for random effects ~(r | f)
      #   r is used to construct the "raw random effects model matrix Xi (n x pi) for term i
      #   f is used to for the grouping factor model matrix Ji (n x qi) for a factor with qi levels.
      #       These levels correspond to levels(flist[[fnmn]]) and are used to index the rows/columns of Ki
      # Following Bates, et al "Fitting Linear Mixed-Effects Models Using lme4",
      # these are combined into the matrix Zi = t(KhatriRao(t(Ji), t(Xi)))
      # 
      # We want to modify Zi to: t(KhatriRao(t(Ji %*% t(R)), t(Xi))), 
      # where t(R) %*% R = K, eg. Cholesky
      #
      # This is equal to Zi %*% kronecker(t(K_sqrt),diag(1,pi))
      #
      # See github.com/variani/lme4qtl/pull/18 contributed by github.com/deruncie
      
      pi <- length(lmod$reTrms$cnms[[i]])
      Zi_t <- lmod$reTrms$Ztlist[[i]] 
      Zi_t <- kronecker(Ki_sqrt, diag(1, pi)) %*% Zi_t # t(Z*)

      # put the new t(Z*) back into the appropriate slot `Ztlist`
      lmod$reTrms$Ztlist[[i]] <- Zi_t
    }
  }
  # update the full Zt matrix (the slot `Zt`) by combining all Zt matrices (the slot `Ztlist`)
  lmod$reTrms[["Zt"]] <- do.call(rBind, lmod$reTrms$Ztlist)
  
  #-------------------------------
  # end of relmatLmer-specific code
  #------------------------------- 
    
  mcout$formula <- lmod$formula
  lmod$formula <- NULL
  
  if(modOnly) return(lmod)
  
  ## create deviance function for covariance parameters (theta)
  devfun <- do.call(mkLmerDevfun, c(lmod,
    list(start = start, verbose = verbose, control = control)))
  # print(with(environment(devfun), pp$theta))
  
  #----
  # vcControl-specific update
  #----  
  if(length(vcControl)) {
    # local update function  
    update_theta_vcControl <- function(theta, vcControl)
    {
      nth <- length(theta)
    
      # vcControl.vareq 
      # Example: c(1, 2, 3)
      vcControl.vareq <- vcControl[["vareq"]]  
      if(length(vcControl.vareq)) {
        for(i in 1:length(vcControl.vareq)) {
          ind <- vcControl.vareq[[i]]
          stopifnot(length(ind) == 3)
          stopifnot(all(ind <= nth))
          
          # t1^2 = t12^2  + t2^2
          theta[ind[1]] <- sqrt(theta[ind[2]]*theta[ind[2]] + theta[ind[3]]*theta[ind[3]])
        }
      }
      
      # vcContol.rho1
      # Example: 3
      vcControl.rho1 <- vcControl[["rho1"]]
      if(length(vcControl.rho1)) {
        for(i in 1:length(vcControl.rho1)) {
          ind <- vcControl.rho1[[i]]
          stopifnot(length(ind) == 1)
          stopifnot(all(ind <= nth))
          
          # t2 = 0          
          theta[ind[1]] <- 0
        }
      }

      # vcContol.rho0
      # Example: 2      
      vcControl.rho0 <- vcControl[["rho0"]]
      if(length(vcControl.rho0)) {
        for(i in 1:length(vcControl.rho0)) {
          ind <- vcControl.rho0[[i]]
          stopifnot(length(ind) == 1)
          stopifnot(all(ind <= nth))
          
          # t12 = 0   
          theta[ind[1]] <- 0
        }
      }
      
      # vcContol.rho
      # Example: 2, 3      
      vcControl.rho <- vcControl[["rho"]]
      if(length(vcControl.rho)) {
        for(i in 1:length(vcControl.rho)) {
          ind <- vcControl.rho[[i]][[1]]

          stopifnot(length(ind) == 2)
          stopifnot(all(ind <= nth))
          
          val <- vcControl.rho[[i]][[2]]
          stopifnot(val != 0)
          stopifnot(val > -1)
          stopifnot(val < 1)
          
          k <- val / sqrt(1 - val*val)

          # t12 = rho / sqrt(1 - rho^2) t2   
          theta[ind[1]] <- k * theta[ind[2]]
        }
      }
      
      # vcContol.th1
      # Example: 1      
      vcControl.th1 <- vcControl[["th1"]]
      if(length(vcControl.th1)) {
        for(i in 1:length(vcControl.th1)) {
          ind <- vcControl.th1[[i]]
          stopifnot(length(ind) == 1)
          stopifnot(all(ind <= nth))
          
          # t1 = 1   
          theta[ind[1]] <- 1
        }
      }

      # vcContol.th0
      # Example: 0    
      vcControl.th0 <- vcControl[["th0"]]
      if(length(vcControl.th0)) {
        for(i in 1:length(vcControl.th0)) {
          ind <- vcControl.th0[[i]]
          stopifnot(length(ind) == 1)
          stopifnot(all(ind <= nth))
          
          # t1 = 0   
          theta[ind[1]] <- 0
        }
      }            

      # vcControl.hom3 
      # Example: c(1, 2, 3, 4, 5, 6)
      vcControl.hom3 <- vcControl[["hom3"]]  
      if(length(vcControl.hom3)) {
        for(i in 1:length(vcControl.hom3)) {
          ind <- vcControl.hom3[[i]]
          stopifnot(length(ind) == 6)
          stopifnot(all(ind <= nth))
          
          # t3 = t2
          theta[ind[3]] <- theta[ind[2]]
          # t4 = sqrt(t1^2 - t2^2)
          stopifnot(abs(theta[ind[1]]) >= abs(theta[ind[2]]))
          theta[ind[4]] <- sqrt(theta[ind[1]]*theta[ind[1]] - 
            theta[ind[2]]*theta[ind[2]])
          # t5 = (t1*t2 - t2^2) / t4
          theta[ind[5]] <- (theta[ind[1]]*theta[ind[2]] - 
            theta[ind[2]]*theta[ind[2]]) / theta[ind[4]]
          # t6 = sqrt(t4^2 - t5^2)
          stopifnot(abs(theta[ind[4]]) >= abs(theta[ind[5]]))
          theta[ind[6]] <- sqrt(theta[ind[4]]*theta[ind[4]] - 
            theta[ind[5]]*theta[ind[5]])
        }
      }
   
         
      return(theta)
    }

    devfun.wrapper <- function(theta) {
      th <- theta
      th <- update_theta_vcControl(th, vcControl)
      
      devfun(th)
    }
    environment(devfun.wrapper) <- environment(devfun)
    assign("devfun", devfun, envir = environment(devfun.wrapper))
    assign("vcControl", vcControl, envir = environment(devfun.wrapper))
    assign("update_theta_vcControl", update_theta_vcControl, envir = environment(devfun.wrapper))
  }
  #----
  # end of vcControl-specific update
  #----  
  
  if(devFunOnly) return(devfun)
  
  ## optimize deviance function over covariance parameters
  if (control$optimizer=="none") {
    #----
    # relmatLmer-specific update
    #----
    if(is.null(start)) {
      opt <- list(par = NA, fval = NA, conv = 1000, message = "no optimization")
    }
    else if("theta" %in% names(start)) {
      opt <- list(par = as.numeric(start$theta), fval = NA, conv = 1000, message = "no optimization + start")
      if(!is.null(start)) {
        if("theta" %in% names(start)) {
          theta <- start$theta
          #devfun(pp$theta <- theta)
          #assign("pp$theta", theta)
          devfun(theta)
        }
      }
    }
  } else {
    if(length(vcControl)) {
      opt <- optimizeLmer(devfun.wrapper,
        optimizer=control$optimizer,
        restart_edge=control$restart_edge,
        boundary.tol=control$boundary.tol,
        control=control$optCtrl,
        verbose=verbose,
        start=start,
        calc.derivs=control$calc.derivs,
        use.last.params=control$use.last.params)    
    } else {
      opt <- optimizeLmer(devfun,
        optimizer=control$optimizer,
        restart_edge=control$restart_edge,
        boundary.tol=control$boundary.tol,
        control=control$optCtrl,
        verbose=verbose,
        start=start,
        calc.derivs=control$calc.derivs,
        use.last.params=control$use.last.params)
        
      cc <- lme4:::checkConv(attr(opt,"derivs"), opt$par,
        ctrl = control$checkConv,
        lbound = environment(devfun)$lower)
    }
  }

  if(length(vcControl)) {
    th <- opt$par
    th <- update_theta_vcControl(th, vcControl)
    devfun(opt$par <- th)

    #mkMerMod(environment(devfun), opt.eq, lmod$reTrms, fr = lmod$fr)
    
    opt <- list(par = th, feval = opt$feval, fval = NA, conv = opt$convergence, message = "no optimization + par from devfun.wrapper")
    mod <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr, mcout) ## prepare output
  } else {
    cc <- lme4:::checkConv(attr(opt,"derivs"), opt$par,
      ctrl = control$checkConv,
      lbound=environment(devfun)$lower)
  
    mod <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr, mcout, lme4conv = cc) ## prepare output
  }
  
  # adhoc solution to store `relmatLmer` variables in an object of `lmerMod` class 
  mod@optinfo$relmat <- list(relfac = relfac)
  
  return(mod)
}

relmatLmer_naive <- function(formula, data = NULL, 
  start = NULL,
  relmat = list()
)
{
  # formula
  control <- lmerControl(check.nobs.vs.rankZ = "ignore", 
    check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
  
  # lmod
  lmod <- lFormula(formula, data, control = control)
  
  #-------------------------------
  # start of relmatLmer-specific code
  #-------------------------------
  stopifnot(is.list(relmat), length(names(relmat)) == length(relmat))
  relnms <- names(relmat)
  relfac <- relmat
  flist <- lmod$reTrms[["flist"]]   ## list of factors

  ind <- (relnms %in% names(flist))
  if(any(ind)) {
    ## random-effects design matrix components
    Ztlist <- lmod$reTrms[["Ztlist"]]
    
    asgn <- attr(flist, "assign")
    for(i in seq_along(relnms[ind])) {
      
      relmati <- relnms[ind][i]
      if(!(relmati %in% names(flist))) {
        stop("a relationship matrix must be (", relmati, ")",
          " associated with only one random effects term (", paste(names(flist), collapse = ", "), ")")
      }
      tn <- which(relmati == names(flist))
      fn <- names(flist)[tn]
      
      zn <- lmod$fr[, fn]
      if(class(rownames(relmat[[i]])) == "character") {
        zn <- as.character(zn)
      }
      
      relmat[[i]] <- Matrix::Matrix(relmat[[i]][zn,zn], sparse = TRUE)
      stopifnot(nrow(relmat[[i]]) == nrow(lmod$fr))
    
      relfac[[i]] <- Matrix::chol(relmat[[i]])
      Ztlist[[i]] <-  relfac[[i]] %*% Ztlist[[i]]
    }
    lmod$reTrms[["Ztlist"]] <- Ztlist
  }
  lmod$reTrms[["Zt"]] <- do.call(rBind, Ztlist)
  #-------------------------------
  # end of relmatLmer-specific code
  #-------------------------------
  
  devfun <- do.call(mkLmerDevfun, c(lmod,
    list(start = start)))#, verbose = verbose, control = control)))
  devfun(opt$par <- start)
  
  opt <- list(par = as.numeric(start$theta), fval = NA,conv = 1000, message="start copied")
  mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
}
