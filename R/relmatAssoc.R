#---------------------
# Function `relmatAssoc`
#---------------------

#' @export
assocLmer <- function(...)
{
  mc <- match.call()
  
  env <- parent.frame(1)
  env$assoc_lmer <- assoc_lmer
    
  mc[[1]] <- quote(assoc_lmer)
  eval(mc, env)
}

#-------------------------
# assoc_lmer functions
#-------------------------

assoc_lmer <- function(formula, data = NULL, REML = TRUE, 
  control = lmerControl(),
  calc.derivs = FALSE, check.scaleX = "stop", optimizer = "bobyqa",
  ...,
  data_snpcov, snps,
  batch_size = 500,
  var_id = "id",
  method = c("Wald", "LR"),
  cores = getOption("cores"),
  verbose = 1)
{
  ### call
  mc <- match.call()
  
  ### args
  missing_snps <- missing(snps)
  cores <- ifelse(is.null(cores), 1, cores)

  method <- match.arg(method)
  REML <- ifelse(method == "LR", FALSE, REML)
  
   # lmerControl
  if(missing(control)) {
    control <- lmerControl(calc.derivs = calc.derivs, check.scaleX = check.scaleX, optimizer = optimizer)
  } else {
    control$checkControl$check.scaleX <- check.scaleX
    control$calc.derivs <- calc.derivs
    control$optimizer <- optimizer
  }
  
  ### parallel
  parallel <- (cores > 1)
  if(parallel) {
    stopifnot(requireNamespace("doParallel", quietly = TRUE))
    doParallel::registerDoParallel(cores = cores)
  }  
  
  ### check requirements to `data` & `data_snpcov`:
  # - both tables have rows in the same order by `var_id`
  # - `var_id` is in the 1st column of `data_snpcov`
  stopifnot(all(data[[var_id]] == data_snpcov[[var_id]]))
  stopifnot(names(data_snpcov)[1] == var_id)

  # drop the 1st column in `data_snpcov`
  data_snpcov <- data_snpcov[-1]
  
  if(missing_snps) {
    snps <- names(data_snpcov)
  }  

  ### prepare batches of `snps`
  beg <- seq(1, length(snps), by = batch_size)
  end <- c(beg[-1] - 1, length(snps))
  
  num_batches <- length(beg)
  
  tab <- ldply(seq(1, num_batches), function(i) {
    if(verbose) {
      cat(" * assoc_lmer:", i, "/", num_batches, "batch\n")
    }

    ind <- seq(beg[i], end[i])
    snps_assoc <- snps[ind]

    data_assoc <- cbind(data, subset(data_snpcov, select = snps_assoc))
    
    # polygenic model (the null)
    mod0 <- relmatLmer(formula, data_assoc, REML = REML, ...)
    start0 <- list(theta = getME(mod0, c("theta")))
    
    n0 <- nobs(mod0)
    p0 <- length(fixef(mod0))
            
    # compute the association models in a loop
    tab_i <- ldply(snps_assoc, function(x) {
      tab_na <- data_frame(snp = x, estimate = as.numeric(NA), se = as.numeric(NA),
        stat = as.numeric(NA), pval = as.numeric(NA))

      mod <- try({
        # update
        #update(mod0, paste(". ~ . +", x), start = start0)
        # relmatLmer
        relmatLmer(update(formula, paste(". ~ . +", x)), data_assoc, REML = REML, start = start0, ...)
      })

      if(class(mod)[1] == "try-error") {
        return(tab_na)
      }
      
      n <- nobs(mod)
      p <- length(fixef(mod))
      
      if(p == p0) { # snp's column was dropped due to rank deficiency
        return(tab_na)
      }

      estimate <- as.numeric(fixef(mod)[x])
      se2 <- vcov(mod)[x, x]

      if(method == "Wald") {
        stat <- estimate * estimate / se2
        pval <- F_test(stat, 1, n - p)
      } 
      else if(method == "LR") {
        stopifnot(!REML)
        
        if(n0 != n) {
          return(tab_na)
        }
        
        logLik0 <- as.numeric(logLik(mod0))
        logLik <- as.numeric(logLik(mod))

        stat <- 2 * (logLik - logLik0) # chi2
        pval <- pchisq(stat, 1, lower.tail = FALSE)
      }
      
      tab <- data_frame(snp = x, estimate = estimate, se = sqrt(se2),
        stat = stat, pval = pval)
      
      return(tab)  
    })
  }, .parallel = parallel)

  ### return
  tab <- as_data_frame(tab)
  
  return(tab)
  

}

assoc_lmer_wald <- function(formula, data = NULL, REML = TRUE,
  ...,
  data_snpcov, snps,
  batch_size = 500,
  cores = getOption("cores"),
  verbose = 1
)
{
  ### call
  mc <- match.call()
  
  ### args
  missing_snps <- missing(snps)
  
  # cores
  if(is.null(cores)) {  
    cores <- 1
  }

  ### parallel
  parallel <- (cores > 1)
  if(parallel) {
    # load required R package doParallel
    stopifnot(requireNamespace("doParallel", quietly = TRUE))
    
    doParallel::registerDoParallel(cores = cores)
  }  
  
  ### variables
  var_id <- "id"
  ids_phen <- data[[var_id]]
  ids_gen <- data_snpcov[[var_id]]

  ids_matched <- all(ids_phen == ids_gen)
  stopifnot(ids_matched)

  mode_combine <- "unknown"
  if(ids_matched) {
    mode_combine <- "cbind"
  }
  
  if(missing_snps) {
    snps <- colnames(data_snpcov)
    
    # remove `id` variable
    if(snps[1] == var_id) {
      snps <- snps[-1]
    } else {
      ind_id <- which(snps == var_id)
      stopifnot(length(ind_id) == 1)
      
      snps <- snps[-ind]
    }
  }

  ### prepare batches of `snps`
  beg <- seq(1, length(snps), by = batch_size)
  end <- c(beg[-1] - 1, length(snps))
  
  num_batches <- length(beg)
  
  tab <- ldply(seq(1, num_batches), function(i) {
    if(verbose) {
      cat(" * assoc_lmer:", i, "/", num_batches, "batch\n")
    }

    ind <- seq(beg[i], end[i])
    snps_i <- snps[ind]

    ### combine data and run polygenic (!) per batch
    data_snpcov_i <- subset(data_snpcov, select = c(var_id, snps_i))
    
    if(mode_combine == "cbind") {
      data_assoc <- cbind(data, data_snpcov_i[-1])
    }
    else {
      stop()
    }      
  
    mod0 <- relmatLmer(formula, data_assoc, REML = REML, ...)

    n0 <- nobs(mod0)
    p0 <- length(fixef(mod0))
            
    ### compute the association models
    tab_i <- ldply(snps_i, function(x) {
      mod <- try({
        update(mod0, paste(". ~ . +", x))
      })
    
      tab <- data.frame(snp = x, pval = as.numeric(NA))
      if(class(mod)[1] != "try-error") {
        pval <- try({
          # function to compute F-test via beta distr.
          ftest <- function(stat, v1, v2) pbeta(v2 / (v2 + v1 * stat), v2/2, v1/2)    
        
          # copmute F-statistics        
          b <- as.numeric(fixef(mod)[x])
          v <- vcov(mod)[x, x]
          stat <- b^2 / v
        
          ftest(stat, 1, n0 - p0 - 1)
        })

        if(class(pval)[1] != "try-error") {
          tab <- data.frame(snp = x, pval = pval)
        }
      }
      
      return(tab)  
    })
  }, .parallel = parallel)

  ### return
  tab <- as_data_frame(tab)
  
  return(tab)
}

assoc_lmer_chisq <- function(formula, data = NULL, REML = FALSE,
  ...,
  data_snpcov, snps,
  format_snpcov = c("012", "01", "-11"),
  batch_size = 500,
  cores = getOption("cores"),
  verbose = 1
)
{
  ### call
  mc <- match.call()
  
  ### args
  missing_snps <- missing(snps)
  
  format_snpcov <- match.arg(format_snpcov)
  
  # cores
  if(is.null(cores)) {  
    cores <- 1
  }

  ### parallel
  parallel <- (cores > 1)
  if(parallel) {
    # load required R package doParallel
    stopifnot(requireNamespace("doParallel", quietly = TRUE))
    
    doParallel::registerDoParallel(cores = cores)
  }  
  
  ### variables
  var_id <- "id"
  ids_phen <- data[[var_id]]
  ids_gen <- data_snpcov[[var_id]]

  ids_matched <- all(ids_phen == ids_gen)
  stopifnot(ids_matched)

  mode_combine <- "unknown"
  if(ids_matched) {
    mode_combine <- "cbind"
  }
  
  if(missing_snps) {
    snps <- colnames(data_snpcov)
    
    # remove `id` variable
    if(snps[1] == var_id) {
      snps <- snps[-1]
    } else {
      ind_id <- which(snps == var_id)
      stopifnot(length(ind_id) == 1)
      
      snps <- snps[-ind]
    }
  }

  ### prepare batches of `snps`
  beg <- seq(1, length(snps), by = batch_size)
  end <- c(beg[-1] - 1, length(snps))
  
  num_batches <- length(beg)
  
  tab <- ldply(seq(1, num_batches), function(i) {
    if(verbose) {
      cat(" * assoc_lmer:", i, "/", num_batches, "batch\n")
    }

    ind <- seq(beg[i], end[i])
    snps_i <- snps[ind]

    ### combine data and run polygenic (!) per batch
    data_snpcov_i <- subset(data_snpcov, select = c(var_id, snps_i))
    
    if(mode_combine == "cbind") {
      data_assoc <- cbind(data, data_snpcov_i[-1])
    }
    else {
      stop()
    }      
  
    mod0 <- relmatLmer(formula, data_assoc, REML = REML, ...)
    start0 <- list(theta = getME(mod0, c("theta")))
    
    ### compute mac/maf per SNP using ids from the polygenic model `mod0`
    ids_model <- as.character(mod0@frame[, var_id])
    ids_ind <- which(data_snpcov_i[[var_id]] %in% ids_model)
  
    maf_snpcov <- function(x, format) 
    {
      f <- switch(format,
        "012" = sum(x) / (2 * length(x)),
        "01" = sum(x) / length(x),
        "-11" = sum(x + 1) / (2 * length(x)),
        stop())
      
      min(f, 1 - f)
    }
    
    mac_snpcov <- function(maf, len, format)
    {
      switch(format,
        "012" = maf * (2 * len),
        "01" = maf * len,
        "-11" = maf * len,
        stop())
    }
    
    maf <- laply(snps_i, function(i) {
      x <- data_snpcov_i[[i]]
      maf_snpcov(x[ids_ind], format_snpcov)
    })
    
    mac <- mac_snpcov(maf, length(ids_model), format_snpcov)
    
    ### compute the association models
    tab_i <- ldply(snps_i, function(x) {
      mod <- try({
        # update
        #update(mod0, paste(". ~ . +", x), start = start0)
        # relmatLmer
        relmatLmer(update(formula, paste(". ~ . +", x)), data_assoc, REML = REML, start = start0, ...)
      })
    
      tab <- data.frame(snp = x, pval = as.numeric(NA))
      if(class(mod)[1] != "try-error") {
        stats <- anova(mod, mod0)
        tab <- as.data.frame(stats)
        pval <- tab["mod", "Pr(>Chisq)"]
        
        tab <- data.frame(snp = x, pval = pval)
      }
      
      return(tab)  
    })
    
    tab_i <- data.frame(tab_i, mac = mac, maf = maf)
  }, .parallel = parallel)

  ### return
  tab <- as_data_frame(tab)
  
  return(tab)
}

