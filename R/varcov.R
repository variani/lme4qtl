
#' @export
varcov <- function(mod, comp, scaled = TRUE, residual = TRUE,
  idvar, ids)
{
  ### arg
  missing_idvar <- missing(idvar)
  missing_ids <- missing(ids)
  all_comp <- missing(comp)
  stopifnot(class(mod)[1] %in% c("lmerMod", "glmerMod"))
  
  # estimate varcov matrix
  n <- lme4::getME(mod, "n")
  
  if(all_comp) { 
    varcov <- tcrossprod(crossprod(lme4::getME(mod, "Zt"), 
      lme4::getME(mod, "Lambda"))) + Diagonal(n)
  } else {
    Lambda <- lme4::getME(mod, "Lambda")
  
    zt_list <- getME(mod, "Ztlist")
    zt_num <- length(zt_list)
    
    zt_names <- names(zt_list)
    zt_names0 <- sapply(strsplit(zt_names, "\\."), function(x) x[1])
    stopifnot(all(comp %in% zt_names0))
    
    zt_ind <- which(zt_names0 %in% comp)
    stopifnot(length(zt_ind) == length(comp))
    
    varcov <- Diagonal(n)
    lambda_start <- 1
    for(i in seq(1, zt_num)) {
      lambda_step <- nrow(zt_list[[i]])
      lambda_ind <- seq(lambda_start, length = lambda_step)
      
      if(i %in% zt_ind) {
        varcov <- varcov + (tcrossprod(crossprod(zt_list[[i]], Lambda[lambda_ind, lambda_ind])))
      }
      
      lambda_start <- lambda_start + lambda_step
    }    
    
    if(!residual) {
      varcov <- varcov - Diagonal(n)
    }
  }
 
  if(!scaled) {
    s2 <- sigma(mod)^2
    varcov <- s2 * varcov
  }
  
  if(!missing_idvar) {
    stopifnot(idvar %in% names(mod@frame))
    fids <- mod@frame[, idvar]
    fids <- as.character(fids)
    
    rownames(varcov) <- fids
    colnames(varcov) <- fids
  }

  if(!missing_ids) {
    stopifnot(!missing_idvar)
    
    stopifnot(all(ids %in% rownames(varcov)))
    ind <- sapply(ids, function(x) which(x == rownames(varcov)))
    
    varcov <- varcov[ind, ind]
  }
  
  return(varcov)
}

