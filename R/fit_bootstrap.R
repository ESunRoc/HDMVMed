#' Fitting the HDMVMediation model 
#' 
#' `bootstrap_model` is the main function for fitting and performing inference for the high-dimensional multivariate outcomes mediation model of Sun, Xiao, and Wu (2025).
#'
#' @param mediators An \eqn{n\times p} numeric matrix of candidate mediators.
#' @param confounders An \eqn{n\times\ell} numeric matrix of confounders.
#' @param trt An \eqn{n\times 1} numeric matrix or vector of treatment/exposure status.
#' @param outcomes An \eqn{n\times q} numeric matrix of outcomes.
#' @param quiet_msglasso A boolean indicator determining whether or not the CV progress of the MSGLasso estimation should be suppressed; defaults to `TRUE`.
#' @param lam1.v A numeric vector of component-wise regularization parameters as in the MSGLasso of Li et al. (2015). By default, this is taken to be `seq(from = 1e-03, to = 0.05, length = 20)`.
#' @param lamG.v A numeric vector of group-wise regularization parameters as in the MSGLasso of Li et al. (2015). By default, this is taken to be `seq(from = 1e-03, to = 0.05, length = 20)`.
#' @param alpha A numeric between 0 and 1 denoting the confidence level for inference; defaults to 0.05.
#' @param nB A numeric determining the number of bootstrap samples taken for inference. By default, 5000 samples will be used.
#' @param msg_folds A numeric determining the number of CV folds to be used in tuning the MSGLasso regularization parameters. By default, this is taken to be 5.
#' @param seed An integer seed for reproducible bootstrap inference.
#' @param outcome_grps A boolean indicator for whether or not there are any groups on the outcome variables. By default, there are no outcome groups so this argument is `FALSE`.
#' @param NumOutGrps An integer denoting the number of outcome groups; default is `NULL` in keeping with the default `outcome_grps = FALSE`
#' @param OutGrpStarts A vector of starting coordinates for the outcome groups. Equivalent to `R.Starts` from the `MSGLasso` package.
#' @param OutGrpEnds A vector of ending coordinates for the outcome groups. Equivalent to `R.Ends` from the `MSGLasso` package.
#' 
#' @returns A \eqn{(p*q+2*q)\times 8} matrix `mod_boot_summ` with columns:
#' * `OrigEst`: the original PIDE/TIDE/DE estimate;
#' * `Mean_boot`: the mean of the bootstrap PIDE/TIDE/DE estimates;
#' * `boot_SE`: the standard deviation of the bootstrap PIDE/TIDE/DE estimates;
#' * `boot_pval`: the unadjusted p-value for the bootstrap PIDE/TIDE/DE estimates;
#' * `bca_lowerCL`: the lower limit of the \eqn{100(1-\texttt{alpha})\%} BCa CI;
#' * `bca_upperCL`: the upper limit of the \eqn{100(1-\texttt{alpha})\%} BCa CI;
#' * `per_lowerCL`: the lower limit of the \eqn{100(1-\texttt{alpha})\%} percentile CI;
#' * `per_upperCL`: the upper limit of the \eqn{100(1-\texttt{alpha})\%} percentile CI.
#' 
#' @references{
#' Li, Y., Nan, B., and Zhu, J. (2015). Multivariate Sparse Group Lasso for the 
#' Multivariate Multiple Linear Regression with an Arbitrary Group Structure. 
#' \emph{Biometrics}, *71*, 354-363.
#' }
#' 
#' @references{
#' Sun, E., Xiao, J., and Wu, T. T. (2025). Causal Mediation Analysis for Multiple 
#' Outcomes and High-dimensional Mediators: Identification, Inference, and Application.
#' \emph{Biometrics}. _Under Review._
#' }
#' 
#' @examples
#' ## Load toy data
#' data(hdmvmed_test_data)
#' 
#' ## Assign mediator/confounder/trt/outcome matrices
#' mediators <- hdmvmed_test_data[,7:106]
#' confounders <- hdmvmed_test_data[,107:111]
#' trt <- hdmvmed_test_data[,1]
#' outcomes <- hdmvmed_test_data[,2:6]
#' 
#' ## Fit the model
#' model_fit <- bootstrap_model(mediators = mediators, confounders = confounders,
#'                              trt = trt, outcomes = outcomes, nB = 10)
#' model_fit
#'
#' @export
bootstrap_model <- function(mediators, confounders, trt, outcomes, quiet_msglasso = TRUE,
                            lam1.v = seq(1e-3, 0.05, length=20), lamG.v = seq(1e-3, 0.05, length=20), 
                            alpha = 0.05, nB = 5e3, msg_folds = 5, seed = 823543, outcome_grps = FALSE,
                            NumOutGrps = NULL, OutGrpStarts = NULL, OutGrpEnds = NULL){
  
  if(msg_folds<=1) stop("You must use at least 2 folds for tuning MSGLasso")
  
  k <- 1                      # number of exposures/treatments
  p <- ncol(mediators)   # number of mediators
  l <- ncol(confounders) # number of confounders
  q <- ncol(outcomes)    # number of responses
  n <- length(trt)       # number of subjects
  
  
  #### Treatment --> Mediators ####
  trtToMedi <- lm(mediators ~ trt + confounders)
  trtToMedi_summ <- broom::tidy(trtToMedi)
  
  mod_Stage1 <- trtToMedi_summ[which(trtToMedi_summ$term=="trt"),]
  
  
  #### Mediators --> Outcomes ####
  P <- p + l + k; Q <- q
  
  G <- p + 2 # Groups on X: p + 2; p mediator singleton groups + 1 confounder group + 1 treatment group
  gmax <- 1 # each variable (resp or pred) belongs to only 1 group
  cmax <- l # a group contains at most l variables (confounders form largest group)
  GarrStarts <- c(0:(p-1), p,  p+l+1); GarrEnds <-   c(0:(p-1), p+l, p+l+2)
  
  if(outcome_grps == T){
    R <- NumOutGrps                                    # Groups on Y: user-specified
    RarrStarts <- OutGrpStarts; RarrEnds <- OutGrpEnds # non-null response groupings
  } else{
    R <- q                                           # Groups on Y: q; all responses are singleton groups
    RarrStarts <- c(0:(q-1)); RarrEnds <- c(0:(q-1)) # singleton response groups
  }
  
  tmp_PQgrps <- FindingPQGrps(P = P, Q = Q, G, R, gmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
  PQgrps <- tmp_PQgrps$PQgrps
  
  tmp_grpWts <- Cal_grpWTs(P = P, Q = Q, G, R, gmax, PQgrps)
  grpWTs <- tmp_grpWts$grpWTs
  
  tmp_GRgrps <- FindingGRGrps(P = P, Q = Q, G, R, cmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
  GRgrps <- tmp_GRgrps$GRgrps
  
  
  Pen_L <- matrix(rep(1, P*Q), P, Q, byrow=T); assign("Pen_L", Pen_L, envir=.GlobalEnv)
  Pen_L[(p+1):P,] <- 0 # don't penalize the confounders, (p+1):(P-1), or treatment, P
  
  Pen_G <- matrix(rep(1,G*R),G,R, byrow=TRUE); assign("Pen_G", Pen_G, envir=.GlobalEnv)
  Pen_G[(G-1):G,] <- 0 # don't penalize confounder group, G-1, or treatment group, G
  
  grp_Norm0 <- matrix(rep(1, G*R), nrow=G, byrow=TRUE)
  
  
  lam1.v <- lam1.v; lamG.v <- lamG.v
  
  if(quiet_msglasso == T){
    capture.output(mod_try.cv <- MSGLasso.cv(X = cbind(mediators,confounders,trt), 
                                             Y = outcomes,
                                             grpWTs, Pen_L, Pen_G,
                                             PQgrps, GRgrps, lam1.v, lamG.v,
                                             # grp_Norm = grp_Norm0,
                                             fold = msg_folds, seed = seed), file = nullfile())
  } else{
    mod_try.cv <- MSGLasso.cv(X = cbind(mediators,confounders,trt), 
                              Y = outcomes,
                              grpWTs, Pen_L, Pen_G,
                              PQgrps, GRgrps, lam1.v, lamG.v,
                              # grp_Norm = grp_Norm0,
                              fold = msg_folds, seed = seed)
  }
  
  MSGLassolam1 <- mod_try.cv$lams.c[which.min(as.vector(mod_try.cv$rss.cv))][[1]]$lam1
  MSGLassolamG <- mod_try.cv$lams.c[which.min(as.vector(mod_try.cv$rss.cv))][[1]]$lam3
  MSGLassolamG.m <- matrix(rep(MSGLassolamG, G*R),G,R,byrow=TRUE) 
  MSGLassolamG.m[(G-1):G,] <- 0
  
  mod_Stage2 <- MSGLasso(X.m = cbind(mediators, confounders, trt),
                         Y.m = outcomes,
                         grpWTs, Pen_L, Pen_G, PQgrps, GRgrps,
                         grp_Norm0, MSGLassolam1, MSGLassolamG.m)
  
  #### Debiase mod_Stage2
  # theta_mod <- theta_calc_parallel(X = mediators)
  # mod_Stage2_debiased <- mod_Stage2$Beta[1:p,] + (1/n)*theta_mod%*%t(mediators)%*%(outcomes - mediators%*%mod_Stage2$Beta[1:p,])
  # mod_Stage2_all <- rbind(mod_Stage2_debiased, mod_Stage2$Beta[(p+1):nrow(mod_Stage2$Beta),])
  # rownames(mod_Stage2_all) <- c(colnames(mediators), colnames(confounders), "Exposure")
  
  
  theta_mod <- theta_calc_parallel(X = cbind(mediators, confounders, trt))
  mod_Stage2_debiased <- mod_Stage2$Beta + (1/n)*theta_mod%*%t(cbind(mediators, confounders, trt))%*%(outcomes - cbind(mediators, confounders, trt)%*%mod_Stage2$Beta)
  mod_Stage2_all <- mod_Stage2_debiased
  rownames(mod_Stage2_all) <- c(colnames(mediators), colnames(confounders), "Exposure")
  
  #### Calculate DE and PIDEs
  mod_DE <- mod_Stage2$Beta[nrow(mod_Stage2$Beta),]
  
  mod_stage2_pide <- mod_Stage2_all[1:p,]
  mod_stage1_pide <- mod_Stage1$estimate
  
  mod_pide_mat <- matrix(nrow = p, ncol = q)
  for(i in 1:q) mod_pide_mat[,i] <- mod_stage1_pide*mod_stage2_pide[,i]
  
  rownames(mod_pide_mat) <- rownames(mod_Stage2_all)[1:p]; colnames(mod_pide_mat) <- colnames(mod_Stage2_all)
  mod_tide_mat <- colSums(mod_pide_mat)
  
  mod_origFit_mat <- matrix(nrow = p*q+2*q, ncol = 1)
  mod_origFit_mat[,1] <- c(vec(mod_pide_mat),unname(mod_tide_mat),unname(mod_DE))
  rownames(mod_origFit_mat) <- c(paste0(rownames(mod_pide_mat),"_ide_resp",rep(1:q, times = rep(p,q))),
                                 paste0("TIDE_resp",1:q),
                                 paste0("DE_resp", 1:q))
  colnames(mod_origFit_mat) <- "Orig_Est"
  
  #### Draw bootstrap samples
  B_draws <- matrix(nrow = n, ncol = nB)
  set.seed(seed)
  for(i in 1:nB) {B_draws[,i] <- sample(1:n, n, replace=T)}
  
  mod_bootRes <- matrix(nrow = p*q+2*q, ncol = nB)
  rownames(mod_bootRes) <- c(paste0(rownames(mod_pide_mat),"_ide_resp",rep(1:q, times = rep(p,q))),
                             paste0("TIDE_resp",1:q),
                             paste0("DE_resp", 1:q))
  colnames(mod_bootRes) <- paste0("BootDraw",1:nB)
  
  
  start_time <- proc.time()
  for(i in 1:nB){
    boot_sample <- B_draws[,i]
    confounders_boot <- confounders[boot_sample,]
    trt_boot <- trt[boot_sample]
    mediators_boot <- mediators[boot_sample,]
    outcomes_boot <- outcomes[boot_sample,]
    
    ## Stage 1 fit
    mod_stage1_boot_fit <- lm(mediators_boot ~ trt_boot + confounders_boot)
    mod_stage1_boot_tidy <- broom::tidy(mod_stage1_boot_fit)
    mod_stage1_boot <- mod_stage1_boot_tidy[which(mod_stage1_boot_tidy$term=="trt_boot"),]
    
    ## Stage 2 fit
    mod_stage2_boot_fit <- MSGLasso(X.m = cbind(mediators_boot, confounders_boot, trt_boot),
                                    Y.m = outcomes_boot,
                                    grpWTs, Pen_L, Pen_G, PQgrps, GRgrps,
                                    grp_Norm0, MSGLassolam1, MSGLassolamG.m)
    mod_stage2_boot_all <- mod_stage2_boot_fit$Beta + (1/n) * theta_mod %*% t(cbind(mediators_boot, confounders_boot, trt_boot)) %*%
      (outcomes_boot - cbind(mediators_boot, confounders_boot, trt_boot) %*% mod_stage2_boot_fit$Beta)
    
    # mod_stage2_boot_debiased <- mod_stage2_boot_fit$Beta[1:p,] + (1/n) * theta_mod %*% t(mediators) %*% (outcomes - mediators %*% mod_stage2_boot_fit$Beta[1:p,])
    # mod_stage2_boot_all <- rbind(mod_stage2_boot_debiased, mod_stage2_boot_fit$Beta[(p+1):nrow(mod_stage2_boot_fit$Beta),])
    
    rownames(mod_stage2_boot_all) <- colnames(cbind(mediators_boot, confounders_boot, trt_boot))
    colnames(mod_stage2_boot_all) <- colnames(outcomes_boot)
    
    ## PIDE and DE calculation
    mod_boot_DE <- unname(mod_stage2_boot_all["trt_boot",])
    
    mod_stage2_pide_boot <- mod_stage2_boot_all[1:p,]
    mod_stage1_pide_boot <- mod_stage1_boot$estimate
    
    mod_pide_mat_boot <- matrix(nrow = p, ncol = q)
    for(j in 1:q) mod_pide_mat_boot[,j] <- mod_stage1_pide_boot*mod_stage2_pide_boot[,j]
    
    mod_tide_mat_boot_sums <- colSums(mod_pide_mat_boot)
    
    mod_bootRes[,i] <- c(vec(mod_pide_mat_boot),mod_tide_mat_boot_sums,mod_boot_DE)
    
    if(i %% 100 == 0){
      cat(paste0("Done with boot sample ", i, "; ", nB-i, " remaining. Time elapsed: ", 
                 round(-1*(start_time[3] - proc.time()[3])/60, 3), " minutes. Apx. ",
                 round(((-1*(start_time[3] - proc.time()[3])/60)/i)*(nB-i), 3), " minutes remaining\n"))
    }
  }
  
  mod_bootRes <- ifelse((is.nan(mod_bootRes) | is.na(mod_bootRes)), rowMeans(mod_bootRes, na.rm=T), mod_bootRes)
  
  ## bootstrap p-values
  mod_boot_pvals <- matrix(nrow=p*q+2*q, ncol = 1)
  for(i in 1:nrow(mod_bootRes)){
    if(mod_origFit_mat[i,] >= 0){
      mod_boot_pvals[i,] <- min(2*sum(mod_bootRes[i,]>=2*mod_origFit_mat[i,],na.rm=T)/nB,1)
    } else{
      mod_boot_pvals[i,] <- min(2*sum(mod_bootRes[i,]<2*mod_origFit_mat[i,],na.rm=T)/nB,1)
    }
  }
  
  ## BCa confidence intervals
  mod_boot_bcaCI <- t(apply(mod_bootRes, 1, coxed::bca, conf.level = 1-alpha))
  mod_boot_bcaCI[is.nan(mod_boot_bcaCI)] <- 0
  
  ## pivotal bootstrap CIs
  mod_boot_perCI <- 2*matrix(data = c(mod_origFit_mat, mod_origFit_mat), ncol = 2) - t(apply(mod_bootRes, 1, quantile, c(1-alpha/2,alpha/2)))
  
  ## Summary table
  mod_boot_summ <- cbind(mod_origFit_mat, apply(mod_bootRes,1,mean), apply(mod_bootRes,1,sd),
                         mod_boot_pvals, mod_boot_bcaCI, mod_boot_perCI)
  colnames(mod_boot_summ) <- c("OrigEst", "Mean_boot", "boot_SE", "boot_pval",
                               "bca_lowerCL", "bca_upperCL",
                               "per_lowerCL", "per_upperCL")
  mod_boot_summ[is.nan(mod_boot_summ)] <- 0
  
  return(mod_boot_summ)
}
