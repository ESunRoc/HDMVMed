#' Fitting the HDMVMediation model
#'
#' `bootstrap_model` is the main function for fitting and performing inference for the high-dimensional multivariate outcomes mediation model of Sun, Xiao, and Wu (2025).
#'
#' @param mediators An \eqn{n\times p} numeric matrix of candidate mediators.
#' @param confounders An \eqn{n\times\ell} numeric matrix of confounders.
#' @param trt An \eqn{n\times 1} numeric matrix or vector of treatment/exposure status.
#' @param outcomes An \eqn{n\times q} numeric matrix of outcomes.
#' @param moderators An optional \eqn{n\times r} numeric matrix of candidate moderators of
#'   the exposure effect; defaults to `NULL` (no moderated-mediation analysis). See details.
#' @param quiet_msglasso A boolean indicator determining whether or not the CV progress of the MSGLasso estimation should be suppressed; defaults to `TRUE`.
#' @param lam1.v A numeric vector of component-wise regularization parameters as in the MSGLasso of Li et al. (2015). By default, this is taken to be `seq(from = 1e-03, to = 0.05, length = 10)`.
#' @param lamG.v A numeric vector of group-wise regularization parameters as in the MSGLasso of Li et al. (2015). By default, this is taken to be `seq(from = 1e-03, to = 0.05, length = 10)`.
#' @param alpha A numeric between 0 and 1 denoting the confidence level for inference; defaults to 0.05.
#' @param nB A numeric determining the number of bootstrap samples taken for inference. By default, 5000 samples will be used.
#' @param msg_folds A numeric determining the number of CV folds to be used in tuning the MSGLasso regularization parameters. By default, this is taken to be 3.
#' @param seed An integer seed for reproducible bootstrap inference.
#' @param theta_parallel A Boolean indicator for whether to compute the debiasing matrices in parallel. By default, this is `TRUE`.
#' @param theta_cores An integer denoting the number of logical cores to use when
#'   `theta_parallel = TRUE`; passed through to [theta_calc_parallel()]. Defaults to `max(1,detectCores()-1)`
#' @param theta_folds An integer denoting the number of CV folds used in the nodewise regressions underlying the debiasing matrix; passed through to [theta_calc()]/[theta_calc_parallel()]. Defaults to 5.
#' @param boot_parallel A Boolean indicator for whether the `nB` bootstrap replicates, the dominant cost, should be computed in parallel via `doParallel`/`foreach`. By default, this is `TRUE`.
#' @param boot_cores An integer denoting the number of logical cores to use when `boot_parallel = TRUE`. Defaults to `max(1,detectCores()-1)`
#' @param penalize_conf A Boolean indicator for whether the confounder block should be unpenalized. By default, this is `FALSE`.
#' @param penalize_moderators A Boolean indicator for whether the treatment x moderator
#'   interaction block should be unpenalized. By default, this is `FALSE`, so that
#'   moderated effects are always estimated and reported rather than possibly shrunk
#'   to zero by the group-Lasso penalty. Ignored if `moderators` is `NULL`.
#' @param outcome_grps A Boolean indicator for whether or not there are any groups on the outcome variables. By default, there are no outcome groups so this argument is `FALSE`.
#' @param NumOutGrps An integer denoting the number of outcome groups; default is `NULL` in keeping with the default `outcome_grps = FALSE`
#' @param OutGrpStarts A vector of starting coordinates for the outcome groups. Equivalent to `R.Starts` from the `MSGLasso` package.
#' @param OutGrpEnds A vector of ending coordinates for the outcome groups. Equivalent to `R.Ends` from the `MSGLasso` package.
#' @param medi_grps A Boolean indicator for whether or not there are any groups on the mediators. By default, there are no mediator groups so this outcome is `FALSE`.
#' @param NumMediGrps An integer denoting the number of mediator groups; default is `NULL` in keeping with the default `medi_grps = FALSE`.
#' @param MediGrpStarts A vector of starting coordinates for the mediators groups. Equivalent to `G.Starts` from the `MSGLasso` package. Starts from 0 and runs to, at most, \eqn{p-1}.
#' @param MediGrpEnds A vector of ending coordinates for the outcome groups. Equivalent to `G.Ends` from the `MSGLasso` package.
#'
#' @returns A list containing `res` and `internals`. `res` is a \eqn{(1+r)(p q+2q)\times 8} matrix (`r` = `ncol(moderators)`,
#'   or simply a \eqn{(p q+2q)\times 8} matrix if `moderators` is `NULL`, as before) with
#'   columns:
#' * `OrigEst`: the original PIDE/TIDE/DE (or moderated-PIDE/TIDE/DE) estimate;
#' * `Mean_boot`: the mean of the bootstrap estimates;
#' * `boot_SE`: the standard deviation of the bootstrap estimates;
#' * `boot_pval`: the unadjusted bootstrap p-value;
#' * `bca_lowerCL`: the lower limit of the \eqn{100(1-\texttt{alpha})\%} BCa CI;
#' * `bca_upperCL`: the upper limit of the \eqn{100(1-\texttt{alpha})\%} BCa CI;
#' * `per_lowerCL`: the lower limit of the \eqn{100(1-\texttt{alpha})\%} percentile CI;
#' * `per_upperCL`: the upper limit of the \eqn{100(1-\texttt{alpha})\%} percentile CI.
#' `internals` is a list containing:
#' * `alpha_hat`, `phi_hat_d`, `xi_hat`, `eta_hat`, `tau_hat`: estimated coefficients
#' * `Sigma_M_hat`, `Sigma_Y_hat`: estimated covariance matrices for the mediators and outcomes;
#' * `Theta`: the nodewise regression debiasing matrix;
#' * `confounders`: the confounders adjusted for;
#' * `n`, `p`, `l`, `q`: the number of observations, mediators, confounders, and outcomes;
#' * `has_moderators`: a Boolean indicating whether moderators were included.
#'
#' Row names identify the estimand as before (`"{mediator}_ide_resp{s}"`,
#' `"TIDE_resp{s}"`, `"DE_resp{s}"`), with moderated-effect rows additionally tagged
#' `"..._MOD_{moderator}_resp{s}"`; [hdmvm_table()] parses this into an extra
#' `Moderator` column.
#'
#'
#' @details
#' The `MediGrpStarts` and `MediGrpEnds` can be structured as follows (`OutGrpStarts` and `OutGrpEnds` are analogous).
#' Say you have \eqn{g} mediator groups, each with \eqn{p_1,\ldots,p_g} mediators. Then
#' `MediGrpStarts=c(0, p_1, ..., p_{g-1})` and `MediGrpEnds=c(p_{1}-1, p_{2}-1, ..., p_g)`. In this way, \eqn{p_g\le p}.
#' Critically, these index starting from 0 to be consistent with internal calls to `MSGLasso`.
#'
#' Each moderator `Z` is assumed to enter every regression as a treatment interaction, i.e.
#' `mediators ~ trt + confounders + moderators + trt:moderators` for eq. (2.1) and
#' `outcomes ~ mediators + trt + confounders + moderators + trt:moderators` for eq. (2.2).
#' The main effect of each moderator is added to the confounder block automatically, so
#' `moderators` should *not* also be duplicated inside `confounders`. When supplied,
#' `bootstrap_model` additionally reports, for every moderator and (mediator,outcome)
#' combination: the index of moderated mediation (the rate of change of the
#' PIDE/TIDE per unit of the moderator) and the moderated direct effect (the rate of
#' change of the DE per unit of the moderator), together with the same bootstrap
#' p-values and BCa/percentile confidence intervals computed for the unmoderated
#' PIDE/TIDE/DE.
#'
#' **Moderated mediation.** For a moderator \eqn{Z_r}, eq. (2.1)-(2.2) of Sun et al. become
#' \eqn{M=A\alpha+U\xi+(A\odot Z_r)\delta^M_r+\epsilon_M} and
#' \eqn{Y=M\phi+A\tau+U\eta+(A\odot Z_r)\delta^Y_r+\epsilon_Y}, so that at a fixed level
#' \eqn{z} of the moderator the PIDE/TIDE/DE from Theorem 2.1 become
#' \eqn{(\alpha+z\delta^M_r)\odot\phi_s}, \eqn{\sum_j(\alpha_{j}+z\delta^M_{r,j})\phi_{js}},
#' and \eqn{\tau_s+z\delta^Y_{r,s}} respectively -- i.e. linear in \eqn{z} with slopes
#' \eqn{\delta^M_r\odot\phi_s} (the index of moderated mediation, by mediator),
#' \eqn{\sum_j\delta^M_{r,j}\phi_{js}} (the index of moderated mediation, summed), and
#' \eqn{\delta^Y_{r,s}} (the moderated direct effect). These slopes (the rate of
#' change of the PIDE/TIDE/DE per unit of \eqn{Z_r}) that `bootstrap_model` reports for
#' each moderator when `moderators` is supplied, using the debiased \eqn{\hat\phi^d} for
#' the mediator-side product exactly as for the unmoderated PIDE/TIDE (Theorem 2.3), and
#' the raw (non-debiased) \eqn{\hat\delta^Y_r} for the moderated DE, exactly as for the
#' unmoderated DE (as, like the treatment and confounder blocks, the interaction block
#' is unpenalized by default; see `penalize_moderators`).
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
#' \emph{Biostatistics}. _Under Review._
#' }
#'
#' @examples
#' \dontrun{
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
#' ## Fit the model with a moderator of the exposure effect (e.g. the first confounder)
#' moderators <- confounders[, 1, drop = FALSE]
#' model_fit_mod <- bootstrap_model(mediators = mediators, confounders = confounders[,-1,drop=FALSE],
#'                              trt = trt, outcomes = outcomes, moderators = moderators, nB = 10)
#' model_fit_mod
#'}
#'
#' @export
bootstrap_model <- function(mediators, confounders, trt, outcomes, moderators = NULL, quiet_msglasso = TRUE,
                            lam1.v = seq(1e-3, 0.05, length=10), lamG.v = seq(1e-3, 0.05, length=10),
                            alpha = 0.05, nB = 5e3, msg_folds = 5, seed = 823543, theta_parallel = TRUE,
                            theta_cores = max(1, parallel::detectCores()-1), theta_folds = 5,
                            boot_parallel = TRUE, boot_cores = max(1, parallel::detectCores()-1),
                            penalize_conf = FALSE, penalize_moderators = FALSE, outcome_grps = FALSE, NumOutGrps = NULL, OutGrpStarts = NULL,
                            OutGrpEnds = NULL, medi_grps = FALSE, NumMediGrps = NULL, MediGrpStarts = NULL, MediGrpEnds = NULL){

  if(msg_folds<=1) stop("You must use at least 2 folds for tuning MSGLasso")

  # Ensure that mediators, confounders, and outcomes are all named
  if(is.null(colnames(mediators))) colnames(mediators) <- paste0("M", 1:ncol(mediators))
  if(is.null(colnames(confounders))) colnames(confounders) <- paste0("C", 1:ncol(confounders))
  if(is.null(colnames(outcomes))) colnames(outcomes) <- paste0("Y", 1:ncol(outcomes))

  #### Moderators: build treatment x moderator interaction terms ####
  # Moderators enter each regression as a treatment interaction. Each moderator's main effect is folded into the
  # confounder block so the interaction coefficient isn't confounded with an omitted main effect of Z; the interaction
  # columns (AxZ) are appended as their own unpenalized-by-default block, placed after the treatment column
  # exactly the way the treatment column itself is appended after the confounder block below.
  has_moderators <- !is.null(moderators)
  if(has_moderators){
    moderators <- as.matrix(moderators)
    if(is.null(colnames(moderators))) colnames(moderators) <- paste0("Z", 1:ncol(moderators))
    r <- ncol(moderators)

    # input validation for moderators:
    # a moderator that is constant (zero variance) has an unidentifiable treatment interaction
    # This column(s) are reported and should be dropped from the input
    const_mods <- colnames(moderators)[apply(moderators, 2, function(z) length(unique(z)) <= 1)]
    if(length(const_mods) > 0){
      stop("The following moderator(s) are constant (zero variance) in this sample: ",
           paste(const_mods, collapse=", "),
           ". A constant moderator has no estimable treatment interaction (trt*Z ",
           "carries no information) and will silently corrupt downstream PIDE/TIDE/DE ",
           "estimates or crash the underlying MSGLasso fit. Remove the constant ",
           "moderator(s) from `moderators` before calling bootstrap_model(), or, if ",
           "this reflects a genuine subgroup with zero observed cases, consider ",
           "whether this moderator is estimable in this sample at all.")
    }

    AxZ <- moderators * as.vector(trt) # elementwise A_i * Z_i,r for every row i, moderator r
    colnames(AxZ) <- paste0("trtx", colnames(moderators))

    confounders <- cbind(confounders, moderators) # moderator main effects join the confounder block
  } else {
    r <- 0
  }

  stage1_design <- if(has_moderators) cbind(1, trt, confounders, AxZ) else cbind(1, trt, confounders)
  stage1_design_names <- if(has_moderators) c("(Intercept)", "trt", colnames(confounders), colnames(AxZ)) else c("(Intercept)", "trt", colnames(confounders))
  qr_design <- qr(stage1_design)
  if(qr_design$rank < ncol(stage1_design)){
    dependent_cols <- stage1_design_names[qr_design$pivot[(qr_design$rank+1):ncol(stage1_design)]]
    stop("The design matrix formed by trt, confounders, and (if used) moderators ",
         "and their treatment interactions is rank-deficient: the following column(s) ",
         "are an exact linear combination of the others and have no identifiable ",
         "coefficient: ", paste(dependent_cols, collapse=", "), ". This may ",
         "happen when two or more moderators (or confounders) are mutually exclusive ",
         "levels of the same underlying categorical variable. Drop one of the implicated moderator/confounder columns (or ",
         "re-code them, e.g. keep only k-1 of k mutually exclusive categories) before ",
         "calling bootstrap_model().")
  }

  k <- 1                 # number of exposures/treatments
  p <- ncol(mediators)   # number of mediators
  l <- ncol(confounders) # number of confounders (includes moderator main effects, if any)
  q <- ncol(outcomes)    # number of responses
  n <- length(trt)       # number of subjects

  X_design <- if(has_moderators) cbind(mediators, confounders, trt, AxZ) else cbind(mediators, confounders, trt)


  #### Treatment --> Mediators ####
  trtToMedi <- if(has_moderators) lm(mediators ~ trt + confounders + AxZ) else lm(mediators ~ trt + confounders)
  trtToMedi_summ <- broom::tidy(trtToMedi)

  mod_Stage1 <- trtToMedi_summ[which(trtToMedi_summ$term=="trt"),]

  # Retain the stage-1 residual covariance and confounder loadings for sensitivity analysis.
  # xi_hat is taken from coef() rather than the tidy frame so that row/column alignment is unambiguous.
  Sigma_M_hat <- crossprod(residuals(trtToMedi))/n
  coef_stage1 <- coef(trtToMedi)
  xi_hat <- coef_stage1[grep("^confounders", rownames(coef_stage1)), , drop = FALSE]
  rownames(xi_hat) <- colnames(confounders)

  if(has_moderators){
    # delta^M: p x r matrix of A x Z_r -> mediator_j coefficients (Stage-1 moderation)
    delta_M <- extract_interaction_estimates(trtToMedi_summ, "AxZ", colnames(AxZ), colnames(mediators))
  }


  #### Mediators --> Outcomes ####
  P <- p + l + k + r; Q <- q

  # Group bookkeeping. FindingPQGrps() (MSGLasso) assigns 0-indexed column p to group
  # g whenever GarrStarts[g] <= p <= GarrEnds[g] (inclusive on both ends; confirmed
  # directly against MSGLasso's C source, Find_PQ_Coord_Grps).
  if(medi_grps){
    n_medi_grps <- NumMediGrps # FIX: was `medi_grps` (a boolean), not the group count
    base_starts <- MediGrpStarts; base_ends <- MediGrpEnds
    cmax <- max(c(MediGrpEnds - MediGrpStarts, l, k, r)) # FIX: also account for the moderator-interaction block width
  } else {
    n_medi_grps <- p
    base_starts <- 0:(p-1); base_ends <- 0:(p-1)
    cmax <- max(c(l, k, r)) # FIX: also account for the moderator-interaction block width
  }
  gmax <- 1 # each variable (resp or pred) belongs to only 1 group; overlapping groups are not currently supported

  conf_grp_row <- n_medi_grps + 1
  trt_grp_row  <- n_medi_grps + 2
  modint_grp_row <- if(has_moderators) n_medi_grps + 3 else NA

  G <- n_medi_grps + 2 + as.integer(has_moderators) # mediator groups + confounder group + treatment group (+ moderator-interaction group)

  GarrStarts <- c(base_starts, p,     p+l,  if(has_moderators) p+l+1)
  GarrEnds   <- c(base_ends,   p+l-1, p+l,  if(has_moderators) p+l+r)


  if(outcome_grps){
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

  # This scoping is required due to a bug in the `MSGLasso` package:
  # In MSGLasso.R (lines 30-37), the transpose-or-not check for its `Pen.L`/`Pen.G` parameters is
  # written as `if(!is.null(ncol(Pen_L)))`/`ncol(Pen_G)`. This naming discrepancy causes an error if
  # `Pen_L`/`Pen_G` don't exist in the global environment of the user.
  Pen_L <<- matrix(rep(1, P*Q), P, Q, byrow=T)
  Pen_G <<- matrix(rep(1,G*R),G,R, byrow=TRUE)

  if(penalize_conf == FALSE){ # if penalize_conf == FALSE,
    Pen_L[(p+1):(p+l+k),] <- 0            # don't penalize the confounders (incl. moderator main effects) or treatment
    Pen_G[c(conf_grp_row, trt_grp_row),] <- 0 # don't penalize confounder group or treatment group
  }
  if(has_moderators && !penalize_moderators){ # if penalize_moderators == FALSE (the default),
    Pen_L[(p+l+k+1):P,] <- 0    # don't penalize the treatment x moderator interaction columns
    Pen_G[modint_grp_row,] <- 0 # don't penalize the treatment x moderator interaction group
  }
  grp_Norm0 <- matrix(rep(1, G*R), nrow=G, byrow=TRUE)


  lam1.v <- lam1.v; lamG.v <- lamG.v

  if(quiet_msglasso == T){
    capture.output(mod_try.cv <- MSGLasso.cv(X = X_design,
                                             Y = outcomes,
                                             grpWTs, Pen_L, Pen_G,
                                             PQgrps, GRgrps, lam1.v, lamG.v,
                                             grp_Norm = grp_Norm0,
                                             fold = msg_folds, seed = seed), file = nullfile())
  } else{
    mod_try.cv <- MSGLasso.cv(X = X_design,
                              Y = outcomes,
                              grpWTs, Pen_L, Pen_G,
                              PQgrps, GRgrps, lam1.v, lamG.v,
                              grp_Norm = grp_Norm0,
                              fold = msg_folds, seed = seed)
  }

  MSGLassolam1 <- mod_try.cv$lams.c[which.min(as.vector(mod_try.cv$rss.cv))][[1]]$lam1
  MSGLassolamG <- mod_try.cv$lams.c[which.min(as.vector(mod_try.cv$rss.cv))][[1]]$lam3
  MSGLassolamG.m <- matrix(rep(MSGLassolamG, G*R),G,R,byrow=TRUE)
  MSGLassolamG.m[c(conf_grp_row, trt_grp_row),] <- 0
  if(has_moderators && !penalize_moderators) MSGLassolamG.m[modint_grp_row,] <- 0

  mod_Stage2 <- MSGLasso(X.m = X_design,
                         Y.m = outcomes,
                         grpWTs, Pen_L, Pen_G, PQgrps, GRgrps,
                         grp_Norm0, MSGLassolam1, MSGLassolamG.m)

  #### Debiase mod_Stage2
  # theta_mod <- theta_calc_parallel(X = mediators)
  # mod_Stage2_debiased <- mod_Stage2$Beta[1:p,] + (1/n)*theta_mod%*%t(mediators)%*%(outcomes - mediators%*%mod_Stage2$Beta[1:p,])
  # mod_Stage2_all <- rbind(mod_Stage2_debiased, mod_Stage2$Beta[(p+1):nrow(mod_Stage2$Beta),])
  # rownames(mod_Stage2_all) <- c(colnames(mediators), colnames(confounders), "Exposure")

  if(theta_parallel){
    theta_mod <- theta_calc_parallel(X = X_design, cores = theta_cores, folds = theta_folds)
  } else {
    theta_mod <- theta_calc(X = X_design)
  }
  mod_Stage2_debiased <- mod_Stage2$Beta + (1/n)*theta_mod%*%t(X_design)%*%(outcomes - X_design%*%mod_Stage2$Beta)
  mod_Stage2_all <- mod_Stage2_debiased
  rownames(mod_Stage2_all) <- c(colnames(mediators), colnames(confounders), "Exposure", if(has_moderators) colnames(AxZ))

  # Retain the stage-2 residual covariance and confounder-to-outcome loadings.
  # Likewise, Sigma_Y_hat uses the raw (penalized, non-debiased) Beta
  Sigma_Y_hat <- crossprod(outcomes - X_design %*% mod_Stage2$Beta)/n
  dimnames(Sigma_Y_hat) <- list(colnames(outcomes), colnames(outcomes))
  eta_hat <- mod_Stage2$Beta[(p+1):(p+l), , drop = FALSE]
  dimnames(eta_hat) <- list(colnames(confounders), colnames(outcomes))


  #### DE and PIDEs
  trt_row_idx <- p + l + 1
  mod_DE <- mod_Stage2$Beta[trt_row_idx, ]

  mod_stage2_pide <- mod_Stage2_all[1:p,]
  mod_stage1_pide <- mod_Stage1$estimate

  mod_pide_mat <- matrix(nrow = p, ncol = q)
  for(i in 1:q) mod_pide_mat[,i] <- mod_stage1_pide*mod_stage2_pide[,i]

  rownames(mod_pide_mat) <- rownames(mod_Stage2_all)[1:p]; colnames(mod_pide_mat) <- colnames(mod_Stage2_all)
  mod_tide_mat <- colSums(mod_pide_mat)

  #### Calculate moderated PIDEs/TIDEs/DEs (index of moderated mediation) ####
  # For moderator Z_r: interaction PIDE = delta_M[,r] %x% phi_hat,
  #                    interaction TIDE = its column sums (index of moderated mediation),
  #                    interaction DE = the raw (non-debiased) Beta row for that moderator's A x Z_r column
  if(has_moderators){
    mod_intpide_list <- lapply(seq_len(r), function(rz) sweep(mod_stage2_pide, 1, delta_M[,rz], `*`))
    mod_inttide_list <- lapply(mod_intpide_list, colSums)
    mod_intDE_mat <- mod_Stage2$Beta[(p+l+2):(p+l+1+r), , drop = FALSE]
    rownames(mod_intDE_mat) <- colnames(AxZ)
  }

  n_main_rows <- p*q + 2*q
  n_rows_total <- n_main_rows * (1 + r)

  mod_origFit_mat <- matrix(nrow = n_rows_total, ncol = 1)
  main_rownames <- c(paste0(rownames(mod_pide_mat),"_ide_resp",rep(1:q, times = rep(p,q))),
                     paste0("TIDE_resp",1:q),
                     paste0("DE_resp", 1:q))
  mod_origFit_mat[1:n_main_rows, 1] <- c(vec(mod_pide_mat),unname(mod_tide_mat),unname(mod_DE))
  all_rownames <- main_rownames

  if(has_moderators){
    for(rz in seq_len(r)){
      modname <- colnames(moderators)[rz]
      mod_rownames <- c(paste0(rownames(mod_pide_mat),"_ide_MOD_",modname,"_resp",rep(1:q, times = rep(p,q))),
                        paste0("TIDE_MOD_",modname,"_resp",1:q),
                        paste0("DE_MOD_",modname,"_resp", 1:q))
      mod_vals <- c(vec(mod_intpide_list[[rz]]), unname(mod_inttide_list[[rz]]), unname(mod_intDE_mat[rz,]))
      block_idx <- (rz*n_main_rows+1):((rz+1)*n_main_rows)
      mod_origFit_mat[block_idx, 1] <- mod_vals
      all_rownames <- c(all_rownames, mod_rownames)
    }
  }
  rownames(mod_origFit_mat) <- all_rownames
  colnames(mod_origFit_mat) <- "Orig_Est"

  #### Draw bootstrap samples
  B_draws <- matrix(nrow = n, ncol = nB)
  set.seed(seed)
  for(i in 1:nB) {B_draws[,i] <- sample(1:n, n, replace=T)}

  # run one bootstrap replicate and return the length-n_rows_total result vector for that replicate
  run_boot_iter <- function(i){
    boot_sample <- B_draws[,i]
    confounders_boot <- confounders[boot_sample,]
    trt_boot <- trt[boot_sample]
    mediators_boot <- mediators[boot_sample,]
    outcomes_boot <- outcomes[boot_sample,]
    if(has_moderators){
      moderators_boot <- moderators[boot_sample, , drop = FALSE]
      AxZ_boot <- moderators_boot * as.vector(trt_boot)
      colnames(AxZ_boot) <- colnames(AxZ)
    }
    X_design_boot <- if(has_moderators) cbind(mediators_boot, confounders_boot, trt_boot, AxZ_boot) else cbind(mediators_boot, confounders_boot, trt_boot)

    ## Stage 1 fit
    mod_stage1_boot_fit <- if(has_moderators) lm(mediators_boot ~ trt_boot + confounders_boot + AxZ_boot) else lm(mediators_boot ~ trt_boot + confounders_boot)
    mod_stage1_boot_tidy <- broom::tidy(mod_stage1_boot_fit)
    mod_stage1_boot <- mod_stage1_boot_tidy[which(mod_stage1_boot_tidy$term=="trt_boot"),]
    if(has_moderators){
      delta_M_boot <- extract_interaction_estimates(mod_stage1_boot_tidy, "AxZ_boot", colnames(AxZ_boot), colnames(mediators_boot))
    }

    ## Stage 2 fit
    mod_stage2_boot_fit <- MSGLasso(X.m = X_design_boot,
                                    Y.m = outcomes_boot,
                                    grpWTs, Pen_L, Pen_G, PQgrps, GRgrps,
                                    grp_Norm0, MSGLassolam1, MSGLassolamG.m)
    mod_stage2_boot_all <- mod_stage2_boot_fit$Beta + (1/n) * theta_mod %*% t(X_design_boot) %*%
      (outcomes_boot - X_design_boot %*% mod_stage2_boot_fit$Beta)

    # mod_stage2_boot_debiased <- mod_stage2_boot_fit$Beta[1:p,] + (1/n) * theta_mod %*% t(mediators) %*% (outcomes - mediators %*% mod_stage2_boot_fit$Beta[1:p,])
    # mod_stage2_boot_all <- rbind(mod_stage2_boot_debiased, mod_stage2_boot_fit$Beta[(p+1):nrow(mod_stage2_boot_fit$Beta),])

    rownames(mod_stage2_boot_all) <- colnames(X_design_boot)
    colnames(mod_stage2_boot_all) <- colnames(outcomes_boot)

    ## PIDE and DE calculation
    mod_boot_DE <- unname(mod_stage2_boot_all["trt_boot",])

    mod_stage2_pide_boot <- mod_stage2_boot_all[1:p,]
    mod_stage1_pide_boot <- mod_stage1_boot$estimate

    mod_pide_mat_boot <- matrix(nrow = p, ncol = q)
    for(j in 1:q) mod_pide_mat_boot[,j] <- mod_stage1_pide_boot*mod_stage2_pide_boot[,j]

    mod_tide_mat_boot_sums <- colSums(mod_pide_mat_boot)

    boot_result <- numeric(n_rows_total)
    boot_result[1:n_main_rows] <- c(vec(mod_pide_mat_boot),mod_tide_mat_boot_sums,mod_boot_DE)

    ## Moderated PIDE/TIDE/DE calculation (mirrors the main fit, per moderator)
    if(has_moderators){
      for(rz in seq_len(r)){
        intpide_boot <- sweep(mod_stage2_pide_boot, 1, delta_M_boot[,rz], `*`)
        inttide_boot <- colSums(intpide_boot)
        intDE_boot <- unname(mod_stage2_boot_all[colnames(AxZ_boot)[rz], ])
        block_idx <- (rz*n_main_rows+1):((rz+1)*n_main_rows)
        boot_result[block_idx] <- c(vec(intpide_boot), inttide_boot, intDE_boot)
      }
    }

    return(boot_result)
  }

  boot_cores_use <- max(1, min(boot_cores, nB))
  use_parallel <- boot_parallel && boot_cores_use > 1

  start_time <- proc.time()

  if(use_parallel){
    chunk_id <- cut(seq_len(nB), boot_cores_use, labels = FALSE)
    chunks <- split(seq_len(nB), chunk_id)

    clust <- parallel::makeCluster(boot_cores_use)
    doParallel::registerDoParallel(clust)
    on.exit(parallel::stopCluster(clust), add = TRUE)

    # Replicate the Pen_L/Pen_G global-environment workaround
    parallel::clusterExport(clust, varlist = c("Pen_L", "Pen_G"))

    cat(paste0("Running ", nB, " bootstrap replicates across ", boot_cores_use,
               " parallel workers...\n"))

    chunk_idx <- NULL # avoid an R CMD check NOTE for the foreach loop variable
    mod_bootRes <- foreach::foreach(chunk_idx = chunks, .combine = "cbind",
                                    .packages = c("MSGLasso", "broom", "HDMVMed")) %dopar% {
                                      sapply(chunk_idx, run_boot_iter)
                                    }

    cat(paste0("Done. Time elapsed: ",
               round((proc.time()-start_time)[3]/60, 3), " minutes.\n"))
  } else {
    mod_bootRes <- matrix(nrow = n_rows_total, ncol = nB)
    for(i in 1:nB){
      mod_bootRes[,i] <- run_boot_iter(i)

      if(i %% 100 == 0){
        cat(paste0("Done with boot sample ", i, "; ", nB-i, " remaining. Time elapsed: ",
                   round(-1*(start_time[3] - proc.time()[3])/60, 3), " minutes. Apx. ",
                   round(((-1*(start_time[3] - proc.time()[3])/60)/i)*(nB-i), 3), " minutes remaining\n"))
      }
    }
  }
  rownames(mod_bootRes) <- all_rownames
  colnames(mod_bootRes) <- paste0("BootDraw",1:nB)

  bad_entries <- is.nan(mod_bootRes) | is.na(mod_bootRes)
  if(any(bad_entries)){
    row_means_boot <- rowMeans(mod_bootRes, na.rm = TRUE)
    mod_bootRes[bad_entries] <- row_means_boot[row(mod_bootRes)[bad_entries]]
  }

  ## bootstrap p-values
  threshold <- as.vector(2*mod_origFit_mat)
  is_nonneg <- as.vector(mod_origFit_mat) >= 0
  ge_counts <- rowSums(mod_bootRes >= threshold, na.rm = TRUE)
  lt_counts <- rowSums(mod_bootRes <  threshold, na.rm = TRUE)
  boot_counts <- ifelse(is_nonneg, ge_counts, lt_counts) # length n_rows_total; ifelse() is fine at this size
  mod_boot_pvals <- matrix(pmin(2*boot_counts/nB, 1), ncol = 1)

  ## BCa confidence intervals
  mod_boot_bcaCI <- t(apply(mod_bootRes, 1, coxed::bca, conf.level = 1-alpha))
  mod_boot_bcaCI[is.nan(mod_boot_bcaCI)] <- 0

  ## pivotal bootstrap CIs
  mod_boot_perCI <- 2*matrix(data = c(mod_origFit_mat, mod_origFit_mat), ncol = 2) - t(apply(mod_bootRes, 1, quantile, c(1-alpha/2,alpha/2)))

  ## Summary table
  mod_boot_rowsd <- sqrt(rowSums((mod_bootRes - rowMeans(mod_bootRes))^2) / (ncol(mod_bootRes) - 1))
  mod_boot_summ <- cbind(mod_origFit_mat, rowMeans(mod_bootRes), mod_boot_rowsd,
                         mod_boot_pvals, mod_boot_bcaCI, mod_boot_perCI)
  colnames(mod_boot_summ) <- c("OrigEst", "Mean_boot", "boot_SE", "boot_pval",
                               "bca_lowerCL", "bca_upperCL",
                               "per_lowerCL", "per_upperCL")
  mod_boot_summ[is.nan(mod_boot_summ)] <- 0

  alpha_hat <- mod_Stage1$estimate
  names(alpha_hat) <- colnames(mediators)

  # MSGLasso's Beta carries no column names, so set them explicitly rather than
  # relying on what propagated through the debiasing arithmetic.
  phi_hat_d <- mod_Stage2_all[1:p, , drop = FALSE]
  dimnames(phi_hat_d) <- list(colnames(mediators), colnames(outcomes))
  tau_hat <- mod_DE; names(tau_hat) <- colnames(outcomes)

  internals <- list(alpha_hat = alpha_hat, Sigma_M_hat = Sigma_M_hat,
                    Sigma_Y_hat = Sigma_Y_hat, phi_hat_d = phi_hat_d,
                    Theta = theta_mod, xi_hat = xi_hat, eta_hat = eta_hat,
                    tau_hat = tau_hat, confounders = confounders,
                    n = n, p = p, l = l, q = q, has_moderators = has_moderators)

  mod_boot <- list("res" = mod_boot_summ,
                   "internals" = internals)

  return(mod_boot)
}
