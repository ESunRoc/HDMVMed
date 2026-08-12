#' Vectorization
#'
#' `vec()` computes the vectorization of a given matrix; the isomorphism \eqn{\text{vec}:\mathbb R^{m\times n}\to\mathbb R^{mn}}.
#'
#' @param Mat A matrix.
#'
#' @returns The vectorization of `Mat`.
#'
#' @examples
#' data <- matrix(1:9, nrow = 3, ncol = 3)
#' vec(data)
#'
#' @export
vec <- function(Mat) return(t(t(as.vector(Mat))))


#' Sample covariance in serial
#'
#' `theta_calc()` computes the sparse sample covariance matrix, \eqn{\Theta\approx(\boldsymbol{X}^\top\boldsymbol{X})^{-1}}, via nodewise regression.
#'
#' @param X A numeric matrix.
#'
#' @returns The sparse sample covariance matrix estimated via the nodewise regression of van de Geer et al. (2014).
#'
#' @details
#' `Theta` is estimated directly for the design matrix `X` as supplied, i.e. `Theta`
#' approximates \eqn{(\boldsymbol{X}^\top\boldsymbol{X})^{-1}} for `X` on its original scale.
#' This is deliberate: [bootstrap_model()] later forms the debiasing correction as
#' `Theta %*% t(X) %*% (Y - X %*% Beta)` using this same (raw, unstandardized) `X`, so
#' `Theta` must be computed on that identical scale for the correction to be valid.
#' Internally, each nodewise LASSO regression is fit with `standardize = TRUE` purely
#' for numerical conditioning (candidate mediators, confounders, and treatment can have
#' very different scales); `glmnet` returns coefficients back-transformed to the
#' original scale of `X`, so no manual re-scaling of `Theta` is required afterwards.
#'
#' @examples
#' \dontrun{
#'## Load the toy data
#'data(hdmvmed_test_data)
#'
#'## Assign mediator/confounder/trt matrices
#'mediators <- hdmvmed_test_data[,7:106]
#'confounders <- hdmvmed_test_data[,107:111]
#'trt <- hdmvmed_test_data[,1]
#'
#'## Calculate column j=3 of Theta
#'theta_calc(X = cbind(mediators, confounders, trt))
#'}
#'
#' @references{
#' van de Geer, S., Bühlmann, P., Ritov, Y., and Dezeure, R. (2014). On Asymptotically Optimal Confidence Regions
#' and Tests for High-Dimensional Models. \emph{The Annals of Statistics}, \bold{72}(3), 1166-1202.
#' }
#'
#' @references{
#' Sun, E., Xiao, J., and Wu, T. T. (2025). Causal Mediation Analysis for Multiple
#' Outcomes and High-dimensional Mediators: Identification, Inference, and Application.
#' \emph{Biometrics}. _Under Review._
#' }
#'
#' @export
theta_calc <- function(X){
  n <- nrow(X)
  p <- ncol(X)

  Theta <- matrix(0, p, p)
  tau_sq <- numeric(p)

  for (j in 1:p) {
    result <- theta_column(X, j)
    # FIX (transpose bug): the j-th nodewise regression populates ROW j of Theta
    # (Theta = T^{-2} C, with C's j-th row built from the j-th nodewise regression;
    # see eq. (2.8)-(2.14) of Sun et al. and eq. (13) of the Supplementary Materials).
    # The previous version wrote this into column j instead, silently returning
    # t(Theta) rather than Theta whenever the nodewise LASSO fits were not exactly
    # symmetric (i.e. essentially always, since Theta is only symmetric in the
    # population and the estimator has no such constraint built in).
    Theta[j, ] <- result$theta_col
    tau_sq[j] <- result$tau_sq
  }

  return(Theta)
}

# Compute row j of Theta \approx solve(t(X) %*% X)
#'
#' Column-wise sample covariance
#'
#' `theta_column()` computes the j-th nodewise regression for use in `theta_calc()`.
#'
#' @param X A numeric matrix.
#' @param j An integer between 1 and `ncol(X)` indicating the column on which to regress \eqn{X_{-j}}.
#'
#' @returns A list containing:
#' * `theta_col`: the j-th nodewise-regression vector, i.e. row j of \eqn{\hat\Theta = \hat T^{-2}\hat C}
#'   (diagonal entry \eqn{1/\hat\tau_j^2}, off-diagonal entries \eqn{-\hat\gamma_{j,k}/\hat\tau_j^2}).
#' * `tau_sq`: the nodewise-regression scale estimate \eqn{\hat\tau_j^2 = n^{-1}\Vert X_j-X_{-j}\hat\gamma_j\Vert_2^2+\lambda_j\Vert\hat\gamma_j\Vert_1}.
#'
#' @references{
#' van de Geer, S., Bühlmann, P., Ritov, Y., and Dezeure, R. (2014). On Asymptotically Optimal Confidence Regions
#' and Tests for High-Dimensional Models. \emph{The Annals of Statistics}, \bold{72}(3), 1166-1202.
#' }
#'
#' @examples
#' \dontrun{
#' ## Load the toy data
#' data(hdmvmed_test_data)
#'
#' ## Assign mediator/confounder/trt matrices
#' mediators <- hdmvmed_test_data[,7:106]
#' confounders <- hdmvmed_test_data[,107:111]
#' trt <- hdmvmed_test_data[,1]
#'
#' ## Calculate column j=3 of Theta
#' theta_column(X = cbind(mediators, confounders, trt), j = 3)
#'}
#'
#' @export
theta_column <- function(X, j){
  p <- ncol(X); n <- nrow(X)

  y_j <- X[, j]
  X_minus_j <- X[, -j, drop = FALSE]

  cv_fit <- glmnet::cv.glmnet(X_minus_j, y_j, intercept = FALSE,
                              standardize = TRUE, nfolds = 5)
  lambda_opt <- cv_fit$lambda.min

  lasso_fit <- glmnet::glmnet(X_minus_j, y_j, lambda = lambda_opt,
                              intercept = FALSE, standardize = TRUE)
  gamma_j <- as.numeric(lasso_fit$beta)

  residuals <- y_j - X_minus_j %*% gamma_j
  tau_sq_j <- mean(residuals^2) + lambda_opt * sum(abs(gamma_j))

  theta_col <- numeric(p)
  theta_col[j] <- 1 / tau_sq_j

  # Fill off-diagonal elements
  k_idx <- 1
  for (k in 1:p) {
    if (k != j) {
      theta_col[k] <- -gamma_j[k_idx] / tau_sq_j
      k_idx <- k_idx + 1
    }
  }

  return(list(theta_col = theta_col, tau_sq = tau_sq_j))
}


#' Sample covariance in parallel
#'
#' `theta_calc_parallel()` computes the sparse sample covariance matrix, \eqn{\Theta\approx(\boldsymbol{X}^\top\boldsymbol{X})^{-1}}, via nodewise regression. The parallel back-end is implemented via `doParallel` and `foreach`.
#'
#' @param X A numeric matrix.
#' @param cores An integer between 1 and `parallel::detectCores()` indicating the number of logical cores to use for parallel computation; defaults to one less than the number of cores on the local system.
#' @param folds An integer between 1 and `nrow(X)` indicating the number of CV folds to use in estimating the nodewise regressions; defaults to 5-fold CV.
#'
#' @returns The sparse sample covariance matrix estimated via the nodewise regression of van de Geer et al. (2014).
#'
#' @details
#' As in [theta_calc()], `Theta` is estimated directly for `X` as supplied (raw scale);
#' `standardize = TRUE` is used only to condition each nodewise LASSO fit, and
#' `glmnet` returns coefficients on the original scale of `X`. This mirrors
#' [theta_calc()] exactly (same tuning-parameter rule, `lambda.min`, and the same
#' \eqn{\hat\tau_j^2} formula) so that the serial and parallel implementations compute
#' the same estimator and differ only in their execution back-end.
#'
#' @examples
#' \dontrun{
#' ## Load the toy data
#' data(hdmvmed_test_data)
#'
#' ## Assign mediator/confounder/trt matrices
#' mediators <- hdmvmed_test_data[,7:106]
#' confounders <- hdmvmed_test_data[,107:111]
#' trt <- hdmvmed_test_data[,1]
#'
#' ## Calculate Theta
#' theta_calc_parallel(X = cbind(mediators, confounders, trt), cores = 4, folds = 5)
#'}
#'
#' @references{
#' van de Geer, S., Bühlmann, P., Ritov, Y., and Dezeure, R. (2014). On Asymptotically Optimal Confidence Regions
#' and Tests for High-Dimensional Models. \emph{The Annals of Statistics}, \bold{72}(3), 1166-1202.
#' }
#'
#' @references{
#' Sun, E., Xiao, J., and Wu, T. T. (2025). Causal Mediation Analysis for Multiple
#' Outcomes and High-dimensional Mediators: Identification, Inference, and Application.
#' \emph{Biometrics}. _Under Review._
#' }
#'
#' @export
theta_calc_parallel <- function(X, cores = parallel::detectCores()-1, folds = 5){
  j <- NULL

  n <- nrow(X)
  res_mat <- matrix(nrow = ncol(X), ncol = ncol(X)+2)

  clust <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(clust)

  res_mat <- foreach::foreach(j = 1:ncol(X), .combine = "rbind", .packages = "glmnet") %dopar% {
    lasso_j <- glmnet::cv.glmnet(x = X[,-j], y = X[,j], nfolds = folds, intercept = F, standardize = TRUE)

    Chat_j <- numeric(ncol(X))
    Chat_j[-j] <- lasso_j$glmnet.fit$beta[,lasso_j$index[1,]]
    Chat_j[j] <- 1

    lambda_min <- lasso_j$lambda.min
    gamma_j <- lasso_j$glmnet.fit$beta[,lasso_j$index[1,]]

    resid_j <- X[,j] - (X[,-j] %*% as.matrix(gamma_j))
    tau2_j <- mean(resid_j^2) + lambda_min * sum(abs(gamma_j))

    c(lambda_min, tau2_j, Chat_j)
  }

  stopCluster(clust)
  colnames(res_mat) <- c("lambda.min", "tau2", paste0("gamma", 1:(ncol(res_mat)-2)))
  rownames(res_mat) <- colnames(X)

  tau2_mat <- res_mat[,"tau2"]
  # lambda_mat <- res_mat[,"lambda.min"]

  That2_inv <- diag(1/tau2_mat)
  Chat_mat <- res_mat[,-which(colnames(res_mat)%in%c("lambda.min", "tau2"))]

  Theta <- That2_inv %*% Chat_mat

  return(unname(Theta))
}

#' Extract moderator interactions
#'
#' Extract treatment x moderator interaction coefficients from a `broom::tidy()` summary of a multivariate `lm()` fit (i.e. a `mediators ~ trt + confounders + AxZ` model), one column per moderator and one row per mediator.
#'
#' @param tidy_df A data frame as returned by `broom::tidy()` on the mlm fit; must have `term` and `response` columns.
#' @param varname A character string for the name of the interaction-matrix argument as it appears in the model formula (e.g. `"AxZ"` or `"AxZ_boot"`); this reproduces the term names R assigns to matrix predictors when combined with each entry of `int_colnames`
#' @param int_colnames A character vector of the column names of the interaction matrix (one per moderator).
#' @param mediator_names A character vector of mediator names; used to both size the output and to align rows by `response` rather than assuming row order, in case `tidy()`'s row order ever changes.
#'
#' @returns A `length(mediator_names)`-by-`length(int_colnames)` matrix.
#'
#' @export
extract_interaction_estimates <- function(tidy_df, varname, int_colnames, mediator_names){
  out <- matrix(nrow = length(mediator_names), ncol = length(int_colnames))
  rownames(out) <- mediator_names
  colnames(out) <- int_colnames
  for(rz in seq_along(int_colnames)){
    term_name <- paste0(varname, int_colnames[rz])
    if(!(term_name %in% tidy_df$term) && length(int_colnames) == 1L) term_name <- varname
    rows <- tidy_df[tidy_df$term == term_name, ]
    out[, rz] <- rows$estimate[match(mediator_names, rows$response)]
  }
  return(out)
}
