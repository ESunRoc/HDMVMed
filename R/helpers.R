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

  X_std <- scale(X)
  X_std[is.na(X_std)] <- 0

  Theta <- matrix(0, p, p)
  tau_sq <- numeric(p)

  for (j in 1:p) {
    result <- theta_column(X_std, j)
    Theta[, j] <- result$theta_col
    tau_sq[j] <- result$tau_sq
  }

  return(Theta)
}

# Compute column j of Theta \approx solve(t(X) %*% X)
#'
#' Column-wise sample covariance
#'
#' `theta_column()` computes the j-th nodewise regression for use in `theta_calc()`.
#'
#' @param X A numeric matrix.
#' @param j An integer between 1 and `ncol(X)` indicating the column on which to regress \eqn{X_{-j}}.
#'
#' @returns A list containing:
#' * `theta_col`:
#' * `tau_sq`:
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
                      standardize = FALSE, nfolds = 5)
  lambda_opt <- cv_fit$lambda.min

  lasso_fit <- glmnet::glmnet(X_minus_j, y_j, lambda = lambda_opt,
                      intercept = FALSE, standardize = FALSE)
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
  X_std <- scale(X); X_std[is.na(X_std)] <- 0

  n <- nrow(X)
  res_mat <- matrix(nrow = ncol(X), ncol = ncol(X)+2)

  clust <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(clust)

  res_mat <- foreach::foreach(j = 1:ncol(X), .combine = "rbind", .packages = "glmnet") %dopar% {
    lasso_j <- glmnet::cv.glmnet(x = X[,-j], y = X[,j], nfolds = folds, intercept = F, standardize = F)

    Chat_j <- numeric(ncol(X))
    Chat_j[-j] <- lasso_j$glmnet.fit$beta[,lasso_j$index[2,]]
    Chat_j[j] <- 1

    lambda_1se <- lasso_j$lambda.1se
    gamma_j <- lasso_j$glmnet.fit$beta[,lasso_j$index[2,]]

    tau2_j <- norm(X[,j]-(X[,-j] %*% as.matrix(gamma_j)), "2")/n + lambda_1se * norm(as.matrix(gamma_j), "1")

    c(lambda_1se, tau2_j, Chat_j)
  }

  stopCluster(clust)
  colnames(res_mat) <- c("lambda.1se", "tau2", paste0("gamma", 1:(ncol(res_mat)-2)))
  rownames(res_mat) <- colnames(X)

  tau2_mat <- res_mat[,"tau2"]
  # lambda_mat <- res_mat[,"lambda.1se"]

  That2_inv <- diag(1/tau2_mat)
  Chat_mat <- res_mat[,-which(colnames(res_mat)%in%c("lambda.1se", "tau2"))]

  Theta <- That2_inv %*% Chat_mat

  return(unname(Theta))
}
