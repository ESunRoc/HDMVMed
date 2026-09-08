#' Sensitivity analysis for unmeasured mediator-outcome confounding
#'
#' `hdmvm_sens()` performs a full sensitivity analysis for a fit [bootstrap_model()].
#' It returns robustness values for every reported estimand, the selection robustness
#' values, and the benchmark values implied by the measured confounders.
#'
#' @param fit The list returned by [bootstrap_model()], with elements `res` and
#'   `internals`.
#' @param mediators,trt Optional overrides for the mediator matrix and exposure vector.
#'   [bootstrap_model()] retains both in `internals`, so these are normally left `NULL`
#'   and taken from the fit; supply them only to compute \eqn{Q} on data other than that
#'   used for the fit.
#' @param ci_type Which reported interval to use, `"bca"` (default) or `"percentile"`.
#' @param level Confidence level of the reported intervals; defaults to 0.95.
#' @param Q_method Method for estimating the design scalar \eqn{Q = \alpha\Sigma_M^{-1}\alpha^\top}.
#'  Options are `"plugin"`, `"crossfit"`, or `"both"`; defaults to `"both"`. See details.
#' @param signal_const A minimum-signal threshold  \eqn{\theta_n = c\sqrt{s\log(pq)/n}};
#'   defaults to 1. This should be varied as part of any report; see [hdmvm_sens_signal_sensitivity()].
#' @param report_only A Boolean. If `TRUE` (the default), the erosion robustness values are
#'   computed over *reported findings* -- pairs that are both selected and whose interval
#'   excludes zero -- rather than over every selected pair. Restricting in this way matters
#'   whenever the penalty is loose: a selected set dominated by pairs with
#'   \eqn{|\hat\phi^d| < \theta_n} drives both the minimum and any low quantile to zero, so
#'   the unrestricted statistic reports the fragility of coefficients no analyst would
#'   present. The unrestricted values are returned alongside as `RV_erode_all`.
#' @param crossfit_folds Number of folds for the cross-fitted \eqn{R^2_{A\sim M}}.
#' @param seed Seed for the cross-fitting split.
#'
#' @returns An object of class `hdmvm_sens` with components:
#'   * `tide`, `de` -- data frames of robustness values per outcome;
#'   * `pide` -- data frame of robustness values per mediator-outcome pair;
#'   * `selection` -- the erosion and spurious-selection robustness values;
#'   * `benchmark` -- sensitivity parameters implied by each measured confounder;
#'   * `Q`, `R2_AM` -- the design scalar and the exposure-prediction \eqn{R^2};
#'   * `meta` -- dimensions, thresholds, and settings.
#'
#' @details
#' **The sensitivity parameters.** A hypothesized unmeasured confounder \eqn{Z},
#' characterized by two bounded numbers: \eqn{R^2_M}, the population \eqn{R^2} from
#' regressing \eqn{Z} on the mediator residuals, shared across outcomes; and
#' \eqn{R^2_Y(k)}, the fraction of outcome-\eqn{k} residual variance attributable to
#' \eqn{Z}. When \eqn{p=q=1} these collapse to the residual correlation of Imai and
#' Yamamoto (2013) by \eqn{\rho = \pm\sqrt{R^2_MR^2_Y}}.
#'
#' **Robustness values.** Reported on the symmetric scale \eqn{R^2_M = R^2_Y = RV}.
#' For an estimand with reported interval \eqn{[L,U]} and worst-case bias
#' \eqn{\kappa\sqrt{R^2_MR^2_Y}}, the robustness value is
#' \eqn{\min\{\min(|L|,|U|)/\kappa,\,1\}}. A value of 1 implies that no
#' confounder of the assumed form can overturn the finding. For individual coefficients,
#' this reduces to \eqn{(|t|-z)/\sqrt n}, free of any matrix quantity.
#'
#' **The design scalar.** \eqn{Q=R^2_{A\sim M}/[\sigma_A^2(1-R^2_{A\sim M})]} governs the
#' worst-case bias in every bound. A design in which the mediators nearly determine
#' the exposure is intrinsically fragile. The plug-in estimator \eqn{\hat\alpha\hat\Theta_M\hat\alpha^\top}
#' is upward biased (conservative) because \eqn{\hat\alpha} is dense; the cross-fitted estimator is
#' downward biased (anti-conservative) through regularization shrinkage. `Q_method = "both"`, the default,
#' takes the larger. The cross-fitted estimator is available without further input, at the cost of one cross-validated
#' lasso fit per fold. Use `Q_method = "plugin"` where that cost matters.
#'
#' Moderated mediation is not (presently) supported: the bias formulas are derived for the
#' unmoderated model.
#'
#' @references{
#' Cinelli, C. and Hazlett, C. (2020). Making Sense of Sensitivity: Extending Omitted
#' Variable Bias. \emph{JRSS-B}, \bold{82}(1), 39-67.
#' }
#' @references{
#' Imai, K., Keele, L., and Yamamoto, T. (2010). Identification, Inference and
#' Sensitivity Analysis for Causal Mediation Effects. \emph{Statistical Science},
#' \bold{25}(1), 51-71.
#' }
#' @references{
#' Sun, E., Xiao, J., and Wu, T. T. (2026). Causal Mediation Analysis for Multiple
#' Outcomes and High-dimensional Mediators. \emph{Biostatistics}, under review.
#' }
#'
#' @examples
#' \dontrun{
#' data(hdmvmed_test_data)
#' fit <- bootstrap_model(mediators = hdmvmed_test_data[, 7:106],
#'                        confounders = hdmvmed_test_data[, 107:111],
#'                        trt = hdmvmed_test_data[, 1],
#'                        outcomes = hdmvmed_test_data[, 2:6], nB = 500)
#' s <- hdmvm_sens(fit) # mediators/trt taken from fit$internals
#' s                    # print method: headline table
#' plot(s, outcome = 1) # contour plot in the (R2_M, R2_Y) plane
#' }
#'
#' @seealso [hdmvm_sens_bias()], [hdmvm_sens_signal_sensitivity()]
#' @export
hdmvm_sens <- function(fit, mediators = NULL, trt = NULL, ci_type = c("bca", "percentile"), level = 0.95,
                       Q_method = c("both", "plugin", "crossfit"), signal_const = 1, sel_quantile = 0.1,
                       report_only = TRUE, crossfit_folds = 5, seed = 823543){

  ci_type  <- match.arg(ci_type); Q_method <- match.arg(Q_method)
  if(is.null(fit$internals)){
    stop("`fit` has no `internals`. `bootstrap_model()` returns it as a named sublist; check that `fit` is the ",
         "unmodified result of the call.")
  }

  int <- fit$internals; res <- fit$res

  if(isTRUE(int$has_moderators)) stop("Sensitivity analysis for moderated mediation is not implemented.")

  n <- int$n; p <- int$p; l <- int$l; q <- int$q; z <- stats::qnorm(1-(1-level)/2)

  mnames <- names(int$alpha_hat); if(is.null(mnames)) mnames <- paste0("M", 1:p)
  onames <- colnames(int$Sigma_Y_hat); if(is.null(onames)) onames <- paste0("Y", 1:q)

  cols <- if(ci_type == "bca") c("bca_lowerCL", "bca_upperCL") else c("per_lowerCL", "per_upperCL")

  ##### Precompute needed matrices #####
  Theta_M   <- int$Theta[seq_len(p), seq_len(p), drop = FALSE]
  theta_jj  <- diag(Theta_M); sigma_Yk  <- sqrt(diag(int$Sigma_Y_hat))
  alpha_hat <- int$alpha_hat

  ## override internal mediators and treatment --- OLD
  if(is.null(mediators)) mediators <- int$mediators
  if(is.null(trt)) trt <- int$trt
  Qest <- .hdmvm_Q(int, Q_method, crossfit_folds, seed, mediators, trt)
  Q <- Qest$Q

  ##### TIDE and DE #####
  ## Worst-case bias per unit sqrt(R2_M*R2_Y): kappa_k = sigma_Yk * sqrt(Q)
  kappa_agg <- sigma_Yk * sqrt(Q)
  tide <- .hdmvm_rv_table(res, paste0("TIDE_resp", seq_len(q)), cols, onames, kappa_agg, "TIDE")
  de   <- .hdmvm_rv_table(res, paste0("DE_resp",   seq_len(q)), cols, onames, kappa_agg, "DE")

  ##### PIDEs #####
  ## kappa_jk = |alpha_j| * sigma_Yk * sqrt(Theta_jj)
  rows <- paste0(rep(mnames, times = q), "_ide_resp", rep(seq_len(q), each = p))
  have <- all(rows %in% rownames(res))
  if(have){
    getm <- function(cn) matrix(res[rows, cn], nrow = p, ncol = q, dimnames = list(mnames, onames))
    est <- getm("OrigEst"); a <- getm(cols[1]); b <- getm(cols[2])
    lo <- pmin(a, b); hi <- pmax(a, b)
    kap <- outer(abs(alpha_hat) * sqrt(theta_jj), sigma_Yk)
    RV  <- pmin(pmin(abs(lo), abs(hi)) / kap, 1)
    RV[!is.finite(RV)] <- NA_real_
    keep <- which(rowSums(abs(est)) > 0)
    if(!length(keep)) keep <- seq_len(p)
    pide <- data.frame(mediator = rep(mnames[keep], times = q),
                       outcome = rep(onames, each = length(keep)),
                       est = as.vector(est[keep, , drop = FALSE]),
                       lower = as.vector(lo[keep, , drop = FALSE]),
                       upper = as.vector(hi[keep, , drop = FALSE]),
                       excl0 = as.vector(!(lo[keep, , drop = FALSE] <= 0 & hi[keep, , drop = FALSE] >= 0)),
                       RV = as.vector(RV[keep, , drop = FALSE]),
                       row.names = NULL, stringsAsFactors = FALSE)
    pide$reading <- ifelse(pide$excl0, "masking", "spurious")
  } else {
    pide <- NULL
  }

  ##### selection robustness #####
  sel <- .hdmvm_selection_rv(int, res, cols, n, p, q, mnames, onames, signal_const, sel_quantile, report_only)

  ##### benchmarking #####
  bench <- .hdmvm_benchmark(int, Theta_M, sigma_Yk, onames)

  out <- list(tide = tide, de = de, pide = pide, selection = sel,
              benchmark = bench,
              Q = Q, Q_detail = Qest, R2_AM = Qest$R2_AM,
              meta = list(n = n, p = p, l = l, q = q, z = z, level = level,
                          Sigma_Y = int$Sigma_Y_hat,
                          ci_type = ci_type, kappa_agg = kappa_agg,
                          theta_n = sel$theta_n, signal_const = signal_const,
                          sigma_Yk = sigma_Yk, theta_jj = theta_jj,
                          alpha_hat = alpha_hat, mnames = mnames, onames = onames))
  class(out) <- "hdmvm_sens"
  out
}


##### Robustness-value table for an aggregate estimand #####
.hdmvm_rv_table <- function(res, rows, cols, onames, kappa, label){
  if(!all(rows %in% rownames(res)))
    stop("`fit$res` does not contain the expected ", label, " rows.")
  est <- res[rows, "OrigEst"]
  a <- res[rows, cols[1]]; b <- res[rows, cols[2]]
  lo <- pmin(a, b); hi <- pmax(a, b)
  RV <- pmin(pmin(abs(lo), abs(hi)) / kappa, 1)
  RV[!is.finite(RV)] <- NA_real_
  excl <- !(lo <= 0 & hi >= 0)
  data.frame(outcome = onames, estimand = label, est = est,
             lower = lo, upper = hi, excl0 = excl, RV = RV,
             reading = ifelse(excl, "masking", "spurious"),
             row.names = NULL, stringsAsFactors = FALSE)
}


##### Design scalar Q #####
## Plug-in:  t(alpha') %*% Theta_M %*% alpha. Upward biased (alpha_hat is dense) => conservative.
## Crossfit: Q = R2_{A~M} / (sigma_A^2 (1 - R2_{A~M})). Downward biased through shrinkage => anti-conservative.
## both:     max(plug-in, crossfit)
.hdmvm_Q <- function(int, method, folds, seed, mediators = NULL, trt = NULL){
  p <- int$p
  Theta_M <- int$Theta[seq_len(p), seq_len(p), drop = FALSE]
  Q_plug <- as.numeric(int$alpha_hat %*% Theta_M %*% int$alpha_hat)

  Q_cf <- NA_real_; R2_AM <- NA_real_
  if(method %in% c("both", "crossfit")){
    ## fallback if glmnet or internals are broken --- OLD
    ok <- requireNamespace("glmnet", quietly = TRUE) &&
      !is.null(mediators) && !is.null(trt)
    if(!ok && method == "crossfit")
      message("Cross-fitted Q needs `mediators` and `trt` in `internals` (and glmnet); ",
              "falling back to the plug-in estimator, which is conservative.")
    if(ok){
      cf <- try(.hdmvm_crossfit_R2(mediators, int$confounders, trt,
                                   folds = folds, seed = seed), silent = TRUE)
      if(!inherits(cf, "try-error")){
        R2_AM <- cf$R2
        if(R2_AM < 1 - 1e-8) Q_cf <- R2_AM / (cf$sigma2_A * (1 - R2_AM))
      }
    }
  }

  Q <- switch(method, plugin = Q_plug,
              crossfit = if(is.finite(Q_cf)) Q_cf else Q_plug,
              both = max(c(Q_plug, Q_cf), na.rm = TRUE))
  list(Q = Q, Q_plugin = Q_plug, Q_crossfit = Q_cf, R2_AM = R2_AM, method = method)
}


##### Cross-fitted linear-projection R^2 of A on M given U #####
## NB: the linear-projection R^2 is required even when A is binary; a logistic
## pseudo-R^2 is a different functional and invalidates the identity Q = R2/(s2(1-R2)).
.hdmvm_crossfit_R2 <- function(mediators, confounders, trt, folds = 5, seed = 823543){
  set.seed(seed)
  M <- as.matrix(mediators); U <- as.matrix(confounders); A <- as.numeric(trt); n <- length(A)

  ## residualize on U so the R^2 is the incremental contribution of M
  A_r <- stats::residuals(stats::lm(A ~ U))
  M_r <- stats::residuals(stats::lm(M ~ U))
  sigma2_A <- mean(A_r^2)

  fold <- sample(rep(seq_len(folds), length.out = n))
  sse  <- 0
  for(v in seq_len(folds)){
    tr <- fold != v; te <- !tr
    cvf <- glmnet::cv.glmnet(M_r[tr, , drop = FALSE], A_r[tr], alpha = 1)
    pr  <- as.numeric(stats::predict(cvf, newx = M_r[te, , drop = FALSE], s = "lambda.min"))
    sse <- sse + sum((A_r[te]-pr)^2)
  }
  list(R2 = max(0, 1-sse/(n*sigma2_A)), sigma2_A = sigma2_A)
}


##### Selection robustness values #####
## theta_n is the minimum-signal threshold. The bound on the entrywise bias is sqrt(n) * se_jk * sqrt(R2_M R2_Y),
## and se_jk = sqrt(sigma2_Yk*Theta_jj / n), so sqrt(n) * se_jk = sigma_Yk sqrt(Theta_jj)
.hdmvm_selection_rv <- function(int, res, cols, n, p, q, mnames, onames, signal_const,
                                sel_quantile = 0.1, report_only = TRUE){
  phi_d <- int$phi_hat_d
  Theta_M <- int$Theta[seq_len(p), seq_len(p), drop = FALSE]
  sigma_Yk <- sqrt(diag(int$Sigma_Y_hat))
  kap <- outer(sqrt(diag(Theta_M)), sigma_Yk)

  ## Selected support. This must come from the pre-debiased phi_hat
  if(is.null(int$phi_hat))
    stop("`internals$phi_hat` is missing. The selected support cannot be read from ",
         "`phi_hat_d`, which is dense by construction. Refit with a version of ",
         "bootstrap_model() that returns the pre-debiased stage-two coefficients.")
  S <- int$phi_hat != 0
  s_hat <- sum(S)
  src <- "penalized_phi"
  theta_n <- signal_const * sqrt(max(s_hat, 1)*log(p*q)/n)

  ## Reported findings: selected AND interval excluding zero. A selected pair whose
  ## interval covers zero is not something an analyst presents, so including it in an
  ## erosion statistic asks how much confounding would unreport a non-finding. With a
  ## loose penalty such pairs dominate the selected set, every one of them contributes an
  ## erosion value of exactly zero, and both the minimum and any low quantile collapse.
  Rep <- S
  pide_rows <- paste0(rep(mnames, times = q), "_ide_resp", rep(seq_len(q), each = p))
  if(report_only && all(pide_rows %in% rownames(res))){
    a <- matrix(res[pide_rows, cols[1]], nrow = p, ncol = q)
    b <- matrix(res[pide_rows, cols[2]], nrow = p, ncol = q)
    lo <- pmin(a, b); hi <- pmax(a, b)
    Rep <- S & !(lo <= 0 & hi >= 0)
    if(!any(Rep)) Rep <- S                          # nothing reported: fall back
  }

  erode_all <- ifelse(S, pmax(abs(phi_d)-theta_n, 0) / kap, NA_real_)
  erode <- ifelse(Rep, pmax(abs(phi_d)-theta_n, 0) / kap, NA_real_)
  spur <- ifelse(!S, pmax(theta_n-abs(phi_d), 0) / kap, NA_real_)

  wm <- function(m){
    if(all(is.na(m))) return(list(rv = NA_real_, which = c(NA, NA)))
    i <- which(m == min(m, na.rm = TRUE), arr.ind = TRUE)[1, ]
    list(rv = min(min(m, na.rm = TRUE), 1),
         which = c(mediator = mnames[i[1]], outcome = onames[i[2]]))
  }
  e <- wm(erode); ea <- wm(erode_all); s <- wm(spur)
  qtl <- function(m){
    v <- m[!is.na(m)]
    if(!length(v)) return(NA_real_)
    min(as.numeric(stats::quantile(v, probs = sel_quantile, names = FALSE)), 1)
  }

  list(RV_erode = e$rv, which_erode = e$which,
       RV_spurious = s$rv, which_spurious = s$which,
       RV_erode_q = qtl(erode), RV_spurious_q = qtl(spur),
       RV_erode_all = ea$rv, RV_erode_all_q = qtl(erode_all),
       sel_quantile = sel_quantile, report_only = report_only,
       n_selected = sum(S), n_reported = sum(Rep),
       theta_n = theta_n, s_hat = s_hat, support_source = src)
}


##### Exact benchmark from the measured confounders (Proposition 5) #####
.hdmvm_benchmark <- function(int, Theta_M, sigma_Yk, onames){
  if(is.null(int$xi_hat) || is.null(int$eta_hat) || is.null(int$confounders)) return(NULL)
  xi <- int$xi_hat; eta <- int$eta_hat
  a <- apply(as.matrix(int$confounders), 2, stats::var)
  l <- nrow(xi)

  w <- vapply(seq_len(l), function(r) as.numeric(xi[r, ] %*% Theta_M %*% xi[r, ]), numeric(1))
  R2M <- a*w/(1+a*w)

  out <- NULL
  for(r in seq_len(l)) for(k in seq_along(sigma_Yk)){
    R2Y <- a[r] * eta[r, k]^2 / (sigma_Yk[k]^2 + a[r] * eta[r, k]^2)
    out <- rbind(out, data.frame(
      confounder = rownames(xi)[r], outcome = onames[k],
      R2_M = R2M[r], R2_Y = R2Y, geo_mean = sqrt(R2M[r] * R2Y),
      row.names = NULL, stringsAsFactors = FALSE))
  }
  out
}


#' Worst-case bias at specified sensitivity parameters
#'
#' `hdmvm_sens_bias()` returns the largest bias any confounder with the stated
#' associations could induce, for each estimand.
#'
#' @param s An `hdmvm_sens` object.
#' @param R2_M,R2_Y Numerics in \eqn{[0,1]}. `R2_Y` may be a scalar (recycled) or a
#'   vector of length \eqn{q}.
#'
#' @returns A list with `tide`, `de` (length \eqn{q}) and, when available, `pide`
#'   (\eqn{p\times q}) worst-case absolute biases.
#'
#' @details These are the worst-case bounds attained at \eqn{\lambda_M\propto\alpha}
#'   for the TIDE and \eqn{\lambda_M\propto e_j} for individual coefficients.
#'   The DE bound equals the TIDE bound with opposite sign, since confounding reallocates effect
#'   between the channels without altering the total effect.
#'
#' @export
hdmvm_sens_bias <- function(s, R2_M, R2_Y){
  stopifnot(inherits(s, "hdmvm_sens"))
  q <- s$meta$q
  if(length(R2_Y) == 1) R2_Y <- rep(R2_Y, q)
  if(any(c(R2_M, R2_Y) < 0 | c(R2_M, R2_Y) > 1)) stop("R2 values must lie in [0,1].")
  agg <- s$meta$kappa_agg * sqrt(R2_M * R2_Y)
  names(agg) <- s$meta$onames
  pide <- outer(abs(s$meta$alpha_hat) * sqrt(s$meta$theta_jj), s$meta$sigma_Yk * sqrt(R2_Y)) * sqrt(R2_M)
  dimnames(pide) <- list(s$meta$mnames, s$meta$onames)
  list(tide = agg, de = -agg, pide = pide)
}


#' Vary the minimum-signal constant
#'
#' The threshold \eqn{\theta_n = c\sqrt{s\log(pq)/n}} specifies an order of magnitude, not a constant.
#' `hdmvm_sens_signal_sensitivity()` recomputes the selection robustness values across a grid of \eqn{c}
#' so the report does not rest on an arbitrary choice.
#'
#' @param fit The list returned by [bootstrap_model()].
#' @param c_grid Values of the constant to try.
#' @param ... Passed to [hdmvm_sens()].
#'
#' @returns A data frame with columns `c`, `theta_n`, `RV_erode`, `RV_spurious`.
#' @export
hdmvm_sens_signal_sensitivity <- function(fit, c_grid = c(0.25, 0.5, 1, 2, 4), ...){
  do.call(rbind, lapply(c_grid, function(cc){
    s <- hdmvm_sens(fit, signal_const = cc, ...)
    data.frame(c = cc, theta_n = s$selection$theta_n, RV_erode = s$selection$RV_erode,
               RV_spurious = s$selection$RV_spurious, row.names = NULL)
  }))
}


#' @export
print.hdmvm_sens <- function(x, ...){
  cat("Sensitivity analysis for unmeasured mediator-outcome confounding\n")
  cat(sprintf("  n = %d, p = %d, q = %d, l = %d;  intervals: %s (%.0f%%)\n",
              x$meta$n, x$meta$p, x$meta$q, x$meta$l, x$meta$ci_type, 100 * x$meta$level))
  cat(sprintf("  Q = %.4g", x$Q))
  if(is.finite(x$Q_detail$Q_crossfit))
    cat(sprintf("  (plug-in %.4g, cross-fitted %.4g; R2_A~M = %.3f)",
                x$Q_detail$Q_plugin, x$Q_detail$Q_crossfit, x$R2_AM))
  cat("\n\n")

  cat("Robustness values (R2_M = R2_Y required to overturn):\n")
  tb <- rbind(x$tide, x$de)
  print(data.frame(Outcome = tb$outcome, Estimand = tb$estimand,
                   Est = round(tb$est, 3),
                   CI = paste0("(", round(tb$lower, 2), ", ", round(tb$upper, 2), ")"),
                   RV = round(tb$RV, 4), Reading = tb$reading),
        row.names = FALSE)

  cat(sprintf("\nSelection (theta_n = %.4f, c = %g, s.hat = %d, support: %s):\n",
              x$selection$theta_n, x$meta$signal_const, x$selection$s_hat,
              x$selection$support_source))
  cat(sprintf("  RV_erode    = %.4f  (%s)\n", x$selection$RV_erode,
              paste(x$selection$which_erode, collapse = " / ")))
  cat(sprintf("  RV_spurious = %.4f  (%s)\n", x$selection$RV_spurious,
              paste(x$selection$which_spurious, collapse = " / ")))
  cat(sprintf("  at the %g quantile: erode = %.4f, spurious = %.4f\n",
              x$selection$sel_quantile, x$selection$RV_erode_q, x$selection$RV_spurious_q))
  cat(sprintf("  erosion computed over %d %s of %d selected pairs (unrestricted: %.4f)\n",
              x$selection$n_reported,
              if(isTRUE(x$selection$report_only)) "reported" else "selected",
              x$selection$n_selected, x$selection$RV_erode_all))

  if(!is.null(x$benchmark)){
    g <- max(x$benchmark$geo_mean, na.rm = TRUE)
    cat(sprintf("\nStrongest measured confounder: sqrt(R2_M R2_Y) = %.4f\n", g))
    cat(sprintf("  An RV of r requires a confounder %.1f x r times as strong as that.\n", 1 / g))
  }
  invisible(x)
}


#' @export
plot.hdmvm_sens <- function(x, outcome = 1, estimand = c("TIDE", "DE"), ...){
  estimand <- match.arg(estimand)
  tb <- if(estimand == "TIDE") x$tide else x$de
  k  <- if(is.character(outcome)) match(outcome, tb$outcome) else outcome
  rv <- tb$RV[k]
  if(!is.finite(rv)) stop("No finite robustness value for that outcome.")

  gr <- seq(0.005, 1, length.out = 300)
  contour_y <- pmin(rv^2 / gr, 1)
  plot(gr, contour_y, type = "l", lwd = 2, xlim = c(0, 1), ylim = c(0, 1),
       xlab = expression(R[M]^2), ylab = expression(R[Y]^2),
       main = sprintf("%s, %s: confounding required to overturn",
                      tb$outcome[k], estimand))
  abline(a = 0, b = 1, lty = 3, col = "grey60")
  points(rv, rv, pch = 19)
  text(rv, rv, labels = sprintf("RV = %.3f", rv), pos = 4, cex = 0.8)
  if(!is.null(x$benchmark)){
    b <- x$benchmark[x$benchmark$outcome == tb$outcome[k], , drop = FALSE]
    points(b$R2_M, b$R2_Y, pch = 17, col = "cornflowerblue")
    text(b$R2_M, b$R2_Y, labels = b$confounder, pos = 4, cex = 0.7, col = "cornflowerblue")
  }
  legend("topright", bty = "n",
         legend = c("tipping contour", "RV (symmetric)", "measured confounders"),
         lty = c(1, NA, NA), pch = c(NA, 19, 17),
         col = c("black", "black", "cornflowerblue"))
  invisible(NULL)
}


#' Check joint attainability of outcome-side sensitivity parameters
#'
#' Each \eqn{R^2_Y(k)} in \eqn{[0,1]} *marginally*, but a vector with entry-wise bounds
#' may not be realizable by any single confounder. The requirement that actually determines
#' feasibility is that the outcome covariance \eqn{\Sigma_{e_Y} = \Sigma_Y - \lambda_Y^\top\lambda_Y} be psd;
#' the marginal bounds are its diagonal consequences and are strictly weaker.
#' `hdmvm_sens_feasible()` evaluates the exact condition.
#'
#' @param s An `hdmvm_sens` object.
#' @param R2_Y Numeric vector of length \eqn{q} (or a scalar, recycled) of proposed outcome-side sensitivity parameters.
#' @param signs Optional \eqn{\pm1} vector giving the signs of the outcome loadings. When
#'   `NULL` (the default) the adversarial sign pattern is used, i.e. the one making the
#'   configuration hardest to realize, which is the conservative choice for a worst-case
#'   report. Searching all \eqn{2^q} patterns is done exactly for \eqn{q\le 20}.
#'
#' @returns A list with `statistic` (the quadratic form \eqn{v^\top C_Y^{-1}v}),
#'   `feasible` (`TRUE` when the statistic is at most 1), `max_scaling` (the largest
#'   factor by which `R2_Y` could be multiplied and remain attainable),
#'   `max_common_R2` (the largest value attainable simultaneously at every outcome), and
#'   `signs` (the pattern used).
#'
#' @details
#' Writing \eqn{C_Y} for the residual outcome correlation matrix and
#' \eqn{v_k = \pm\sqrt{R^2_Y(k)}}, feasibility is exactly \eqn{v^\top C_Y^{-1}v\le1}.
#' For uncorrelated outcomes this reads \eqn{\sum_kR^2_Y(k)\le1}: a single confounder has a fixed
#' "budget" to spread across outcomes and cannot be maximally damaging everywhere at once.
#' Correlation relaxes the budget, and as outcomes become perfectly correlated the condition returns
#' to the single-outcome bound, since perfectly correlated outcomes are effectively one outcome.
#'
#' Per-outcome robustness values remain individually valid whatever this function
#' reports as each is a correct statement about confounding of that outcome. What an
#' infeasible verdict means is that the *vector* describes a scenario no single confounder
#' can produce, and so is jointly over-pessimistic if presented as one narrative.
#'
#' @examples
#' \dontrun{
#' s <- hdmvm_sens(fit)
#' ## are the reported per-outcome robustness values jointly attainable?
#' hdmvm_sens_feasible(s, s$tide$RV^2)
#' }
#'
#' @seealso [hdmvm_sens()]
#' @export
hdmvm_sens_feasible <- function(s, R2_Y, signs = NULL){
  stopifnot(inherits(s, "hdmvm_sens"))
  q <- s$meta$q
  if(length(R2_Y) == 1) R2_Y <- rep(R2_Y, q)
  if(length(R2_Y) != q) stop("`R2_Y` must have length 1 or q.")
  if(any(R2_Y < 0 | R2_Y > 1, na.rm = TRUE)) stop("`R2_Y` entries must lie in [0,1].")
  R2_Y[is.na(R2_Y)] <- 0

  SY <- s$meta$Sigma_Y; d  <- sqrt(diag(SY))
  CY <- SY/outer(d,d); Ci <- solve(CY)

  quad <- function(sg) as.numeric(t(sg * sqrt(R2_Y)) %*% Ci %*% (sg * sqrt(R2_Y)))

  if(is.null(signs)){
    ## adversarial: the sign pattern that makes the configuration hardest to attain
    if(q <= 20){
      pats <- as.matrix(expand.grid(rep(list(c(-1, 1)), q)))
      vals <- apply(pats, 1, quad)
      i <- which.max(vals); signs <- as.numeric(pats[i, ]); st <- vals[i]
    } else {
      signs <- rep(1, q); st <- quad(signs)
      warning("q > 20: using the all-positive sign pattern rather than searching.")
    }
  } else {
    signs <- sign(signs); st <- quad(signs)
  }

  ones <- rep(1, q)
  list(statistic = st, feasible = st <= 1 + 1e-9,
       max_scaling = if(st > 0) 1 / st else Inf,
       max_common_R2 = min(1 / as.numeric(t(ones) %*% Ci %*% ones), 1),
       signs = signs)
}


#' Convert between marginal and partial outcome-side sensitivity parameters
#'
#' The outcome-side parameter of this framework, \eqn{R^2_Y(k)}, is the _marginal_
#' association of the confounder with the outcome residual --- the convention of Imai,
#' Keele, and Yamamoto (2010). Cinelli and Hazlett (2020) instead use the _partial_
#' association given the included regressors, which here means conditioning on the
#' mediators. The two differ whenever \eqn{R^2_M>0}, since the confounder is correlated
#' with the mediators. `hdmvm_sens_r2_convert()` provides a map between the two.
#'
#' @param R2 Numeric in \eqn{[0,1]}: the value to convert.
#' @param R2_M Numeric in \eqn{[0,1)}: the mediator-side sensitivity parameter.
#' @param to Either `"partial"` (the default; treats `R2` as marginal) or `"marginal"`
#'   (treats `R2` as partial).
#'
#' @returns A numeric of the same length as `R2`.
#'
#' @details
#' The relationship is
#' \deqn{R^{2,\mathrm{par}}_Y(k) = \frac{R^2_Y(k)\,(1-R^2_M)}{1-R^2_Y(k)\,R^2_M},}
#' which is monotone in \eqn{R^2_Y(k)} and satisfies
#' \eqn{R^{2,\mathrm{par}}_Y(k)\le R^2_Y(k)}, with equality iff
#' \eqn{R^2_M = 0}. Either parameterization may be elicited provided the choice is stated.
#'
#' @examples
#' hdmvm_sens_r2_convert(0.20, R2_M = 0.30)                  # marginal -> partial
#' hdmvm_sens_r2_convert(0.20, R2_M = 0.30, to = "marginal") # and back
#'
#' @seealso [hdmvm_sens()]
#' @export
hdmvm_sens_r2_convert <- function(R2, R2_M, to = c("partial", "marginal")){
  to <- match.arg(to)
  if(any(R2 < 0 | R2 > 1)) stop("`R2` must lie in [0,1].")
  if(R2_M < 0 || R2_M >= 1) stop("`R2_M` must lie in [0,1).")
  if(to == "partial"){
    R2*(1-R2_M)/(1-R2*R2_M)
  } else {
    R2/(1-R2_M+R2*R2_M)
  }
}
