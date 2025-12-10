#' Summary table for HDMVMediation models
#'
#' `hdmvm_table()` converts the matrix output of the main `bootstrap_model()` function into a searchable/filterable html table using the `DT` package.
#'
#' @param mod_boot_summ A numeric matrix; the output from `bootstrap_model()`
#' @param p An integer representing the number of candidate mediators
#' @param outcomes A character vector with the names of each outcome
#'
#' @returns An html widget of class `datatables`.
#'
#' @examples
#'\dontrun{
#'## Load toy data
#'data(hdmvmed_test_data)
#'
#'## Assign mediator/confounder/trt/outcome matrices
#'mediators <- hdmvmed_test_data[,7:106]
#'confounders <- hdmvmed_test_data[,107:111]
#'trt <- hdmvmed_test_data[,1]
#'outcomes <- hdmvmed_test_data[,2:6]
#'
#'## Fit the model
#'model_fit <- bootstrap_model(mediators = mediators, confounders = confounders,
#'                              trt = trt, outcomes = outcomes, nB = 10)
#'
#'## Generate summary table
#'hdmvm_table(model_fit, p = ncol(mediators), outcomes = colnames(outcomes))
#'}
#'
#' @export
hdmvm_table <- function(mod_boot_summ, p, outcomes){
  mod_boot_summ <- as.data.frame(mod_boot_summ)

  bca_inter <- paste0("(",round(mod_boot_summ$bca_lowerCL,4),", ",round(mod_boot_summ$bca_upperCL,4),")")
  per_inter <- paste0("(",round(mod_boot_summ$per_lowerCL,4),", ",round(mod_boot_summ$per_upperCL,4),")")

  mod_boot_outcome <- rep(outcomes, each = p+2)

  mod_boot_table_df <- data.frame("Outcome" = mod_boot_outcome,
                                  "Estimand" = gsub("(_resp[1-9])", "", gsub("(_ide_resp[1-9])", "", rownames(mod_boot_summ))),
                                  "Orig. Est." = round(mod_boot_summ$OrigEst,4),
                                  "Mean(boot)" = round(mod_boot_summ$Mean_boot,4),
                                  "sd(boot)" = round(mod_boot_summ$boot_SE,4),
                                  "pval" = signif(mod_boot_summ$boot_pval,4),
                                  "BC_a CI" = bca_inter,
                                  "PBCI" = per_inter)
  # mod_boot_table_df_print <- mod_boot_table_df[which(mod_boot_table_df$Orig..Est.!=0),]
  mod_boot_table_df_print <- mod_boot_table_df

  DT::datatable(mod_boot_table_df_print, rownames = F,
                colnames = c("Outcome", "Estimand", "Orig. Est.",
                             "Mean(boot)", "sd(boot)", "p-value",
                             "BCa CI", "PBCI"))
  # clipr::write_clip(knitr::kable(mod_boot_table_df_print, digits = 4, format = "latex", align = "c"))
}
