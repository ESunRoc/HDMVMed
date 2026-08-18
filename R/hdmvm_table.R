#' Summary table for HDMVMediation models
#'
#' `hdmvm_table()` converts the matrix output of the main [bootstrap_model()] function into a searchable/filterable html table using the `DT` package.
#'
#' @param mod_boot A list; the output from [bootstrap_model()]
#' @param p An integer representing the number of candidate mediators
#' @param outcomes A character vector with the names of each outcome
#' @param p.adj.method A character string for the p-value correction method; see \link[stats]{p.adjust}. Defaults to `"BH"`.
#' @param DT_table A boolean indicating whether to return a DT table
#'
#' @returns Either an html widget of class `datatables` (if `DT_table = TRUE`) or the equivalent dataframe (if `DT_table = FALSE`). The latter is used when constructing a DAG with the function `to_DAG()`. Includes a `Moderator` column identifying, for each row, which moderator (if any) the reported effect is indexed by; unmoderated PIDE/TIDE/DE rows show `"-"`.
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
hdmvm_table <- function(mod_boot, p, outcomes, p.adj.method = "BH", DT_table = TRUE){
  p.adj.method <- match.arg(p.adj.method)

  mod_boot$res <- as.data.frame(mod_boot$res)

  bca_inter <- paste0("(",round(mod_boot$res$bca_lowerCL,4),", ",round(mod_boot$res$bca_upperCL,4),")")
  per_inter <- paste0("(",round(mod_boot$res$per_lowerCL,4),", ",round(mod_boot$res$per_upperCL,4),")")

  row_ids <- rownames(mod_boot$res)

  # Response/outcome index: every row (main or moderated) ends "..._resp<k>".
  resp_idx <- as.integer(sub(".*_resp(\\d+)$", "\\1", row_ids))
  mod_boot$res_outcome <- outcomes[resp_idx]

  # Moderator tag: bootstrap_model() names moderated-effect rows
  # "..._MOD_<moderator>_resp<k>"; unmoderated rows have no such tag.
  stripped_resp <- sub("_resp\\d+$", "", row_ids)              # drop "_resp<k>"
  has_mod       <- grepl("_MOD_", stripped_resp)
  moderator     <- ifelse(has_mod, sub(".*_MOD_", "", stripped_resp), "-")

  # Estimand: mediator name for PIDE rows, or "TIDE"/"DE"; drop the moderator tag
  # (if any) and then the "_ide" suffix used only on PIDE rows.
  base     <- sub("_MOD_.*$", "", stripped_resp)                # e.g. "M1_ide", "TIDE", "DE"
  estimand <- sub("_ide$", "", base)                            # "M1_ide" -> "M1"; "TIDE"/"DE" unaffected


  pvals_adjusted <- p.adjust(mod_boot$res$boot_pval, method = p.adj.method)
  mod_boot$res_table_df <- data.frame("Outcome" = mod_boot$res_outcome,
                                  "Estimand" = estimand,
                                  "Moderator" = moderator,
                                  "Orig. Est." = round(mod_boot$res$OrigEst,4),
                                  "Mean(boot)" = round(mod_boot$res$Mean_boot,4),
                                  "sd(boot)" = round(mod_boot$res$boot_SE,4),
                                  "pval" = signif(mod_boot$res$boot_pval,4),
                                  "pval_adj" = signif(pvals_adjusted, 4),
                                  "BC_a CI" = bca_inter,
                                  "PBCI" = per_inter,
                                  check.names = FALSE)
  # mod_boot$res_table_df_print <- mod_boot$res_table_df[which(mod_boot$res_table_df$Orig..Est.!=0),]
  mod_boot$res_table_df_print <- mod_boot$res_table_df

  if(DT_table){
    DT::datatable(mod_boot$res_table_df_print, rownames = F,
                  colnames = c("Outcome", "Mediator/Effect", "Moderator", "Orig. Est.",
                               "Mean(boot)", "sd(boot)", "p-value",
                               "BCa CI", "PBCI"))
  } else{
    return(mod_boot$res_table_df_print)
  }
}
