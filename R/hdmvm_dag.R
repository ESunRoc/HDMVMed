#' Generate a tikz-cd DAG
#'
#' `hdmvm_dag` generates a DAG using the syntax of the `tikz-cd` LaTeX package.
#'
#' @param mod_boot_summ A data frame; the output from [hdmvm_table()] with `DT_table = FALSE`.
#' @param p An integer denoting the number of candidate mediators i.e., `ncol(mediators)`.
#' @param q An integer denoting the number of outcomes i.e., `ncol(outcomes)`.
#' @param trt_name A character vector denoting the name of the treatment/exposure. By default, "trt" is used as a placeholder.
#' @param alpha A numeric between 0 and 1 denoting the nominal significance threshold for inclusion in the DAG. By default, 0.1 is used.
#'
#' @return A text string that, when passed through `cat()`, returns the `tikz-cd` code for the selected DAG.
#'
#' @examples
#'\dontrun{
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
#'
#' ## Generate summary table
#' model_df <- hdmvm_table(model_fit, p = ncol(mediators),
#'                         outcomes = colnames(outcomes), DT_table = FALSE)
#'
#' ## Generate tikz-cd code
#' cat(hdmvm_dag(model_df, p = ncol(mediators), q = ncol(outcomes)))
#'}
#'
#' @export
hdmvm_dag <- function(mod_boot_summ, p, q, trt_name = "trt", alpha = 0.1){
  dag_matrix <- as.matrix(mod_boot_summ[which(mod_boot_summ$pval[1:(q*p)]<=alpha),c("Outcome","Estimand")])

  # extract node names
  out_nodes <- sort(unique(dag_matrix[, 1]))
  med_nodes <- sort(unique(dag_matrix[, 2]))

  # create the node names (for aliasing in tikz-cd)
  med_map <- setNames(med_nodes, med_nodes)
  out_map <- setNames(out_nodes, out_nodes)

  #### set-up diagram ####
  tikz_cd_code <- c(
    "\\begin{tikzcd}",
    "    % Node Declarations"
  )

  #### Generate outcome nodes ####
  out_declarations <- character(length(out_nodes))
  for(i in seq_along(out_nodes)){
    num_empty_cells <- 2 # to align to the right

    out_label <- out_map[as.character(out_nodes[i])]
    out_declarations[i] <- paste0(
      paste(rep("  ", num_empty_cells), collapse = "&&"),
      "&&",
      paste0(" |[alias=", out_label, "]|{", out_label, "} \\\\")
    )
  }

  #### Generate mediator and treatment nodes ####
  max_rows <- max(length(med_nodes), length(out_nodes))
  row_ids <- 1:max_rows

  for(i in 1:max_rows){
    row_content <- character()

    # treatment node --- only in appx. middle row
    if(i == round(median(row_ids))){
      row_content <- c(
        row_content,
        paste0(" |[alias=", trt_name, "]|{", trt_name, "} ")
      )
    } else {
      row_content <- c(row_content, " ")
    }

    # mediator nodes
    if(i <= length(med_nodes)){
      med_label <- med_map[as.character(med_nodes[i])]
      row_content <- c(row_content,
                       paste0(" && |[alias=", med_label, "]|{", med_label, "} "))
    } else {
      row_content <- c(row_content, " && ") # Empty cell for M
    }

    # outcome nodes
    if(i <= length(out_nodes)){
      out_label <- out_map[as.character(out_nodes[i])]
      row_content <- c(
        row_content,
        paste0(" && |[alias=", out_label, "]|{", out_label, "} ")
      )
    } else {
      row_content <- c(row_content, " &&& ") # Empty cell for Y
    }

    # append line
    tikz_cd_code <- c(tikz_cd_code,
                      paste(" ", paste(row_content, collapse = ""), "\\\\"))
  }

  #### Generate arrows ####
  tikz_cd_code <- c(tikz_cd_code, "    % Arrow Declarations")

  # treatment --> mediators
  for(med_index in med_nodes){
    med_label <- med_map[as.character(med_index)]
    tikz_cd_code <- c(tikz_cd_code,
                      paste0(" \\arrow[from=", trt_name, ", to=", med_label, "]"))
  }

  # mediators --> outcomes
  for (i in 1:nrow(dag_matrix)) {
    out_index <- dag_matrix[i, 1]
    med_index <- dag_matrix[i, 2]

    out_label <- out_map[as.character(out_index)]
    med_label <- med_map[as.character(med_index)]

    arrow_code <- paste0(" \\arrow[from=", med_label, ", to=", out_label, "]")

    # check for duplicates
    if (!arrow_code %in% tikz_cd_code) {
      tikz_cd_code <- c(tikz_cd_code, arrow_code)
    }
  }

  #### generate final figure ####
  tikz_cd_code <- c(tikz_cd_code, "\\end{tikzcd}")

  return(paste(tikz_cd_code, collapse = "\n"))
}
