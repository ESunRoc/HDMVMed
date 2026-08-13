#' Render a DAG directly from the `hdmvm_table()` summary
#'
#' `hdmvm_plot_dag()` automates the construction of a rendered DAG from the
#' searchable summary table produced by [hdmvm_table()], rather than from the
#' raw [bootstrap_model()] output consumed by [hdmvm_tex_dag()]. It selects all
#' mediators significant at the chosen `alpha`, all outcomes with at least
#' one such significant mediator, adds significant direct-effect arrows, and
#' overlays significant moderation effects, all read directly from the
#' `Outcome`/`Estimand`/`Moderator` columns of `mod_boot_summ` rather than by
#' assuming a fixed row order. The diagram is drawn with [igraph::plot.igraph()]
#' as a side effect (in the currently active graphics device).
#'
#' @param mod_boot_summ A data frame; the output from [hdmvm_table()] with `DT_table = FALSE`.
#' @param p An integer denoting the number of candidate mediators i.e., `ncol(mediators)`.
#' @param q An integer denoting the number of outcomes i.e., `ncol(outcomes)`.
#' @param trt_name A character vector denoting the name of the treatment/exposure. By default, "trt" is used as a placeholder.
#' @param alpha A numeric between 0 and 1 denoting the nominal significance threshold for inclusion in the DAG. By default, 0.1 is used.
#' @param use_adj A Boolean indicating whether to use the adjusted p-values; see [hdmvm_table()].
#' @param include_estimates A Boolean indicating whether the effect estimates (the `` `Orig. Est.` `` column of [hdmvm_table()]) should be printed as edge labels. By default, `FALSE`.
#' @param color_arrows A Boolean indicating whether arrows should be colored according to the sign of their associated estimate (positive vs. negative). By default, `FALSE`.
#' @param est_sign_colors A length-2 character vector of colors used for positive and negative estimates (in that order) when `color_arrows = TRUE`. Defaults to `c("blue", "red")`.
#' @param ... Additional arguments passed through to [igraph::plot.igraph()] (e.g. `vertex.label.cex`, `edge.arrow.size`, `main`) for further customizing the rendered figure.
#'
#' @details
#' Node/edge selection follows these rules:
#' \itemize{
#' \item A mediator is included iff it has at least one significant
#' *unmoderated* PIDE (`Moderator == "-"`), and an outcome is included iff it
#' has at least one such significant mediator (mirroring [hdmvm_tex_dag()]'s
#' convention and the package's definition of a "relevant" outcome).
#' \item Mediator -> outcome edges (solid) are drawn for every significant
#' unmoderated PIDE pair, and additionally for a pair that is significant
#' only through a *moderated* PIDE, provided both endpoints are already in
#' the diagram.
#' \item Direct-effect (treatment -> outcome) edges (dashed) are drawn for
#' any included outcome with a significant unmoderated or moderated DE.
#' \item Because [hdmvm_table()] only ever reports the *product* estimand
#' (PIDE = the stage-1 x stage-2 path product; see [bootstrap_model()]),
#' never the individual treatment -> mediator or mediator -> outcome paths
#' in isolation, only mediator -> outcome, treatment -> outcome, and
#' moderator-overlay edges can be honestly labeled/colored by
#' `include_estimates`/`color_arrows`; treatment -> mediator edges are
#' always drawn plain gray, since a single mediator's structural edge to
#' treatment may be shared by several outcomes with different PIDEs.
#' \item Since a moderator conceptually modifies an *edge* rather than a
#' node, and a graph edge in `igraph` must connect two actual vertices,
#' moderation is drawn as a dotted edge from the moderator to the node
#' whose incoming pathway it modifies: `moderator -> mediator` for a
#' significant moderated PIDE (labeled with the specific outcome it
#' applies to, since the same moderator/mediator pair can be significant
#' for one outcome and not another), or `moderator -> outcome` for a
#' significant moderated DE. `igraph` supports multiple parallel edges, so
#' a moderator significant for several outcomes via the same mediator
#' draws one labeled edge per outcome rather than collapsing them.
#' \item Node classes are visually encoded (treatment = gray square;
#' mediators = light-blue circles; outcomes = light-green circles;
#' moderators = gold squares), and nodes are placed in a deterministic
#' 4-column layout (treatment, mediators, outcomes, moderators, each column
#' vertically centered). A legend explaining the three edge line
#' styles is drawn automatically.
#' }
#'
#' @return Invisibly returns the underlying `igraph` graph object (with
#' vertex attribute `type` and edge attributes `lty`, `color`, `label`), so
#' it can be re-plotted, further customized, or exported (e.g. via
#' [igraph::write_graph()]), in addition to being drawn as a side effect.
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
#' ## Fit the model (optionally with a moderator)
#' model_fit <- bootstrap_model(mediators = mediators, confounders = confounders,
#'                              trt = trt, outcomes = outcomes, nB = 10)
#'
#' ## Generate summary table
#' model_df <- hdmvm_table(model_fit, p = ncol(mediators),
#'                         outcomes = colnames(outcomes), DT_table = FALSE)
#'
#' ## Render the DAG directly from the table, with labeled, colored arrows
#' to_DAG(model_df, p = ncol(mediators), q = ncol(outcomes),
#'        include_estimates = TRUE, color_arrows = TRUE)
#'}
#'
#' @export
hdmvm_plot_dag <- function(mod_boot_summ, p, q, trt_name = "trt", alpha = 0.1, use_adj = TRUE,
                           include_estimates = FALSE, color_arrows = FALSE,
                           est_sign_colors = c("blue", "red"), ...){

  if(!requireNamespace("igraph", quietly = TRUE)) stop("The `igraph` package is required by to_DAG(); please install it with install.packages('igraph').")
  if(!is.data.frame(mod_boot_summ)) stop("`mod_boot_summ` must be a data frame; the output from hdmvm_table() with DT_table = FALSE.")
  if(color_arrows && length(est_sign_colors) != 2) stop("`est_sign_colors` must be a length-2 character vector: c(<positive color>, <negative color>).")
  if(!all(c("Outcome","Estimand","Moderator","pval","pval_adj") %in% colnames(mod_boot_summ))) stop("`mod_boot_summ` does not look like hdmvm_table() output (DT_table = FALSE): missing expected columns.")

  est_col <- "Orig. Est."
  if(!est_col %in% colnames(mod_boot_summ)) stop("`mod_boot_summ` is missing the '", est_col, "' column.")

  pcol <- if(use_adj) "pval_adj" else "pval"

  fmt_est <- function(x) formatC(round(as.numeric(x), 3), format = "f", digits = 3)
  arrow_color <- function(x) ifelse(as.numeric(x) >= 0, est_sign_colors[1], est_sign_colors[2])

  is_sig         <- !is.na(mod_boot_summ[[pcol]]) & mod_boot_summ[[pcol]] <= alpha
  is_unmoderated <- mod_boot_summ$Moderator == "-"
  is_pide        <- !(mod_boot_summ$Estimand %in% c("TIDE", "DE"))
  is_DE          <- mod_boot_summ$Estimand == "DE"

  # Node selection: mediators / relevant outcomes (unmoderated PIDE significance only)
  sig_pide_main <- mod_boot_summ[is_sig & is_unmoderated & is_pide, c("Outcome", "Estimand", est_col), drop = FALSE]
  med_nodes <- sort(unique(sig_pide_main$Estimand))
  out_nodes <- sort(unique(sig_pide_main$Outcome))

  if(length(med_nodes) == 0 || length(out_nodes) == 0){
    warning("No mediators are significant at the chosen alpha; returning an (empty) treatment-only diagram.")
  }

  # Mediator -> outcome edges: unmoderated-significant pairs, plus moderated-only-significant pairs
  # for (mediator, outcome) combinations already in the diagram
  sig_pide_mod <- mod_boot_summ[is_sig & !is_unmoderated & is_pide &
                                  mod_boot_summ$Estimand %in% med_nodes &
                                  mod_boot_summ$Outcome %in% out_nodes,
                                c("Outcome", "Estimand", "Moderator", est_col), drop = FALSE]

  my_edges <- unique(rbind(
    if(nrow(sig_pide_main) > 0) sig_pide_main[, c("Outcome","Estimand")] else NULL,
    if(nrow(sig_pide_mod)  > 0) sig_pide_mod[, c("Outcome","Estimand")]  else NULL
  ))
  if(is.null(my_edges)) my_edges <- data.frame(Outcome = character(0), Estimand = character(0))

  # attach a labeling/coloring estimate: prefer the unmoderated PIDE; fall back to the
  # first significant moderated PIDE for that (mediator, outcome) pair when unmoderated is absent
  my_edges[[est_col]] <- sig_pide_main[[est_col]][match(paste(my_edges$Outcome, my_edges$Estimand),
                                                        paste(sig_pide_main$Outcome, sig_pide_main$Estimand))]
  if(nrow(sig_pide_mod) > 0){
    need_fallback <- is.na(my_edges[[est_col]])
    if(any(need_fallback)){
      fb <- sig_pide_mod[match(paste(my_edges$Outcome[need_fallback], my_edges$Estimand[need_fallback]),
                               paste(sig_pide_mod$Outcome, sig_pide_mod$Estimand)), est_col]
      my_edges[[est_col]][need_fallback] <- fb
    }
  }

  # treatment -> outcome direct-effect edges: significant (unmoderated or moderated) DE, restricted to already-relevant outcomes
  sig_DE_main <- mod_boot_summ[is_sig & is_unmoderated & is_DE & mod_boot_summ$Outcome %in% out_nodes,
                               c("Outcome", est_col), drop = FALSE]
  sig_DE_mod  <- mod_boot_summ[is_sig & !is_unmoderated & is_DE & mod_boot_summ$Outcome %in% out_nodes,
                               c("Outcome", "Moderator", est_col), drop = FALSE]
  DE_outcomes <- sort(unique(c(sig_DE_main$Outcome, sig_DE_mod$Outcome)))

  DE_edges <- data.frame(Outcome = DE_outcomes, stringsAsFactors = FALSE)
  DE_edges[[est_col]] <- sig_DE_main[[est_col]][match(DE_edges$Outcome, sig_DE_main$Outcome)]
  if(nrow(sig_DE_mod) > 0){
    need_fallback <- is.na(DE_edges[[est_col]])
    if(any(need_fallback)){
      fb <- sig_DE_mod[match(DE_edges$Outcome[need_fallback], sig_DE_mod$Outcome), est_col]
      DE_edges[[est_col]][need_fallback] <- fb
    }
  }

  # Moderator overlay: real Z -> M / Z -> Y edges
  mod_pide_overlay <- sig_pide_mod   # columns: Outcome, Estimand (mediator), Moderator, est_col
  mod_DE_overlay   <- sig_DE_mod     # columns: Outcome, Moderator, est_col
  moderators_used  <- sort(unique(c(mod_pide_overlay$Moderator, mod_DE_overlay$Moderator)))

  # Build the igraph vertex/edge lists
  vertices <- data.frame(
    name  = c(trt_name, med_nodes, out_nodes, moderators_used),
    type  = c("trt", rep("mediator", length(med_nodes)), rep("outcome", length(out_nodes)),
              rep("moderator", length(moderators_used))),
    stringsAsFactors = FALSE
  )

  edges <- data.frame(from = character(0), to = character(0), lty = character(0),
                      color = character(0), label = character(0), stringsAsFactors = FALSE)

  # trt -> M : structural, plain
  if(length(med_nodes) > 0){
    edges <- rbind(edges, data.frame(from = trt_name, to = med_nodes, lty = "solid",
                                     color = "gray40", label = "", stringsAsFactors = FALSE))
  }

  # M -> Y : mediated pathway
  if(nrow(my_edges) > 0){
    lbl <- if(include_estimates) fmt_est(my_edges[[est_col]]) else ""
    col <- if(color_arrows) ifelse(is.na(my_edges[[est_col]]), "gray40", arrow_color(my_edges[[est_col]])) else "gray20"
    edges <- rbind(edges, data.frame(from = my_edges$Estimand, to = my_edges$Outcome, lty = "solid",
                                     color = col, label = ifelse(is.na(lbl), "", lbl), stringsAsFactors = FALSE))
  }

  # trt -> Y : direct effect
  if(nrow(DE_edges) > 0){
    lbl <- if(include_estimates) fmt_est(DE_edges[[est_col]]) else ""
    col <- if(color_arrows) ifelse(is.na(DE_edges[[est_col]]), "gray40", arrow_color(DE_edges[[est_col]])) else "gray20"
    edges <- rbind(edges, data.frame(from = trt_name, to = DE_edges$Outcome, lty = "dashed",
                                     color = col, label = ifelse(is.na(lbl), "", lbl), stringsAsFactors = FALSE))
  }

  # Z -> M : moderated PIDE overlay
  if(nrow(mod_pide_overlay) > 0){
    lbl <- if(include_estimates){
      paste0(mod_pide_overlay$Outcome, " (", fmt_est(mod_pide_overlay[[est_col]]), ")")
    } else {
      as.character(mod_pide_overlay$Outcome)
    }
    col <- if(color_arrows) arrow_color(mod_pide_overlay[[est_col]]) else "gray50"
    edges <- rbind(edges, data.frame(from = mod_pide_overlay$Moderator, to = mod_pide_overlay$Estimand,
                                     lty = "dotted", color = col, label = lbl, stringsAsFactors = FALSE))
  }

  # Z -> Y : moderated DE overlay
  if(nrow(mod_DE_overlay) > 0){
    lbl <- if(include_estimates) fmt_est(mod_DE_overlay[[est_col]]) else ""
    col <- if(color_arrows) arrow_color(mod_DE_overlay[[est_col]]) else "gray50"
    edges <- rbind(edges, data.frame(from = mod_DE_overlay$Moderator, to = mod_DE_overlay$Outcome,
                                     lty = "dotted", color = col, label = lbl, stringsAsFactors = FALSE))
  }

  g <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = vertices)

  # Layout: deterministic 4-column layout (trt | mediators | outcomes | moderators)
  col_x <- c(trt = 0, mediator = 1, outcome = 2, moderator = 3)
  center_y <- function(k) if(k <= 1) 0 else seq((k - 1) / 2, -(k - 1) / 2, length.out = k)

  y_lookup <- c(stats::setNames(0, trt_name),
                stats::setNames(center_y(length(med_nodes)), med_nodes),
                stats::setNames(center_y(length(out_nodes)), out_nodes),
                stats::setNames(center_y(length(moderators_used)), moderators_used))

  layout_mat <- cbind(x = col_x[igraph::V(g)$type], y = y_lookup[igraph::V(g)$name])

  # Visual encoding
  vshape <- c(trt = "square", mediator = "circle", outcome = "circle", moderator = "square")
  vcolor <- c(trt = "gray70", mediator = "lightblue", outcome = "lightgreen", moderator = "gold")

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  igraph::plot.igraph(g, layout = layout_mat,
                      vertex.shape = vshape[igraph::V(g)$type],
                      vertex.color = vcolor[igraph::V(g)$type],
                      vertex.label = igraph::V(g)$name,
                      vertex.label.color = "black",
                      vertex.size = 28,
                      edge.lty = igraph::E(g)$lty,
                      edge.color = igraph::E(g)$color,
                      edge.label = igraph::E(g)$label,
                      edge.label.cex = 0.75,
                      edge.arrow.size = 0.6,
                      edge.curved   = 0.15,
                      ...)
  graphics::legend("topleft", bty = "n", cex = 0.7,
                   legend = c("mediated (M -> Y)", "direct effect (trt -> Y)", "moderation"),
                   lty = c("solid", "dashed", "dotted"), col = "gray30", seg.len = 2)

  return(invisible(g))
}
