#' Correlation Analysis
#'
#' Conducts correlation analysis for each comparison in a
#' provided SummarizedExperiment. Returns the top compounds for
#' each correlation pair, statistical results, and a correlation
#' heatmap. Also supports inputs of two different matrices for
#' calculating pairwise correlations.
#'
#' @param ptitle Heatmap title.
#' @param exp1 SummarizedExperiment object.
#' @param asy1 Assay to use for the first matrix.
#' @param exp2 (optional) Second SummarizedExperiment object or
#' matrix if conducting correlation analysis between two matrices.
#' @param asy2 Assay to use for the second matrix, if applicable.
#' @param var_cl1 Class or sample group variable for adding an
#' annotation to the y-axis (exp1) of the correlation heatmap.
#' @param var_cl2 Class or sample group variable for adding an
#' annotation to the x-axis (exp2) of the correlation heatmap.
#' @param mode_cor Should correlation coefficients be calculated
#' based on "samples" (columns) or "features" (rows)?
#' @param filt_sig Show only significant correlations on the heatmap?
#' @param snc Show heatmap column names?
#' @param snr Show heatmap row names?
#' @param hh Heatmap height.
#' @param hw Heatmap width.
#' @param fsc Heatmap column font size.
#' @param fsr Heatmap row font size.
#' @param hmlinewidth Heatmap box line width.
#' @return A correlation analysis result list.
#' @examples
#' # ms_stat_corr(
#' #   exp1 = d1,
#' #   asy1 = "raw"
#' # )
#'
#' @export
ms_stat_corr <- function(
  ptitle = "Correlation Matrix",
  exp1,
  asy1 = "norm",
  exp2 = NULL,
  asy2 = "norm",
  var_cl1 = NULL,
  var_cl2 = NULL,
  mode_cor = "features",
  filt_sig = TRUE,
  snc = TRUE,
  snr = TRUE,
  hh = 20,
  hw = 20,
  fsc = 8,
  fsr = 8,
  hmlinewidth = 0.0
) {
  #---- Setup ----
  if (mode_cor == "features") {
    ex1 <- t(assay(exp1, asy1))
    ex1 <- ex1[, colSums(ex1) > 0]
    if (!is.null(exp2)) {
      if ("SummarizedExperiment" %in% class(exp2)) {
        ex2 <- t(assay(exp2, asy2))
      }
      if ("matrix" %in% class(exp2)) {
        ex2 <- exp2
      }
      ex2 <- ex2[, colSums(ex2) > 0]
    }
  }
  if (mode_cor == "samples") {
    ex1 <- assay(exp1, asy1)
    ex1 <- ex1[, colSums(ex1) > 0]
    if (!is.null(exp2)) {
      if ("SummarizedExperiment" %in% class(exp2)) {
        ex2 <- assay(exp2, asy2)
      }
      if ("matrix" %in% class(exp2)) {
        ex2 <- exp2
      }
      ex2 <- ex2[, colSums(ex2) > 0]
    }
  }
  #---- Run Correlation Analysis  ----
  if (!exists("ex2")) {
    cor1 <- cor(ex1, method = "spearman", use = "complete.obs")
  }
  if (exists("ex2")) {
    cor1 <- cor(ex1, ex2, method = "spearman", use = "complete.obs")
  }
  #---- Return top values per compound/sample ----
  cor1top <- dplyr::bind_rows(lapply(
    seq.int(1, nrow(cor1), 1),
    function(i) {
      dcor <- cor1[i, ]
      t1 <- sort(dcor, decreasing = TRUE)[1:10]
      t2 <- sort(dcor, decreasing = FALSE)[1:10]
      out1 <- data.frame(
        "Compound" = rep(rownames(cor1)[i], 10),
        "Pos.Top" = names(t1),
        "Pos.Rho" = t1,
        "Neg.Top" = names(t2),
        "Neg.Rho" = t2
      )
      return(out1) # nolint
    }
  ))
  #---- Calculate FDR corrected p-values ----
  # Convert rho to t-statistic, then to p-value
  rho_to_pval <- function(rho, n) {
    t_stat <- rho * sqrt((n - 2) / (1 - rho^2))
    p_val  <- 2 * pt(-abs(t_stat), df = n - 2)
    return(p_val) # nolint
  }
  cross_pval <- rho_to_pval(cor1, n = nrow(cor1))
  # FDR correction
  cross_padj <- matrix(
    p.adjust(as.vector(cross_pval), method = "BH"),
    nrow     = nrow(cross_pval),
    ncol     = ncol(cross_pval),
    dimnames = dimnames(cross_pval)
  )
  #---- Plot Heatmap ----
  # Filter non-significant correlations
  if (filt_sig == TRUE) {
    total_pairs <- prod(dim(cor1))
    cor1[cross_padj >= 0.05] <- NA
    sig_pairs  <- sum(!is.na(cor1))
    cat("Total pairs:       ", total_pairs, "\n")
    cat("Significant pairs: ", sig_pairs, "\n")
    cat("Retained after filtering (%):      ", round(sig_pairs / total_pairs * 100, 2), "%\n") # nolint
    cor1[is.na(cor1)] <- 0
  }
  # Include heatmap annotation
  if (!is.null(var_cl1)) {
    if (mode_cor == "features") {
      cl1 <- as.character(
        rowData(exp1)[[var_cl1]][rownames(exp1) %in% colnames(ex1)]
      )
    }
    if (mode_cor == "samples") {
      cl1 <- as.character(
        colData(exp1)[[var_cl1]][colnames(exp1) %in% colnames(ex1)]
      )
    }
    fun_hm_bar_row <- list(
      "var_row1" = setNames(
        col_univ()[1:length(unique(cl1))], # nolint
        factor(
          as.character(unique(cl1)),
          levels = gtools::mixedsort(
            as.character(unique(cl1))
          )
        ) # nolint
      )
    )
    ### Annotations
    hm_anno_row <- list(
      "row1" = ComplexHeatmap::rowAnnotation( # nolint
        `var_row1` = factor(
          as.character(cl1),
          levels = unique(gtools::mixedsort(cl1))
        ), # nolint
        col = fun_hm_bar_row,
        show_annotation_name = FALSE
      )
    )
  }
  if (!is.null(var_cl2)) {
    if (mode_cor == "features") {
      cl2 <- as.character(
        rowData(exp2)[[var_cl2]][rownames(exp2) %in% colnames(ex2)]
      )
    }
    if (mode_cor == "samples") {
      cl2 <- as.character(
        colData(exp2)[[var_cl2]][colnames(exp2) %in% colnames(ex2)]
      )
    }
    # Heatmap annotation colors
    fun_hm_bar <- list(
      cl_var = setNames(
        col_univ()[1:length(unique(cl2))], # nolint
        gtools::mixedsort(unique(cl2))
      )
    )
    # Annotations
    hm_anno_col <- list(
      hm_col <- ComplexHeatmap::HeatmapAnnotation( # nolint
        `Group` = cl2, # nolint
        col = fun_hm_bar,
        show_annotation_name = FALSE,
        show_legend = TRUE
      )
    )
  }
  # Heatmap colors
  fun_hm_col <- circlize::colorRamp2(
    c(-1, - 0.5, 0, 0.5, 1),
    colors = col_grad(scm = 7) # nolint
  )
  h3 <- ComplexHeatmap::Heatmap(
    cor1,
    top_annotation = if (!is.null(var_cl2)) hm_anno_col[[1]] else NULL,
    left_annotation = if (!is.null(var_cl1)) hm_anno_row[[1]] else NULL,
    col = fun_hm_col,
    name = "Correlation",
    show_column_names = snc,
    show_row_names = snr,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    heatmap_width = ggplot2::unit(hw, "cm"),
    heatmap_height = ggplot2::unit(hh, "cm"),
    column_title = ptitle,
    row_names_gp = grid::gpar(fontsize = fsr),
    column_names_gp = grid::gpar(fontsize = fsc),
    rect_gp = grid::gpar(col = "black", lwd = hmlinewidth)
  )
  #---- Output all as list ----
  out1 <- list(
    "cor_mat" = cor1,
    "cor_top" = cor1top,
    "cor_pval" = cross_pval,
    "cor_padj" = cross_padj,
    "cor_heatmap" = h3
  )
  return(out1) # nolint
}
