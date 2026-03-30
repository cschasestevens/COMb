#' Dimension Reduction
#'
#' Performs dimension reduction for a selected dataset using one
#' of the following methods: PCA, UMAP, PLSDA.
#'
#' @param exp1 SummarizedExperiment object.
#' @param asy Assay to use for dimension reduction.
#' @param md_var Name of column data variable for plot overlay.
#' @param dim_type Dimension reduction technique to use for plotting.
#' Currently available options are "PCA", "UMAP", or "PLSDA".
#' @param ret_umap Only used if dim_type = "UMAP"; Use Python umap-learn
#' for calculating UMAP (TRUE or FALSE).
#' @param dim1 Number of dimensions to plot data (Either "2D" or "3D").
#' Only 2D plots are used if plotting multiple panels.
#' @param p_lab Label groups on plot panels (Only available if "2D").
#' @param show_axes Should the x and y axis information be shown on the plot?
#' @param namefile If dim1 is "3D", provide the file name for saving the plot.
#' @param flag_outliers If TRUE, flag samples that are +/- 3 standard
#' deviations from the individual group means. Generally used only when
#' evaluating results from ms_data_norm.
#' @param outlier_calc If flag_outliers is TRUE, provide the column name
#' containing the TIC for determining outliers.
#' @param sid If flag_outliers is TRUE, provide the name of the sample ID
#' column for labeling outliers.
#' @param scale_data If TRUE, data will be pareto scaled
#' prior to dimension reduction.
#' @return A dimension reduction plot of the specified type.
#' @examples
#'
#' # ms_dim_rd(
#' #   exp1 = data1,
#' #   md_var = "Group"
#' # )
#'
#' @export
ms_dim_rd <- function( # nolint
  exp1,
  asy = "norm.pareto",
  md_var,
  dim_type = "PCA",
  ret_umap = FALSE,
  dim1 = "2D",
  p_lab = TRUE,
  show_axes = TRUE,
  namefile = NULL,
  flag_outliers = FALSE,
  outlier_calc = NULL,
  sid = NULL,
  scale_data = TRUE
) {
  # Load data
  d1 <- exp1
  if (scale_data == FALSE) {
    SummarizedExperiment::assay(d1, "scaled") <- assay(d1, asy) # nolint
  }
  if (scale_data == TRUE) {
    print("Scaling data for dimension reduction...")
    # log2-transformation
    assay(d1, "log2") <- log2(assay(d1, asy)) # nolint
    # pareto scaling
    assay(d1, "scaled") <- apply( # nolint
      apply(
        assay(d1, "log2"),
        2,
        function(x) x - mean(x)
      ),
      2,
      function(y) round(y / sqrt(sd(y)), digits = 3)
    )
  }
  #---- PCA ----
  if (dim_type == "PCA") {
    # Calculate PCA
    p1 <- prcomp(t(assay(d1, "scaled"))) # nolint
    vrc <- summary(p1)$importance[2, ]
    ## Format input
    p2 <- setNames(
      as.data.frame(p1[["x"]][, 1:3]),
      c(
        unlist(
          lapply(
            seq.int(1, 3, 1),
            function(x) {
              paste0(
                paste("PC", x, sep = ""), "(",
                formatC(100 * vrc[[x]], format = "f", digits = 2), "%)"
              )
            }
          )
        )
      )
    )
  }
  #---- UMAP ----
  if (dim_type == "UMAP") {
    if(ret_umap == TRUE) { # nolint
      p1 <- umap::umap(
        t(assay(d1, "scaled")), # nolint
        method = "umap-learn",
        n_epochs = 500,
        n_components = 3,
        verbose = TRUE,
        preserve.seed = FALSE
      )
    }
    if(ret_umap == FALSE) { # nolint
      p1 <- umap::umap(
        t(assay(d1, "scaled")), # nolint
        n_epochs = 500,
        n_components = 3,
        verbose = TRUE,
        preserve.seed = FALSE
      )
    }
    p2 <- setNames(
      as.data.frame(p1[["layout"]]),
      c("UMAP.1", "UMAP.2", "UMAP.3")
    )
  }
  #---- PLSDA ----
  if (dim_type == "PLSDA") {

  }
  #---- Outlier flag ----
  if (flag_outliers == TRUE) {
    # Determine samples +/- 3 SD from mean mTIC of each group
    ## Calculate normalized group TIC mean and sd
    colData(d1) <- cbind( # nolint
      colData(d1),
      data.frame(
        "TIC.group.mean" = unlist(
          lapply(
            seq.int(1, ncol(d1), 1),
            function(i) {
              mean(
                colData(d1)[ # nolint
                  grepl(
                    colData(d1)[i, ][[md_var]], # nolint
                    colData(d1)[[md_var]]
                  ),
                ][[outlier_calc]]
              )
            }
          )
        ),
        "TIC.group.sd" = unlist(
          lapply(
            seq.int(1, ncol(d1), 1),
            function(i) {
              sd(
                colData(d1)[ # nolint
                  grepl(
                    colData(d1)[i, ][[md_var]], # nolint
                    colData(d1)[[md_var]]
                  ),
                ][[outlier_calc]]
              )
            }
          )
        )
      )
    )
    ## Calculate number of standard deviations from mean
    colData(d1)[["sd.number"]] <- abs( # nolint
      colData(d1)[[outlier_calc]] - colData(d1)[["TIC.group.mean"]]
    ) /
      colData(d1)[["TIC.group.sd"]]
    ## flag samples +/- 3 SD from mean
    colData(d1)[["Outlier.flag"]] <- unlist( # nolint
      lapply(
        seq.int(1, ncol(d1), 1),
        function(i) {
          ifelse(
            colData(d1)[i, ][[outlier_calc]] > # nolint
              colData(d1)[i, ][["TIC.group.mean"]] +
                3 * colData(d1)[i, ][["TIC.group.sd"]] ||
              colData(d1)[i, ][[outlier_calc]] <
                colData(d1)[i, ][["TIC.group.mean"]] -
                  3 * colData(d1)[i, ][["TIC.group.sd"]],
            1, 0
          )
        }
      )
    )
  }
  #---- Plot ----
  if(dim1 == "2D") { # nolint
    if(p_lab == TRUE) { # nolint
      # Calculate median label positions
      p3lab <- setNames(
        aggregate(
          p2[, 1:2],
          list(colData(d1)[[md_var]]), # nolint
          FUN = median
        ),
        c(md_var, names(p2)[1:2]
        )
      )
      # Plot
      p3 <- ggplot2::ggplot(
        data = p2,
        ggplot2::aes(
          x = .data[[names(p2)[[1]]]], # nolint
          y = .data[[names(p2)[[2]]]] # nolint
        )
      ) +
        ggplot2::geom_point(
          ggplot2::aes(
            color = colData(d1)[[md_var]] # nolint
          ),
          shape = 16,
          size = 2,
          alpha = 0.6
        ) +
        ggrepel::geom_text_repel(
          data = p3lab,
          ggplot2::aes(
            label = p3lab[[md_var]]
          ),
          size = 4,
          bg.color = "white"
        ) +
        ggplot2::scale_color_manual(
          paste(""),
          values = col_univ() # nolint
        ) +
        ms_theme() # nolint
      # If show_axes is FALSE
      if (show_axes == FALSE) {
        p3 <- p3 +
          ggplot2::theme(
            panel.grid.major.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            plot.margin = ggplot2::unit(
              c(0.1, 0.1, 0.1, 0.1), "cm"
            ),
            legend.position = "none"
          )
      }
    }
    if(p_lab == FALSE) { # nolint
      # Plot
      p3 <- ggplot2::ggplot(
        data = p2,
        ggplot2::aes(
          x = .data[[names(p2)[[1]]]], # nolint
          y = .data[[names(p2)[[2]]]] # nolint
        )
      ) +
        ggplot2::geom_point(
          ggplot2::aes(
            color = colData(d1)[[md_var]] # nolint
          ),
          shape = 16,
          size = 2,
          alpha = 0.6
        ) +
        ggplot2::scale_color_manual(
          paste(""),
          values = col_univ() # nolint
        ) +
        ms_theme() # nolint
      # If show_axes is FALSE
      if (show_axes == FALSE) {
        p3 <- p3 +
          ggplot2::theme(
            panel.grid.major.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            plot.margin = ggplot2::unit(
              c(0.1, 0.1, 0.1, 0.1), "cm"
            ),
            legend.position = "none"
          )
      }
    }
    p3 <- p3 + ggplot2::labs(title = dim_type)
    if (flag_outliers == TRUE) {
      # Generate loadings plot
      p3 <- p3 +
        ggrepel::geom_label_repel(
          ggplot2::aes(
            label = ifelse(
              colData(d1)[["Outlier.flag"]] == 1, # nolint
              colData(d1)[[sid]],
              ""
            ),
            alpha = ifelse(
              colData(d1)[["Outlier.flag"]] == 1,
              1, 0
            ),
            linewidth = ifelse(
              colData(d1)[["Outlier.flag"]] == 1,
              0.15, NA
            )
          ),
          show.legend = FALSE
        )
      # Generate standard deviation plot
      d2 <- as.data.frame(colData(d1)) # nolint
      d2 <- d2[order(d2[[md_var]]), ] # nolint
      d2[["ID"]] <- seq.int(1, nrow(d2), 1)
      p4 <- ggplot2::ggplot(
        data = d2,
        ggplot2::aes(x = .data[["ID"]]) # nolint
      ) +
        # mean line
        ggplot2::geom_point(
          ggplot2::aes(
            y = .data[["sd.number"]],
            color = .data[[md_var]] # nolint
          )
        ) +
        ggplot2::geom_hline(
          yintercept = 3,
          linetype = "dashed"
        ) +
        # outlier labels
        ggrepel::geom_label_repel(
          ggplot2::aes(
            y = .data[["sd.number"]],
            label = ifelse(
              .data[["sd.number"]] > 3, # nolint
              d2[[sid]],
              ""
            ),
            alpha = ifelse(
              .data[["sd.number"]] > 3,
              1, 0
            ),
            linewidth = ifelse(
              .data[["sd.number"]] > 3,
              0.15, NA
            )
          ),
          show.legend = FALSE
        ) +
        ggplot2::labs(
          title = "SD Plot",
          x = "Sample ID",
          y = "Standard Deviations From Group Mean"
        ) +
        ggplot2::scale_color_manual(
          paste(""),
          values = col_univ() # nolint
        ) +
        ms_theme() # nolint
      p3 <- ggpubr::ggarrange(
        p3, p4, nrow = 1, ncol = 2
      )
    }
  }
  if(dim1 == "3D") { # nolint
    p3 <- plotly::plot_ly(
      p2,
      x = ~.data[[names(p2)[[1]]]], # nolint
      y = ~.data[[names(p2)[[2]]]], # nolint
      z = ~.data[[names(p2)[[3]]]], # nolint
      color = ~as.factor(colData(d1)[[md_var]]), # nolint
      colors = col_univ() # nolint
    ) %>% # nolint
      plotly::add_markers(marker = list(size = 4)) %>%
      plotly::layout(
        autosize = FALSE,
        width = 800,
        height = 600,
        margin = list(
          l = 50,
          r = 50,
          b = 25,
          t = 25,
          pad = 1
        )
      )
    htmlwidgets::saveWidget(
      p3,
      file = namefile # nolint
    )
  }
  return(p3)
}
