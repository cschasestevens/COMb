#' Dimension Reduction
#'
#' Performs dimension reduction for a selected dataset using one
#' of the following methods: PCA, UMAP, PLSDA.
#'
#' @param mat1 Pareto-scaled data matrix returned by ms_data_check().
#' @param md Data frame containing sample information.
#' @param md_var Vector of either a single variable or multiple variables
#' to plot. For single dimension reduction plots, only a single variable
#' should be listed.
#' @param dim_type Dimension reduction technique to use for plotting.
#' Currently available options are "PCA", "UMAP", or "PLSDA".
#' @param ret_umap Only used if dim_type = "UMAP"; Use Python umap-learn
#' for calculating UMAP (TRUE or FALSE).
#' @param dim1 Number of dimensions to plot data (Either "2D" or "3D").
#' Only 2D plots are used if plotting multiple panels.
#' @param p_type Plot either a single variable or a list of metadata
#' variables. Enter either "single" or "list" to control output.
#' @param p_lab Label groups on plot panels (Only available if "2D").
#' @return A single plot or series of dimension reduction plots.
#' @examples
#'
#' # ms_dim_rd(
#' #   mat1 = d[["data"]][["data.scale"]],
#' #   md = d[["data"]][["meta"]],
#' #   md_var = c("Organ"),
#' #   dim_type = "PCA",
#' #   ret_umap = TRUE,
#' #   dim1 = "3D",
#' #   p_type = "single",
#' #   p_lab = FALSE
#' # )
#'
#' @export
ms_dim_rd <- function( # nolint
  mat1,
  md,
  md_var,
  dim_type,
  ret_umap = FALSE,
  dim1,
  p_type,
  p_lab = TRUE
) {
  d1 <- mat1
  md1 <- md
  if(dim_type == "PCA") { # nolint
    p1 <- prcomp(d1)
    vrc <- summary(p1)$importance[2, ]
    p2 <- setNames(
      p1[["x"]][, 1:3],
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
    p2 <- cbind(
      p2,
      md1
    )
    # Plot
    if(p_type == "single") { # nolint
      if(dim1 == "2D") { # nolint
        if(p_lab == TRUE) { # nolint
          p3 <- ggplot2::ggplot(
            p2,
            ggplot2::aes(
              x = .data[["PC1"]], # nolint
              y = .data[["PC2"]], # nolint
              color = .data[[md_var]], # nolint
              label = .data[[md_var]] # nolint
            )
          ) +
            ggplot2::geom_point(
              shape = 16,
              size = 2,
              alpha = 0.6
            ) +
            ggrepel::geom_text_repel(
              data = setNames(
                aggregate(
                  p2[, 1:2],
                  list(p2[[md_var]]),
                  FUN = median
                ),
                c(md_var, names(
                  p2[, 1:2]
                )
                )
              ),
              size = 4,
              bg.color = "white"
            ) +
            ggplot2::scale_color_manual(
              paste(""),
              values = col_univ() # nolint
            ) +
            ms_theme1() + # nolint
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
        if(p_lab == FALSE) { # nolint
          p3 <- ggplot2::ggplot(
            p2,
            ggplot2::aes(
              x = .data[["PC1"]], # nolint
              y = .data[["PC2"]], # nolint
              color = .data[[md_var]], # nolint
              label = .data[[md_var]] # nolint
            )
          ) +
            ggplot2::geom_point(
              shape = 16,
              size = 2,
              alpha = 0.6
            ) +
            ggplot2::labs(
              y = names(p2[2]),
              x = names(p2[1])
            ) +
            ggplot2::scale_color_manual(
              paste(""),
              values = col_univ() # nolint
            ) +
            ms_theme1() + # nolint
            ggplot2::theme(
              plot.margin = ggplot2::unit(
                c(0.1, 0.1, 0.1, 0.1), "cm"
              )
            )
        }
      }
      if(dim1 == "3D") { # nolint
        p3 <- plotly::plot_ly(
          p2,
          x = ~.data[["PC1"]], # nolint
          y = ~.data[["PC2"]], # nolint
          z = ~.data[["PC3"]], # nolint
          color = ~as.factor(.data[[md_var]]), # nolint
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
          file = paste("analysis/plot.3D.pca.", md_var, ".html", sep = "") # nolint
        )
      }
    }
    if(p_type == "list") { # nolint
      d2_list <- setNames(
        lapply(
          c(md_var), # nolint
          function(x) {
            p2[, c(
              "PC1",
              "PC2",
              x
            )
            ]
          }
        ),
        c(md_var)
      )
      d2_plot <- lapply(
        c(md_var), # nolint
        function(x) {
          p <-  ggplot2::ggplot(
            d2_list[[x]],
            ggplot2::aes(
              x = .data[["PC1"]], # nolint
              y = .data[["PC2"]], # nolint
              color = as.factor(.data[[x]]), # nolint
              label = as.factor(.data[[x]]) # nolint
            )
          ) +
            ggplot2::labs(
              y = names(p2[2]),
              x = names(p2[1])
            ) +
            ggplot2::scale_color_manual(
              paste(""),
              values = col_univ() # nolint
            ) +
            # Add points
            ggplot2::geom_point(
              shape = 16,
              size = 2,
              alpha = 0.6
            ) +
            ms_theme1() + # nolint
            ggplot2::theme(
              panel.grid.major.y = ggplot2::element_blank(),
              axis.text.x = ggplot2::element_blank(),
              axis.text.y = ggplot2::element_blank(),
              axis.title.x = ggplot2::element_blank(),
              axis.title.y = ggplot2::element_blank(),
              axis.ticks = ggplot2::element_blank(),
              plot.margin = ggplot2::unit(
                c(0.1, 0.1, 0.1, 0.1),
                "cm"
              ),
              legend.position = c(
                0.9,
                0.85
              )
            )
          return(p)
        }
      )
      # Combine output
      p3 <- ggpubr::ggarrange(
        plotlist = d2_plot,
        labels = c(
          names(
            d2_plot
          )
        ),
        ncol = ifelse(
          length(
            d2_plot
          ) <= 3,
          length(
            d2_plot
          ),
          3
        ),
        nrow = ifelse(
          length(
            d2_plot
          ) > 3,
          ceiling(
            length(
              d2_plot
            ) /
              3
          ),
          1
        )
      )
    }
  }
  if(dim_type == "UMAP") { # nolint
    if(ret_umap == TRUE) { # nolint
      p1 <- umap::umap(
        d1,
        method = "umap-learn",
        n_epochs = 500,
        n_components = 3,
        verbose = TRUE,
        preserve.seed = FALSE
      )
    }
    if(ret_umap == FALSE) { # nolint
      p1 <- umap::umap(
        d1,
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
    p2 <- cbind(
      p2,
      md1
    )
    if(p_type == "single") { # nolint
      if(dim1 == "2D") { # nolint
        if(p_lab == TRUE) { # nolint
          p3 <- ggplot2::ggplot(
            p2,
            ggplot2::aes(
              x=`UMAP.1`, # nolint
              y=`UMAP.2`, # nolint
              color = .data[[md_var]], # nolint
              label = .data[[md_var]] # nolint
            )
          ) +
            ggplot2::geom_point(
              shape = 16,
              size = 2,
              alpha = 0.6
            ) +
            ggrepel::geom_text_repel(
              data = setNames(
                aggregate(
                  p2[, c("UMAP.1", "UMAP.2")],
                  list(p2[[md_var]]),
                  FUN = median
                ),
                c(md_var, names(
                  p2[, c("UMAP.1", "UMAP.2")]
                )
                )
              ),
              size = 4,
              bg.color = "white"
            ) +
            ggplot2::scale_color_manual(
              paste(""),
              values = col_univ() # nolint
            ) +
            ms_theme1() + # nolint
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
        if(p_lab == FALSE) { # nolint
          p3 <- ggplot2::ggplot(
            p2,
            ggplot2::aes(
              x=`UMAP.1`, # nolint
              y=`UMAP.2`, # nolint
              color = .data[[md_var]], # nolint
              label = .data[[md_var]] # nolint
            )
          ) +
            ggplot2::geom_point(
              shape = 16,
              size = 2,
              alpha = 0.6
            ) +
            ggplot2::scale_color_manual(
              paste(""),
              values = col_univ() # nolint
            ) +
            ms_theme1() + # nolint
            ggplot2::theme(
              plot.margin = ggplot2::unit(
                c(0.1, 0.1, 0.1, 0.1), "cm"
              )
            )
        }
      }
      if(dim1 == "3D") { # nolint
        p3 <- plotly::plot_ly(
          p2,
          x = ~`UMAP.1`, # nolint
          y = ~`UMAP.2`, # nolint
          z = ~`UMAP.3`, # nolint
          color = ~as.factor(.data[[md_var]]), # nolint
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
          file = paste("analysis/plot.3D.umap.", md_var, ".html", sep = "") # nolint
        )
      }
    }
    if(p_type == "list") { # nolint
      d2_list <- setNames(
        lapply(
          c(md_var), # nolint
          function(x) {
            p2[, c(
              x,
              "UMAP.1",
              "UMAP.2"
            )
            ]
          }
        ),
        c(md_var)
      )
      d2_plot <- lapply(
        c(md_var), # nolint
        function(x) {
          p <-  ggplot2::ggplot(
            d2_list[[x]],
            ggplot2::aes(
              x=`UMAP.1`, # nolint
              y=`UMAP.2`, # nolint
              color = as.factor(.data[[x]]), # nolint
              label = as.factor(.data[[x]]) # nolint
            )
          ) +
            ggplot2::scale_color_manual(
              paste(""),
              values = col_univ() # nolint
            ) +
            # Add points
            ggplot2::geom_point(
              shape = 16,
              size = 2,
              alpha = 0.6
            ) +
            ms_theme1() + # nolint
            ggplot2::theme(
              panel.grid.major.y = ggplot2::element_blank(),
              axis.text.x = ggplot2::element_blank(),
              axis.text.y = ggplot2::element_blank(),
              axis.title.x = ggplot2::element_blank(),
              axis.title.y = ggplot2::element_blank(),
              axis.ticks = ggplot2::element_blank(),
              plot.margin = ggplot2::unit(
                c(0.1, 0.1, 0.1, 0.1),
                "cm"
              ),
              legend.position = c(
                0.9,
                0.85
              )
            )
          return(p)
        }
      )
      # Combine output
      p3 <- ggpubr::ggarrange(
        plotlist = d2_plot,
        labels = c(
          names(
            d2_plot
          )
        ),
        ncol = ifelse(
          length(
            d2_plot
          ) <= 3,
          length(
            d2_plot
          ),
          3
        ),
        nrow = ifelse(
          length(
            d2_plot
          ) > 3,
          ceiling(
            length(
              d2_plot
            ) /
              3
          ),
          1
        )
      )
    }
  }
  if(dim_type == "PLSDA") { # nolint

  }
  return(p3)
}
