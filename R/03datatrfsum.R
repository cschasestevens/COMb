#' Data Transformation, Scaling and Quality Check
#'
#' Performs data imputation, log2-transformation, scaling, and
#' quality checks on an untargeted dataset.
#'
#' @param ld Data list generated from ms_input().
#' @param cl_var Variable for annotating sample correlation heatmap.
#' @param samp_id Variable containing sample IDs.
#' @param sc_d Transform and scale data?
#' @param sc_meth Data scaling method (performed after log2 transformation);
#' Currently available methods are "Median" or "Pareto" (default).
#' @param hm_w Correlation heatmap width.
#' @param hm_h Correlation heatmap height.
#' @return List containing formatted input data and metadata
#' for downstream analysis.
#' @examples
#'
#' # ms_data_check(d)
#'
#' @export
ms_data_check <- function(
  ld,
  cl_var,
  samp_id,
  sc_d = TRUE,
  sc_meth = "Pareto",
  hm_w = 20,
  hm_h = 20
) {
  # input data from ms_input()
  ld1 <- ld
  if(sc_d == FALSE) { # nolint
    ld1[["data.scale"]] <- ld1[["data"]]
  }
  if(sc_d == TRUE) { # nolint
    # data imputation
    ## replace zeroes or NA with 1/10th of
    ## the lowest non-zero value present
    ld1[["imputed"]] <- ld1[["data"]]
    ld1[["imputed"]] <- setNames(
      as.data.frame(
        lapply(
          seq.int(1, ncol(ld1[["imputed"]]), 1),
          function(x) {
            ld1[["imputed"]][[x]][is.na(ld1[["imputed"]][[x]])] <- 0
            d1 <- ld1[["imputed"]][[x]]
            d1 <- ifelse(
              d1 == 0,
              round(0.1 * min(d1[d1 > 0]), digits = 0),
              d1
            )
            return(d1)
          }
        )
      ),
      c(names(ld1[["imputed"]]))
    )
    # Log2-transformation
    ld1[["data.log2"]] <- setNames(
      as.data.frame(
        lapply(
          seq.int(
            1,
            ncol(ld1[["imputed"]]),
            1
          ),
          function(x) {
            log2(ld1[["imputed"]][[x]])
          }
        )
      ),
      c(names(ld1[["imputed"]]))
    )
    # Median centering
    if(sc_meth == "Median") { # nolint
      ld1[["data.scale"]] <- setNames(
        dplyr::bind_rows(
          lapply(
            seq.int(1, nrow(ld1[["data.log2"]]), 1),
            function(x) {
              d <- as.data.frame(t(ld1[["data.log2"]]))[[x]]
              d <- as.data.frame(t(round(d - median(d), digits = 2)))
              return(d)
            }
          )
        ),
        names(ld1[["data.log2"]])
      )
    }
    # Pareto scaling
    if(sc_meth == "Pareto") { # nolint
      ld1[["data.scale"]] <- setNames(
        as.data.frame(
          apply(
            apply(
              as.matrix(ld1[["data.log2"]]),
              2,
              function(x) x - mean(x)
            ),
            2,
            function(y) round(y / sqrt(sd(y)), digits = 2)
          )
        ),
        c(names(ld1[["data.log2"]]))
      )
    }
  }
  # Data distributions
  dst <- function(df, n1) {
    dsty <- as.data.frame(colMeans(df))
    names(dsty) <- c("Value")
    p <- ggplot2::ggplot(
      dsty,
      ggplot2::aes(x = Value) # nolint 
    ) +
      # Density Plot
      ggplot2::geom_density(
        color = "darkslategrey",
        fill = col_univ()[[1]] # nolint
      ) +
      # Plot Theme
      ggplot2::labs(
        title = paste(n1, "Data Distribution"),
        y = "Density"
      ) +
      ms_theme1() + # nolint
      ggplot2::theme(
        plot.margin = ggplot2::unit(c(rep(0.5, 4)), "cm")
      )
    return(p)
  }

  ## Sample correlation heatmap
  # Input
  h1 <- ld1[["data.scale"]]
  h1md <- data.frame(
    "Group" = as.factor(ld1[["meta"]][[cl_var]]),
    "ID" = ld1[["meta"]][[samp_id]]
  )
  h2 <- t(
    magrittr::set_rownames(
      as.matrix(h1),
      h1md[["ID"]]
    )
  )
  ld1[["data.samp.corr"]] <- round(
    cor(
      h2,
      method = "spearman"
    ),
    digits = 2
  )
  # Heatmap colors
  fun_hm_col <- circlize::colorRamp2(
    c(-1, 0, 0.5, 1),
    colors = col_grad()[c(1, 3, 6, 12)] # nolint
  )
  set.seed(1234)
  fun_hm_bar <- list(
    cl_var = setNames(
      col_univ()[1:length(levels(h1md[["Group"]]))], # nolint
      as.character(levels(h1md[["Group"]]))
    )
  )
  # Annotations
  hm_anno_list <- list(
    hm_col <- ComplexHeatmap::HeatmapAnnotation( # nolint
      `Group` = h1md[["Group"]],
      col = fun_hm_bar,
      show_annotation_name = FALSE,
      show_legend = TRUE
    )
  )
  # Plot
  h_out <- ComplexHeatmap::Heatmap(
    ld1[["data.samp.corr"]],
    col = fun_hm_col,
    name = "Correlation",
    top_annotation = hm_anno_list[[1]],
    show_column_names = TRUE,
    show_row_names = TRUE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    heatmap_width = ggplot2::unit(hm_w, "cm"),
    heatmap_height = ggplot2::unit(hm_h, "cm"),
    column_title = "Sample Correlation"
  )
  if(sc_d == FALSE) { # nolint
    lp1 <- list(
      "plot.dist" = dst(
        ld1[["data.scale"]],
        paste(sc_meth, "Scaled", sep = " ")
      ),
      "plot.cor" = h_out
    )
  }
  if(sc_d == TRUE) { # nolint
    lp1 <- list(
      "plot.dist" = ggpubr::ggarrange(
        dst(ld1[["data"]], "Input"),
        dst(ld1[["data.log2"]], "Log2"),
        dst(ld1[["data.scale"]], paste(sc_meth, "Scaled", sep = " ")),
        nrow = 1,
        ncol = 3
      ),
      "plot.cor" = h_out
    )
  }
  return(
    list(
      "data" = ld1,
      "plots" = lp1
    )
  )
}
