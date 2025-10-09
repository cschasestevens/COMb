# QC (MA plot and normalization plot)
ggplot2::ggplot() +
  ggplot2::geom_point(
    data = fit,
    ggplot2::aes(
      x = .data[["x"]],
      y = .data[["y"]],
      color = .data[["type"]]
    ),
    shape = 16
  ) +
  ggplot2::geom_point(
    data = fit,
    ggplot2::aes(
      x = .data[["x"]],
      y = .data[["int"]]
    ),
    color = "purple",
    shape = 16
  ) +
  ggplot2::geom_hline(
    yintercept = mean(fit[["y"]]),
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    data = fit,
    ggplot2::aes(
      x = .data[["x"]],
      y = .data[["norm"]]
    ),
    color = "orange"
  )

ggplot2::ggplot(data = fit) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = .data[["A"]],
      y = .data[["M"]],
      color = .data[["type"]]
    ),
    shape = 16
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = .data[["A"]],
      y = .data[["MAnorm"]]
    ),
    shape = 16,
    color = "purple"
  ) +
  ggplot2::geom_smooth(
    ggplot2::aes(
      x = .data[["A"]],
      y = .data[["M"]]
    ),
    color = "red"
  ) +
  ggplot2::geom_smooth(
    ggplot2::aes(
      x = .data[["A"]],
      y = .data[["MAnorm"]]
    ),
    color = "blue"
  )

#' Data Transformation, Scaling and Quality Check
#'
#' Performs data imputation, log2-transformation, scaling, and
#' quality checks on an untargeted dataset.
#'
#' @param dat A data matrix containing compounds as rownames and sample IDs
#' as column names.
#' @param md A data frame containing metadata accompanying the data matrix.
#' @param md_var Grouping variable for data checks; this is generally
#' a variable containing the group IDs for each sample.
#' @param sid Sample ID variable.
#' @param sc_meth Data scaling method (performed after log2 transformation);
#' Either "Median" or "Pareto" (default).
#' @param dimtype Dimension reduction to calculate for scaled data; either "PCA"
#' or "UMAP" (default).
#' @param plab Should plot labels be shown on the loadings plot?
#' @param hw Correlation heatmap width.
#' @param hh Correlation heatmap height.
#' @param fsr Heatmap row names fontsize.
#' @param fsc Heatmap column names fontsize.
#' @param trans Transpose input data?
#' @return List containing formatted input data and metadata
#' for downstream analysis.
#' @examples
#'
#' # ms_data_check(d, md)
#'
#' @import dplyr
#' @import reshape2
#' @import magrittr
#' @import ComplexHeatmap
#' @import circlize
#' @import ggplot2
#' @import grid
#' @import ggpubr
#'
#' @export
ms_data_check <- function(
  dat,
  md,
  md_var = "Group",
  sid = "Label",
  sc_meth = "Pareto",
  dimtype = "UMAP",
  plab = TRUE,
  hw = 20,
  hh = 20,
  fsr = 8,
  fsc = 8,
  trans = FALSE
) {
  # input data and metadata
  ld1 <- dat
  md1 <- md
  mdv <- md_var
  ## Input data
  dimr <- list(
    "input" = ld1,
    "meta" = md1
  )
  dimr[["imputed"]] <- dimr[["input"]]
  ### replace NA with 0
  dimr[["imputed"]][is.na(dimr[["imputed"]])] <- 0
  ### remove missing compounds
  if(trans == FALSE) { # nolint
    dimr[["imputed"]] <- dimr[["imputed"]][colSums(dimr[["imputed"]]) > 0]
  }
  if(trans == TRUE) { # nolint
    dimr[["imputed"]] <- dimr[["imputed"]][rowSums(dimr[["imputed"]]) > 0, ]
  }
  ### replace zeroes with 1/10th of the lowest non-zero
  ### value per compound
  dimr[["imputed"]] <- magrittr::set_rownames(dplyr::bind_cols(
    lapply(
      seq.int(1, ncol(dimr[["imputed"]]), 1),
      function(x) {
        d1 <- dimr[["imputed"]][[x]]
        d1 <- ifelse(
          d1 == 0,
          0.1 * min(d1[d1 > 0]),
          d1
        )
        d1 <- setNames(as.data.frame(d1), names(dimr[["imputed"]])[[x]])
        return(d1) # nolint
      }
    )
  ), rownames(dimr[["imputed"]]))
  ## Log2-transformation
  dimr[["data.log2"]] <- magrittr::set_rownames(setNames(
    as.data.frame(
      lapply(
        seq.int(1, ncol(dimr[["imputed"]]), 1),
        function(x) {
          log2(dimr[["imputed"]][[x]])
        }
      )
    ),
    c(names(dimr[["imputed"]]))
  ), row.names(dimr[["imputed"]]))
  ## Median centering
  if(sc_meth == "Median") { # nolint
    dimr[["data.scale"]] <- setNames(
      dplyr::bind_rows(
        lapply(
          seq.int(1, nrow(dimr[["data.log2"]]), 1),
          function(x) {
            d <- as.data.frame(t(dimr[["data.log2"]]))[[x]]
            d <- as.data.frame(t(round(d - median(d), digits = 2)))
            return(d) # nolint
          }
        )
      ),
      names(dimr[["data.log2"]])
    )
  }
  ## Pareto scaling
  if(sc_meth == "Pareto") { # nolint
    dimr[["data.scale"]] <- setNames(
      as.data.frame(
        apply(
          apply(
            as.matrix(dimr[["data.log2"]]),
            2,
            function(x) x - mean(x)
          ),
          2,
          function(y) round(y / sqrt(sd(y)), digits = 2)
        )
      ),
      c(names(dimr[["data.log2"]]))
    )
    head(dimr[["data.scale"]])
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
      ms_theme() + # nolint
      ggplot2::theme(
        plot.margin = ggplot2::unit(c(rep(0.5, 4)), "cm")
      )
    return(p) # nolint
  }
  ## Combined distribution plot
  dimr[["pdist"]] <- ggpubr::ggarrange(
    dst(dimr[["input"]], "Input"),
    dst(dimr[["data.log2"]], "Log2"),
    dst(dimr[["data.scale"]], paste(sc_meth, "Scaled", sep = " ")),
    nrow = 1,
    ncol = 3
  )
  # Sample correlation heatmap
  ## Input data
  h2 <- dimr[["data.scale"]]
  h1md <- data.frame(
    "Group" = as.factor(sort(dimr[["meta"]][[mdv]])),
    "ID" = dimr[["meta"]][[sid]]
  )
  h2 <- t(
    magrittr::set_rownames(
      as.matrix(h2),
      h1md[["ID"]]
    )
  )
  dimr[["heatsampin"]] <- round(
    cor(
      h2,
      method = "spearman"
    ),
    digits = 2
  )
  # Heatmap colors
  fun_hm_col <- circlize::colorRamp2(
    c(-1, 0, 0.5, 1),
    colors = col_grad(scm = 3) # nolint
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
  dimr[["heatsamp"]] <- ComplexHeatmap::Heatmap(
    dimr[["heatsampin"]],
    col = fun_hm_col,
    name = "Correlation",
    top_annotation = hm_anno_list[[1]],
    show_column_names = TRUE,
    show_row_names = TRUE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    heatmap_width = ggplot2::unit(hw, "cm"),
    heatmap_height = ggplot2::unit(hh, "cm"),
    column_title = "Sample Correlation",
    row_names_gp = grid::gpar(fontsize = fsr),
    column_names_gp = grid::gpar(fontsize = fsc)
  )
  # Dimension reduction plot
  dimr[["dimr"]] <- ms_dim_rd(
    mat1 = dimr[["data.scale"]],
    md = dimr[["meta"]],
    md_var = mdv,
    sid = sid,
    dim_type = dimtype,
    p_lab = plab,
    data_scale = FALSE
  )
  return(dimr)
}

#' Data Transformation, Scaling and Quality Check (Targeted)
#'
#' Performs data imputation, log2-transformation, scaling, and
#' quality checks on a targeted dataset.
#'
#' @param dat A data matrix containing compounds as rownames and sample IDs
#' as column names.
#' @param md A data frame containing metadata accompanying the data matrix.
#' @param md_var Grouping variable for data checks; this is generally
#' a variable containing the group IDs for each sample.
#' @param sid Sample ID variable.
#' @param sc_meth Data scaling method (performed after log2 transformation);
#' Either "Median" or "Pareto" (default).
#' @param dimtype Dimension reduction to calculate for scaled data; either "PCA"
#' or "UMAP" (default).
#' @param plab Should plot labels be shown on the loadings plot?
#' @param hw Presence/absence heatmap width.
#' @param hh Presence/absence heatmap height.
#' @param hw2 Correlation heatmap width.
#' @param hh2 Correction heatmap height.
#' @param fsr Heatmap row names fontsize.
#' @param fsc Heatmap column names fontsize.
#' @param trans Transpose input data?
#' @return List containing formatted input data and metadata
#' for downstream analysis.
#' @examples
#'
#' # ms_data_check(d)
#'
#' @import dplyr
#' @import reshape2
#' @import magrittr
#' @import ComplexHeatmap
#' @import circlize
#' @import ggplot2
#' @import grid
#' @import ggpubr
#'
#' @export
tms_data_check <- function(
  dat,
  md,
  md_var = "Group",
  sid = "SampleID",
  sc_meth = "Pareto",
  dimtype = "UMAP",
  plab = TRUE,
  hw = 20,
  hh = 20,
  hw2 = 20,
  hh2 = 20,
  fsr = 8,
  fsc = 8,
  trans = FALSE
) {
  # input data and metadata
  ld1 <- dat
  md1 <- md
  mdv <- md_var
  # Check presence/absence for each sample group
  if(trans == FALSE) {dchk <- ld1} # nolint
  if(trans == TRUE) { # nolint
    dchk <- setNames(as.data.frame(t(ld1)), rownames(ld1))
  }
  dchk2 <- dplyr::bind_rows(lapply(
    seq.int(1, ncol(dchk), 1),
    function(x) {
      dcnt <- setNames(aggregate(
        dchk[[x]],
        list(md1[[mdv]]),
        function(y) {
          d1 <- y > 0 & !is.na(y)
          d1 <- length(d1[d1 == TRUE])
        }
      ), c("Group", "detected"))
      dtot <- setNames(aggregate(
        dchk[[x]],
        list(md1[[mdv]]),
        function(y) length(y)
      ), c("Group", "total"))
      dout <- dplyr::left_join(dcnt, dtot, by = c("Group"))
      dout[["prop"]] <- round(
        dout[["detected"]] / dout[["total"]],
        digits = 3
      )
      dout[["Name"]] <- names(dchk)[[x]]
      dout <- dplyr::select(dout, "Name", everything()) # nolint
      return(dout) # nolint
    }
  ))
  dchk2 <- reshape2::dcast(
    dchk2,
    Name ~ Group,
    value.var = "prop"
  )
  ### Presence/absence matrix
  dchk2 <- magrittr::set_rownames(
    dchk2, dchk2[["Name"]]
  )[-1]
  # Presence/absence heatmap
  h1 <- as.matrix(dchk2)
  # Heatmap colors
  fun_hm_col <- circlize::colorRamp2(
    c(0, 0.5, 1),
    colors = col_grad(scm = 4) # nolint
  )
  h_out <- ComplexHeatmap::Heatmap(
    h1,
    col = fun_hm_col,
    name = "Proportion",
    show_column_names = TRUE,
    show_row_names = TRUE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    heatmap_width = ggplot2::unit(hw, "cm"),
    heatmap_height = ggplot2::unit(hh, "cm"),
    column_title = paste("Presence/absence"),
    row_names_gp = grid::gpar(fontsize = fsr),
    column_names_gp = grid::gpar(fontsize = fsc)
  )
  # Imputed data sets for dimension reduction
  ## Input data
  dimr <- list(
    "input" = ld1,
    "meta" = md1,
    "heatpres" = h_out,
    "heatpresin" = dchk2
  )
  dimr[["imputed"]] <- dimr[["input"]]
  ### replace NA with 0
  dimr[["imputed"]][is.na(dimr[["imputed"]])] <- 0
  ### remove missing compounds
  if(trans == FALSE) { # nolint
    dimr[["imputed"]] <- dimr[["imputed"]][colSums(dimr[["imputed"]]) > 0]
  }
  if(trans == TRUE) { # nolint
    dimr[["imputed"]] <- dimr[["imputed"]][rowSums(dimr[["imputed"]]) > 0, ]
  }
  ### replace zeroes with 1/10th of the lowest non-zero
  ### value per compound
  dimr[["imputed"]] <- magrittr::set_rownames(dplyr::bind_cols(
    lapply(
      seq.int(1, ncol(dimr[["imputed"]]), 1),
      function(x) {
        d1 <- dimr[["imputed"]][[x]]
        d1 <- ifelse(
          d1 == 0,
          0.1 * min(d1[d1 > 0]),
          d1
        )
        d1 <- setNames(as.data.frame(d1), names(dimr[["imputed"]])[[x]])
        return(d1) # nolint
      }
    )
  ), rownames(dimr[["imputed"]]))
  ## Log2-transformation
  dimr[["data.log2"]] <- magrittr::set_rownames(setNames(
    as.data.frame(
      lapply(
        seq.int(1, ncol(dimr[["imputed"]]), 1),
        function(x) {
          log2(dimr[["imputed"]][[x]])
        }
      )
    ),
    c(names(dimr[["imputed"]]))
  ), row.names(dimr[["imputed"]]))
  ## Median centering
  if(sc_meth == "Median") { # nolint
    dimr[["data.scale"]] <- setNames(
      dplyr::bind_rows(
        lapply(
          seq.int(1, nrow(dimr[["data.log2"]]), 1),
          function(x) {
            d <- as.data.frame(t(dimr[["data.log2"]]))[[x]]
            d <- as.data.frame(t(round(d - median(d), digits = 2)))
            return(d) # nolint
          }
        )
      ),
      names(dimr[["data.log2"]])
    )
  }
  ## Pareto scaling
  if(sc_meth == "Pareto") { # nolint
    dimr[["data.scale"]] <- setNames(
      as.data.frame(
        apply(
          apply(
            as.matrix(dimr[["data.log2"]]),
            2,
            function(x) x - mean(x)
          ),
          2,
          function(y) round(y / sqrt(sd(y)), digits = 2)
        )
      ),
      c(names(dimr[["data.log2"]]))
    )
    head(dimr[["data.scale"]])
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
      ms_theme() + # nolint
      ggplot2::theme(
        plot.margin = ggplot2::unit(c(rep(0.5, 4)), "cm")
      )
    return(p) # nolint
  }
  ## Combined distribution plot
  dimr[["pdist"]] <- ggpubr::ggarrange(
    dst(dimr[["input"]], "Input"),
    dst(dimr[["data.log2"]], "Log2"),
    dst(dimr[["data.scale"]], paste(sc_meth, "Scaled", sep = " ")),
    nrow = 1,
    ncol = 3
  )
  # Sample correlation heatmap
  ## Input data
  h2 <- dimr[["data.scale"]]
  h1md <- data.frame(
    "Group" = as.factor(sort(dimr[["meta"]][[mdv]])),
    "ID" = dimr[["meta"]][[sid]]
  )
  h2 <- t(
    magrittr::set_rownames(
      as.matrix(h2),
      h1md[["ID"]]
    )
  )
  dimr[["heatsampin"]] <- round(
    cor(
      h2,
      method = "spearman"
    ),
    digits = 2
  )
  # Heatmap colors
  fun_hm_col <- circlize::colorRamp2(
    c(-1, 0, 0.5, 1),
    colors = col_grad(scm = 3) # nolint
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
  dimr[["heatsamp"]] <- ComplexHeatmap::Heatmap(
    dimr[["heatsampin"]],
    col = fun_hm_col,
    name = "Correlation",
    top_annotation = hm_anno_list[[1]],
    show_column_names = TRUE,
    show_row_names = TRUE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    heatmap_width = ggplot2::unit(20, "cm"),
    heatmap_height = ggplot2::unit(20, "cm"),
    column_title = "Sample Correlation",
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8)
  )
  # Dimension reduction plot
  dimr[["dimr"]] <- ms_dim_rd(
    mat1 = dimr[["data.scale"]],
    md = dimr[["meta"]],
    md_var = md_var,
    dim_type = dimtype,
    p_lab = plab,
    data_scale = FALSE
  )
  return(dimr)
}

#' Blank Subtraction
#'
#' Removes compounds from a dataset based on intensity in method blanks.
#'
#' @param dat A data matrix containing compounds as rownames and sample IDs
#' as column names.
#' @param dat_bl A data matrix of method blank intensities for each compound
#' in the same format as the input data.
#' @param md A data frame containing metadata accompanying the data matrix.
#' Must have a column indicating sample groups.
#' @param md_var Grouping variable for data checks; this is generally
#' a variable containing the group IDs for each sample.
#' @param sid Sample ID variable.
#' @param threshold Signal intensity threshold for removing a compound from
#' the dataset. The default is to remove compounds that are < 3-fold the
#' signal intensity of the method blanks.
#' @param sub_blanks Should compounds present in blanks be removed?
#' If FALSE (default), only missing compounds are removed from the dataset.
#' @return A filtered data frame retaining compounds that exceed the
#' blank intensity threshold.
#' @examples
#'
#' # ms_blank_sub(
#' #   dat = ld[["wash"]][["input"]],
#' #   dat_bl = ld[["qc"]][["input"]],
#' #   md = ld[["wash"]][["meta"]]
#' # )
#'
#' @import dplyr
#' @import magrittr
#' @export
ms_blank_sub <- function(
  dat,
  dat_bl,
  md,
  md_var = "Group",
  sid = "SampleID",
  threshold = 3,
  sub_blanks = FALSE
) {
  # load data
  dexp <- dat
  dblnk <- dat_bl
  meta <- md
  # Calculate group averages for each compound
  ## blanks
  dblnk <- as.data.frame(t(dblnk[grepl("MB|blank|Blank|Method", names(dblnk))]))
  dblnk
  dblnk2 <- data.frame(
    "Name" = names(dblnk),
    "mean_blank" = unlist(
      lapply(
        seq.int(1, ncol(dblnk), 1),
        function(i) mean(dblnk[[i]])
      )
    )
  )
  ## experimental samples
  dexp <- as.data.frame(t(dexp))
  dexp <- dexp[sort(meta[[sid]]), ]
  dexp2 <- as.data.frame(t(dplyr::bind_cols(lapply(
    seq.int(1, ncol(dexp), 1),
    function(i) {
      ld1 <- aggregate(
        dexp[[i]],
        list(meta[[md_var]]),
        function(x) mean(x)
      )
      ld1 <- magrittr::set_rownames(
        setNames(
          as.data.frame(ld1[[2]]),
          names(dexp)[[i]]
        ),
        ld1[[1]]
      )
      return(ld1) # nolint
    }
  ))))
  dexp2 <- setNames(
    dplyr::select(
      dplyr::mutate(
        dexp2,
        "Name" = row.names(dexp2)
      ),
      "Name", everything() # nolint
    ),
    c("Name", paste("mean_", names(dexp2), sep = ""))
  ) # nolint
  # Merge and filter
  dfilt <- dplyr::left_join(
    dexp2,
    dblnk2,
    by = "Name"
  )
  ## flag for removal
  if(sub_blanks == TRUE) { # nolint
    lfilt <- data.frame(
      "Name" = dfilt[["Name"]],
      dplyr::bind_cols(
        lapply(
          seq.int(2, ncol(dexp2), 1),
          function(x) {
            setNames(
              as.data.frame(
                # remove compounds present in blanks
                dfilt[[x]] > (dfilt[["mean_blank"]] * threshold) &
                  # remove missing compounds
                  dfilt[[x]] > 0
              ),
              names(dfilt)[[x]]
            )
          }
        )
      )
    )
  }
  if(sub_blanks == FALSE) { # nolint
    lfilt <- data.frame(
      "Name" = dfilt[["Name"]],
      dplyr::bind_cols(
        lapply(
          seq.int(2, ncol(dexp2), 1),
          function(x) {
            setNames(
              as.data.frame(
                # remove missing compounds only
                dfilt[[x]] > 0
              ),
              names(dfilt)[[x]]
            )
          }
        )
      )
    )
  }
  ## remove compounds that do not exceed threshold
  ## in all groups (except any internal standards)
  ## or compounds that are missing in all samples
  lfilt[
    grepl("PUHA|CUDA|d3-|d4-|d5-|d6-|d7-|d8-|d11-", lfilt[["Name"]]),
    2:ncol(lfilt)
  ] <- TRUE
  lfilt <- data.frame(
    "Name" = lfilt[["Name"]],
    "Remove" = unlist(
      lapply(
        seq.int(1, nrow(lfilt), 1),
        function(i) {
          ifelse(
            length(unique(as.vector(lfilt[i, 2:ncol(lfilt)]))[
              unique(as.vector(lfilt[i, 2:ncol(lfilt)])) == FALSE
            ]) == 1,
            TRUE,
            FALSE
          )
        }
      )
    )
  )
  print("Removed the following compounds from the dataset:")
  print(dfilt[lfilt[["Remove"]], "Name"])
  print(
    paste(
      length(dfilt[lfilt[["Remove"]], "Name"]),
      " of ",
      nrow(dfilt) -
        length(lfilt[
          grepl("PUHA|CUDA|d3-|d4-|d5-|d6-|d7-|d8-|d11-", lfilt[["Name"]]),
          "Name"
        ]),
      " compounds removed in total.",
      sep = ""
    )
  )
  dexp <- dexp[lfilt[lfilt[["Remove"]] == FALSE, "Name"]]
  out <- list(
    "filt" = dexp,
    "flag" = lfilt
  )
  return(out) # nolint
}

#' Check Compound Presence Per Group
#'
#' Checks compound presence per treatment group for the selected dataset
#' and removes compounds that are not present in at least 50% of samples
#' in one group.
#'
#' @param dat A list generated by COMb functions containing compounds
#' intensities, metadata, and annotation information.
#' @param hw Presence/absence heatmap width.
#' @param hh Presence/absence heatmap height.
#' @param fsc Heatmap column names fontsize.
#' @param trans Transpose input data?
#' @return List containing formatted input data and metadata
#' for downstream analysis.
#' @examples
#'
#' # ms_data_check_pres(
#' #   dat = d_raw[[1]]
#' # )
#'
#' @import dplyr
#' @import reshape2
#' @import magrittr
#' @import ComplexHeatmap
#' @import circlize
#' @import grid
#'
#' @export
ms_data_check_pres <- function(
  dat,
  md_var = "Group",
  md_anno = "Label",
  hw = 20,
  hh = 20,
  fsc = 8,
  trans = TRUE
) {
  # input data and metadata
  ld1 <- dat[["data"]]
  md1 <- dat[["meta"]]
  md2 <- dat[["anno"]]
  mdv <- md_var
  # Check presence/absence for each sample group
  if(trans == FALSE) {dchk <- ld1} # nolint
  if(trans == TRUE) { # nolint
    dchk <- setNames(as.data.frame(t(ld1)), rownames(ld1))
  }
  dchk2 <- dplyr::bind_rows(lapply(
    seq.int(1, ncol(dchk), 1),
    function(x) {
      dcnt <- setNames(aggregate(
        dchk[[x]],
        list(md1[[mdv]]),
        function(y) {
          d1 <- y > 0 & !is.na(y)
          d1 <- length(d1[d1 == TRUE])
        }
      ), c("Group", "detected"))
      dtot <- setNames(aggregate(
        dchk[[x]],
        list(md1[[mdv]]),
        function(y) length(y)
      ), c("Group", "total"))
      dout <- dplyr::left_join(dcnt, dtot, by = c("Group"))
      dout[["prop"]] <- round(
        dout[["detected"]] / dout[["total"]],
        digits = 3
      )
      dout[["Name"]] <- names(dchk)[[x]]
      dout <- dplyr::select(dout, "Name", everything()) # nolint
      return(dout) # nolint
    }
  ))
  dchk2 <- reshape2::dcast(
    dchk2,
    Name ~ Group,
    value.var = "prop"
  )
  ### Presence/absence matrix
  dchk2 <- magrittr::set_rownames(
    dchk2, make.names(md2[[md_anno]], unique = TRUE)
  )[-1]
  # Presence/absence heatmap
  h1 <- as.matrix(dchk2)
  # Heatmap colors
  fun_hm_col <- circlize::colorRamp2(
    c(0, 0.5, 1),
    colors = col_grad(scm = 4) # nolint
  )
  h_out <- ComplexHeatmap::Heatmap(
    h1,
    col = fun_hm_col,
    name = "Proportion",
    show_column_names = TRUE,
    show_row_names = FALSE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    heatmap_width = ggplot2::unit(hw, "cm"),
    heatmap_height = ggplot2::unit(hh, "cm"),
    column_title = paste("Presence/absence"),
    column_names_gp = grid::gpar(fontsize = fsc)
  )
  # Imputed data sets for dimension reduction
  ## Input data
  dimr <- list(
    "data" = dchk,
    "meta" = md1,
    "anno" = md2,
    "heatpres" = h_out,
    "heatpresin" = dchk2
  )
  ### replace NA with 0
  dimr[["data"]][is.na(dimr[["data"]])] <- 0
  ### remove missing compounds
  drem <- unlist(apply(dchk2, 1, function(y) max(y) > 0.5))
  if(trans == FALSE) { # nolint
    dimr[["data"]] <- dimr[["data"]][, drem]
  }
  if(trans == TRUE) { # nolint
    dimr[["data"]] <- dimr[["data"]][, drem]
    dimr[["anno"]] <- dimr[["anno"]][drem, ]
    print(paste(
      "Removed ", ncol(dchk) - ncol(dimr[["data"]]),
      " low abundance compounds.", sep = ""
    ))
  }
  return(dimr)
}
