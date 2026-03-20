#' QC Summary Report
#'
#' Creates a markdown-style report from a ms_qc QC summary.
#'
#' @param data_qc Input object from ms_qc.
#' @param study_name Study name.
#' @param author Report author name.
#' @param an_plat Name of analysis platform.
#' @param an_date Analysis date(s).
#' @param meth_dp Data processing software/methods.
#' @param meth_nm Normalization method.
#' @param outdir Output directory.
#' @param outname Report filename.
#' @return Group RSD values for each feature.
#' @examples
#'
#' # ms_qc_report(
#' #   qcdata,
#' #   outdir = "path/to/directory/",
#' #   outname = "qcreport"
#' # )
#'
#' @export
ms_qc_report <- function(
  data_qc,
  study_name = NULL,
  author = NULL,
  an_plat = NULL,
  an_date = NULL,
  meth_dp = NULL,
  meth_nm = NULL,
  outdir,
  outname
) {
  # Load data and define parameters
  ## Input QC summary list
  d1 <- data_qc
  ## Study name
  if (!is.null(study_name)) {sn1 <- metadata(d1[["input data"]])[[study_name]]} # nolint
  if (is.null(study_name)) {sn1 <- "N/A"} # nolint
  ## Report author
  if (!is.null(author)) {sa1 <- metadata(d1[["input data"]])[[author]]} # nolint
  if (is.null(author)) {sa1 <- "N/A"} # nolint
  ## analysis platform
  if (!is.null(an_plat)) {ap1 <- metadata(d1[["input data"]])[[an_plat]]} # nolint
  if (is.null(an_plat)) {ap1 <- "N/A"} # nolint
  ## analysis date
  if (!is.null(an_date)) {ad1 <- metadata(d1[["input data"]])[[an_date]]} # nolint
  if (is.null(an_date)) {ad1 <- "N/A"} # nolint
  ## modes/batches included in report
  sd1 <- unlist(lapply(
    seq.int(1, length(names(d1[["transformed data"]])), 1),
    function(i) {
      lapply(
        seq.int(1, length(names(d1[["transformed data"]][[i]])), 1),
        function(j) {
          paste(
            names(d1[["transformed data"]])[[i]],
            names(d1[["transformed data"]][[i]])[[j]]
          )
        }
      )
    }
  ))
  ## sample numbers
  sfn1 <- paste(
    d1[["data summary"]][[1]][[1]], " Total Samples + QCs and ",
    d1[["data summary"]][[1]][[2]], " Unique Features; ",
    paste(
      d1[["data summary"]][[2]][[2]],
      d1[["data summary"]][[2]][[1]],
      collapse = ", "
    )
  )
  # RSD summary table
  rsd01 <- dplyr::bind_rows(
    lapply(
      seq.int(1, length(d1[["rsd"]]), 1),
      function(i) {
        dplyr::bind_rows(
          lapply(
            seq.int(1, length(d1[["rsd"]][[i]]), 1),
            function(j) {
              dplyr::bind_rows(
                lapply(
                  seq.int(1, length(d1[["rsd"]][[i]][[j]])),
                  function(k) {
                    cbind(
                      data.frame(
                        "Set" = paste(
                          names(d1[["rsd"]])[[i]],
                          names(d1[["rsd"]][[i]])[[j]],
                          names(d1[["rsd"]][[i]][[j]])[[k]],
                          sep = "."
                        )
                      ),
                      d1[["rsd"]][[i]][[j]][[k]][["mRSD"]]
                    )
                  }
                )
              )
            }
          )
        )
      }
    )
  )
  ## data processing information
  if (!is.null(meth_dp)) {dp1 <- metadata(d1[["input data"]])[[meth_dp]]} # nolint
  if (is.null(meth_dp)) {dp1 <- "N/A"} # nolint
  ## normalization method
  if (!is.null(metadata(d1[["input data"]])[[meth_nm]])) {mn1 <- metadata(d1[["input data"]])[[meth_nm]]} # nolint
  if (is.null(metadata(d1[["input data"]])[[meth_nm]])) {mn1 <- "N/A"} # nolint
  # Generate report
  # Save all params as list
  rep_inline <- list(
    "sn" = sn1, "sa" = sa1, "ap" = ap1,
    "ad" = ad1, "sd" = sd1, "sfn" = sfn1,
    "dp" = dp1, "mn" = mn1, "rsd" = rsd01
  )
  rmarkdown::render(
    input = system.file("rmd", "qc_report.Rmd", package = "COMb"),
    output_file = paste(outname, ".html", sep = ""),
    output_dir = outdir,
    params = list(rep1 = rep_inline, data1 = d1),
    envir = new.env(parent = baseenv())
  )
}

#' COMb QC Report
#'
#' Performs quality control for a SummarizedExperiment that can be
#' summarized as a markdown report of key quality control metrics,
#' to assess normalization performance, including: dataset information,
#' data distribution, technical variance, sample correlation,
#' and outlier detection. The following minimum information should be
#' specified in the object metadata and included in the column
#' and/or feature data:
#' 'Group column name', 'Acquisition order column name',
#' 'Sample column name', 'QC sample name', 'Sample ID column name'.
#' Additionally, include 'Blank sample name' if the data contain
#' analytical method blank samples and 'Batch column name' if the
#' data are intended to be analyzed as separate batches
#' (ex. two distinct analytical runs or two separate
#' sample types, regions, matrices, etc.). By specifying the
#' 'Batch column name', this function will automatically split the
#' data by batch and ionization mode. If 'Batch column name' is not
#' specified, the function will assume all data belong to the same
#' analytical run and will only be split by ionization mode.
#'
#' @param exp1 SummarizedExperiment containing one or more assays for
#' evaluating QC metrics.
#' @return A markdown report containing all quality control information.
#' @examples
#'
#' # ms_qc(data1)
#'
#' @export
ms_qc <- function(
  exp1
) {
  #---- Load data ----
  print("1. Loading SummarizedExperiment")
  qcdata <- exp1
  #---- Check metadata ----
  print("2. Checking metadata for prerequisite information")
  if (
    length(unique(c(
      "Feature column name",
      "Group column name",
      "Acquisition order column name",
      "Sample column name",
      "QC sample name",
      "Sample ID column name"
    ) %in%
      names(metadata(qcdata)))) != 1 && # nolint
      any(
        unique(c(
          "Feature column name",
          "Group column name",
          "Acquisition order column name",
          "Sample column name",
          "QC sample name",
          "Sample ID column name"
        ) %in%
          names(metadata(qcdata))) == FALSE
      )
  ) {
    print("Must include 'Feature column name', 'Group column name', 'Acquisition order column name', 'Sample column name', 'QC sample name', and 'Sample ID column name' in metadata to run data normalization!") # nolint
  }
  #---- Validate sample and feature data ----
  print("3. Validating information in sample and feature data")
  # Extract metadata from SummarizedExperiment
  md1 <- list(
    "fc" = metadata(qcdata)[["Feature column name"]], # nolint
    "gc" = metadata(qcdata)[["Group column name"]], # nolint
    "ac" = metadata(qcdata)[["Acquisition order column name"]], # nolint
    "sc" = metadata(qcdata)[["Sample column name"]], # nolint
    "qc" = metadata(qcdata)[["QC sample name"]], # nolint
    "sid" = metadata(qcdata)[["Sample ID column name"]] # nolint
  )
  if ("Class column name" %in% names(metadata(qcdata))) { # nolint
    print("Class column name present in metadata; will split feature counts by compound class...") # nolint
    md1 <- c(md1, "cc" = metadata(qcdata)[["Class column name"]]) # nolint
  }
  if (
    "Blank sample name" %in% names(metadata(qcdata)) && # nolint
      any(grepl(metadata(qcdata)[["Blank sample name"]], colData(qcdata)[[metadata(qcdata)[["Sample column name"]]]])) == FALSE # nolint
  ) { # nolint
    print("Data contain blank samples but metadata information does not match sample names!") # nolint
  }
  if (
    "Blank sample name" %in% names(metadata(qcdata)) && # nolint
      any(grepl(metadata(qcdata)[["Blank sample name"]], colData(qcdata)[[metadata(qcdata)[["Sample column name"]]]])) # nolint
  ) { # nolint
    md1 <- c(md1, "bc" = metadata(qcdata)[["Blank sample name"]]) # nolint
  }
  if (
    "Batch column name" %in% names(metadata(qcdata)) && # nolint
      any(grepl(metadata(qcdata)[["Batch column name"]], names(colData(qcdata)))) # nolint
  ) { # nolint
    print("Batch column name present in metadata; will split report by batch...") # nolint
    md1 <- c(md1, "batc" = metadata(qcdata)[["Batch column name"]]) # nolint
  }
  if (
    "Batch column name" %in% names(metadata(qcdata)) == FALSE # nolint
  ) { # nolint
    print("Batch column name not present in metadata; initializing metadata column for report...") # nolint
    md1 <- c(md1, "batc" = "Batch") # nolint
    colData(qcdata)[["Batch"]] <- "batch1" # nolint
  }
  if (
    "Internal standard prefix" %in% names(metadata(qcdata)) && # nolint
      any(grepl(metadata(qcdata)[["Internal standard prefix"]], rowData(qcdata)[[md1[["fc"]]]])) # nolint
  ) { # nolint
    print("Internal standard prefix present in metadata; will include additional information reporting on internal standards...") # nolint
    md1 <- c(md1, "is" = metadata(qcdata)[["Internal standard prefix"]]) # nolint
  }
  outs <- list(
    "input data" = qcdata,
    "metadata" = md1
  )
  #---- General dataset information ----
  print("4. General dataset characteristics")
  # Data structure
  gen1 <- c(ncol(qcdata), nrow(qcdata))
  # Samples
  gen2 <- dplyr::count(
    setNames(as.data.frame(colData(qcdata)), names(colData(qcdata))), # nolint
    .data[[md1[["sc"]]]]
  )
  gen3 <- dplyr::count(
    setNames(as.data.frame(colData(qcdata)), names(colData(qcdata))), # nolint
    .data[[md1[["gc"]]]]
  )
  # Features
  if ("cc" %in% names(md1) == FALSE) { # nolint
    gen4 <- dplyr::count(
      setNames(as.data.frame(rowData(qcdata)), names(rowData(qcdata))), # nolint
      .data[["Mode"]]
    )
  }
  if ("cc" %in% names(md1)) { # nolint
    gen4 <- dplyr::count(
      setNames(as.data.frame(rowData(qcdata)), names(rowData(qcdata))), # nolint
      .data[["Mode"]], .data[["Class"]]
    )
  }
  if ("is" %in% names(md1)) { # nolint
    gen5 <- rowData(qcdata)[grepl(md1[["is"]], rowData(qcdata)[[md1[["fc"]]]]), ][[md1[["fc"]]]] # nolint
  }
  outs[["data summary"]] <- list(
    gen1, gen2, gen3, gen4, gen5
  )
  #---- Data transformation, scaling, and subsetting ----
  print("5. Transforming and subsetting data for report")
  # For imputed and normalized assays, calculate log2 and pareto distributions
  ## log2 assays
  if ("imputed" %in% names(assays(qcdata)) == FALSE) { # nolint
    print("Must have 'imputed' assay in SummarizedExperiment; log2 transformation and Pareto scaling require non-zero and no missing values!") # nolint
  }
  assay(qcdata, "imputed.log2") <- log2(assay(qcdata, "imputed")) # nolint
  assay(qcdata, "norm.log2") <- log2(assay(qcdata, "norm")) # nolint
  ## pareto scale
  ### split into separate SummarizedExperiments if multiple modes/batches
  qc_inputs <- setNames(lapply(
    # Split by batch
    unique(colData(qcdata)[[md1[["batc"]]]]), # nolint
    function(i) {
      setNames(lapply(
        # Split by mode
        unique(rowData(qcdata)[["Mode"]]), # nolint
        function(j) {
          # Subset data by mode/batch
          qc1 <- qcdata[
            grepl(j, rowData(qcdata)[["Mode"]]), # nolint
            grepl(i, colData(qcdata)[[md1[["batc"]]]]) # nolint
          ]
          # remove iSTD and blanks if present
          qc2 <- qc1
          if ("is" %in% names(md1)) {
            qc2 <- qc2[!grepl(md1[["is"]], rowData(qc2)[[md1[["fc"]]]]), ] # nolint
          }
          if ("bc" %in% names(md1)) {
            qc2 <- qc2[, !grepl(md1[["bc"]], colData(qc2)[[md1[["sc"]]]])] # nolint
          }
          # pareto scale imputed and norm log2
          assay(qc2, "imputed.pareto") <- apply( # nolint
            apply(
              assay(qc2, "imputed.log2"), # nolint
              2,
              function(x) x - mean(x)
            ),
            2,
            function(y) round(y / sqrt(sd(y)), digits = 3)
          )
          assay(qc2, "norm.pareto") <- apply( # nolint
            apply(
              assay(qc2, "norm.log2"), # nolint
              2,
              function(x) x - mean(x)
            ),
            2,
            function(y) round(y / sqrt(sd(y)), digits = 3)
          )
          return(qc2) # nolint
        }
      ), unique(rowData(qcdata)[["Mode"]]))
    }
  ), unique(colData(qcdata)[[md1[["batc"]]]]))
  outs[["transformed data"]] <- qc_inputs
  #---- Presence/absence ----
  print("6. Checking compound presence/absence")
  pres1 <- setNames(
    lapply(
      seq.int(1, length(qc_inputs), 1),
      function(h) {
        p1 <- setNames(lapply(
          seq.int(1, length(qc_inputs[[h]]), 1),
          function(k) {
            print(
              paste(
                "Checking presence/absence for ",
                names(qc_inputs)[[h]], " ",
                names(qc_inputs[[h]])[[k]],
                ":", sep = ""
              )
            )
            ms_check_pres(
              qc_inputs[[h]][[k]],
              title1 = paste(
                names(qc_inputs)[[h]], " ",
                names(qc_inputs[[h]])[[k]],
                ":", sep = ""
              )
            )
          }
        ), names(qc_inputs[[h]]))
        return(p1) # nolint
      }
    ), names(qc_inputs)
  )
  outs[["compound presence"]] <- pres1
  #---- Data distributions ----
  print("7. Visualize data distributions")
  # View distribution
  pdist <- function(se1, d1, d2, c1) {
    df1 <- setNames(as.data.frame(colMeans(assay(se1, d1))), "Value") # nolint
    df2 <- setNames(as.data.frame(colMeans(assay(se1, d2))), "Value") # nolint
    p1 <- ggplot2::ggplot() +
      # Imputed
      ggplot2::geom_density(
        data = df1,
        ggplot2::aes(x = df1[["Value"]]), # nolint
        color = "white",
        fill = col_univ()[[c1]], # nolint
        alpha = 0.75,
        linetype = "dashed"
      ) +
      # Norm
      ggplot2::geom_density(
        data = df2,
        ggplot2::aes(x = df2[["Value"]]), # nolint
        color = "darkslategrey",
        fill = col_univ()[[c1]], # nolint
        alpha = 0.75
      ) +
      ggplot2::guides(x = "axis") +
      # Plot Theme
      ggplot2::labs(
        title = gsub("imputed.", "", d1),
        x = "Value",
        y = "Density"
      ) +
      ms_theme() + # nolint
      ggplot2::theme(
        plot.margin = ggplot2::unit(c(rep(0.5, 4)), "cm")
      )
    return(p1) # nolint
  }
  ## Plot as grid
  dist1 <- setNames(
    lapply(
      seq.int(1, length(qc_inputs), 1),
      function(h) {
        p1 <- setNames(lapply(
          seq.int(1, length(qc_inputs[[h]]), 1),
          function(k) {
            ggpubr::ggarrange(
              # intensity
              pdist(qc_inputs[[h]][[k]], "imputed", "norm", 1),
              # log2
              pdist(qc_inputs[[h]][[k]], "imputed.log2", "norm.log2", 2),
              # pareto
              pdist(qc_inputs[[h]][[k]], "imputed.pareto", "norm.pareto", 3),
              nrow = 1,
              ncol = 3
            )
          }
        ), names(qc_inputs[[h]]))
        p2 <- ggpubr::ggarrange(
          plotlist = p1, nrow = length(p1),
          ncol = 1, labels = names(p1),
          label.x = 0.05
        )
        return(p2) # nolint
      }
    ), names(qc_inputs)
  )
  outs[["data distribution"]] <- dist1
  #---- PCA and Outlier Detection ----
  print("8. Calculating PCA and flagging outlier samples from normalized assays") # nolint
  dimr1 <- setNames(
    lapply(
      seq.int(1, length(qc_inputs), 1),
      function(h) {
        p1 <- setNames(lapply(
          seq.int(1, length(qc_inputs[[h]]), 1),
          function(k) {
            # select TIC columns
            tic1 <- c(
              paste("TIC", names(qc_inputs[[h]])[[k]], "imputed", sep = "."), # nolint
              paste("TIC", names(qc_inputs[[h]])[[k]], "norm", sep = ".")
            )
            ggpubr::ggarrange( # nolint
              # imputed
              ms_dim_rd(
                qc_inputs[[h]][[k]], "imputed",
                md1[["gc"]], p_lab = FALSE,
                flag_outliers = TRUE,
                outlier_calc = tic1[[1]],
                sid = md1[["sid"]]
              ),
              # normalized
              ms_dim_rd(
                qc_inputs[[h]][[k]], "norm",
                md1[["gc"]], p_lab = FALSE,
                flag_outliers = TRUE,
                outlier_calc = tic1[[2]],
                sid = md1[["sid"]]
              ),
              nrow = 1,
              ncol = 2,
              labels = c("imputed", "norm"),
              label.x = 0.05
            )
          }
        ), paste(names(qc_inputs[[h]])))
        p2 <- ggpubr::ggarrange(
          plotlist = p1,
          nrow = length(p1), ncol = 1,
          labels = names(p1)
        )
        return(p2) # nolint
      }
    ), names(qc_inputs)
  )
  outs[["dimension reduction"]] <- dimr1
  #---- Signal Intensity Drift ----
  print("9. Determining intensity drifts across analytical runs")
  drift1 <- setNames(
    lapply(
      seq.int(1, length(qc_inputs), 1),
      function(h) {
        p1 <- setNames(lapply(
          seq.int(1, length(qc_inputs[[h]]), 1),
          function(k) {
            # select TIC columns
            tc1 <- c(
              paste("TIC", names(qc_inputs[[h]])[[k]], "imputed", sep = "."),
              paste("TIC", names(qc_inputs[[h]])[[k]], "norm", sep = ".")
            )
            ms_plot_tic(
              qc_inputs[[h]][[k]],
              var_x = md1[["ac"]],
              var_g = md1[["gc"]],
              tic1 = tc1[[1]],
              tic2 = tc1[[2]]
            )
          }
        ), names(qc_inputs[[h]]))
        p2 <- ggpubr::ggarrange(
          plotlist = p1,
          nrow = length(p1), ncol = 1,
          labels = names(p1),
          label.x = 0.05
        )
        return(p2) # nolint
      }
    ), names(qc_inputs)
  )
  outs[["signal drift"]] <- drift1
  #---- RSD and Technical Variance ----
  print("10. Calculating group and median QC RSD")
  grsd <- setNames(
    lapply(
      seq.int(1, length(qc_inputs), 1),
      function(h) {
        p1 <- setNames(lapply(
          seq.int(1, length(qc_inputs[[h]]), 1),
          function(k) {
            setNames(
              lapply(
                c("imputed", "norm"),
                function(l) {
                  rsd1 <- ms_group_rsd(
                    exp1 = qc_inputs[[h]][[k]],
                    asy = l,
                    var_g = md1[["gc"]],
                    var_feat = md1[["fc"]]
                  )
                  # calculate median per group
                  rsd2 <- setNames(
                    as.data.frame(lapply(
                      seq.int(1, ncol(rsd1), 1),
                      function(m) median(rsd1[[m]])
                    )),
                    names(rsd1)
                  )
                  # combine as list
                  rsd1 <- list(
                    "RSD" = rsd1,
                    "mRSD" = rsd2
                  )
                  return(rsd1) # nolint
                }
              ),
              c("imputed", "norm")
            )
          }
        ), names(qc_inputs[[h]]))
        return(p1) # nolint
      }
    ), names(qc_inputs)
  )
  outs[["rsd"]] <- grsd
  #---- Sample Correlation ----
  print("11. Determining correlation between study samples")
  cor1 <- setNames(
    lapply(
      seq.int(1, length(qc_inputs), 1),
      function(h) {
        p1 <- setNames(lapply(
          seq.int(1, length(qc_inputs[[h]]), 1),
          function(k) {
            setNames(
              lapply(
                c("imputed", "norm"),
                function(l) {
                  ms_samp_cor(
                    exp1 = qc_inputs[[h]][[k]],
                    asy = l,
                    var_g = md1[["gc"]],
                    sid = md1[["sid"]],
                    pname = paste(
                      names(qc_inputs)[[h]],
                      names(qc_inputs[[h]])[[k]],
                      l
                    )
                  )
                }
              ),
              c("imputed", "norm")
            )
          }
        ), names(qc_inputs[[h]]))
        return(p1) # nolint
      }
    ), names(qc_inputs)
  )
  outs[["correlation"]] <- cor1
  #---- Report Output ----
  print("QC Complete! Run ms_qc_report to generate a summary report for the selected dataset.") # nolint
  return(outs) # nolint
}

#' Normalization QC Scatter Plot
#'
#' Scatter plot visualization to evaluate normalization performance.
#'
#' @param exp1 SummarizedExperiment containing sample metadata.
#' @param var_x X-axis variable (usually sample ID).
#' @param var_g Grouping variable.
#' @param tic1 Pre-normalization TIC column.
#' @param tic2 Post-normalization TIC column.
#' @return A scatter plot visualizing overall sample intensities.
#' @examples
#'
#' # ptic <- ms_plot_tic(
#' #   exp1 = data1,
#' #   var_x = "sampleID",
#' #   var_g = "Group",
#' #   tic1 = "TIC.imputed",
#' #   tic2 = "TIC.norm"
#' # )
#'
#' @export
ms_plot_tic <- function(
  exp1,
  var_x,
  var_g = "Group",
  tic1,
  tic2
) {
  ## TIC plot
  ggplot2::ggplot(
    data = as.data.frame(colData(exp1)) # nolint
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data[[var_x]], # nolint
        y = .data[[tic1]],
        color = factor(
          .data[[var_g]],
          levels = gtools::mixedsort(unique(.data[[var_g]]))
        )
      ),
      shape = 16,
      size = 4,
      alpha = 0.5
    ) +
    ggplot2::geom_smooth(
      color = "firebrick1",
      alpha = 0.5,
      ggplot2::aes(x = .data[[var_x]], y = .data[[tic1]]),
      linetype = "dashed"
    ) +
    ggplot2::scale_color_manual(values = col_univ()) + # nolint
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data[[var_x]],
        y = .data[[tic2]],
        color = factor(
          .data[[var_g]],
          levels = gtools::mixedsort(unique(.data[[var_g]]))
        )
      ),
      shape = 18,
      size = 4,
      alpha = 0.8
    ) +
    ggplot2::geom_smooth(
      color = "dodgerblue1",
      alpha = 0.8,
      ggplot2::aes(x = .data[[var_x]], y = .data[[tic2]])
    ) +
    ggplot2::geom_label(
      x = max(as.data.frame(colData(exp1))[[var_x]]), # nolint
      y = median(as.data.frame(colData(exp1))[[tic1]]), # nolint
      label = "Pre-norm",
      color = "firebrick1"
    ) +
    ggplot2::geom_label(
      x = max(as.data.frame(colData(exp1))[[var_x]]), # nolint
      y = median(as.data.frame(colData(exp1))[[tic2]]), # nolint
      label = "Post-norm",
      color = "dodgerblue1"
    ) +
    ggplot2::labs(
      y = "TIC",
      x = "Sample ID",
      color = "Group"
    ) +
    ms_theme() # nolint
}

#' Normalization QC RSD
#'
#' Calculates relative standard deviation for each feature per
#' group within a SummarizedExperiment object.
#'
#' @param exp1 SummarizedExperiment containing sample metadata.
#' @param asy Assay for calculating RSD values.
#' @param var_g Grouping variable.
#' @param var_feat Feature name column from rowData.
#' @return Group RSD values for each feature.
#' @examples
#'
#' # ptic <- ms_group_rsd(
#' #   exp1 = data1
#' # )
#'
#' @export
ms_group_rsd <- function(
  exp1,
  asy = "imputed",
  var_g = "Group",
  var_feat = "Label"
) {
  # Load data
  d1 <- exp1
  if (Sys.info()[["sysname"]] != "Windows") {
    # Calculate values
    rsd1 <- setNames(as.data.frame(t(dplyr::bind_cols(parallel::mclapply(
      mc.cores = ceiling(parallel::detectCores() / 4),
      seq.int(1, nrow(d1), 1),
      function(i) {
        c1 <- aggregate(
          assay(d1, asy)[i, ], # nolint
          list(colData(d1)[[var_g]]), # nolint
          function(j) round((sd(j) / mean(j)) * 100, digits = 3)
        )
        c1 <- setNames(
          as.data.frame(c1[gtools::mixedorder(c1[["Group.1"]]), 2]),
          rowData(d1)[i, ][[var_feat]] # nolint
        )
        return(c1) # nolint
      }
    )))), gtools::mixedsort(unique(colData(d1)[[var_g]])))
  }
  if (Sys.info()[["sysname"]] == "Windows") {
    # Calculate values
    rsd1 <- setNames(as.data.frame(t(dplyr::bind_cols(lapply(
      seq.int(1, nrow(d1), 1),
      function(i) {
        c1 <- aggregate(
          assay(d1, asy)[i, ], # nolint
          list(colData(d1)[[var_g]]), # nolint
          function(j) round((sd(j) / mean(j)) * 100, digits = 3)
        )
        c1 <- setNames(
          as.data.frame(c1[gtools::mixedorder(c1[["Group.1"]]), 2]),
          rowData(d1)[i, ][[var_feat]] # nolint
        )
        return(c1) # nolint
      }
    )))), gtools::mixedsort(unique(colData(d1)[[var_g]])))
  }
  return(rsd1) # nolint
}

#' Normalization QC Sample Correlation
#'
#' Calculates correlations between study samples
#' within a SummarizedExperiment object.
#'
#' @param exp1 SummarizedExperiment containing sample metadata.
#' @param asy Assay for calculating RSD values.
#' @param var_g Grouping variable.
#' @param sid Sample ID column name.
#' @param hh Heatmap height/width.
#' @param fsr Row/column font size.
#' @param pname Plot name.
#' @return Sample correlation heatmap.
#' @examples
#'
#' # ms_samp_cor(
#' #   exp1 = data1
#' # )
#'
#' @export
ms_samp_cor <- function(
  exp1,
  asy = "norm",
  var_g = "Group",
  sid = "Label",
  hh = 18,
  fsr = 6,
  pname = NULL
) {
  # Load data
  h1 <- exp1
  # Calculate correlation matrix
  h2 <- cor(assay(h1, asy), method = "spearman") # nolint
  # Add row/column names
  rownames(h2) <- colData(h1)[[sid]] # nolint
  colnames(h2) <- colData(h1)[[sid]] # nolint
  # Heatmap colors
  fun_hm_col <- circlize::colorRamp2(
    c(-1, - 0.5, 0, 0.5, 1),
    colors = col_grad(scm = 7) # nolint
  )
  # Heatmap annotation colors
  fun_hm_bar <- list(
    cl_var = setNames(
      col_univ()[1:length(unique(colData(h1)[[var_g]]))], # nolint
      gtools::mixedsort(unique(colData(h1)[[var_g]]))
    )
  )
  # Annotations
  hm_anno_list <- list(
    hm_col <- ComplexHeatmap::HeatmapAnnotation( # nolint
      `Group` = colData(h1)[[var_g]], # nolint
      col = fun_hm_bar,
      show_annotation_name = FALSE,
      show_legend = TRUE
    )
  )
  # Plot
  h3 <- ComplexHeatmap::Heatmap(
    h2,
    col = fun_hm_col,
    name = "Correlation",
    top_annotation = hm_anno_list[[1]],
    show_column_names = TRUE,
    show_row_names = TRUE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    heatmap_width = ggplot2::unit(hh, "cm"),
    heatmap_height = ggplot2::unit(hh, "cm"),
    column_title = pname,
    row_names_gp = grid::gpar(fontsize = fsr),
    column_names_gp = grid::gpar(fontsize = fsr)
  )
  return(h3) # nolint
}

#' QC Check Group Presence
#'
#' Checks compound presence per treatment group for the selected dataset
#' and flags compounds that are not present in at least 50% of samples
#' in one group.
#'
#' @param exp1 A SummarizedExperiment.
#' @param asy Assay name to calculate group presence/absence.
#' @param var_g Group variable name.
#' @param var_feat Feature name column from rowData.
#' @param hw Presence/absence heatmap width.
#' @param hh Presence/absence heatmap height.
#' @param fsc Heatmap column names fontsize.
#' @param fsr Heatmap row names fontsize.
#' @param clc Cluster heatmap columns?
#' @param clr Cluster heatmap rows?
#' @param showcols Show heatmap column names?
#' @param showrows Show heatmap row names?
#' @param trans Transpose heatmap?
#' @param title1 Heatmap title.
#' @return Heatmap displaying compound presence/absence.
#' @examples
#'
#' # ms_check_pres(data1)
#'
#'
#' @export
ms_check_pres <- function(
  exp1,
  asy = "raw",
  var_g = "Group",
  var_feat = "Label",
  hw = 20,
  hh = 20,
  fsc = 8,
  fsr = 8,
  clc = TRUE,
  clr = TRUE,
  showcols = TRUE,
  showrows = TRUE,
  trans = TRUE,
  title1 = NULL
) {
  # input data
  ld1 <- exp1
  # Calculate proportion present per group
  # Calculate pre-norm values
  dchk <- setNames(as.data.frame(t(dplyr::bind_cols(lapply(
    seq.int(1, nrow(ld1), 1),
    function(i) {
      # Calculate presence
      c1 <- setNames(aggregate(
        assay(ld1, asy)[i, ], # nolint
        list(colData(ld1)[[var_g]]), # nolint
        function(j) {
          d1 <- j > 0 & !is.na(j) # nolint
          d1 <- length(d1[d1 == TRUE])
        }
      ), c("Group", "n.pres"))
      c1 <- c1[gtools::mixedorder(c1[["Group"]]), ]
      # Calculate group totals
      c2 <- setNames(
        dplyr::count(
          as.data.frame(colData(ld1)[[var_g]]), colData(ld1)[[var_g]] # nolint
        ), c("Group", "n.tot")
      )
      c2 <- c2[gtools::mixedorder(c2[["Group"]]), ]
      # Return proportion per group
      c1 <- setNames(
        as.data.frame(round(c1[["n.pres"]] / c2[["n.tot"]], digits = 3)),
        rowData(ld1)[i, ][[var_feat]] # nolint
      )
      return(c1) # nolint
    }
  )))), gtools::mixedsort(unique(colData(ld1)[[var_g]])))
  # Heatmap displaying proportions for each compound
  h1 <- as.matrix(dchk)
  # Heatmap colors
  fun_hm_col <- circlize::colorRamp2(
    c(0, 0.5, 1),
    colors = col_grad(scm = 4) # nolint
  )
  if (trans == FALSE) {
    h_out <- ComplexHeatmap::Heatmap(
      h1,
      col = fun_hm_col,
      name = "Proportion",
      show_column_names = showcols,
      show_row_names = showrows,
      cluster_columns = clc,
      cluster_rows = clr,
      heatmap_width = ggplot2::unit(hw, "cm"),
      heatmap_height = ggplot2::unit(hh, "cm"),
      column_title = paste(title1),
      column_names_gp = grid::gpar(fontsize = fsc),
      row_names_gp = grid::gpar(fontsize = fsr)
    )
  }
  if (trans == TRUE) {
    h_out <- ComplexHeatmap::Heatmap(
      t(h1),
      col = fun_hm_col,
      name = "Proportion",
      show_column_names = showcols,
      show_row_names = showrows,
      cluster_columns = clc,
      cluster_rows = clr,
      heatmap_width = ggplot2::unit(hw, "cm"),
      heatmap_height = ggplot2::unit(hh, "cm"),
      column_title = paste(title1),
      column_names_gp = grid::gpar(fontsize = fsc),
      row_names_gp = grid::gpar(fontsize = fsr)
    )
  }
  # return vector of compounds with less than 50% presence
  # in at least one group
  ### replace NA with 0
  h1[is.na(h1)] <- 0
  ### flag missing or low presence compounds
  h2 <- unlist(apply(h1, 1, function(y) max(y) < 0.5))
  h2 <- rownames(h1[h2, ])
  print("The following compounds do not exceed 50% presence in at least one group:") # nolint
  print(h2)
  print("Consider removing these compounds from the dataset")
  return(h_out)
}
