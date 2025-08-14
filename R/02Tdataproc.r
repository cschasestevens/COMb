#' Chromatogram Viewer (Targeted)
#'
#' Plots chromatograms for selected samples and compounds from
#' input x and y coordinates.
#'
#' @param dat Input data list.
#' @param intype Input data type: either "Skyline," a list of 2 vectors
#' containing values separated by commas (default) or "xy," which is
#' a list of two variables containing vectors of x and y coordinates.
#' @param comp Compound name.
#' @param sid Sample ID.
#' @param sidcol Sample ID column name.
#' @param xcol X column name.
#' @param ycol Y column name.
#' @param rtcol Retention time column containing the specific retention
#' times of identified peaks.
#' @param lab Column to use for labeling samples.
#' @param pkcol Column generated from tms_peakdetect providing the
#' upper and lower RT limits of the detected peak.
#' @param pkview Highlight detected peaks if present in input data.
#' Requires result from tms_peakdetect.
#' @return A chromatogram based on X and Y coordinates.
#' @examples
#'
#' # tms_chromatogram <- function(
#' #   dat = pks,
#' #   comp = "15-HETE",
#' #   sid = "Sample1"
#' # )
#'
#' @import ggplot2
#' @import stringr
#' @export
tms_chromatogram <- function( # nolint
  dat,
  intype = "Skyline",
  comp,
  sid,
  sidcol = "SampleID",
  xcol = "RT_all",
  ycol = "Intensity",
  rtcol = "pkRT",
  lab = "SampleID",
  pkview = TRUE,
  pkcol = "pkBoundary"
) {
  # Load data
  ld <- dat
  # Subset data for plot input
  ld2 <- ld[
    ld[["Name"]] == comp &
      ld[[sidcol]] == sid,
  ]
  if(intype == "Skyline") { # nolint
    ld3 <- data.frame(
      "X" = as.numeric(unlist(stringr::str_split(ld2[[xcol]], ","))),
      "Y" = as.numeric(unlist(stringr::str_split(ld2[[ycol]], ",")))
    )
    if(pkview == TRUE) { # nolint
      ld4 <- ld3[
        ld3[["X"]] >= as.numeric(
          unlist(stringr::str_split(ld2[[pkcol]], ","))
        )[[1]] &
          ld3[["X"]] <= as.numeric(
            unlist(stringr::str_split(ld2[[pkcol]], ","))
          )[[3]],
      ]
      # Create plot
      p1 <- ggplot2::ggplot(data = ld3) +
        ggplot2::geom_area(
          data = ld4,
          ggplot2::aes(
            x = ld4[["X"]],
            y = ld4[["Y"]]
          ),
          color = "grey50",
          alpha = 0.5
        ) +
        ggplot2::geom_line(
          ggplot2::aes(
            x = .data[["X"]], # nolint
            y = .data[["Y"]]
          )
        ) +
        ggplot2::geom_vline(
          xintercept = ld2[[rtcol]],
          linetype = "dashed"
        ) +
        ggplot2::geom_segment(
          data = ld4,
          ggplot2::aes(
            x = min(ld4[["X"]]),
            xend = min(ld4[["X"]]),
            y = 0,
            yend = ld4[ld4[["X"]] == min(ld4[["X"]]), "Y"]
          ),
          linetype = "dotdash",
          color = "grey25"
        ) +
        ggplot2::geom_segment(
          data = ld4,
          ggplot2::aes(
            x = max(ld4[["X"]]),
            xend = max(ld4[["X"]]),
            y = 0,
            yend = ld4[ld4[["X"]] == max(ld4[["X"]]), "Y"]
          ),
          linetype = "dotdash",
          color = "grey25"
        ) +
        ggplot2::scale_color_manual(values = col_univ()) +
        ms_theme() +
        ggplot2::labs(
          x = "Retention Time (min.)",
          y = "Intensity",
          title = paste(ld2[["Name"]], ld2[[sidcol]], sep = ": ")
        )
    }
    if(pkview == FALSE) { # nolint
      # Create plot
      p1 <- ggplot2::ggplot(data = ld3) +
        ggplot2::geom_line(
          ggplot2::aes(
            x = .data[["X"]], # nolint
            y = .data[["Y"]]
          )
        ) +
        ggplot2::geom_vline(
          xintercept = ld2[[rtcol]],
          linetype = "dashed"
        ) +
        ggplot2::scale_color_manual(values = col_univ()) +
        ms_theme() +
        ggplot2::labs(
          x = "Retention Time (min.)",
          y = "Intensity",
          title = paste(ld2[["Name"]], ld2[[sidcol]], sep = ": ")
        )
    }
  }
  return(p1)
}

#' Peak Integration (Targeted)
#'
#' Integrates detected peaks from a tms_peakdetect result.
#'
#' @param dat Input data list.
#' @param intype Input data type: either "Skyline," a list of 2 vectors
#' containing values separated by commas (default) or "xy," which is
#' a list of two variables containing vectors of x and y coordinates.
#' @param xcol X column name.
#' @param ycol Y column name.
#' @param comp Compound name to plot.
#' @param sid Name of a specific sample for peak integration.
#' @param sidcol Sample ID column.
#' @param pkcol Peak boundary column, generated from tms_peakdetect.
#' @param phcol Peak height column.
#' @param snr_type Summary statistic for calculating signal to noise
#' ratio (either "median" or "mean").
#' @param snr Signal to noise ratio for filtering peaks.
#' Use 3 for LOD and 10 for LOQ.
#' @return Integrated raw and normalized peak areas and heights.
#' @examples
#'
#' # d1 <- tms_integrate(
#' #   dat = pks,
#' #   comp = "15-HETE",
#' #   sid = "Stevens_879654_1_cal_8"
#' # )
#'
#' @import dplyr
#' @import stringr
#' @import DescTools
#' @export
tms_integrate <- function( # nolint
  dat,
  comp,
  sid,
  sidcol = "SampleID",
  intype = "Skyline",
  xcol = "RT_all",
  ycol = "Intensity",
  pkcol = "pkBoundary",
  phcol = "pkHeight",
  snr_type = "median",
  snr = 3
) {
  # Load data
  ld <- dat
  # Subset data
  ld2 <- ld[
    ld[["Name"]] == comp &
      ld[[sidcol]] == sid,
  ]
  # Integrate detected peak
  if(intype == "Skyline") { # nolint
    # Convert data
    ld3 <- data.frame(
      "X" = as.numeric(unlist(stringr::str_split(ld2[[xcol]], ","))),
      "Y" = as.numeric(unlist(stringr::str_split(ld2[[ycol]], ",")))
    )
    ld4 <- ld3[
      ld3[["X"]] >= as.numeric(
        unlist(stringr::str_split(ld2[[pkcol]], ","))
      )[[1]] &
        ld3[["X"]] <= as.numeric(
          unlist(stringr::str_split(ld2[[pkcol]], ","))
        )[[3]],
    ]
    # Determine signal/noise threshold
    if(snr_type == "median") {ld3_thresh <- median(ld3[["Y"]]) * snr} # nolint
    if(snr_type == "mean") {ld3_thresh <- mean(ld3[["Y"]]) * snr} # nolint
    # Calculate raw and background-subtracted AUC and PH
    # Calculate AUC for detected peak
    if(ld3_thresh > max(ld4[["Y"]])) { # nolint
      print("SNR < 3 for selected peak.")
      print("Returning NA for AUC calculation...")
      auc <- data.frame(
        "Name" = comp,
        "SampleID" = sid,
        "AUC_raw" = NA,
        "AUC_cor" = NA,
        "PH_raw" = NA,
        "PH_cor" = NA
      )
    }
    if(ld3_thresh < max(ld4[["Y"]])) { # nolint
      ## AUC - background
      auc <- data.frame(
        "Name" = comp,
        "SampleID" = sid,
        "AUC_raw" = DescTools::AUC(
          x = ld4[["X"]],
          y = ld4[["Y"]],
          method = "spline"
        ),
        "AUC_cor" = DescTools::AUC( # total area
          x = ld4[["X"]],
          y = ld4[["Y"]],
          method = "spline"
        ) - DescTools::AUC( # SNR area
          x = ld4[["X"]],
          y = rep(ld3_thresh, nrow(ld4)),
          method = "spline"
        ),
        "PH_raw" = ld2[[phcol]],
        "PH_cor" = ld2[[phcol]] - ld3_thresh
      )
    }
  }
  # Return the appended data frame with raw and normalized values
  out <- dplyr::left_join(
    ld2,
    auc,
    by = c("Name", sidcol)
  )
  return(out)
}

#' Retention Time Shift (Targeted)
#'
#' Calculates retention time shifts of standards
#' used for data acquisition.
#'
#' @param dat Input standard list.
#' @param mzrt Input mzrt list containing reference precursor m/z and
#' retention time columns.
#' @param intype Input data type: either "Skyline," a list of 2 vectors
#' containing values separated by commas (default) or "xy," which is
#' a list of two variables containing vectors of x and y coordinates.
#' @param xcol X column name.
#' @param ycol Y column name.
#' @param comp Compound name to plot.
#' @param rttol Retention time tolerance window (in minutes).
#' @param spn Number of points used for peak detection; must be an odd number.
#' @param snr_type Summary statistic to calculate signal-to-noise ratio
#' (either "median" [default] or "mean").
#' @param lab Column containing sample labels.
#' @param sampid Name of a specific sample for peak integration.
#' @param mcc Cores to use if not using Windows.
#' @return A data frame containing the calculated difference between
#' the expected and measured retention times for a list of reference
#' standards.
#' @examples
#'
#' # msT_integrate <- function(
#' #   dat = d1[["merged"]][d1[["merged"]][["Group"]] == "Std_curve", ],
#' #   comp = "15-HETE",
#' #   ref = "Std_curve",
#' #   grp = "Std_curve",
#' #   sampid = "Stevens_879654_1_cal_8"
#' # )
#'
#' @import ggplot2
#' @import dplyr
#' @import stringr
#' @import ggpmisc
#' @import DescTools
#' @export
tms_rtshift <- function( # nolint
  dat,
  mzrt,
  intype = "Skyline",
  xcol = "RT_all",
  ycol = "Intensity",
  mcc = 8,
  comp,
  rttol = 0.2,
  spn = 7,
  lab = "SampleID",
  sampid,
  snr_type = "median"
) {
  # Load data and set parameters
  ld <- dat
  xc <- xcol
  yc <- ycol
  lb <- lab
  cmp <- comp
  # Detect peaks from Skyline datasets
  if(intype == "Skyline") { # nolint
    # Convert RT and intensities to X and Y coordinates
    if(Sys.info()[["sysname"]] != "Windows") { # nolint
      ## Extract standards from data
      ld1 <- ld[grepl("standard|STD", ld[["synthesis_pathway"]]), ]
      ## Format
      ld1 <- dplyr::bind_rows(parallel::mclapply(
        mc.cores = mcc,
        seq.int(1, nrow(ld), 1),
        function(x) {
          ld2 <- data.frame(
            "X" = as.numeric(
              unlist(stringr::str_split(ld[x, xc], ","))
            ),
            "Y" = as.numeric(
              unlist(stringr::str_split(ld[x, yc], ","))
            ),
            "SampleID" = ld[x, lb],
            "Group" = ld[x, "Group"],
            "Name" = ld[x, "Name"]
          )
          return(ld2) # nolint
        }
      ))
    }
    if(Sys.info()[["sysname"]] == "Windows") { # nolint
      ## Extract standards from data
      ld1 <- ld[grepl("standard|STD", ld[["synthesis_pathway"]]), ]
      ## Format
      ld1 <- dplyr::bind_rows(lapply(
        seq.int(1, nrow(ld), 1),
        function(x) {
          ld2 <- data.frame(
            "X" = as.numeric(
              unlist(stringr::str_split(ld[x, xc], ","))
            ),
            "Y" = as.numeric(
              unlist(stringr::str_split(ld[x, yc], ","))
            ),
            "SampleID" = ld[x, lb],
            "Group" = ld[x, "Group"],
            "Name" = ld[x, "Name"]
          )
          return(ld2) # nolint
        }
      ))
    }
    # Select metadata
    ld2 <- dplyr::select(
      ld[grepl("standard|STD", ld[["synthesis_pathway"]]), ],
      -c(xc, yc)
    )
    # Merge meta with reference list
    ld2 <- dplyr::select(dplyr::left_join(
      ld2,
      mzrt,
      by = c("Name", "synthesis_pathway")
    ), c("ref_mz", "ref_RT"), everything()) # nolint
    # Convert sample IDs to factor
    ld1[["SampleID"]] <- factor(
      ld1[["SampleID"]],
      levels = unique(ld1[["SampleID"]])
    )
    ld2[["SampleID"]] <- factor(
      ld2[["SampleID"]],
      levels = unique(ld2[["SampleID"]])
    )
    # Subset data for plot input
    ## Select compound and sample group
    ld3 <- ld1[
      ld1[["Name"]] == cmp &
        ld1[["SampleID"]] == sampid,
    ]
    ## Subset metadata for plotting reference RT
    ld3c <- ld2[
      ld2[["Name"]] == cmp &
        ld2[["SampleID"]] == sampid,
    ]
    # Determine RT shift between reference and measured RT
    ## Peak detection
    ### Local maximum within reference RT window
    #### Subset data
    ld3_peak <- ld3
    ld3_peak <- ld3_peak[
      ld3_peak[["X"]] > (unique(ld3c[["ref_RT"]]) - rttol) &
        ld3_peak[["X"]] < (unique(ld3c[["ref_RT"]]) + rttol),
    ]
    #### set global threshold to find peaks > 3x signal-to-noise
    if(snr_type == "median") {ld3_threshold <- median(ld3_peak[["Y"]]) * 3} # nolint
    if(snr_type == "mean") {ld3_threshold <- mean(ld3_peak[["Y"]]) * 3} # nolint
    ld3_threshold <- ld3_threshold / max(ld3_peak[["Y"]])
    if(ld3_threshold < -1 || ld3_threshold > 1) { # nolint
      ld3_max <- ld3_peak[
        ggpmisc::find_peaks(
          ld3_peak[["Y"]],
          span = spn
        ),
      ]
    }
    if(ld3_threshold > -1 && ld3_threshold < 1) { # nolint
      ld3_max <- ld3_peak[
        ggpmisc::find_peaks(
          ld3_peak[["Y"]],
          span = spn,
          global.threshold = ld3_threshold
        ),
      ]
    }
    if(nrow(ld3_max) == 1) { # nolint
      ld3_max <- dplyr::mutate(
        ld3_max,
        "Xdiff" = abs(ld3_max[["X"]] - unique(ld3c[["ref_RT"]]))
      )
    }
    if(nrow(ld3_max) > 1) { # nolint
      ld3_max <- dplyr::mutate(
        ld3_max,
        "Xdiff" = abs(ld3_max[["X"]] - unique(ld3c[["ref_RT"]]))
      )
      # Select maximum peak that is closest to reference
      ld3_max <- ld3_max[ld3_max[["Xdiff"]] == min(ld3_max[["Xdiff"]]), ]
    }
    ### Local minimums relative to detected peak
    #### Subset
    ld3_peak <- ld3
    ld3_min_1 <- ld3_peak[
      ld3_peak[["X"]] < ld3_max[["X"]] &
        ld3_peak[["X"]] > ld3_max[["X"]] - rttol,
    ]
    ld3_min_2 <- ld3_peak[
      ld3_peak[["X"]] > ld3_max[["X"]] &
        ld3_peak[["X"]] < ld3_max[["X"]] + rttol,
    ]
    #### trailing tail
    min1 <- ggpmisc::find_valleys(ld3_min_1[["Y"]], span = spn)
    if(length(min1[min1 == TRUE]) == 0) { # nolint
      ld3_min_1 <- ld3_min_1[ld3_min_1[["Y"]] == min(ld3_min_1[["Y"]]), ]
    }
    if(length(min1[min1 == TRUE]) > 0) { # nolint
      ld3_min_1 <- ld3_min_1[min1, ]
    }
    #### leading tail
    min2 <- ggpmisc::find_valleys(ld3_min_2[["Y"]], span = spn)
    if(length(min2[min2 == TRUE]) == 0) { # nolint
      ld3_min_2 <- ld3_min_2[ld3_min_2[["Y"]] == min(ld3_min_2[["Y"]]), ]
    }
    if(length(min2[min2 == TRUE]) > 0) { # nolint
      ld3_min_2 <- ld3_min_2[min2, ]
    }
    #### filter if more than one value present
    if(nrow(ld3_min_1) == 1) { # nolint
      ld3_min_1 <- dplyr::mutate(
        ld3_min_1,
        "Xdiff" = abs(ld3_min_1[["X"]] - ld3_max[["X"]])
      )
      ld3_min_1 <- ld3_min_1[ld3_min_1[["X"]] < ld3_max[["X"]], ]
    }
    if(nrow(ld3_min_1) > 1) { # nolint
      ld3_min_1 <- dplyr::mutate(
        ld3_min_1,
        "Xdiff" = abs(ld3_min_1[["X"]] - ld3_max[["X"]])
      )
      ld3_min_1 <- ld3_min_1[ld3_min_1[["X"]] < ld3_max[["X"]], ]
      ld3_min_1 <- ld3_min_1[
        ld3_min_1[["Xdiff"]] == min(ld3_min_1[["Xdiff"]]),
      ]
    }
    if(nrow(ld3_min_2) == 1) { # nolint
      ld3_min_2 <- dplyr::mutate(
        ld3_min_2,
        "Xdiff" = abs(ld3_min_2[["X"]] - ld3_max[["X"]])
      )
      ld3_min_2 <- ld3_min_2[ld3_min_2[["X"]] > ld3_max[["X"]], ]
    }
    if(nrow(ld3_min_2) > 1) { # nolint
      ld3_min_2 <- dplyr::mutate(
        ld3_min_2,
        "Xdiff" = abs(ld3_min_2[["X"]] - ld3_max[["X"]])
      )
      ld3_min_2 <- ld3_min_2[ld3_min_2[["X"]] > ld3_max[["X"]], ]
      ld3_min_2 <- ld3_min_2[
        ld3_min_2[["Xdiff"]] == min(ld3_min_2[["Xdiff"]]),
      ]
    }
    ### Combine peak points
    ld3_peak <- dplyr::bind_rows(
      ld3_min_1,
      ld3_max,
      ld3_min_2
    )
    ### Calculate RT shift relative to reference
    ld_shift <- round(
      ld3_max[["X"]] - unique(ld3c[["ref_RT"]]),
      digits = 2
    )
    ld_shift2 <- data.frame(
      "x" = c(unique(ld3c[["ref_RT"]]), ld3_max[["X"]]),
      "y" = rep(ld3_max[["Y"]], 2)
    )
    ## Peak validation
    p1 <- ggplot2::ggplot(
      data = ld3[ld3[["SampleID"]] == sampid, ],
      ggplot2::aes(
        x = .data[["X"]], # nolint
        y = .data[["Y"]],
        color = .data[["SampleID"]]
      )
    ) +
      ggplot2::geom_area(
        data = ld3[
          ld3[["X"]] >= min(ld3_peak[["X"]]) &
            ld3[["X"]] <= max(ld3_peak[["X"]]) &
            ld3[["SampleID"]] == sampid,
        ],
        ggplot2::aes(
          x = ld3[
            ld3[["X"]] >= min(ld3_peak[["X"]]) &
              ld3[["X"]] <= max(ld3_peak[["X"]]) &
              ld3[["SampleID"]] == sampid,
            "X"
          ],
          y = ld3[
            ld3[["X"]] >= min(ld3_peak[["X"]]) &
              ld3[["X"]] <= max(ld3_peak[["X"]]) &
              ld3[["SampleID"]] == sampid,
            "Y"
          ]
        ),
        color = "grey50",
        alpha = 0.5
      ) +
      ggplot2::geom_line(show.legend = FALSE) +
      # Reference RT
      ggplot2::geom_vline(
        xintercept = unique(ld3c[["ref_RT"]]),
        linetype = "dashed",
        color = "grey25"
      ) +
      # Detected peak RT
      ggplot2::geom_vline(
        xintercept = ld3_max[["X"]],
        linetype = "dashed",
        color = "firebrick"
      ) +
      # arrow showing RT shift
      ggplot2::geom_segment(
        data = ld_shift2,
        ggplot2::aes(
          x = ld_shift2[["x"]][[1]],
          xend = ld_shift2[["x"]][[2]],
          y = ld_shift2[["y"]][[1]],
          yend = ld_shift2[["y"]][[2]]
        ),
        linetype = "solid",
        color = "black",
        arrow = grid::arrow(type = "closed")
      ) +
      # label for RT shift
      ggplot2::annotate(
        "label",
        x = c(
          ld3_max[["X"]] - ((ld3_max[["X"]] - unique(ld3c[["ref_RT"]])) / 2)
        ),
        y = c(ld3_max[["Y"]]),
        label = ifelse(
          ld3_max[["X"]] > unique(ld3c[["ref_RT"]]),
          paste("+", round(ld3_max[["Xdiff"]], digits = 2), " min.", sep = ""),
          paste("-", round(ld3_max[["Xdiff"]], digits = 2), " min.", sep = "")
        ),
        color = "grey10"
      ) +
      # lower peak area boundary
      ggplot2::geom_segment(
        data = ld3_peak,
        ggplot2::aes(
          x = min(ld3_peak[["X"]]),
          xend = min(ld3_peak[["X"]]),
          y = 0,
          yend = ld3_peak[ld3_peak[["X"]] == min(ld3_peak[["X"]]), "Y"]
        ),
        linetype = "dotdash",
        color = "grey25"
      ) +
      # upper peak area boundary
      ggplot2::geom_segment(
        data = ld3_peak,
        ggplot2::aes(
          x = max(ld3_peak[["X"]]),
          xend = max(ld3_peak[["X"]]),
          y = 0,
          yend = ld3_peak[ld3_peak[["X"]] == max(ld3_peak[["X"]]), "Y"]
        ),
        linetype = "dotdash",
        color = "grey25"
      ) +
      ggplot2::scale_color_manual(values = col_univ()) +
      ms_theme() +
      ggplot2::labs(
        x = "Retention Time (min.)",
        y = "Intensity",
        title = sampid
      )
    # Return plot and calculated RT shift for correcting peak detection
    ld3_max[["RTshift"]] <- ld_shift
    out <- list("peak" = p1, "RTshift" = ld3_max)
  }
  return(out)
}

#' Retention Time Shift Summary (Targeted)
#'
#' Mzrt scatterplot of internal standards highlighting
#' shift in expected vs. measured retention times.
#'
#' @param dat Input standard list.
#' @param ref Reference list containing retention times.
#' @param colcomp Variable containing standard names.
#' @param colsid Variable containing sample IDs.
#' @param colrtref Variable containing reference standard retention times.
#' @param colrtexp Variable containing measured retention times.
#' @return mzrt scatter plot visualizing differences between
#' the expected and measured retention times of all reference
#' standards.
#' @examples
#'
#' # d1 <- tms_mzrt(rt1, rt2)
#'
#' @import ggplot2
#' @export
tms_mzrt <- function( # nolint
  dat,
  ref,
  colcomp = "Name",
  colsid = "SampleID",
  colrtref = "ref_RT",
  colrtexp = "X"
) {
  # Load data and set parameters
  ld <- dat
  ld2 <- ref
  ## mzrt plot
  p1 <- ggplot2::ggplot(
    data = ld
  ) +
    ggplot2::geom_point(
      data = ld2,
      ggplot2::aes(
        x = ld2[["Name"]],
        y = ld2[["ref_RT"]]
      ),
      shape = 18,
      size = 2,
      color = "red",
      alpha = 0.9
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data[["Name"]], # nolint
        y = .data[["X"]],
        color = .data[["SampleID"]]
      ),
      shape = 16,
      size = 2,
      alpha = 0.5
    ) +
    ggplot2::scale_color_manual(values = col_univ()) +
    ms_theme() +
    ggplot2::labs(
      x = "Compound Name",
      y = "Retention Time (min.)"
    )
  return(p1) # nolint
}

#' Retention Time Correction (Targeted)
#'
#' Applies retention time corrections based on observed
#' differences between expected vs. measured retention times
#' of reference standards in a targeted assay.
#'
#' @param dat Input standard list.
#' @param ref Reference list containing retention times.
#' @param colrt Retention time column.
#' @param colshift Variable containing measured retention time shifts.
#' @param type Regression type for predicting retention times ("lm" or "poly").
#' @return mzrt scatter plot visualizing differences between
#' the expected and measured retention times of all reference
#' standards.
#' @examples
#'
#' # d1 <- tms_rtcor(ref = rtc, dat = mzlist)
#'
#' @import ggplot2
#' @export
tms_rtcor <- function(
  dat,
  ref,
  colrt = "ref_RT",
  colshift = "med_shift",
  type = "lm"
) {
  # load data
  ## reference standards
  ld1 <- ref
  ld1[["exp_RT"]] <- ld1[[colrt]] + ld1[["med_shift"]]
  ## all compounds
  ld2 <- dat
  # Use reference standards for training regression model
  if(type == "lm") {train <- lm(data = ld1, exp_RT ~ ref_RT)} # nolint
  ## Predict retention times based on model
  ld2[["exp_RT"]] <- ld2[[colrt]] * (train[[1]][[2]]) +
    train[[1]][[1]]
  ## Plot predictions
  p1 <- ggplot2::ggplot(
    data = ld1
  ) +
    ggplot2::geom_line(
      ggplot2::aes(
        x = .data[[colrt]], # nolint
        y = .data[[colrt]]
      ),
      linetype = "dashed"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        x = .data[[colrt]], # nolint
        y = .data[["exp_RT"]]
      )
    ) +
    ggplot2::geom_point(
      data = ld2,
      ggplot2::aes(
        x = ld2[[colrt]],
        y = ld2[[colrt]]
      ),
      shape = 16,
      color = col_univ()[[1]],
      size = 2
    ) +
    ggplot2::geom_point(
      data = ld2,
      ggplot2::aes(
        x = ld2[[colrt]],
        y = ld2[["exp_RT"]]
      ),
      shape = 16,
      color = col_univ()[[2]],
      size = 2
    ) +
    ms_theme() +
    ggplot2::labs(
      x = "Reference RT",
      y = "Predicted RT"
    )
  # return data frame with expected retention times for each compound
  out <- list("data" = ld2, "plot" = p1)
  return(out)
}

#' Peak Detection (Targeted)
#'
#' Calculates local minima and maxima within a specified
#' retention time window to determine peaks
#'
#' @param dat Input standard list.
#' @param ref Input mzrt list containing reference precursor m/z and
#' retention time columns.
#' @param intype Input data type: either "Skyline," a list of 2 vectors
#' containing values separated by commas (default) or "xy," which is
#' a list of two variables containing vectors of x and y coordinates.
#' @param xcol X column name.
#' @param ycol Y column name.
#' @param comp Compound name to plot.
#' @param lab Column containing sample labels.
#' @param complab Column containing compound names.
#' @param rttol Retention time tolerance window (in minutes).
#' @param rtcol Retention time column.
#' @param snr_type Summary statistic to calculate signal-to-noise ratio
#' (either "median" [default] or "mean").
#' @param snr Signal to noise ratio for filtering peaks.
#' Use 3 for LOD and 10 for LOQ.
#' @param spn Number of points used for peak detection; must be an odd number.
#' @param mcc Cores to use if not using Windows.
#' @return mzrt scatter plot visualizing differences between
#' the expected and measured retention times of all reference
#' standards.
#' @examples
#'
#' # d1 <- tms_peakdetect(
#' #   dat = d1[["merged"]],
#' #   ref = rtcor[["data"]],
#' #   comp = "15-HETE"
#' # )
#'
#' @import dplyr
#' @import stringr
#' @import ggpmisc
#' @import parallel
#' @export
tms_peakdetect <- function(
  dat,
  ref,
  intype = "Skyline",
  xcol = "RT_all",
  ycol = "Intensity",
  mcc = 8,
  comp,
  rttol = 0.2,
  spn = 7,
  lab = "SampleID",
  complab = "Name",
  rtcol = "exp_RT",
  snr_type = "median",
  snr = 3
) {
  # Load data and set parameters
  ld <- dat
  mzrt <- ref
  xc <- xcol
  yc <- ycol
  lb <- lab
  cmp <- comp
  # Detect peaks from Skyline datasets
  if(intype == "Skyline") { # nolint
    # Convert RT and intensities to X and Y coordinates
    ## Filter data
    ld1 <- ld[ld[["Name"]] == cmp, ]
    ## Format
    ld1 <- dplyr::bind_rows(lapply(
      seq.int(1, nrow(ld1), 1),
      function(x) {
        ld2 <- data.frame(
          "X" = as.numeric(
            unlist(stringr::str_split(ld1[x, xc], ","))
          ),
          "Y" = as.numeric(
            unlist(stringr::str_split(ld1[x, yc], ","))
          ),
          "SampleID" = ld1[x, lb],
          "Name" = ld1[x, complab]
        )
        return(ld2) # nolint
      }
    ))
    # Select metadata
    ld2 <- dplyr::select(
      ld,
      -c(xc, yc)
    )
    # Merge meta with reference list
    ld2 <- dplyr::select(dplyr::left_join(
      ld2,
      mzrt,
      by = c("Name", "synthesis_pathway")
    ), c("ref_mz", rtcol), everything()) # nolint
    # Convert sample IDs to factor
    ld1[["SampleID"]] <- factor(
      ld1[["SampleID"]],
      levels = unique(ld1[["SampleID"]])
    )
    ld2[["SampleID"]] <- factor(
      ld2[["SampleID"]],
      levels = unique(ld2[["SampleID"]])
    )
    ld2 <- ld2[ld2[["Name"]] == cmp, ]
    # Detect peaks for all samples in data
    s1 <- unique(ld1[["SampleID"]])
    if(Sys.info()[["sysname"]] != "Windows") { # nolint
      pks <- dplyr::bind_rows(
        parallel::mclapply(
          mc.cores = mcc,
          seq.int(1, length(s1), 1),
          function(i) {
            tryCatch(
              {
                ## Select sample
                ld3 <- ld1[ld1[["SampleID"]] == s1[[i]], ]
                ## Select metadata
                ld3c <- ld2[ld2[["SampleID"]] == s1[[i]], ]
                # Peak detection
                ## Local maximum within reference RT window
                ### Subset
                ld3_peak <- ld3
                ld3_peak <- ld3_peak[
                  ld3_peak[["X"]] > (unique(ld3c[["exp_RT"]]) - rttol) &
                    ld3_peak[["X"]] < (unique(ld3c[["exp_RT"]]) + rttol),
                ]
                ### set global threshold to find peaks > signal-to-noise ratio
                ### use snr = 3 for LOD and snr = 10 for LOQ
                if(snr_type == "median") {ld3_threshold <- median(ld3_peak[["Y"]]) * snr} # nolint
                if(snr_type == "mean") {ld3_threshold <- mean(ld3_peak[["Y"]]) * snr} # nolint
                ld3_threshold <- ld3_threshold / max(ld3_peak[["Y"]])
                if(ld3_threshold < -1 || ld3_threshold > 1) { # nolint
                  ld3_max <- ld3_peak[
                    ggpmisc::find_peaks(
                      ld3_peak[["Y"]],
                      span = spn
                    ),
                  ]
                }
                if(ld3_threshold > -1 && ld3_threshold < 1) { # nolint
                  ld3_max <- ld3_peak[
                    ggpmisc::find_peaks(
                      ld3_peak[["Y"]],
                      span = spn,
                      global.threshold = ld3_threshold
                    ),
                  ]
                }
                if(nrow(ld3_max) == 1) { # nolint
                  ld3_max <- dplyr::mutate(
                    ld3_max,
                    "Xdiff" = abs(ld3_max[["X"]] - unique(ld3c[["exp_RT"]]))
                  )
                }
                if(nrow(ld3_max) > 1) { # nolint
                  ld3_max <- dplyr::mutate(
                    ld3_max,
                    "Xdiff" = abs(ld3_max[["X"]] - unique(ld3c[["exp_RT"]]))
                  )
                  # Select maximum peak that is closest to reference
                  ld3_max <- ld3_max[
                    ld3_max[["Xdiff"]] == min(ld3_max[["Xdiff"]]),
                  ]
                }
                ## Local minima relative to detected peak
                ### Subset
                ld3_peak <- ld3
                ld3_min_1 <- ld3_peak[
                  ld3_peak[["X"]] < ld3_max[["X"]] &
                    ld3_peak[["X"]] > ld3_max[["X"]] - rttol,
                ]
                ld3_min_2 <- ld3_peak[
                  ld3_peak[["X"]] > ld3_max[["X"]] &
                    ld3_peak[["X"]] < ld3_max[["X"]] + rttol,
                ]
                ### trailing edge
                min1 <- ggpmisc::find_valleys(ld3_min_1[["Y"]], span = spn)
                if(length(min1[min1 == TRUE]) == 0) { # nolint
                  ld3_min_1 <- ld3_min_1[
                    ld3_min_1[["Y"]] == min(ld3_min_1[["Y"]]),
                  ]
                }
                if(length(min1[min1 == TRUE]) > 0) { # nolint
                  ld3_min_1 <- ld3_min_1[min1, ]
                }
                #### leading edge
                min2 <- ggpmisc::find_valleys(ld3_min_2[["Y"]], span = spn)
                if(length(min2[min2 == TRUE]) == 0) { # nolint
                  ld3_min_2 <- ld3_min_2[
                    ld3_min_2[["Y"]] == min(ld3_min_2[["Y"]]),
                  ]
                }
                if(length(min2[min2 == TRUE]) > 0) { # nolint
                  ld3_min_2 <- ld3_min_2[min2, ]
                }
                #### filter if more than one value present
                if(nrow(ld3_min_1) == 1) { # nolint
                  ld3_min_1 <- dplyr::mutate(
                    ld3_min_1,
                    "Xdiff" = abs(ld3_min_1[["X"]] - ld3_max[["X"]])
                  )
                  ld3_min_1 <- ld3_min_1[ld3_min_1[["X"]] < ld3_max[["X"]], ]
                }
                if(nrow(ld3_min_1) > 1) { # nolint
                  ld3_min_1 <- dplyr::mutate(
                    ld3_min_1,
                    "Xdiff" = abs(ld3_min_1[["X"]] - ld3_max[["X"]])
                  )
                  ld3_min_1 <- ld3_min_1[ld3_min_1[["X"]] < ld3_max[["X"]], ]
                  ld3_min_1 <- ld3_min_1[
                    ld3_min_1[["Xdiff"]] == min(ld3_min_1[["Xdiff"]]),
                  ]
                }
                if(nrow(ld3_min_2) == 1) { # nolint
                  ld3_min_2 <- dplyr::mutate(
                    ld3_min_2,
                    "Xdiff" = abs(ld3_min_2[["X"]] - ld3_max[["X"]])
                  )
                  ld3_min_2 <- ld3_min_2[ld3_min_2[["X"]] > ld3_max[["X"]], ]
                }
                if(nrow(ld3_min_2) > 1) { # nolint
                  ld3_min_2 <- dplyr::mutate(
                    ld3_min_2,
                    "Xdiff" = abs(ld3_min_2[["X"]] - ld3_max[["X"]])
                  )
                  ld3_min_2 <- ld3_min_2[ld3_min_2[["X"]] > ld3_max[["X"]], ]
                  ld3_min_2 <- ld3_min_2[
                    ld3_min_2[["Xdiff"]] == min(ld3_min_2[["Xdiff"]]),
                  ]
                }
                ## Define peak boundaries
                ld3_peak <- dplyr::bind_rows(
                  ld3_min_1,
                  ld3_max,
                  ld3_min_2
                )
                ## Recombine detected peak with metadata
                ld3_peak <- data.frame(
                  "pkRT" = ld3_peak[["X"]][[2]],
                  "pkHeight" = ld3_peak[["Y"]][[2]],
                  "pkBoundary" = paste(ld3_peak[["X"]], collapse = ","),
                  "SampleID" = ld3_peak[["SampleID"]][[2]],
                  "Name" = ld3_peak[["Name"]][[2]],
                  "pkRTdiff" = ld3_peak[["Xdiff"]][[2]]
                )
                pks2 <- dplyr::select(dplyr::left_join(
                  ld3c, # nolint
                  ld3_peak,
                  by = c("Name", "SampleID")
                ), c("Name", "SampleID", everything())) # nolint
              },
              error = function(e) {
                print(paste("No peaks detected in", s1[[i]]))
              }
            )
            if(!exists("pks2")) { # nolint
              pks2 <- data.frame(
                ld2[ld2[["SampleID"]] == s1[[i]], ],
                "pkRT" = NA,
                "pkHeight" = NA,
                "pkBoundary" = NA,
                "pkRTdiff" = NA
              )
              colnames(pks2) <- c(
                names(ld2),
                "pkRT", "pkHeight", "pkBoundary", "pkRTdiff"
              )
            }
            return(pks2) # nolint
          }
        )
      )
    }
  }
  return(pks) # nolint
}
