#' Chromatogram Viewer (Targeted)
#'
#' Plots chromatograms for selected samples and compounds from
#' input x and y coordinates.
#'
#' @param dat Input data list.
#' @param intype Input data type: either "Skyline," a list of 2 vectors
#' containing values separated by commas (default) or "xy," which is
#' a list of two variables containing vectors of x and y coordinates.
#' @param xcol X column name.
#' @param ycol Y column name.
#' @param comp Compound name to plot.
#' @param ref Reference group name for determining peak retention time.
#' @param rtcol Retention time column containing the specific retention
#' times of identified peaks.
#' @param grp Sample treatment/QC group to plot.
#' @param lab Column to use for labelling samples.
#' @param mcc Cores to use if not using Windows.
#' @return A chromatogram based on X and Y input values.
#' @examples
#'
#' # msT_chromatogram <- function( # nolint
#' #   dat = d1[["merged"]],
#' #   comp = "15-HETE",
#' #   ref = "Std_curve"
#' # )
#'
#' @import ggplot2
#' @export
msT_chromatogram <- function( # nolint
  dat,
  intype = "Skyline",
  xcol = "RT_all",
  ycol = "Intensity",
  mcc = 8,
  comp,
  ref,
  grp,
  rtcol = "RT",
  lab = "SampleID"
) {
  # Load data
  ld <- dat
  if(intype == "Skyline") { # nolint
    # Convert RT and intensities to X and Y coordinates
    if(Sys.info()[["sysname"]] != "Windows") { # nolint
      ld1 <- dplyr::bind_rows(parallel::mclapply(
        mc.cores = mcc,
        seq.int(1, nrow(ld), 1),
        function(x) {
          ld2 <- data.frame(
            "X" = as.numeric(
              unlist(stringr::str_split(ld[x, xcol], ","))
            ),
            "Y" = as.numeric(
              unlist(stringr::str_split(ld[x, ycol], ","))
            ),
            "SampleID" = ld[x, lab],
            "Group" = ld[x, "Group"],
            "Name" = ld[x, "Name"]
          )
          return(ld2) # nolint
        }
      ))
    }
    if(Sys.info()[["sysname"]] == "Windows") { # nolint
      ld1 <- dplyr::bind_rows(lapply(
        seq.int(1, nrow(ld), 1),
        function(x) {
          ld2 <- data.frame(
            "X" = as.numeric(
              unlist(stringr::str_split(ld[x, xcol], ","))
            ),
            "Y" = as.numeric(
              unlist(stringr::str_split(ld[x, xcol], ","))
            ),
            "SampleID" = ld[x, lab],
            "Group" = ld[x, "Group"],
            "Name" = ld[x, "Name"]
          )
          return(ld2) # nolint
        }
      ))
    }
    # Select metadata
    ld2 <- dplyr::select(ld, -c(xcol, ycol))
    # Convert sample IDs to factor
    ld1[["SampleID"]] <- factor(
      ld1[["SampleID"]],
      levels = unique(ld1[["SampleID"]])
    )
    # Determine median retention time of reference peak
    # for selected compound
    ld2_peak <- median(
      ld2[ld2[["Name"]] == comp & ld2[["Group"]] == ref, rtcol]
    )
    # Subset data for plot input
    ld3 <- ld1[
      ld1[["Name"]] == comp &
        ld1[["Group"]] == grp,
    ]
    # Create separate data frame for storing sample labels
    ld3b <- levels(ld3[["SampleID"]])[
      levels(ld3[["SampleID"]]) %in% as.character(ld3[["SampleID"]])
    ]
    ld4 <- dplyr::bind_rows(
      lapply(
        seq.int(1, length(ld3b), 1),
        function(x) {
          tryCatch(
            {
              data.frame(
                "SampleID" = ld3b[[x]],
                "X" = median(ld3[ld3[["SampleID"]] == ld3b[[x]], "X"]),
                "Y" = ld3[
                  ld3[["SampleID"]] == ld3b[[x]] &
                    ld3[["X"]] == median(
                      ld3[ld3[["SampleID"]] == ld3b[[x]], "X"]
                    ),
                  "Y"
                ]
              )
            },
            error = function(e) {
              data.frame("SampleID" = ld3b[[x]], "X" = NA, "Y" = NA)
            }
          )
        }
      )
    )
    ld4[["SampleID"]] <- factor(
      ld4[["SampleID"]],
      levels = levels(ld3[["SampleID"]])
    )
    # Create plot
    p1 <- ggplot2::ggplot(
      data = ld3,
      ggplot2::aes(
        x = .data[["X"]], # nolint
        y = .data[["Y"]],
        color = .data[["SampleID"]]
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_text(
        data = ld4,
        ggplot2::aes(
          x = ld4[["X"]],
          y = ld4[["Y"]],
          label = ld4[["SampleID"]]
        ),
        nudge_x = -0.05
      ) +
      ggplot2::geom_vline(
        xintercept = ld2_peak,
        linetype = "dashed"
      ) +
      ggplot2::scale_color_manual(values = col_univ()) +
      ms_theme() +
      ggplot2::labs(
        x = "Retention Time (min.)",
        y = "Intensity"
      )
  }
  return(p1)
}

#' Peak Integration (Targeted)
#'
#' Detects and integrates chromatogram peaks using reference standards.
#'
#' @param dat Input data list.
#' @param intype Input data type: either "Skyline," a list of 2 vectors
#' containing values separated by commas (default) or "xy," which is
#' a list of two variables containing vectors of x and y coordinates.
#' @param xcol X column name.
#' @param ycol Y column name.
#' @param comp Compound name to plot.
#' @param ref Reference group name for determining peak retention time.
#' @param rtcol Retention time column containing the specific retention
#' times of identified peaks.
#' @param rttol Retention time tolerance window (in minutes).
#' @param spn Number of points used for peak detection; must be an odd number.
#' @param grp Sample treatment/QC group to plot.
#' @param lab Column containing sample labels.
#' @param sampid Name of a specific sample for peak integration.
#' @param mcc Cores to use if not using Windows.
#' @return A chromatogram based on X and Y input values.
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
tms_integrate <- function( # nolint
  dat,
  intype = "Skyline",
  xcol = "RT_all",
  ycol = "Intensity",
  mcc = 8,
  comp,
  ref,
  grp,
  rtcol = "RT",
  rttol = 0.15,
  spn = 7,
  lab = "SampleID",
  sampid
) {
  # Load data
  ld <- dat
  xc <- xcol
  yc <- ycol
  lb <- lab
  cmp <- comp
  ref1 <- ref
  rtc <- rtcol
  grp1 <- grp
  # Peak detection and integration
  if(intype == "Skyline") { # nolint
    # Convert RT and intensities to X and Y coordinates
    if(Sys.info()[["sysname"]] != "Windows") { # nolint
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
    ld2 <- dplyr::select(ld, -c(xc, yc))
    # Convert sample IDs to factor
    ld1[["SampleID"]] <- factor(
      ld1[["SampleID"]],
      levels = unique(ld1[["SampleID"]])
    )
    # Determine median retention time of reference peak
    # for selected compound
    ld2_peak <- median(
      ld2[ld2[["Name"]] == cmp & ld2[["Group"]] == ref1, rtc]
    )
    # Subset data for plot input
    ld3 <- ld1[
      ld1[["Name"]] == cmp &
        ld1[["Group"]] == grp1,
    ]
    # Create separate data frame for storing sample labels
    ld3b <- levels(ld3[["SampleID"]])[
      levels(ld3[["SampleID"]]) %in% as.character(ld3[["SampleID"]])
    ]
    ld4 <- dplyr::bind_rows(
      lapply(
        seq.int(1, length(ld3b), 1),
        function(x) {
          tryCatch(
            {
              data.frame(
                "SampleID" = ld3b[[x]],
                "X" = median(ld3[ld3[["SampleID"]] == ld3b[[x]], "X"]),
                "Y" = ld3[
                  ld3[["SampleID"]] == ld3b[[x]] &
                    ld3[["X"]] == median(
                      ld3[ld3[["SampleID"]] == ld3b[[x]], "X"]
                    ),
                  "Y"
                ]
              )
            },
            error = function(e) {
              data.frame("SampleID" = ld3b[[x]], "X" = NA, "Y" = NA)
            }
          )
        }
      )
    )
    ld4[["SampleID"]] <- factor(
      ld4[["SampleID"]],
      levels = levels(ld3[["SampleID"]])
    )
    # Calculate x coordinates of local minimum and maximum based on reference RT
    ## Subset data
    ld3_peak <- ld3[ld3[["SampleID"]] == sampid, ]
    ld3_peak <- ld3_peak[
      ld3_peak[["X"]] > (ld2_peak - rttol) &
        ld3_peak[["X"]] < (ld2_peak + rttol),
    ]
    ## Local Maximum closest to corresponding to reference peak
    ## set global threshold to find peaks that are above 3x signal-to-noise
    ld3_threshold <- range(ld3_peak[["Y"]])[[1]] * 3
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
    if(nrow(ld3_max) > 1) { # nolint
      ld3_max <- dplyr::mutate(
        ld3_max,
        "Xdiff" = abs(ld3_max[["X"]] - ld2_peak)
      )
      ld3_max <- dplyr::select(
        ld3_max[ld3_max[["Xdiff"]] == min(ld3_max[["Xdiff"]]), ],
        -c("Xdiff")
      )
    }
    ## Local Minimum at either tail of the detected peak
    ld3_min_1 <- ld3_peak[ld3_peak[["X"]] < ld2_peak, ]
    ld3_min_2 <- ld3_peak[ld3_peak[["X"]] > ld2_peak, ]
    min1 <- ggpmisc::find_valleys(ld3_min_1[["Y"]], span = spn)
    if(length(min1[min1 == TRUE]) == 0) { # nolint
      ld3_min_1 <- ld3_min_1[ld3_min_1[["Y"]] == min(ld3_min_1[["Y"]]), ]
    }
    if(length(min1[min1 == TRUE]) > 0) { # nolint
      ld3_min_1 <- ld3_min_1[min1, ]
    }
    min2 <- ggpmisc::find_valleys(ld3_min_2[["Y"]], span = spn)
    if(length(min2[min2 == TRUE]) == 0) { # nolint
      ld3_min_2 <- ld3_min_2[ld3_min_2[["Y"]] == min(ld3_min_2[["Y"]]), ]
    }
    if(length(min2[min2 == TRUE]) > 0) { # nolint
      ld3_min_2 <- ld3_min_2[min2, ]
    }
    if(nrow(ld3_min_1) > 1) { # nolint
      ld3_min_1 <- dplyr::mutate(
        ld3_min_1,
        "Xdiff" = abs(ld3_min_1[["X"]] - ld2_peak)
      )
      ld3_min_1 <- dplyr::select(
        ld3_min_1[ld3_min_1[["Xdiff"]] == min(ld3_min_1[["Xdiff"]]), ],
        -c("Xdiff")
      )
    }
    if(nrow(ld3_min_2) > 1) { # nolint
      ld3_min_2 <- dplyr::mutate(
        ld3_min_2,
        "Xdiff" = abs(ld3_min_2[["X"]] - ld2_peak)
      )
      ld3_min_2 <- dplyr::select(
        ld3_min_2[ld3_min_2[["Xdiff"]] == min(ld3_min_2[["Xdiff"]]), ],
        -c("Xdiff")
      )
    }
    ## Calculate Area under curve
    ld3_peak <- dplyr::bind_rows(
      ld3_min_1,
      ld3_max,
      ld3_min_2
    )
    # Check area prior to integration
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
      ggplot2::geom_text(
        data = ld4[ld4[["SampleID"]] == sampid, ],
        ggplot2::aes(
          x = ld4[ld4[["SampleID"]] == sampid, "X"],
          y = ld4[ld4[["SampleID"]] == sampid, "Y"],
          label = ld4[ld4[["SampleID"]] == sampid, "SampleID"]
        ),
        color = "grey20",
        nudge_x = -0.15,
        nudge_y = 0.15,
        show.legend = FALSE
      ) +
      ggplot2::geom_vline(
        xintercept = ld2_peak,
        linetype = "dashed"
      ) +
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
        y = "Intensity"
      )
    # Calculate AUC for detected peak
    if(min(ld3_peak[["Y"]]) * 3 > max(ld3_peak[["Y"]])) { # nolint
      print("No peak detected for selected sample.")
      print("SNR < 3 or compound is not present in sample.")
      print("Returning NA for AUC calculation")
      auc <- NA
    }
    if(min(ld3_peak[["Y"]]) * 3 < max(ld3_peak[["Y"]])) { # nolint
      auc <- DescTools::AUC(
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
      )
    }
  }
  out <- list("peak" = p1, "AUC" = auc)
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
tms_mzrt <- function( # nolint
  dat
) {
  # Load data and set parameters
  ld <- dat
  ## mzrt plot
  p1 <- ggplot2::ggplot(
    data = ld
  ) +
    ggplot2::geom_point(
      data = rt3,
      ggplot2::aes(
        x = rt3[["Name"]],
        y = rt3[["ref_RT"]]
      ),
      shape = 18,
      size = 2,
      color = "red",
      alpha = 0.9
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data[["Name"]],
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
      x = "Retention Time (min.)",
      y = "Compound Name"
    )
  return(p1)
}
