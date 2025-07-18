#' Data Normalization
#'
#' Includes multiple methods for data normalization. Currently supported methods
#' are total ion count (TIC) normalization.
#' Planned methods include iSTD, SERRF, and LOESS.
#'
#' @param ld1 List data object generated from ms_input().
#' @param mtd Normalization method (either "none" or "tic").
#' @param parl (optional) Logical indicating if normalization
#' should be performed in parallel.
#' @param core_perc (optional) Percentage of cores to use for parallel
#' normalization.
#' @return A list containing normalized compound intensities and metadata.
#' @examples
#'
#' # d_norm <- ms_data_norm(
#' #   dm = d_norm[["Data"]],
#' #   mtd = "tic"
#' # )
#'
#' @export
ms_data_norm <- function( # nolint
  ld1,
  mtd,
  parl = FALSE,
  core_perc = NULL
) {
  # Load objects
  d <- ld1[["data"]]
  mpx <- ld1[["meta"]]
  mft <- ld1[["anno"]]
  # No normalization
  if(mtd == "none") { # nolint
    print("Performing no normalization...")
    if(Sys.info()[["sysname"]] == "Windows" | parl == FALSE) { # nolint
      d <- as.data.frame(d)
      # data imputation
      ## replace zeroes or NA with 1/10th of
      ## the lowest non-zero value present
      d <- setNames(
        as.data.frame(
          lapply(
            seq.int(1, ncol(d), 1),
            function(x) {
              d[[x]][is.na(d[[x]])] <- 0
              d1 <- d[[x]]
              d1 <- ifelse(
                d1 == 0,
                round(0.1 * min(d1[d1 > 0]), digits = 0),
                d1
              )
              return(d1)
            }
          )
        ),
        c(names(d))
      )
      mpx[["TIC.prenorm"]] <- unlist(lapply(
        seq.int(1, ncol(t(as.matrix(d))), 1),
        function(x) sum(d[x, ])
      ))
      mpx[["TIC.postnorm"]] <- mpx[["TIC.prenorm"]]
      ## Calculate Group RSD pre- and post-normalization
      mpx <- dplyr::left_join(
        mpx,
        data.frame(
          "Group" = aggregate(
            d[[1]],
            list(mpx[["Group"]]),
            function(z) (sd(z) / mean(z)) * 100
          )[[1]],
          "RSD.median.pre" = round(unlist(lapply(
            seq.int(1, length(unique(mpx[["Group"]])), 1),
            function(a) {
              median(
                as.data.frame(t(dplyr::bind_cols(
                  lapply(
                    seq.int(1, ncol(d), 1),
                    function(y) {
                      setNames(as.data.frame(aggregate(
                        d[[y]],
                        list(mpx[["Group"]]),
                        function(z) (sd(z) / mean(z)) * 100
                      )[[2]]), paste("X", y, sep = ""))
                    }
                  )
                )))[[a]]
              )
            }
          )), digits = 2)
        ),
        by = "Group"
      )
      mpx[["RSD.median.post"]] <- mpx[["RSD.median.pre"]]
      mpx[["ID"]] <- seq.int(1, nrow(mpx), 1)
    }
    if(parl == TRUE) { # nolint
      d <- as.data.frame(d)
      # data imputation
      ## replace zeroes or NA with 1/10th of
      ## the lowest non-zero value present
      d <- setNames(
        as.data.frame(
          lapply(
            seq.int(1, ncol(d), 1),
            function(x) {
              d[[x]][is.na(d[[x]])] <- 0
              d1 <- d[[x]]
              d1 <- ifelse(
                d1 == 0,
                round(0.1 * min(d1[d1 > 0]), digits = 0),
                d1
              )
              return(d1)
            }
          )
        ),
        c(names(d))
      )
      mpx[["TIC.prenorm"]] <- unlist(lapply(
        seq.int(1, ncol(t(as.matrix(d))), 1),
        function(x) sum(d[x, ])
      ))
      mpx[["TIC.postnorm"]] <- mpx[["TIC.prenorm"]]
      ## Calculate Group RSD pre- and post-normalization
      mpx <- dplyr::left_join(
        mpx,
        data.frame(
          "Group" = aggregate(
            d[[1]],
            list(mpx[["Group"]]),
            function(z) (sd(z) / mean(z)) * 100
          )[[1]],
          "RSD.median.pre" = round(unlist(lapply(
            seq.int(1, length(unique(mpx[["Group"]])), 1),
            function(a) {
              median(
                as.data.frame(t(dplyr::bind_cols(
                  parallel::mclapply(
                    mc.cores = core_perc,
                    seq.int(1, ncol(d), 1),
                    function(y) {
                      setNames(as.data.frame(aggregate(
                        d[[y]],
                        list(mpx[["Group"]]),
                        function(z) (sd(z) / mean(z)) * 100
                      )[[2]]), paste("X", y, sep = ""))
                    }
                  )
                )))[[a]]
              )
            }
          )), digits = 2)
        ),
        by = "Group"
      )
      mpx[["RSD.median.post"]] <- mpx[["RSD.median.pre"]]
      mpx[["ID"]] <- seq.int(1, nrow(mpx), 1)
    }
    d1 <- list(
      "data" = d,
      "meta" = mpx,
      "anno" = mft,
      "norm.method" = "none"
    )
  }
  # TIC normalization
  if(mtd == "tic") { # nolint
    print("Performing TIC normalization...")
    if(Sys.info()[["sysname"]] == "Windows" | parl == FALSE) { # nolint
      print(
        "Detected OS is Windows or parl is FALSE;
        Defaulting to sequential processing..."
      )
      d <- as.data.frame(d)
      d <- setNames(
        as.data.frame(
          lapply(
            seq.int(1, ncol(d), 1),
            function(x) {
              d[[x]][is.na(d[[x]])] <- 0
              d1 <- d[[x]]
              d1 <- ifelse(
                d1 == 0,
                round(0.1 * min(d1[d1 > 0]), digits = 0),
                d1
              )
              return(d1)
            }
          )
        ),
        c(names(d))
      )
      dtic <- data.frame(
        "TIC.prenorm" = unlist(
          lapply(
            seq.int(1, ncol(t(as.matrix(d))), 1),
            function(x) sum(d[x, ])
          )
        )
      )
      ## Calculate Group RSD pre- and post-normalization
      mpx <- dplyr::left_join(
        mpx,
        data.frame(
          "Group" = aggregate(
            d[[1]],
            list(mpx[["Group"]]),
            function(z) (sd(z) / mean(z)) * 100
          )[[1]],
          "RSD.median.pre" = round(unlist(lapply(
            seq.int(1, length(unique(mpx[["Group"]])), 1),
            function(a) {
              median(
                as.data.frame(t(dplyr::bind_cols(
                  lapply(
                    seq.int(1, ncol(d), 1),
                    function(y) {
                      setNames(as.data.frame(aggregate(
                        d[[y]],
                        list(mpx[["Group"]]),
                        function(z) (sd(z) / mean(z)) * 100
                      )[[2]]), paste("X", y, sep = ""))
                    }
                  )
                )))[[a]]
              )
            }
          )), digits = 2)
        ),
        by = "Group"
      )
      dtic <- dplyr::select(
        dplyr::left_join(
          data.frame("Group" = mpx[["Group"]], dtic),
          setNames(
            aggregate(
              dtic[["TIC.prenorm"]],
              list(
                mpx[["Group"]]
              ),
              function(x) mean(x)
            ),
            c("Group", "TIC.avg")
          ),
          by = "Group"
        ),
        c("TIC.avg", dplyr::everything())
      )
      d <- setNames(as.data.frame(lapply(
        seq.int(1, ncol(d), 1),
        function(x) {
          (d[, x] / dtic[["TIC.prenorm"]]) *
            dtic[["TIC.avg"]]
        }
      )), names(d))
      ## Calculate Group RSD pre- and post-normalization
      mpx <- dplyr::left_join(
        mpx,
        data.frame(
          "Group" = aggregate(
            d[[1]],
            list(mpx[["Group"]]),
            function(z) (sd(z) / mean(z)) * 100
          )[[1]],
          "RSD.median.post" = round(unlist(lapply(
            seq.int(1, length(unique(mpx[["Group"]])), 1),
            function(a) {
              median(
                as.data.frame(t(dplyr::bind_cols(
                  lapply(
                    seq.int(1, ncol(d), 1),
                    function(y) {
                      setNames(as.data.frame(aggregate(
                        d[[y]],
                        list(mpx[["Group"]]),
                        function(z) (sd(z) / mean(z)) * 100
                      )[[2]]), paste("X", y, sep = ""))
                    }
                  )
                )))[[a]]
              )
            }
          )), digits = 2)
        ),
        by = "Group"
      )
      mpx[["TIC.prenorm"]] <- dtic[["TIC.prenorm"]]
      mpx[["TIC.postnorm"]] <- unlist(lapply(
        seq.int(1, ncol(t(as.matrix(d))), 1),
        function(x) sum(d[x, ])
      ))
      mpx[["ID"]] <- seq.int(1, nrow(mpx), 1)
    }
    if(Sys.info()[["sysname"]] != "Windows" && # nolint
        parl == TRUE
    ) {
      d <- as.data.frame(d)
      d <- setNames(
        as.data.frame(
          lapply(
            seq.int(1, ncol(d), 1),
            function(x) {
              d[[x]][is.na(d[[x]])] <- 0
              d1 <- d[[x]]
              d1 <- ifelse(
                d1 == 0,
                round(0.1 * min(d1[d1 > 0]), digits = 0),
                d1
              )
              return(d1)
            }
          )
        ),
        c(names(d))
      )
      dtic <- data.frame(
        "TIC.prenorm" = unlist(
          lapply(
            seq.int(1, ncol(t(as.matrix(d))), 1),
            function(x) sum(d[x, ])
          )
        )
      )
      ## Calculate Group RSD pre- and post-normalization
      mpx <- dplyr::left_join(
        mpx,
        data.frame(
          "Group" = aggregate(
            d[[1]],
            list(mpx[["Group"]]),
            function(z) (sd(z) / mean(z)) * 100
          )[[1]],
          "RSD.median.pre" = round(unlist(lapply(
            seq.int(1, length(unique(mpx[["Group"]])), 1),
            function(a) {
              median(
                as.data.frame(t(dplyr::bind_cols(
                  parallel::mclapply(
                    mc.cores = core_perc,
                    seq.int(1, ncol(d), 1),
                    function(y) {
                      setNames(as.data.frame(aggregate(
                        d[[y]],
                        list(mpx[["Group"]]),
                        function(z) (sd(z) / mean(z)) * 100
                      )[[2]]), paste("X", y, sep = ""))
                    }
                  )
                )))[[a]]
              )
            }
          )), digits = 2)
        ),
        by = "Group"
      )
      dtic <- dplyr::select(
        dplyr::left_join(
          data.frame("Group" = mpx[["Group"]], dtic),
          setNames(
            aggregate(
              dtic[["TIC.prenorm"]],
              list(
                mpx[["Group"]]
              ),
              function(x) mean(x)
            ),
            c("Group", "TIC.avg")
          ),
          by = "Group"
        ),
        c("TIC.avg", dplyr::everything())
      )
      d <- setNames(as.data.frame(lapply(
        seq.int(1, ncol(d), 1),
        function(x) {
          (d[, x] / dtic[["TIC.prenorm"]]) *
            dtic[["TIC.avg"]]
        }
      )), names(d))
      ## Calculate Group RSD pre- and post-normalization
      mpx <- dplyr::left_join(
        mpx,
        data.frame(
          "Group" = aggregate(
            d[[1]],
            list(mpx[["Group"]]),
            function(z) (sd(z) / mean(z)) * 100
          )[[1]],
          "RSD.median.post" = round(unlist(lapply(
            seq.int(1, length(unique(mpx[["Group"]])), 1),
            function(a) {
              median(
                as.data.frame(t(dplyr::bind_cols(
                  parallel::mclapply(
                    mc.cores = core_perc,
                    seq.int(1, ncol(d), 1),
                    function(y) {
                      setNames(as.data.frame(aggregate(
                        d[[y]],
                        list(mpx[["Group"]]),
                        function(z) (sd(z) / mean(z)) * 100
                      )[[2]]), paste("X", y, sep = ""))
                    }
                  )
                )))[[a]]
              )
            }
          )), digits = 2)
        ),
        by = "Group"
      )
      mpx[["TIC.prenorm"]] <- dtic[["TIC.prenorm"]]
      mpx[["TIC.postnorm"]] <- unlist(lapply(
        seq.int(1, ncol(t(as.matrix(d))), 1),
        function(x) sum(d[x, ])
      ))
      mpx[["ID"]] <- seq.int(1, nrow(mpx), 1)
    }
    d1 <- list(
      "data" = d,
      "meta" = mpx,
      "anno" = mft,
      "norm.method" = "tic"
    )
  }
  return(d1)
}

#' Normalization QC
#'
#' Scatter plot visualization to evaluate normalization performance.
#'
#' @param df input data from list data object.
#' @param var_x X-axis variable (usually sample ID).
#' @param var_g Grouping variable.
#' @return A scatter plot visualizing overall sample intensities.
#' @examples
#'
#' # ptic <- msi_plot_tic(
#' #   df = ld[["meta"]],
#' #   var_x = "sampleID",
#' #   var_g = "Group"
#' # )
#'
#' @export
ms_plot_tic <- function(
  df,
  var_x,
  var_g = "Group"
) {
  ## TIC plot
  ggplot2::ggplot(
    df
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data[[var_x]], # nolint
        y = .data[["TIC.prenorm"]],
        color = as.factor(.data[[var_g]])
      ),
      shape = 16,
      size = 1,
      alpha = 0.5
    ) +
    ggplot2::geom_smooth(
      color = "firebrick1",
      alpha = 0.5,
      ggplot2::aes(x = .data[[var_x]], y = .data[["TIC.prenorm"]])
    ) +
    ggplot2::labs(y = "TIC",
      x = "Sample ID"
    ) +
    ggplot2::scale_color_manual(values = col_univ()) + # nolint
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data[[var_x]],
        y = df[["TIC.postnorm"]],
        color = as.factor(.data[[var_g]])
      ),
      shape = 16,
      size = 1,
      alpha = 0.8
    ) +
    ggplot2::geom_smooth(
      color = "dodgerblue1",
      alpha = 0.8,
      ggplot2::aes(x = .data[[var_x]], y = df[["TIC.postnorm"]])
    ) +
    ggplot2::labs(y = "TIC",
      x = "Sample ID"
    ) +
    ms_theme1() # nolint
}
