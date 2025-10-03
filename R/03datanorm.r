#' Data Normalization
#'
#' Includes multiple methods for data normalization. Currently supported methods
#' are total ion count (TIC) normalization.
#' Planned methods include iSTD, SERRF, and LOESS.
#'
#' @param ld1 List data object generated from ms_input().
#' @param mtd Normalization method (either "none", "tic", or "LOESS").
#' @param bl_rem Should blank samples be removed from the input data?
#' @param bl_col If bl_rem is TRUE, provides the column name including
#' sample groups or identities, including blanks.
#' @param bl_nm Character string for detecting blank samples.
#' @param col_id Name of the sample ID column.
#' @param col_order If mtd is "LOESS," provide the column name specifying
#' the sample injection order.
#' @param col_grp If mtd is "LOESS," provide the column name specifying
#' the sample group information for calculating span.
#' @param col_nm If mtd is "LOESS," provide the column name specifying
#' the feature and sample names.
#' @param spn LOESS span.
#' @param qc_nm QC sample name within the group column (e.g. "Pool").
#' @param msg Show group RSDs in the console output.
#' @return A list containing normalized compound intensities and metadata.
#' @examples
#'
#' # d_norm <- ms_data_norm(ld1 = d_norm[["Data"]])
#'
#' @import dplyr
#' @import parallel
#' @import magrittr
#' @export
ms_data_norm <- function( # nolint
  ld1,
  mtd = "LOESS",
  bl_rem = TRUE,
  bl_col = "File type",
  bl_nm = "Blank",
  col_id = "ID",
  col_order = "Order",
  col_grp = "Group",
  col_nm = "Label",
  qc_nm = "Pool",
  msg = FALSE
) {
  # Load objects
  d <- ld1[["data"]]
  mpx <- ld1[["meta"]]
  mft <- ld1[["anno"]]
  # Remove blanks
  if (bl_rem == TRUE && mtd == "LOESS") {
    mpx <- mpx[mpx[[bl_col]] != bl_nm, ]
    d <- d[names(d) %in% mpx[[col_id]]]
  }
  # No normalization
  if(mtd == "none") { # nolint
    print("Performing no normalization...")
    if (Sys.info()[["sysname"]] == "Windows") { # nolint
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
              return(d1) # nolint
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
    if (Sys.info()[["sysname"]] != "Windows") { # nolint
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
              return(d1) # nolint
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
                    mc.cores = ceiling(parallel::detectCores() / 2),
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
    if (Sys.info()[["sysname"]] == "Windows") { # nolint
      print(
        "Detected OS is Windows;
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
              return(d1) # nolint
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
    if (Sys.info()[["sysname"]] != "Windows") {
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
              return(d1) # nolint
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
                    mc.cores = ceiling(parallel::detectCores() / 2),
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
                    mc.cores = ceiling(parallel::detectCores() / 2),
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
  # LOESS normalization
  if (mtd == "LOESS") {
    print("Performing LOESS normalization...")
    if (Sys.info()[["sysname"]] != "Windows") {
      # LOESS Normalization
      d1 <- dplyr::bind_cols(
        parallel::mclapply(
          mc.cores = ceiling(parallel::detectCores() / 2),
          seq.int(1, nrow(d), 1),
          function(i) {
            iter <- i
            print(
              paste("Normalizing compound ", iter, " of ", nrow(d), sep = "")
            )
            #---- Train LOESS ----
            dl <- data.frame(
              "type" = mpx[[bl_col]],
              "x" = mpx[[col_order]],
              "y" = unlist(d[
                iter,
              ])
            )
            dl <- dl[order(dl[["x"]]), ]
            if (sum(dl[["y"]]) == 0) {
              print("Compound intensities are all 0; returning original data")
              fit[["norm"]] <- dl[["y"]]
            }
            if (sum(dl[["y"]]) > 0) {
              if (length(dl[dl[["y"]] <= 0, ][["y"]]) > 0) {
                dl[dl[["y"]] <= 0, ][["y"]] <- round(
                  0.1 * min(dl[dl[["y"]] > 0, ][["y"]]),
                  digits = 2
                )
              }
              # Fit model
              fit <- dplyr::mutate(
                dl,
                "int" = predict(
                  loess(
                    dl[["y"]] ~ dl[["x"]],
                    data = dl,
                    ## set span equal to proportion of largest
                    ## group size relative to total samples
                    span = max(
                      dplyr::count(mpx, mpx[[col_grp]])[[2]]
                    ) / nrow(mpx)
                  )
                )
              )
              # Replace negative values with 1/10th minimum value
              # of normalized intensity
              if (length(fit[fit[["int"]] <= 0, ][["int"]]) > 0) {
                fit[fit[["int"]] <= 0, ][["int"]] <- round(
                  0.1 * min(fit[fit[["int"]] > 0, ][["int"]]),
                  digits = 2
                )
              }
              #---- Normalize based on fitted model using MA method ----
              ## Minus
              fit[["M"]] <- log2(fit[["y"]] / fit[["int"]])
              ## Average
              fit[["A"]] <- log2(fit[["y"]] * fit[["int"]]) / 2
              # Corrected MA
              fit[["MAnorm"]] <- loess(
                fit[["M"]] ~ fit[["A"]],
                span = max(
                  dplyr::count(mpx, mpx[[col_grp]])[[2]]
                ) / nrow(mpx)
              )[["residuals"]]
              if (length(fit[is.na(fit[["MAnorm"]]), ][["MAnorm"]]) > 0) {
                fit[["MAnorm2"]] <- fit[["MAnorm"]]
                fit[is.na(fit[["MAnorm2"]]), ][["MAnorm2"]] <- 0
              }
              if (length(fit[is.na(fit[["MAnorm"]]), ][["MAnorm"]]) == 0) {
                fit[["MAnorm2"]] <- fit[["MAnorm"]]
              }
              ## calculate normalized values
              fit[["norm"]] <- (
                2 ^ ((log2(fit[["y"]]) - fit[["MAnorm2"]]) / 2) /
                  2 ^ ((log2(fit[["int"]]) - fit[["MAnorm2"]]) / 2)
              ) *
                mean(fit[["y"]])
              if (
                length(
                  fit[
                    fit[["norm"]] <= 0 |
                      is.na(fit[["MAnorm"]]),
                  ][["norm"]]
                ) > 0
              ) {
                fit[
                  fit[["norm"]] <= 0 |
                    is.na(fit[["MAnorm"]]),
                ][["norm"]] <- round(
                  0.1 * min(fit[fit[["norm"]] > 0, ][["norm"]]),
                  digits = 0
                )
              }
              ## pre and post-normalization QC RSD
              ### pre
              (sd(fit[
                fit[["x"]] %in%
                  mpx[mpx[[col_grp]] == qc_nm, ][[col_order]]
                ,
              ][["y"]]) / mean(fit[
                fit[["x"]] %in%
                  mpx[mpx[[col_grp]] == qc_nm, ][[col_order]]
                ,
              ][["y"]])) * 100
              ### post
              (sd(fit[
                fit[["x"]] %in%
                  mpx[mpx[[col_grp]] == qc_nm, ][[col_order]]
                ,
              ][["norm"]]) / mean(fit[
                fit[["x"]] %in%
                  mpx[mpx[[col_grp]] == qc_nm, ][[col_order]]
                ,
              ][["norm"]])) * 100
              # Output normalized data
              if (msg == TRUE) {
                print("Group RSD changed from:")
                print(aggregate(
                  x = fit[["y"]],
                  list(mpx[order(mpx[[col_order]]), ][[col_grp]]),
                  function(x) {
                    d1 <- data.frame(
                      "RSD" = round((sd(x) / mean(x)) * 100, digits = 2),
                      "mean" = round(mean(x), digits = 0),
                      "sd" = round(sd(x), digits = 0)
                    )
                    return(d1) # nolint
                  }
                ))
                print("to:")
                print(aggregate(
                  x = fit[["norm"]],
                  list(mpx[order(mpx[[col_order]]), ][[col_grp]]),
                  function(x) {
                    d1 <- data.frame(
                      "RSD" = round((sd(x) / mean(x)) * 100, digits = 2),
                      "mean" = round(mean(x), digits = 0),
                      "sd" = round(sd(x), digits = 0)
                    )
                    return(d1) # nolint
                  }
                ))
                print("After LOESS normalization")
              }
            }
            dl2 <- as.data.frame(magrittr::set_rownames(
              setNames(
                as.data.frame(fit[["norm"]]),
                mft[iter, ][[col_nm]]
              ),
              mpx[order(mpx[[col_order]]), ][[col_nm]]
            ))
            return(dl2) # nolint
          }
        )
      )
    }
    if (Sys.info()[["sysname"]] == "Windows") {
      # LOESS Normalization
      d1 <- dplyr::bind_cols(
        lapply(
          seq.int(1, nrow(d), 1),
          function(i) {
            iter <- i
            print(
              paste("Normalizing compound ", iter, " of ", nrow(d), sep = "")
            )
            #---- Train LOESS ----
            dl <- data.frame(
              "type" = mpx[[bl_col]],
              "x" = mpx[[col_order]],
              "y" = unlist(d[
                iter,
              ])
            )
            dl <- dl[order(dl[["x"]]), ]
            if (sum(dl[["y"]]) == 0) {
              print("Compound intensities are all 0; returning original data")
              fit[["norm"]] <- dl[["y"]]
            }
            if (sum(dl[["y"]]) > 0) {
              if (length(dl[dl[["y"]] <= 0, ][["y"]]) > 0) {
                dl[dl[["y"]] <= 0, ][["y"]] <- round(
                  0.1 * min(dl[dl[["y"]] > 0, ][["y"]]),
                  digits = 2
                )
              }
              # Fit model
              fit <- dplyr::mutate(
                dl,
                "int" = predict(
                  loess(
                    dl[["y"]] ~ dl[["x"]],
                    data = dl,
                    ## set span equal to proportion of largest
                    ## group size relative to total samples
                    span = max(
                      dplyr::count(mpx, mpx[[col_grp]])[[2]]
                    ) / nrow(mpx)
                  )
                )
              )
              # Replace negative values with 1/10th minimum value
              # of normalized intensity
              if (length(fit[fit[["int"]] <= 0, ][["int"]]) > 0) {
                fit[fit[["int"]] <= 0, ][["int"]] <- round(
                  0.1 * min(fit[fit[["int"]] > 0, ][["int"]]),
                  digits = 2
                )
              }
              #---- Normalize based on fitted model using MA method ----
              ## Minus
              fit[["M"]] <- log2(fit[["y"]] / fit[["int"]])
              ## Average
              fit[["A"]] <- log2(fit[["y"]] * fit[["int"]]) / 2
              # Corrected MA
              fit[["MAnorm"]] <- loess(
                fit[["M"]] ~ fit[["A"]],
                span = max(
                  dplyr::count(mpx, mpx[[col_grp]])[[2]]
                ) / nrow(mpx)
              )[["residuals"]]
              if (length(fit[is.na(fit[["MAnorm"]]), ][["MAnorm"]]) > 0) {
                fit[["MAnorm2"]] <- fit[["MAnorm"]]
                fit[is.na(fit[["MAnorm2"]]), ][["MAnorm2"]] <- 0
              }
              if (length(fit[is.na(fit[["MAnorm"]]), ][["MAnorm"]]) == 0) {
                fit[["MAnorm2"]] <- fit[["MAnorm"]]
              }
              ## calculate normalized values
              fit[["norm"]] <- (
                2 ^ ((log2(fit[["y"]]) - fit[["MAnorm2"]]) / 2) /
                  2 ^ ((log2(fit[["int"]]) - fit[["MAnorm2"]]) / 2)
              ) *
                mean(fit[["y"]])
              if (
                length(
                  fit[
                    fit[["norm"]] <= 0 |
                      is.na(fit[["MAnorm"]]),
                  ][["norm"]]
                ) > 0
              ) {
                fit[
                  fit[["norm"]] <= 0 |
                    is.na(fit[["MAnorm"]]),
                ][["norm"]] <- round(
                  0.1 * min(fit[fit[["norm"]] > 0, ][["norm"]]),
                  digits = 0
                )
              }
              ## pre and post-normalization QC RSD
              ### pre
              (sd(fit[
                fit[["x"]] %in%
                  mpx[mpx[[col_grp]] == qc_nm, ][[col_order]]
                ,
              ][["y"]]) / mean(fit[
                fit[["x"]] %in%
                  mpx[mpx[[col_grp]] == qc_nm, ][[col_order]]
                ,
              ][["y"]])) * 100
              ### post
              (sd(fit[
                fit[["x"]] %in%
                  mpx[mpx[[col_grp]] == qc_nm, ][[col_order]]
                ,
              ][["norm"]]) / mean(fit[
                fit[["x"]] %in%
                  mpx[mpx[[col_grp]] == qc_nm, ][[col_order]]
                ,
              ][["norm"]])) * 100
              # Output normalized data
              if (msg == TRUE) {
                print("Group RSD changed from:")
                print(aggregate(
                  x = fit[["y"]],
                  list(mpx[order(mpx[[col_order]]), ][[col_grp]]),
                  function(x) {
                    d1 <- data.frame(
                      "RSD" = round((sd(x) / mean(x)) * 100, digits = 2),
                      "mean" = round(mean(x), digits = 0),
                      "sd" = round(sd(x), digits = 0)
                    )
                    return(d1) # nolint
                  }
                ))
                print("to:")
                print(aggregate(
                  x = fit[["norm"]],
                  list(mpx[order(mpx[[col_order]]), ][[col_grp]]),
                  function(x) {
                    d1 <- data.frame(
                      "RSD" = round((sd(x) / mean(x)) * 100, digits = 2),
                      "mean" = round(mean(x), digits = 0),
                      "sd" = round(sd(x), digits = 0)
                    )
                    return(d1) # nolint
                  }
                ))
                print("After LOESS normalization")
              }
            }
            dl2 <- as.data.frame(magrittr::set_rownames(
              setNames(
                as.data.frame(fit[["norm"]]),
                mft[iter, ][[col_nm]]
              ),
              mpx[order(mpx[[col_order]]), ][[col_nm]]
            ))
            return(dl2) # nolint
          }
        )
      )
    }
    d1 <- list(
      "data" = d1,
      "meta" = mpx[order(mpx[["Order"]]), ],
      "anno" = mft,
      "norm.method" = "LOESS"
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
#' @import ggplot2
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
    ms_theme() # nolint
}
