#' Pre-processing of Exported MS-DIAL Alignment Results
#'
#' Performs processing of alignment results exported from MS-DIAL
#' (v5.5 and higher). Adds a column flagging features for removal
#' based on presence in blanks, low abundance, or duplicate IDs.
#' Optionally annotates unknown features from matching to a reference
#' mzrt list.
#'
#' @param dattype Data type (either "lipidomics" or "metabolomics"). If the data
#' type is "lipidomics," features will be flagged for the presence of
#' characteristic MS/MS fragments for specific lipid classes.
#' @param ionmode Ionization mode (either "pos" or "neg").
#' @param dir1 Directory containing the exported alignment file.
#' @param dat File name of an alignment result directly exported from MS-DIAL.
#' @param ref2 (optional) The path to a text file containing compound names,
#' m/z, and retention time information for annotating unknowns.
#' @param ref1 Name of a text file containing sample metadata. Must ensure that
#' the order of samples listed in the alignment file matches the order of the
#' input sample list.
#' @param tolmz m/z tolerance window (in mDa) for annotating unknown peaks
#' from a reference mzrt list.
#' @param tolmz2 m/z tolerance window (in min) for flagging adducts and isotopes
#' of the same peak.
#' @param tolrt RT tolerance window (in min) for annotating unknown peaks
#' from a reference mzrt list.
#' @param tolrt2 RT tolerance window (in min) for flagging adducts and isotopes
#' of the same peak.
#' @param gcol Sample group column name.
#' @param gblnk Method blank group name for calculating sample/blank intensity
#' ratio.
#' @param blank_filt Numerical value indicating the threshold for removing
#' compounds based on intensity in method blanks.
#' Defaults to 5x sample/blank average intensity per group.
#' @param min_fill Minimum fill percentage (expressed as a proportion) for
#' including features in the preprocessed output.
#' @return A chromatogram based on X and Y coordinates.
#' @examples
#'
#' # preproc <- ms_preproc(
#' #   dat = "20250919_Chase_LSAE_CSH_pos_alignment",
#' #   ref1 = "ref/sample_list.csv",
#' #   ref2 = "ref/pos_mzrt_list.txt"
#' # )
#'
#' @import dplyr
#' @import magrittr
#' @export
ms_preproc <- function(
  dattype = "lipidomics",
  ionmode = "pos",
  dir1 = "data/01alignment/",
  dat,
  ref1,
  ref2 = NULL,
  tolmz = 0.01,
  tolmz2 = 0.01,
  tolrt = 0.15,
  tolrt2 = 0.01,
  gcol = "Group",
  gblnk = "Blank",
  blank_filt = 5,
  min_fill = 0.01
) {
  # Load alignment file and set parameters
  f1 <- dat
  r1 <- ref1
  r2 <- ref2
  d1 <- dir1
  t1 <- dattype
  t2 <- ionmode
  g1 <- gcol
  filt1 <- blank_filt
  ld <- list(
    "data" = read.table(
      paste(
        d1,
        f1,
        ".txt",
        sep = ""
      ),
      header = FALSE,
      sep = "\t",
      skip = 5,
      comment.char = "",
      quote = ""
    )[, 36:ncol(read.table(
      paste(
        d1,
        f1,
        ".txt",
        sep = ""
      ),
      header = FALSE,
      sep = "\t",
      skip = 5,
      comment.char = "",
      quote = ""
    ))],
    "feat" = read.table(
      paste(
        d1,
        f1,
        ".txt",
        sep = ""
      ),
      header = TRUE,
      sep = "\t",
      skip = 4,
      comment.char = "",
      quote = ""
    )[, 1:35],
    "samp" = read.table(
      paste(
        d1,
        f1,
        ".txt",
        sep = ""
      ),
      header = FALSE,
      sep = "\t",
      comment.char = "",
      quote = ""
    )[1:5, 35:ncol(read.table(
      paste(
        d1,
        f1,
        ".txt",
        sep = ""
      ),
      header = FALSE,
      sep = "\t",
      comment.char = "",
      quote = ""
    ))]
  )
  # Format alignment file inputs
  ## data
  ld[["data"]] <- as.data.frame(t(ld[["data"]]))
  ## features (leave as-is)
  ## samples
  ld[["samp"]] <- setNames(as.data.frame(t(ld[["samp"]])), ld[["samp"]][, 1])
  # Remove redundant columns and rows from inputs
  ## remove non-sample columns and merge sample metadata
  ld[["samp"]] <- ld[["samp"]][
    unlist(
      lapply(
        seq.int(1, nrow(ld[["samp"]]), 1),
        function(x) {
          !anyNA(ld[["samp"]][x, ])
        }
      )
    ),
  ][-1, ]
  if (grepl(".csv", r1)) {
    dsamp <- read.table(r1, sep = ",", header = TRUE)
  }
  if (grepl(".txt", r1)) {
    dsamp <- read.table(r1, sep = "\t", header = TRUE)
  }
  if ("Class" %in% names(dsamp)) {
    dsamp <- dplyr::select(dsamp, -c("Class"))
  }
  if (nrow(dsamp) != nrow(ld[["samp"]])) {
    print("Number of samples in alignment file and sample list do not match!")
  }
  ld[["samp"]] <- cbind(ld[["samp"]], dsamp)
  ## remove redundant columns from data and set column names and row names
  ld[["data"]] <- setNames(magrittr::set_rownames(
    ld[["data"]][rownames(ld[["samp"]]), ], ld[["samp"]][["ID"]]
  ), paste("Alignment.ID_", ld[["feat"]][["Alignment.ID"]], sep = ""))
  #----Annotate features based on reference mzrt list----
  if (!is.null(r2)) {
    dfeat <- read.table(
      r2,
      sep = "\t",
      header = TRUE
    )
    ## For each entry in the reference list,
    ## search unknowns in the feature list
    ## for matching alignment ID(s).
    if (Sys.info()[["sysname"]] != "Windows") {
      dfeat2 <- parallel::mclapply(
        mc.cores = ceiling(parallel::detectCores() / 2),
        seq.int(1, nrow(dfeat), 1),
        function(i) {
          d1 <- ld[["feat"]][
            ld[["feat"]][["Average.Mz"]] <
              dfeat[i, ][[2]] + tolmz &
              ld[["feat"]][["Average.Mz"]] >
                dfeat[i, ][[2]] - tolmz &
              ld[["feat"]][["Average.Rt.min."]] <
                dfeat[i, ][[3]] + tolrt &
              ld[["feat"]][["Average.Rt.min."]] >
                dfeat[i, ][[3]] - tolrt,
          ]
          if (nrow(d1) > 0) {
            d1[["Metabolite.name"]] <- dfeat[i, ][[1]]
            d1[["Reference.m.z"]] <- dfeat[i, ][[2]]
            d1[["Reference.RT"]] <- dfeat[i, ][[3]]
            d1[["RT.matched"]] <- rep(TRUE, nrow(d1))
            d1[["m.z.matched"]] <- rep(TRUE, nrow(d1))
          }
          return(d1) # nolint
        }
      )
    }
    if (Sys.info()[["sysname"]] == "Windows") {
      dfeat2 <- lapply(
        seq.int(1, nrow(dfeat), 1),
        function(i) {
          d1 <- ld[["feat"]][
            ld[["feat"]][["Average.Mz"]] <
              dfeat[i, ][[2]] + tolmz &
              ld[["feat"]][["Average.Mz"]] >
                dfeat[i, ][[2]] - tolmz &
              ld[["feat"]][["Average.Rt.min."]] <
                dfeat[i, ][[3]] + tolrt &
              ld[["feat"]][["Average.Rt.min."]] >
                dfeat[i, ][[3]] - tolrt,
          ]
          if (nrow(d1) > 0) {
            d1[["Metabolite.name"]] <- dfeat[i, ][[1]]
            d1[["Reference.m.z"]] <- dfeat[i, ][[2]]
            d1[["Reference.RT"]] <- dfeat[i, ][[3]]
            d1[["RT.matched"]] <- rep(TRUE, nrow(d1))
            d1[["m.z.matched"]] <- rep(TRUE, nrow(d1))
          }
          return(d1) # nolint
        }
      )
    }
    # Filter IDs without matches in the feature list
    dfeat2 <- dplyr::bind_rows(
      dfeat2[unlist(lapply(dfeat2, function(i) nrow(i) > 0))]
    )
    # Annotate alignment IDs in the feature list
    ld[["feat"]][dfeat2[["Alignment.ID"]], ] <- dfeat2
  }
  #----Complete pre-curation check for each alignment ID----
  d_check <- ld[["feat"]][
    , c(
      "Alignment.ID", "Metabolite.name",
      "Average.Mz", "Average.Rt.min."
    )
  ]
  # Identify standards
  d_check[["iSTD"]] <- grepl("iSTD", d_check[["Metabolite.name"]])
  # Total present in samples and per group
  d_check[["Prop.total"]] <- unlist(
    lapply(
      seq.int(1, ncol(ld[["data"]]), 1),
      function(i) {
        round(length(ld[["data"]][[i]][ld[["data"]][[i]] > 0]) /
            length(ld[["data"]][[i]]), digits = 2
        )
      }
    )
  )
  if (Sys.info()[["sysname"]] != "Windows") {
    d_check[["Prop.max.per.group"]] <- unlist(
      parallel::mclapply(
        mc.cores = ceiling(parallel::detectCores() / 2),
        seq.int(1, ncol(ld[["data"]]), 1),
        function(i) {
          d1 <- aggregate(
            ld[["data"]][[i]],
            list(ld[["samp"]][[g1]]),
            function(x) round(length(x[x > 0]) / length(x), digits = 3)
          )
          d1 <- max(d1[[2]])
          return(d1) # nolint
        }
      )
    )
  }
  if (Sys.info()[["sysname"]] == "Windows") {
    d_check[["Prop.max.per.group"]] <- unlist(
      lapply(
        seq.int(1, ncol(ld[["data"]]), 1),
        function(i) {
          d1 <- aggregate(
            ld[["data"]][[i]],
            list(ld[["samp"]][[g1]]),
            function(x) round(length(x[x > 0]) / length(x), digits = 3)
          )
          d1 <- max(d1[[2]])
          return(d1) # nolint
        }
      )
    )
  }
  # Sample / Blank intensity ratio
  if (Sys.info()[["sysname"]] != "Windows") {
    d_check[["Blank.ratio.max"]] <- unlist(
      parallel::mclapply(
        mc.cores = ceiling(parallel::detectCores() / 2),
        seq.int(1, ncol(ld[["data"]]), 1),
        function(i) {
          d1 <- aggregate(
            ld[["data"]][[i]],
            list(ld[["samp"]][[g1]]),
            function(x) mean(x)
          )
          if (d1[d1[[1]] == gblnk, ][[2]] == 0) {
            d1[d1[[1]] == gblnk, 2] <- round(
              min(d1[d1[[1]] != gblnk & d1[[2]] != 0, ][[2]]) / 10,
              digits = 3
            )
          }
          d1[["rat"]] <- round(
            d1[[2]] / d1[d1[[1]] == gblnk, ][[2]],
            digits = 3
          )
          d1 <- d1[d1[[1]] != gblnk, ]
          d1 <- max(d1[["rat"]])
          return(d1) # nolint
        }
      )
    )
  }
  if (Sys.info()[["sysname"]] == "Windows") {
    d_check[["Blank.ratio.max"]] <- unlist(
      lapply(
        seq.int(1, ncol(ld[["data"]]), 1),
        function(i) {
          d1 <- aggregate(
            ld[["data"]][[i]],
            list(ld[["samp"]][[g1]]),
            function(x) mean(x)
          )
          if (d1[d1[[1]] == gblnk, ][[2]] == 0) {
            d1[d1[[1]] == gblnk, 2] <- round(
              min(d1[d1[[1]] != gblnk & d1[[2]] != 0, ][[2]]) / 10,
              digits = 3
            )
          }
          d1[["rat"]] <- round(
            d1[[2]] / d1[d1[[1]] == gblnk, ][[2]],
            digits = 3
          )
          d1 <- d1[d1[[1]] != gblnk, ]
          d1 <- max(d1[["rat"]])
          return(d1) # nolint
        }
      )
    )
  }
  # Identify adducts and isotopes based on mzrt and peak height correlation
  ## Adapted from MS-FLO
  if (Sys.info()[["sysname"]] != "Windows") {
    d_check2 <- dplyr::bind_rows(parallel::mclapply(
      mc.cores = ceiling(parallel::detectCores() / 2),
      seq.int(1, nrow(ld[["feat"]]), 1),
      function(i) {
        # Return features within RT window
        d1 <- ld[["feat"]][
          ld[["feat"]][["Average.Rt.min."]] <
            ld[["feat"]][i, ][["Average.Rt.min."]] + tolrt2 &
            ld[["feat"]][["Average.Rt.min."]] >
              ld[["feat"]][i, ][["Average.Rt.min."]] - tolrt2,
        ]
        # Calculate relative correlation with selected feature
        d1 <- dplyr::select(
          dplyr::mutate(
            d1,
            "cor" = unlist(
              lapply(
                seq.int(1, nrow(d1), 1),
                function(z) {
                  round(cor.test(
                    ld[["data"]][[
                      paste(
                        "Alignment.ID_",
                        ld[["feat"]][i, ][["Alignment.ID"]],
                        sep = ""
                      )
                    ]],
                    ld[["data"]][[
                      paste(
                        "Alignment.ID_",
                        d1[["Alignment.ID"]][[z]],
                        sep = ""
                      )
                    ]],
                    method = "pearson"
                  )[["estimate"]], digits = 2)
                }
              )
            )
          ),
          "cor", everything()
        )
        ## filter alignment IDs with low correlation
        d1 <- d1[
          d1[["cor"]] > 0.95 &
            d1[["Alignment.ID"]] %in%
              ld[["feat"]][i, ][["Alignment.ID"]] == FALSE,
        ]
        ## Return original alignment ID if no matches present
        if (nrow(d1) == 0) {
          print(
            paste(
              "No features correlated with alignment ID ",
              ld[["feat"]][i, ][["Alignment.ID"]],
              ".", sep = ""
            )
          )
          d2 <- ld[["feat"]][i, ]
          d2 <- dplyr::select(
            dplyr::mutate(
              d2,
              "Alignment.ID.orig" = d2[["Alignment.ID"]],
              "adduct.match.ID" = NA,
              "adduct.name" = NA
            ),
            "Alignment.ID.orig",
            "adduct.match.ID", "adduct.name",
            everything()
          )
        }
        # Return features with matches to other adducts
        if (t2 == "pos" && nrow(d1) > 0) {
          print(
            paste(
              "Potential adducts found for alignment ID ",
              ld[["feat"]][i, ][["Alignment.ID"]],
              ". Performing adduct and isotope detection...",
              sep = ""
            )
          )
          ## adduct list
          adduct_pos <- data.frame(
            "Adduct" = c(
              "[M+H]+",
              "[M+NH4]+",
              "[M+Na]+",
              "[M+K]+",
              "[M+H-H2O]+",
              "[M+H-2H2O]+",
              "[2M+H]+",
              "[2M+NH4]+",
              "[2M+Na]+",
              "[2M+K]+"
            ),
            "mz" = c(
              1.00782503207,
              18.03437413,
              22.9897692809,
              38.96370668,
              -17.00273964793,
              -35.01330432793,
              1.00782503207,
              18.03437413,
              22.9897692809,
              38.96370668
            )
          )
          adduct_pos <- dplyr::bind_rows(
            adduct_pos,
            data.frame(
              "Adduct" = paste(adduct_pos[["Adduct"]], "C13", sep = "_"),
              "mz" = adduct_pos[["mz"]] + 1.0034
            )
          )
          # Find potential adduct matches within m/z threshold
          d2 <- dplyr::bind_rows(
            lapply(
              seq.int(1, nrow(d1), 1),
              function(y) {
                mz3 <- tryCatch(
                  {
                    d2 <- d1[y, ]
                    ## Calculate neutral mass
                    mz1 <- data.frame(
                      "Alignment.ID" = d2[["Alignment.ID"]],
                      "Metabolite.name" = d2[["Metabolite.name"]],
                      "mz.neutral" = ifelse(
                        grepl("2M", d2[["Adduct.type"]]),
                        (d2[["Average.Mz"]] -
                          adduct_pos[
                            adduct_pos[["Adduct"]] %in% d2[["Adduct.type"]],
                          ][["mz"]]
                        ) * 2,
                        d2[["Average.Mz"]] -
                          adduct_pos[
                            adduct_pos[["Adduct"]] %in% d2[["Adduct.type"]],
                          ][["mz"]]
                      )
                    )
                    ## Calculate adduct masses
                    mz2 <- data.frame(
                      adduct_pos,
                      "mz.adduct" = adduct_pos[["mz"]] + mz1[["mz.neutral"]]
                    )
                    mz2[grepl("2M", mz2[["Adduct"]]), "mz.adduct"] <-
                      adduct_pos[
                        adduct_pos[["Adduct"]] %in%
                        mz2[grepl("2M", mz2[["Adduct"]]), ][["Adduct"]],
                      ][["mz"]] + mz1[["mz.neutral"]] * 2
                    ## Return list of m/z matches
                    mz3 <- dplyr::bind_rows(
                      lapply(
                        seq.int(1, nrow(mz2), 1),
                        function(j) {
                          ad1 <- mz2[j, ]
                          ad2 <- d1[
                            d1[["Average.Mz"]] < ad1[["mz.adduct"]] + tolmz2 &
                              d1[["Average.Mz"]] > ad1[["mz.adduct"]] - tolmz2,
                          ]
                          ad3 <- dplyr::select(
                            dplyr::mutate(
                              ad2,
                              "adduct.match.ID" = paste(
                                rep(d2[["Alignment.ID"]], nrow(ad2)),
                                ad2[["Alignment.ID"]], sep = "_"
                              ),
                              "adduct.name" = mz2[j, ][["Adduct"]]
                            ),
                            "adduct.match.ID", "adduct.name", everything()
                          )
                          return(ad3) # nolint
                        }
                      )
                    )
                    mz3 <- dplyr::select(
                      dplyr::mutate(
                        mz3,
                        "Alignment.ID.orig" = d2[["Alignment.ID"]]
                      ),
                      "Alignment.ID.orig", everything()
                    )
                  },
                  error = function(e) {
                    print(
                      paste(
                        "error in adduct/isotope detection for alignment ID",
                        d1[y, ][["Alignment.ID"]]
                      )
                    )
                    mz3 <- dplyr::select(
                      dplyr::mutate(
                        d1[y, ],
                        "Alignment.ID.orig" = d1[y, ][["Alignment.ID"]],
                        "adduct.match.ID" = NA,
                        "adduct.name" = NA
                      ),
                      "Alignment.ID.orig",
                      "adduct.match.ID", "adduct.name",
                      everything()
                    )
                    mz3 <- mz3[!duplicated(dplyr::select(mz3, -c("cor"))), ]
                    return(mz3) # nolint
                  }
                )
                return(mz3) # nolint
              }
            )
          )
          d2a <- d2[, 1:10]
          d2a <- !duplicated(d2a)
          d2 <- d2[d2a, ]
        }
        if (t2 == "neg" && nrow(d1) > 0) {
          print(
            paste(
              "Potential adducts found for alignment ID ",
              ld[["feat"]][i, ][["Alignment.ID"]],
              ". Performing adduct and isotope detection...",
              sep = ""
            )
          )
          ## adduct list
          adduct_neg <- data.frame(
            "Adduct" = c(
              "[M-H]-",
              "[M-H2O-H]-",
              "[M+HCOO]-",
              "[M+CH3COO]-",
              "[2M-H]-",
              "[2M+FA-H]-",
              "[2M+Hac-H]-"
            ),
            "mz" = c(
              -1.00782503207,
              -19.0184,
              44.9977,
              59.0133,
              -1.00782503207,
              44.9977,
              59.0133
            )
          )
          adduct_neg <- dplyr::bind_rows(
            adduct_neg,
            data.frame(
              "Adduct" = paste(adduct_neg[["Adduct"]], "C13", sep = "_"),
              "mz" = adduct_neg[["mz"]] + 1.0034
            )
          )
          # Find potential adduct matches within m/z threshold
          d2 <- dplyr::bind_rows(
            lapply(
              seq.int(1, nrow(d1), 1),
              function(y) {
                mz3 <- tryCatch(
                  {
                    d2 <- d1[y, ]
                    ## Calculate neutral mass
                    mz1 <- data.frame(
                      "Alignment.ID" = d2[["Alignment.ID"]],
                      "Metabolite.name" = d2[["Metabolite.name"]],
                      "mz.neutral" = ifelse(
                        grepl("2M", d2[["Adduct.type"]]),
                        (d2[["Average.Mz"]] -
                          adduct_neg[
                            adduct_neg[["Adduct"]] %in% d2[["Adduct.type"]],
                          ][["mz"]]
                        ) * 2,
                        d2[["Average.Mz"]] -
                          adduct_neg[
                            adduct_neg[["Adduct"]] %in% d2[["Adduct.type"]],
                          ][["mz"]]
                      )
                    )
                    ## Calculate adduct masses
                    mz2 <- data.frame(
                      adduct_neg,
                      "mz.adduct" = adduct_neg[["mz"]] + mz1[["mz.neutral"]]
                    )
                    mz2[grepl("2M", mz2[["Adduct"]]), "mz.adduct"] <-
                      adduct_neg[
                        adduct_neg[["Adduct"]] %in%
                        mz2[grepl("2M", mz2[["Adduct"]]), ][["Adduct"]],
                      ][["mz"]] + mz1[["mz.neutral"]] * 2
                    ## Return list of m/z matches
                    mz3 <- dplyr::bind_rows(
                      lapply(
                        seq.int(1, nrow(mz2), 1),
                        function(j) {
                          ad1 <- mz2[j, ]
                          ad2 <- d1[
                            d1[["Average.Mz"]] < ad1[["mz.adduct"]] + tolmz2 &
                              d1[["Average.Mz"]] > ad1[["mz.adduct"]] - tolmz2,
                          ]
                          ad3 <- dplyr::select(
                            dplyr::mutate(
                              ad2,
                              "adduct.match.ID" = paste(
                                rep(d2[["Alignment.ID"]], nrow(ad2)),
                                ad2[["Alignment.ID"]], sep = "_"
                              ),
                              "adduct.name" = mz2[j, ][["Adduct"]]
                            ),
                            "adduct.match.ID", "adduct.name", everything()
                          )
                          return(ad3) # nolint
                        }
                      )
                    )
                    mz3 <- dplyr::select(
                      dplyr::mutate(
                        mz3,
                        "Alignment.ID.orig" = d2[["Alignment.ID"]]
                      ),
                      "Alignment.ID.orig", everything()
                    )
                  },
                  error = function(e) {
                    print(
                      paste(
                        "error in adduct/isotope detection for alignment ID",
                        d1[y, ][["Alignment.ID"]]
                      )
                    )
                    mz3 <- dplyr::select(
                      dplyr::mutate(
                        d1[y, ],
                        "Alignment.ID.orig" = d1[y, ][["Alignment.ID"]],
                        "adduct.match.ID" = NA,
                        "adduct.name" = NA
                      ),
                      "Alignment.ID.orig",
                      "adduct.match.ID", "adduct.name",
                      everything()
                    )
                    mz3 <- mz3[!duplicated(dplyr::select(mz3, -c("cor"))), ]
                    return(mz3) # nolint
                  }
                )
                return(mz3) # nolint
              }
            )
          )
          d2a <- d2[, 1:10]
          d2a <- !duplicated(d2a)
          d2 <- d2[d2a, ]
        }
        return(d2) # nolint
      }
    ))
    ## Remove duplicate IDs
    d_check2 <- d_check2[!duplicated(dplyr::select(d_check2, -c("cor"))), ]
  }
  if (Sys.info()[["sysname"]] == "Windows") {
    d_check2 <- dplyr::bind_rows(lapply(
      seq.int(1, nrow(ld[["feat"]]), 1),
      function(i) {
        # Return features within RT window
        d1 <- ld[["feat"]][
          ld[["feat"]][["Average.Rt.min."]] <
            ld[["feat"]][i, ][["Average.Rt.min."]] + tolrt2 &
            ld[["feat"]][["Average.Rt.min."]] >
              ld[["feat"]][i, ][["Average.Rt.min."]] - tolrt2,
        ]
        # Calculate relative correlation with selected feature
        d1 <- dplyr::select(
          dplyr::mutate(
            d1,
            "cor" = unlist(
              lapply(
                seq.int(1, nrow(d1), 1),
                function(z) {
                  round(cor.test(
                    ld[["data"]][[
                      paste(
                        "Alignment.ID_",
                        ld[["feat"]][i, ][["Alignment.ID"]],
                        sep = ""
                      )
                    ]],
                    ld[["data"]][[
                      paste(
                        "Alignment.ID_",
                        d1[["Alignment.ID"]][[z]],
                        sep = ""
                      )
                    ]],
                    method = "pearson"
                  )[["estimate"]], digits = 2)
                }
              )
            )
          ),
          "cor", everything()
        )
        ## filter alignment IDs with low correlation
        d1 <- d1[
          d1[["cor"]] > 0.95 &
            d1[["Alignment.ID"]] %in%
              ld[["feat"]][i, ][["Alignment.ID"]] == FALSE,
        ]
        ## Return original alignment ID if no matches present
        if (nrow(d1) == 0) {
          print(
            paste(
              "No features correlated with alignment ID ",
              ld[["feat"]][i, ][["Alignment.ID"]],
              ".", sep = ""
            )
          )
          d2 <- ld[["feat"]][i, ]
          d2 <- dplyr::select(
            dplyr::mutate(
              d2,
              "Alignment.ID.orig" = d2[["Alignment.ID"]],
              "adduct.match.ID" = NA,
              "adduct.name" = NA
            ),
            "Alignment.ID.orig",
            "adduct.match.ID", "adduct.name",
            everything()
          )
        }
        # Return features with matches to other adducts
        if (t2 == "pos" && nrow(d1) > 0) {
          print(
            paste(
              "Potential adducts found for alignment ID ",
              ld[["feat"]][i, ][["Alignment.ID"]],
              ". Performing adduct and isotope detection...",
              sep = ""
            )
          )
          ## adduct list
          adduct_pos <- data.frame(
            "Adduct" = c(
              "[M+H]+",
              "[M+NH4]+",
              "[M+Na]+",
              "[M+K]+",
              "[M+H-H2O]+",
              "[M+H-2H2O]+",
              "[2M+H]+",
              "[2M+NH4]+",
              "[2M+Na]+",
              "[2M+K]+"
            ),
            "mz" = c(
              1.00782503207,
              18.03437413,
              22.9897692809,
              38.96370668,
              -17.00273964793,
              -35.01330432793,
              1.00782503207,
              18.03437413,
              22.9897692809,
              38.96370668
            )
          )
          adduct_pos <- dplyr::bind_rows(
            adduct_pos,
            data.frame(
              "Adduct" = paste(adduct_pos[["Adduct"]], "C13", sep = "_"),
              "mz" = adduct_pos[["mz"]] + 1.0034
            )
          )
          # Find potential adduct matches within m/z threshold
          d2 <- dplyr::bind_rows(
            lapply(
              seq.int(1, nrow(d1), 1),
              function(y) {
                mz3 <- tryCatch(
                  {
                    d2 <- d1[y, ]
                    ## Calculate neutral mass
                    mz1 <- data.frame(
                      "Alignment.ID" = d2[["Alignment.ID"]],
                      "Metabolite.name" = d2[["Metabolite.name"]],
                      "mz.neutral" = ifelse(
                        grepl("2M", d2[["Adduct.type"]]),
                        (d2[["Average.Mz"]] -
                          adduct_pos[
                            adduct_pos[["Adduct"]] %in% d2[["Adduct.type"]],
                          ][["mz"]]
                        ) * 2,
                        d2[["Average.Mz"]] -
                          adduct_pos[
                            adduct_pos[["Adduct"]] %in% d2[["Adduct.type"]],
                          ][["mz"]]
                      )
                    )
                    ## Calculate adduct masses
                    mz2 <- data.frame(
                      adduct_pos,
                      "mz.adduct" = adduct_pos[["mz"]] + mz1[["mz.neutral"]]
                    )
                    mz2[grepl("2M", mz2[["Adduct"]]), "mz.adduct"] <-
                      adduct_pos[
                        adduct_pos[["Adduct"]] %in%
                        mz2[grepl("2M", mz2[["Adduct"]]), ][["Adduct"]],
                      ][["mz"]] + mz1[["mz.neutral"]] * 2
                    ## Return list of m/z matches
                    mz3 <- dplyr::bind_rows(
                      lapply(
                        seq.int(1, nrow(mz2), 1),
                        function(j) {
                          ad1 <- mz2[j, ]
                          ad2 <- d1[
                            d1[["Average.Mz"]] < ad1[["mz.adduct"]] + tolmz2 &
                              d1[["Average.Mz"]] > ad1[["mz.adduct"]] - tolmz2,
                          ]
                          ad3 <- dplyr::select(
                            dplyr::mutate(
                              ad2,
                              "adduct.match.ID" = paste(
                                rep(d2[["Alignment.ID"]], nrow(ad2)),
                                ad2[["Alignment.ID"]], sep = "_"
                              ),
                              "adduct.name" = mz2[j, ][["Adduct"]]
                            ),
                            "adduct.match.ID", "adduct.name", everything()
                          )
                          return(ad3) # nolint
                        }
                      )
                    )
                    mz3 <- dplyr::select(
                      dplyr::mutate(
                        mz3,
                        "Alignment.ID.orig" = d2[["Alignment.ID"]]
                      ),
                      "Alignment.ID.orig", everything()
                    )
                  },
                  error = function(e) {
                    print(
                      paste(
                        "error in adduct/isotope detection for alignment ID",
                        d1[y, ][["Alignment.ID"]]
                      )
                    )
                    mz3 <- dplyr::select(
                      dplyr::mutate(
                        d1[y, ],
                        "Alignment.ID.orig" = d1[y, ][["Alignment.ID"]],
                        "adduct.match.ID" = NA,
                        "adduct.name" = NA
                      ),
                      "Alignment.ID.orig",
                      "adduct.match.ID", "adduct.name",
                      everything()
                    )
                    mz3 <- mz3[!duplicated(dplyr::select(mz3, -c("cor"))), ]
                    return(mz3) # nolint
                  }
                )
                return(mz3) # nolint
              }
            )
          )
          d2a <- d2[, 1:10]
          d2a <- !duplicated(d2a)
          d2 <- d2[d2a, ]
        }
        if (t2 == "neg" && nrow(d1) > 0) {
          print(
            paste(
              "Potential adducts found for alignment ID ",
              ld[["feat"]][i, ][["Alignment.ID"]],
              ". Performing adduct and isotope detection...",
              sep = ""
            )
          )
          ## adduct list
          adduct_neg <- data.frame(
            "Adduct" = c(
              "[M-H]-",
              "[M-H2O-H]-",
              "[M+HCOO]-",
              "[M+CH3COO]-",
              "[2M-H]-",
              "[2M+FA-H]-",
              "[2M+Hac-H]-"
            ),
            "mz" = c(
              -1.00782503207,
              -19.0184,
              44.9977,
              59.0133,
              -1.00782503207,
              44.9977,
              59.0133
            )
          )
          adduct_neg <- dplyr::bind_rows(
            adduct_neg,
            data.frame(
              "Adduct" = paste(adduct_neg[["Adduct"]], "C13", sep = "_"),
              "mz" = adduct_neg[["mz"]] + 1.0034
            )
          )
          # Find potential adduct matches within m/z threshold
          d2 <- dplyr::bind_rows(
            lapply(
              seq.int(1, nrow(d1), 1),
              function(y) {
                mz3 <- tryCatch(
                  {
                    d2 <- d1[y, ]
                    ## Calculate neutral mass
                    mz1 <- data.frame(
                      "Alignment.ID" = d2[["Alignment.ID"]],
                      "Metabolite.name" = d2[["Metabolite.name"]],
                      "mz.neutral" = ifelse(
                        grepl("2M", d2[["Adduct.type"]]),
                        (d2[["Average.Mz"]] -
                          adduct_neg[
                            adduct_neg[["Adduct"]] %in% d2[["Adduct.type"]],
                          ][["mz"]]
                        ) * 2,
                        d2[["Average.Mz"]] -
                          adduct_neg[
                            adduct_neg[["Adduct"]] %in% d2[["Adduct.type"]],
                          ][["mz"]]
                      )
                    )
                    ## Calculate adduct masses
                    mz2 <- data.frame(
                      adduct_neg,
                      "mz.adduct" = adduct_neg[["mz"]] + mz1[["mz.neutral"]]
                    )
                    mz2[grepl("2M", mz2[["Adduct"]]), "mz.adduct"] <-
                      adduct_neg[
                        adduct_neg[["Adduct"]] %in%
                        mz2[grepl("2M", mz2[["Adduct"]]), ][["Adduct"]],
                      ][["mz"]] + mz1[["mz.neutral"]] * 2
                    ## Return list of m/z matches
                    mz3 <- dplyr::bind_rows(
                      lapply(
                        seq.int(1, nrow(mz2), 1),
                        function(j) {
                          ad1 <- mz2[j, ]
                          ad2 <- d1[
                            d1[["Average.Mz"]] < ad1[["mz.adduct"]] + tolmz2 &
                              d1[["Average.Mz"]] > ad1[["mz.adduct"]] - tolmz2,
                          ]
                          ad3 <- dplyr::select(
                            dplyr::mutate(
                              ad2,
                              "adduct.match.ID" = paste(
                                rep(d2[["Alignment.ID"]], nrow(ad2)),
                                ad2[["Alignment.ID"]], sep = "_"
                              ),
                              "adduct.name" = mz2[j, ][["Adduct"]]
                            ),
                            "adduct.match.ID", "adduct.name", everything()
                          )
                          return(ad3) # nolint
                        }
                      )
                    )
                    mz3 <- dplyr::select(
                      dplyr::mutate(
                        mz3,
                        "Alignment.ID.orig" = d2[["Alignment.ID"]]
                      ),
                      "Alignment.ID.orig", everything()
                    )
                  },
                  error = function(e) {
                    print(
                      paste(
                        "error in adduct/isotope detection for alignment ID",
                        d1[y, ][["Alignment.ID"]]
                      )
                    )
                    mz3 <- dplyr::select(
                      dplyr::mutate(
                        d1[y, ],
                        "Alignment.ID.orig" = d1[y, ][["Alignment.ID"]],
                        "adduct.match.ID" = NA,
                        "adduct.name" = NA
                      ),
                      "Alignment.ID.orig",
                      "adduct.match.ID", "adduct.name",
                      everything()
                    )
                    mz3 <- mz3[!duplicated(dplyr::select(mz3, -c("cor"))), ]
                    return(mz3) # nolint
                  }
                )
                return(mz3) # nolint
              }
            )
          )
          d2a <- d2[, 1:10]
          d2a <- !duplicated(d2a)
          d2 <- d2[d2a, ]
        }
        return(d2) # nolint
      }
    ))
    ## Remove duplicate IDs
    d_check2 <- d_check2[!duplicated(dplyr::select(d_check2, -c("cor"))), ]
  }
  # Add adduct flag to summary list
  d_check3 <- dplyr::count(
    dplyr::group_by(
      d_check2,
      d_check2[["Alignment.ID.orig"]]
    ),
    !is.na(d_check2[["adduct.match.ID"]])
  )
  d_check3 <- d_check3[d_check3[[2]] == TRUE, ][[1]]
  d_check[["check.adducts"]] <- ifelse(
    d_check[["Alignment.ID"]] %in% d_check3,
    TRUE,
    FALSE
  )
  # Identify characteristic fragments in lipidomics data
  if (t1 == "lipidomics") {
    if (t2 == "pos") {

    }
    if (t2 == "neg") {

    }
  }
  # Remove alignment IDs based on flagged criteria
  ## includes blank removal, sample presence, and fill %
  ## and retains internal standards
  d_remove <- d_check[
    d_check[["Blank.ratio.max"]] < filt1 &
      d_check[["iSTD"]] == FALSE &
      d_check[["Prop.max.per.group"]] < 0.50,
  ][["Alignment.ID"]]
  d_out <- d_check2[
    d_check2[["Alignment.ID"]] %in% d_remove == FALSE &
      d_check2[["Fill.."]] > min_fill,
  ]
  d_out_int <- ld[["data"]][
    names(ld[["data"]]) %in%
      paste(
        "Alignment.ID_", unique(d_out[["Alignment.ID"]]), sep = ""
      )
  ]
  ## Return output data and summary
  list_output <- list(
    "data_int" = d_out_int,
    "data_samples" = ld[["samp"]],
    "data_preproc_feat" = d_out,
    "data_summary" = d_check
  )
  return(list_output)
}

#' Post-processing of Curated Lipidomics/Metabolomics Data
#'
#' Combines flagged adducts in a curated feature list by alignment
#' ID and filters input data for downstream normalization and analysis.
#' The curated feature list should be free of duplicate compounds
#' across both ionization modes (if applicable).
#'
#' @param dat Curated feature list, which is formatted as either a .csv
#' or tab-delimited .txt file. A column named "Mode" indicating the
#' ionization mode should be present in the file.
#' @param dat_in A RDS file generated by ms_preproc containing peak intensity
#' and sample information. If two file names are provided (e.g. pos and neg
#' ionization modes), then the data are merged into a single
#' SummarizedExperiment.
#' @param col_id The name of the column containing the original alignment IDs
#' for each feature in an exported MS-DIAL feature list. MUST match alignment
#' IDs in the peak intensity matrix to function properly.
#' @param col_ad Alignment ID column name; combines all features with duplicate
#' IDs, assuming that duplicate IDs correspond to adducts of the same feature.
#' @param samp_id Sample ID column name.
#' @param add_metadata Add additional metadata information to the output object.
#' These data should be provided as named list elements and are useful for
#' simplifying the use of COMb functions. See examples of common metadata below.
#' @return A SummarizedExperiment containing the data matrix
#' with the final feature and sample metadata for the specified ionization mode.
#' @examples
#'
#' # postproc <- ms_postproc(
#' #   dat = "data/03postprocessing/Chase_LSAE_curated.csv",
#' #   ionmode = "pos",
#' #   dat_in = "data/02preprocessing/Chase_LSAE_preprocessed_pos.rds"
#' # )
#'
#' @import dplyr
#' @import magrittr
#' @export
ms_postproc <- function(
  dat,
  dat_in,
  col_id = "Alignment.ID.original",
  col_ad = "Alignment.ID.matched",
  samp_id = "ID",
  add_metadata = NULL
) {
  # Set parameters
  ## Curated feature list
  d1 <- dat
  if (grepl(".csv", dat)) {
    d1 <- read.csv(d1, quote = "")
  }
  if (grepl(".txt", dat)) {
    d1 <- read.table(d1, header = TRUE, sep = "\t", quote = "")
  }
  ## Ionization mode
  m1 <- unique(d1[["Mode"]])
  ## Original alignment ID
  c1 <- col_id
  ## Adduct alignment ID
  c2 <- col_ad
  ## Input RDS
  d2 <- dat_in
  # If both modes provided
  if (length(d2) == 2) {
    # Load each data set
    d2 <- setNames(lapply(d2, function(x) readRDS(x)), m1)
    # Perform peak filtering and adduct combine
    d_out <- setNames(lapply(
      seq.int(1, 2, 1),
      function(y) {
        # Filter peak intensity matrix and combine adducts
        ## Filter peak data
        d2filt <- d2[[y]][["data_int"]][
          names(d2[[y]][["data_int"]]) %in%
            paste(
              "Alignment.ID_",
              d1[d1[["Mode"]] == names(d2)[[y]], ][[c1]],
              sep = ""
            )
        ]
        ## Transpose, combine with feature list, and aggregate by adduct ID
        d2filt2 <- setNames(data.frame(
          "Feature.name" = names(d2filt),
          "Alignment.ID.original" = as.numeric(
            gsub("Alignment.ID_", "", names(d2filt))
          )
        ), c("Feature.name", c1))
        d2filt2 <- dplyr::left_join(
          d2filt2,
          d1[d1[["Mode"]] == names(d2)[[y]], ],
          by = c1
        )
        d2filt <- as.data.frame(t(d2filt))
        d2filt <- dplyr::bind_cols(
          c(
            setNames(as.data.frame(aggregate(
              d2filt[[1]],
              list(d2filt2[[c2]]),
              function(x) sum(x)
            )[, 1]), c1),
            lapply(
              seq.int(1, ncol(d2filt), 1),
              function(i) {
                setNames(as.data.frame(aggregate(
                  d2filt[[i]],
                  list(d2filt2[[c2]]),
                  function(x) sum(x)
                )[[2]]), i)
              }
            )
          )
        )
        ## Merge data matrix with feature and sample metadata
        d2filt2 <- d2filt2[d2filt2[[c1]] %in% d2filt[[c1]], ]
        d2filt <- d2filt[d2filt[[c1]] %in% d2filt2[[c1]], ]
        ### Combine as list
        d3 <- list(
          # data matrix with features as rows and samples as columns
          "data" = d2filt[-1],
          # feature metadata
          "anno" = d2filt2,
          # sample metadata
          "meta" = d2[[y]][["data_samples"]]
        )
        return(d3) # nolint
      }
    ), names(d2))
    # Combine lists before converting to SummarizedExperiment
    d_out <- list(
      "data" = rbind(
        d_out[[1]][["data"]],
        d_out[[2]][["data"]]
      ),
      "meta" = d_out[[1]][["meta"]],
      "anno" = rbind(
        d_out[[1]][["anno"]],
        d_out[[2]][["anno"]]
      )
    )
    rownames(d_out[["data"]]) <- NULL
    colnames(d_out[["data"]]) <- NULL
    # Convert to SummarizedExperiment
    d4 <- SummarizedExperiment::SummarizedExperiment(
      assays = list("raw" = as.matrix(d_out[["data"]])),
      rowData = d_out[["anno"]],
      colData = d_out[["meta"]],
      metadata = list(
        "Ionization mode" = names(d2)
      )
    )
    rownames(d4) <- d_out[["anno"]][["Feature.name"]]
    colnames(d4) <- d_out[["meta"]][[samp_id]]
    if (!is.null(add_metadata)) {
      metadata(d4) <- c( # nolint
        metadata(d4),
        add_metadata
      )
    }
    # Export data
    print(paste("Final data for both ionization modes contains:", sep = ""))
    print(paste(nrow(d4), " unique compounds, including internal standards, and", sep = "")) # nolint
    print(paste(ncol(d4), " samples, including QCs.", sep = ""))
  }
  # If single ionization mode provided
  if (length(d2) == 1) {
    # Load RDS
    d2 <- readRDS(d2)
    # Filter peak intensity matrix and combine adducts
    ## Filter peak data
    d2filt <- d2[["data_int"]][
      names(d2[["data_int"]]) %in%
        paste("Alignment.ID_", d1[d1[["Mode"]] == m1, ][[c1]], sep = "")
    ]
    ## Transpose, combine with feature list, and aggregate by adduct ID
    d2filt2 <- setNames(data.frame(
      "Feature.name" = names(d2filt),
      "Alignment.ID.original" = as.numeric(
        gsub("Alignment.ID_", "", names(d2filt))
      )
    ), c("Feature.name", c1))
    d2filt2 <- dplyr::left_join(
      d2filt2,
      d1[d1[["Mode"]] == m1, ],
      by = c1
    )
    d2filt <- as.data.frame(t(d2filt))
    d2filt <- dplyr::bind_cols(
      c(
        setNames(as.data.frame(aggregate(
          d2filt[[1]],
          list(d2filt2[[c2]]),
          function(x) sum(x)
        )[, 1]), c1),
        lapply(
          seq.int(1, ncol(d2filt), 1),
          function(i) {
            setNames(as.data.frame(aggregate(
              d2filt[[i]],
              list(d2filt2[[c2]]),
              function(x) sum(x)
            )[[2]]), i)
          }
        )
      )
    )
    ## Merge data matrix with feature and sample metadata
    d2filt2 <- d2filt2[d2filt2[[c1]] %in% d2filt[[c1]], ]
    d2filt <- d2filt[d2filt[[c1]] %in% d2filt2[[c1]], ]
    ### Combine as list
    d3 <- list(
      # data matrix with features as rows and samples as columns
      "data" = d2filt[-1],
      # feature metadata
      "anno" = d2filt2,
      # sample metadata
      "meta" = d2[["data_samples"]]
    )
    rownames(d3[["data"]]) <- NULL
    colnames(d3[["data"]]) <- NULL
    # Convert to SummarizedExperiment
    d4 <- SummarizedExperiment::SummarizedExperiment(
      assays = list("raw" = as.matrix(d3[["data"]])),
      rowData = d3[["anno"]],
      colData = d3[["meta"]],
      metadata = list(
        "Ionization mode" = unique(d3[["anno"]][["Mode"]])
      )
    )
    rownames(d4) <- d3[["anno"]][["Feature.name"]]
    colnames(d4) <- d3[["meta"]][[samp_id]]
    if (!is.null(add_metadata)) {
      metadata(d4) <- c( # nolint
        metadata(d4),
        add_metadata
      )
    }
    # Export data
    print(paste("Final data for ", m1, " ionization mode contains:", sep = ""))
    print(paste(nrow(d3[[1]]), " unique compounds, including internal standards, and", sep = "")) # nolint
    print(paste(ncol(d3[[1]]), " samples, including QCs.", sep = ""))
  }
  return(d4)
}
