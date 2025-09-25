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
