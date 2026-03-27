#' Group-wise Fold Change Function
#'
#' Calculates the log2 fold changes for all
#' possible combinations of selected metadata variables.
#'
#' @param exp1 A SummarizedExperiment.
#' @param asy Assay to use for fold change calculation.
#' @param mat_type Input data type, either "norm" or "scaled." If
#' "scaled," fold changes are calculated as the difference between
#' two group averages rather than by dividing normalized compound
#' average intensities per group.
#' @param col_grp Vector of metadata variables for splitting data.
#' @param col_class Compound class variable name; only required if
#' fc_method is "class" or "HiVE".
#' @param col_lab Compound name column.
#' @param fc_method Fold change calculation mode; either "standard",
#' which is the most common, "class", or "HiVE"
#' @param comp_grp Only applicable if fc_method is "HiVE"; Should
#' group fold changes be calculated? If FALSE, returns ratios
#' without aggregation by treatment group.
#' @return A data frame containing the group-wise
#' fold changes for each compound or class.
#' @examples
#'
#' # ms_stat_fc(exp1 = data2)
#'
#' @export
ms_stat_fc <- function( # nolint
  exp1,
  asy = "norm",
  mat_type = "norm",
  col_grp = "Group",
  col_class = "Class",
  col_lab = "Label",
  fc_method = "standard",
  comp_grp = TRUE
) {
  #---- Setup ----
  # Load summarizedexperiment and associated data
  ## Input matrix
  d1 <- as.data.frame(t(assay(exp1, asy)))
  mtype <- mat_type
  ## column data
  an1 <- as.data.frame(rowData(exp1))
  ## row data
  md1 <- as.data.frame(colData(exp1))
  ## metadata variable names
  md_var1 <- col_grp
  cc1 <- col_class
  cn1 <- col_lab
  #---- Functions ----
  # Define fold change function
  if (mtype == "norm") {
    fun.comb <- function(x, y) {x / y} # nolint
  }
  if (mtype == "scaled") {
    fun.comb <- function(x, y) {x - y} # nolint
  }
  ## group combinations
  fun.grpcomb <- function() { # nolint
    # Group Combinations
    var_comb <- dplyr::bind_rows(lapply(
      seq.int(0, length(md_var1) - 1, 1),
      function(x) {
        cb1 <- combn(md_var1, x + 1)
        cb1 <- as.data.frame(t(
          dplyr::bind_cols(
            lapply(
              as.data.frame(cb1),
              function(y) paste(y, collapse = ":")
            )
          )
        ))
        return(cb1) # nolint
      }
    ))
    if(nrow(var_comb) == 1) { # nolint
      cb1 <- paste(unlist(strsplit(var_comb[1, 1], ":")), sep = ", ")
      # Change variable to factor
      cb2 <- setNames(
        as.data.frame(
          as.character(md1[[md_var1]])
        ),
        c(cb1)
      )
      cb2 <- as.data.frame(unique(cb2))
      # Combine columns
      cb2[["comb"]] <- unlist(
        lapply(
          as.data.frame(t(cb2)),
          function(z) paste(z, collapse = ":")
        )
      )
      # Determine combinations
      fold_comb <- as.data.frame(
        t(
          unique(
            combn(
              sort(cb2[["comb"]], decreasing = TRUE),
              2
            )
          )
        )
      )
    }
    if(nrow(var_comb) > 1) { # nolint
      fold_comb <- dplyr::bind_rows(
        lapply(
          seq.int(1, nrow(var_comb), 1),
          function(x) {
            cb1 <- paste(unlist(strsplit(var_comb[x, 1], ":")), sep = ", ")
            # Change variable to factor
            cb2 <- setNames(
              as.data.frame(
                lapply(
                  cb1,
                  function(y) as.character(md1[, y])
                )
              ),
              c(cb1)
            )
            cb2 <- unique(cb2)
            # Combine columns
            cb2[["comb"]] <- unlist(
              lapply(
                as.data.frame(t(cb2)),
                function(z) paste(z, collapse = ":")
              )
            )
            # Determine combinations
            cb2 <- as.data.frame(
              t(
                unique(
                  combn(
                    sort(cb2[["comb"]], decreasing = TRUE),
                    2
                  )
                )
              )
            )
            return(cb2) # nolint
          }
        )
      )
    }
    return(list("var_comb" = var_comb, "fold_comb" = fold_comb)) # nolint
  }
  # Define mean group function
  fun.fcmean <- function() { # nolint
    if (Sys.info()[["sysname"]] != "Windows") {
      fold_mean <- dplyr::bind_rows(
        setNames(
          parallel::mclapply(
            mc.cores = ceiling(parallel::detectCores() / 4),
            seq.int(1, nrow(fold_comb[[1]]), 1),
            function(x) {
              cm1 <- paste(
                unlist(strsplit(fold_comb[[1]][x, 1], ":")), sep = ", "
              )
              cm2 <- aggregate(
                d1,
                lapply(cm1, function(z) as.character(md1[[z]])),
                function(y) mean(y)
              )
              cm2 <- setNames(
                data.frame(
                  "Group" = unlist(
                    lapply(
                      as.data.frame(t(cm2[, 1:length(cm1)])), # nolint
                      function(z) paste(z, collapse = ":")
                    )
                  ),
                  cm2[, (length(cm1) + 1):ncol(cm2)]
                ),
                c("Group", names(cm2[, (length(cm1) + 1):ncol(cm2)]))
              )
              return(cm2) # nolint
            }
          ),
          fold_comb[[1]][[1]]
        )
      )
    }
    if (Sys.info()[["sysname"]] == "Windows") {
      fold_mean <- dplyr::bind_rows(
        setNames(
          lapply(
            seq.int(1, nrow(fold_comb[[1]]), 1), # nolint
            function(x) {
              cm1 <- paste(unlist(strsplit(fold_comb[[1]][x, 1], ":")), sep = ", ") # nolint
              cm2 <- aggregate(
                d1,
                lapply(cm1, function(z) as.character(md1[[z]])),
                function(y) mean(y)
              )
              cm2 <- setNames(
                data.frame(
                  "Group" = unlist(
                    lapply(
                      as.data.frame(t(cm2[, 1:length(cm1)])), # nolint
                      function(z) paste(z, collapse = ":")
                    )
                  ),
                  cm2[, (length(cm1) + 1):ncol(cm2)]
                ),
                c("Group", names(cm2[, (length(cm1) + 1):ncol(cm2)]))
              )
              return(cm2) # nolint
            }
          ),
          fold_comb[[1]][[1]]
        )
      )
    }
    return(fold_mean) # nolint
  }
  ## Fold change calculation
  fun.fccalc <- function() { # nolint
    # Group Fold Changes
    fold_all <- setNames(
      reshape2::melt(
        dplyr::mutate(dplyr::bind_cols(
          lapply(
            seq.int(1, nrow(fold_comb[[2]]), 1),
            function(x) {
              fm1 <- t(
                fold_mean[
                  fold_mean[[1]] == fold_comb[[2]][x, 1], 2:ncol(fold_mean)
                ]
              )
              fm2 <- t(
                fold_mean[
                  fold_mean[[1]] == fold_comb[[2]][x, 2], 2:ncol(fold_mean)
                ]
              )
              fc1 <- setNames(
                as.data.frame(fun.comb(fm1, fm2)),
                paste(fold_comb[[2]][x, 1], fold_comb[[2]][x, 2], sep = "-")
              )
              if(mtype == "norm") { # nolint
                fc1 <- log2(fc1)
              }
              if(mtype == "scaled") { # nolint
                fc1 <- fc1
              }
              return(fc1) # nolint
            }
          )
        ), "Name" = names(d1)), id.vars = "Name"
      ),
      c("Name", "Comparison.fc", "Log2FC")
    )
    fold_all[["Comparison.fc"]] <- as.character(fold_all[["Comparison.fc"]])
    return(fold_all) # nolint
  }
  #---- Standard fold change calculation ----
  if (fc_method == "standard") {
    print("Calculating fold changes using standard method...")
    d1 <- setNames(d1, an1[[cn1]])
  }
  #---- Compound class fold changes ----
  if (fc_method == "class") { # nolint
    print("Calculating fold changes for each compound class...") # nolint
    d1 <- aggregate(
      t(as.matrix(d1)),
      list(an1[[cc1]]),
      function(x) mean(x)
    )
    d1 <- setNames(as.data.frame(t(d1[, 2:ncol(d1)])), d1[[1]])
  }
  #---- HiVE lipid ratio fold changes ----
  if (fc_method == "HiVE") { # nolint
    print("Calculating fold changes for HiVE comparisons...") # nolint
    d1 <- aggregate(
      t(as.matrix(d1)),
      list(an1[[cc1]]),
      function(x) mean(x)
    )
    d1 <- setNames(as.data.frame(t(d1[, 2:ncol(d1)])), d1[[1]])
    ## Retrieve ratios from HiVE edge list
    rat1 <- data.frame(
      "associated.ratio" = HiVE::nedge[["associated.ratio"]]
    )
    rat1 <- dplyr::mutate(
      rat1,
      "class1" = unlist(lapply(
        seq.int(1, nrow(rat1), 1),
        function(i) {
          unlist(strsplit(rat1[[1]][[i]], "_"))[[1]]
        }
      )),
      "class2" = unlist(lapply(
        seq.int(1, nrow(rat1), 1),
        function(i) {
          unlist(strsplit(rat1[[1]][[i]], "_"))[[2]]
        }
      ))
    )
    rat1 <- rat1[
      rat1[["class1"]] %in% names(d1) |
        rat1[["class2"]] %in% names(d1),
    ]
    rat2 <- setNames(lapply(
      seq.int(1, nrow(rat1), 1),
      function(x) {
        tryCatch(
          {
            fun.comb(d1[, rat1[x, "class1"]], d1[, rat1[x, "class2"]])
          },
          error = function(e) {
            print(
              paste(
                "Skipping ",
                rat1[["class1"]][[x]], "_",
                rat1[["class2"]][[x]],
                ": Not present in dataset",
                sep = ""
              )
            )
          }
        )
      }
    ), rat1[["associated.ratio"]])
    d1 <- dplyr::bind_cols(rat2[lengths(rat2) > 1])
  }
  #---- Calculation and Output ----
  if (comp_grp == FALSE) {
    fold_all <- data.frame("Group" = md1[[md_var1]], d1)
  }
  if (comp_grp == TRUE) {
    # Group Combinations
    fold_comb <- fun.grpcomb()
    # Metabolite Group Means (excluding missing samples)
    fold_mean <- fun.fcmean()
    # Group fold changes
    fold_all <- fun.fccalc()
  }
  return(fold_all) # nolint
}
