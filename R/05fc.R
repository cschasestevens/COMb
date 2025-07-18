#' Group-wise Fold Change Function
#'
#' Calculates the log2 fold changes for all
#' possible combinations of selected metadata variables.
#'
#' @param mat1 Data matrix assigned by ms_input().
#' @param mat_type Input data type, either "norm" or "scaled." If
#' "scaled," fold changes are calculated as the difference between
#' two group averages rather than by dividing normalized compound
#' average intensities per group.
#' @param md Dataframe containing study metadata.
#' @param md_var Vector of metadata variables to calculate fold
#' changes for.
#' @param an Dataframe containing annotation metadata.
#' @param fc_class Logical indicating whether class-wide
#' fold changes should be calculated. Requires a valid
#' annotation list contained within ms_input() result.
#' @param grp_class Name of the class annotation column
#' to use. Ignored when fc_class is set to FALSE.
#' @param lipid_fc (Lipidomics Only) Logical indicating
#' whether related lipid class fold changes should be calculated.
#' If TRUE, calculates lipid class averages and calculates related class
#' ratios prior to treatment group comparisons.
#' @param comp_grp Only applicable if lipid_fc is TRUE; Should treatment
#' group fold changes be calculated? If FALSE, returns lipid class ratios
#' without aggregating and comparing by treatment group.
#' @return A data frame containing the group-wise
#' fold changes for each compound.
#' @examples
#'
#' # ms_stat_fc(
#' #   mat1 = d[["data"]][["imputed"]],
#' #   md = d[["data"]][["meta"]],
#' #   md_var = c("Organ", "Sex", "Dose"),
#' #   an = d[["data"]][["anno"]],
#' #   fc_class = TRUE,
#' #   grp_class = "Subclass"
#' # )
#'
#' @export
ms_stat_fc <- function( # nolint
  mat1,
  mat_type = "norm",
  md,
  md_var,
  an = NULL,
  fc_class = FALSE,
  grp_class = NULL,
  lipid_fc = FALSE,
  comp_grp = TRUE
) {
  if(mat_type == "norm") { # nolint
    fun.comb <- function(x, y) {x / y} # nolint
  }
  if(mat_type == "scaled") { # nolint
    fun.comb <- function(x, y) {x - y} # nolint
  }
  if(fc_class == FALSE && lipid_fc == FALSE) { # nolint
    d1 <- mat1
    md1 <- md
    if(is.null(an)) { # nolint
      an1 <- data.frame("Name" = names(mat1))
    }
    if(!is.null(an)) { # nolint
      an1 <- an
    }
    md_var1 <- md_var
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
        return(cb1)
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
      cb2
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
      fold_comb
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
            return(cb2)
          }
        )
      )
    }
    # Metabolite Group Means (excluding missing samples)
    fold_mean <- dplyr::bind_rows(
      setNames(
        parallel::mclapply(
          mc.cores = ceiling(parallel::detectCores() * 0.75),
          seq.int(1, nrow(var_comb), 1),
          function(x) {
            cm1 <- paste(unlist(strsplit(var_comb[x, 1], ":")), sep = ", ")
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
            return(cm2)
          }
        ),
        var_comb[[1]]
      )
    )
    # Group Fold Changes
    fold_all <- setNames(
      reshape2::melt(
        dplyr::mutate(dplyr::bind_cols(
          lapply(
            seq.int(1, nrow(fold_comb), 1),
            function(x) {
              fm1 <- t(
                fold_mean[fold_mean[[1]] == fold_comb[x, 1], 2:ncol(fold_mean)]
              )
              fm2 <- t(
                fold_mean[fold_mean[[1]] == fold_comb[x, 2], 2:ncol(fold_mean)]
              )
              fc1 <- setNames(
                as.data.frame(fun.comb(fm1, fm2)),
                paste(fold_comb[x, 1], fold_comb[x, 2], sep = "-")
              )
              if(mat_type == "norm") { # nolint
                fc1 <- log2(fc1)
              }
              if(mat_type == "scaled") { # nolint
                fc1 <- fc1
              }
              return(fc1)
            }
          )
        ), "Name" = an1[["Name"]]), id.vars = "Name"
      ),
      c("Name", "Comparison.fc", "Log2FC")
    )
    fold_all[["Comparison.fc"]] <- as.character(fold_all[["Comparison.fc"]])
  }

  if(fc_class == TRUE && lipid_fc == FALSE) { # nolint
    d1 <- mat1
    an1 <- an
    d1 <- aggregate(
      t(as.matrix(d1)),
      list(an1[[grp_class]]),
      function(x) mean(x)
    )
    d1 <- setNames(as.data.frame(t(d1[, 2:ncol(d1)])), d1[[1]])
    md1 <- md

    # Group Combinations
    var_comb <- dplyr::bind_rows(lapply(
      seq.int(0, length(md_var) - 1, 1),
      function(x) {
        cb1 <- combn(md_var, x + 1)
        cb1 <- as.data.frame(t(
          dplyr::bind_cols(
            lapply(
              as.data.frame(cb1),
              function(y) paste(y, collapse = ":")
            )
          )
        ))
        return(cb1)
      }
    ))

    if(nrow(var_comb) == 1) { # nolint
      cb1 <- paste(unlist(strsplit(var_comb[1, 1], ":")), sep = ", ")
      # Change variable to factor
      cb2 <- setNames(
        as.data.frame(
          as.character(md1[[md_var]])
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
      cb2
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
      fold_comb
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
            return(cb2)
          }
        )
      )
    }
    # Metabolite Group Means (excluding missing samples)
    fold_mean <- dplyr::bind_rows(
      setNames(
        parallel::mclapply(
          mc.cores = ceiling(parallel::detectCores() * 0.75),
          seq.int(1, nrow(var_comb), 1),
          function(x) {
            cm1 <- paste(unlist(strsplit(var_comb[x, 1], ":")), sep = ", ")
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
            return(cm2)
          }
        ),
        var_comb[[1]]
      )
    )

    # Group Fold Changes
    fold_all <- setNames(
      reshape2::melt(
        dplyr::mutate(dplyr::bind_cols(
          lapply(
            seq.int(1, nrow(fold_comb), 1),
            function(x) {
              fm1 <- t(
                fold_mean[fold_mean[[1]] == fold_comb[x, 1], 2:ncol(fold_mean)]
              )
              fm2 <- t(
                fold_mean[fold_mean[[1]] == fold_comb[x, 2], 2:ncol(fold_mean)]
              )
              fc1 <- setNames(
                as.data.frame(fun.comb(fm1, fm2)),
                paste(fold_comb[x, 1], fold_comb[x, 2], sep = "-")
              )
              if(mat_type == "norm") { # nolint
                fc1 <- log2(fc1)
              }
              if(mat_type == "scaled") { # nolint
                fc1 <- fc1
              }
              return(fc1)
            }
          )
        ), "Name" = names(d1)), id.vars = "Name"
      ),
      c("Name", "Comparison.fc", "Log2FC")
    )
    fold_all[["Comparison.fc"]] <- as.character(fold_all[["Comparison.fc"]])
  }
  # Calculate related lipid class fold changes for each group
  if(fc_class == FALSE && lipid_fc == TRUE) { # nolint
    # aggregate by class then pathway ratio
    d1 <- mat1
    an1 <- an
    d1 <- aggregate(
      t(as.matrix(d1)),
      list(an1[[grp_class]]),
      function(x) mean(x)
    )
    d1 <- setNames(as.data.frame(t(d1[, 2:ncol(d1)])), d1[[1]])
    md1 <- md
    ## pathway ratios
    rat1 <- DBI::dbGetQuery(
      db1, # nolint
      'select * from lipid_class_ratios' # nolint
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
              "Specified class combination not present in dataset; skipping to next combination" # nolint
            )
          }
        )
      }
    ), paste(rat1[["class1"]], rat1[["class2"]], sep = "_"))
    d1 <- dplyr::bind_cols(rat2[lengths(rat2) > 1])
    if(comp_grp == FALSE) { # nolint
      fold_all <- data.frame(
        "Group" = md1[[md_var]],
        d1
      )
    }
    if(comp_grp == TRUE) { # nolint
      # Group Combinations
      var_comb <- dplyr::bind_rows(lapply(
        seq.int(0, length(md_var) - 1, 1),
        function(x) {
          cb1 <- combn(md_var, x + 1)
          cb1 <- as.data.frame(t(
            dplyr::bind_cols(
              lapply(
                as.data.frame(cb1),
                function(y) paste(y, collapse = ":")
              )
            )
          ))
          return(cb1)
        }
      ))
      if(nrow(var_comb) == 1) { # nolint
        cb1 <- paste(unlist(strsplit(var_comb[1, 1], ":")), sep = ", ")
        # Change variable to factor
        cb2 <- setNames(
          as.data.frame(
            as.character(md1[[md_var]])
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
        fold_comb
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
              return(cb2)
            }
          )
        )
      }
      # Metabolite Group Means (excluding missing samples)
      fold_mean <- dplyr::bind_rows(
        setNames(
          parallel::mclapply(
            mc.cores = ceiling(parallel::detectCores() * 0.75),
            seq.int(1, nrow(var_comb), 1),
            function(x) {
              cm1 <- paste(unlist(strsplit(var_comb[x, 1], ":")), sep = ", ")
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
              return(cm2)
            }
          ),
          var_comb[[1]]
        )
      )
      # Group Fold Changes
      fold_all <- setNames(
        reshape2::melt(
          dplyr::mutate(dplyr::bind_cols(
            lapply(
              seq.int(1, nrow(fold_comb), 1),
              function(x) {
                fm1 <- t(
                  fold_mean[
                    fold_mean[[1]] == fold_comb[x, 1], 2:ncol(fold_mean)
                  ]
                )
                fm2 <- t(
                  fold_mean[
                    fold_mean[[1]] == fold_comb[x, 2], 2:ncol(fold_mean)
                  ]
                )
                fc1 <- setNames(
                  as.data.frame(fun.comb(fm1, fm2)),
                  paste(fold_comb[x, 1], fold_comb[x, 2], sep = "-")
                )
                if(mat_type == "norm") { # nolint
                  fc1 <- log2(fc1)
                }
                if(mat_type == "scaled") { # nolint
                  fc1 <- fc1
                }
                return(fc1)
              }
            )
          ), "Name" = names(d1)), id.vars = "Name"
        ),
        c("Name", "Comparison.fc", "Log2FC")
      )
      fold_all[["Comparison.fc"]] <- as.character(fold_all[["Comparison.fc"]])
    }
  }
  return(fold_all)
}
