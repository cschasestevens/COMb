#' Multivariate ANOVA
#'
#' Conducts an ANOVA for each metadata variable with Tukey's post-hoc
#' analysis for each compound and combines with fold changes
#' for each group comparison.
#'
#' @param exp1 A SummarizedExperiment object.
#' @param asy Input assay for statistical analysis.
#' @param mat_type Input data type, either "norm" or "scaled." If
#' "scaled," fold changes are calculated as the difference between
#' two group averages rather than by dividing normalized compound
#' average intensities per group.
#' @param col_grp Grouping variables for testing.
#' @param col_class Compound class variable name; only required if
#' fc_method is "class" or "HiVE".
#' @param col_lab Compound name column.
#' @param fc_method Fold change calculation mode; either "standard",
#' which is the most common, "class", or "HiVE".
#' @param name_group If statistics comparing two groups is desired,
#' provide a vector containing each group name for comparison.
#' @return A data frame containing the combined ANOVA and fold change results
#' for the specified variable(s).
#' @examples
#'
#' # ms_stat_uni(exp1 = data2)
#'
#' @export
ms_stat_uni <- function( # nolint
  exp1,
  asy = "norm",
  mat_type = "norm",
  col_grp = "Group",
  col_class = "Class",
  col_lab = "Label",
  fc_method = "standard",
  name_group = NULL
) { # nolint
  #---- Setup ----
  # Load summarizedexperiment and associated data
  ## Input matrix
  d1 <- as.data.frame(t(assay(exp1, asy))) # nolint
  mtype <- mat_type
  ## row data
  an1 <- as.data.frame(rowData(exp1)) # nolint
  ## column data
  md1 <- as.data.frame(colData(exp1)) # nolint
  ## metadata variable names
  md_var1 <- col_grp
  cc1 <- col_class
  cn1 <- col_lab
  grpname <- name_group
  #---- Functions ----
  fun.check.norm <- function() { # nolint
    # Check for normality
    d1_check <- as.data.frame(
      apply(
        as.matrix(d1),
        2,
        function(x) {
          tryCatch(
            {
              mvnormtest::mshapiro.test(t(as.matrix(x)))[[2]]
            },
            error = function(e) {
              print(
                "Normality test unsuccessful; NAs likely present in input data" # nolint
              )
            }
          )
        }
      )
    )
    print(
      paste(
        length(d1_check[d1_check[[1]] > (0.05 / nrow(d1_check)), ]),
        "of",
        nrow(d1_check),
        "compounds have normally distributed intensities."
      )
    )
    return(d1_check) # nolint
  }
  # One-way ANOVA
  fun.stat.anova <- function() { # nolint
    fun.test.in <- function(y) { # nolint
      d2 <- tryCatch(
        {
          if (length(md_var1) == 3) {
            d2 <- dplyr::bind_rows(
              lapply(
                seq.int(1, 1, 1),
                function(x) {
                  as.data.frame(
                    as.matrix(
                      TukeyHSD(
                        aov(
                          d1[[y]] ~
                            as.factor(md1[[md_var1[[1]]]]) *
                              as.factor(md1[[md_var1[[2]]]]) *
                              as.factor(md1[[md_var1[[3]]]])
                        )
                      )[[x]]
                    )
                  )
                }
              )
            )
          }
          if (length(md_var1) == 2) {
            d2 <- dplyr::bind_rows(
              lapply(
                seq.int(1, 1, 1),
                function(x) {
                  as.data.frame(
                    as.matrix(
                      TukeyHSD(
                        aov(
                          d1[[y]] ~
                            as.factor(md1[[md_var1[[1]]]]) *
                              as.factor(md1[[md_var1[[2]]]])
                        )
                      )[[x]]
                    )
                  )
                }
              )
            )
          }
          if (length(md_var1) == 1) {
            d2 <- dplyr::bind_rows(
              lapply(
                seq.int(1, 1, 1),
                function(x) {
                  as.data.frame(
                    as.matrix(
                      TukeyHSD(
                        aov(
                          d1[[y]] ~
                            as.factor(md1[[md_var1[[1]]]])
                        )
                      )[[x]]
                    )
                  )
                }
              )
            )
          }
        },
        error = function(e) {
          print("ANOVA failed for selected comparison...")
          d2 <- data.frame("p adj" = NA)
          names(d2) <- c("p adj")
          return(d2) # nolint
        }
      )
      return(d2) # nolint
    }
    if (Sys.info()[["sysname"]] != "Windows") {
      d_mult <- dplyr::bind_rows(
        setNames(
          parallel::mclapply(
            mc.cores = ceiling(parallel::detectCores() / 4),
            seq.int(1, ncol(d1), 1),
            function(y) {
              d2 <- fun.test.in(y)
              d2 <- dplyr::mutate(
                d2,
                "Name" = names(d1)[[y]],
                "Comparison" = rownames(d2),
                "p.value" = d2[["p adj"]]
              )
              d2 <- d2[, c("Name", "Comparison", "p.value")]
              d2[["FDR"]] <- p.adjust(d2[["p.value"]], method = "fdr")
              return(d2) # nolint
            }
          ),
          names(d1)
        )
      )
    }
    if (Sys.info()[["sysname"]] == "Windows") {
      d_mult <- dplyr::bind_rows(
        setNames(
          lapply(
            seq.int(1, ncol(d1), 1),
            function(y) {
              d2 <- fun.test.in(y)
              d2 <- dplyr::mutate(
                d2,
                "Name" = names(d1)[[y]],
                "Comparison" = rownames(d2),
                "p.value" = d2[["p adj"]]
              )
              d2 <- d2[, c("Name", "Comparison", "p.value")]
              d2[["FDR"]] <- p.adjust(d2[["p.value"]], method = "fdr")
              return(d2) # nolint
            }
          ),
          names(d1)
        )
      )
    }
    d_mult <- d_mult[d_mult[["Comparison"]] != "1", ]
    return(d_mult) # nolint
  }
  # Wilcox Rank Sum Test
  fun.stat.wilcox <- function() { # nolint
    fun.test.in <- function(y) { # nolint
      d2 <- tryCatch(
        {
          g1 <- d1[grepl(grpname[[1]], md1[[md_var1]]), ][[y]]
          g2 <- d1[grepl(grpname[[2]], md1[[md_var1]]), ][[y]]
          d2 <- setNames(
            as.data.frame(
              wilcox.test(g1, g2)[["p.value"]]
            ),
            "p.value"
          )
        },
        error = function(e) {
          print("Wilcox Test failed for selected comparison...")
          d2 <- data.frame("p.value" = NA)
          names(d2) <- c("p.value")
          return(d2) # nolint
        }
      )
      d2 <- dplyr::mutate(
        d2,
        "Name" = names(d1)[[y]],
        "Comparison" = paste(
          grpname[[1]],
          grpname[[2]], sep = "-"
        )
      )
      d2 <- d2[, c("Name", "Comparison", "p.value")]
      d2[["FDR"]] <- p.adjust(d2[["p.value"]], method = "fdr")
      return(d2) # nolint
    }
    if (Sys.info()[["sysname"]] != "Windows") {
      d_mult <- dplyr::bind_rows(
        setNames(
          parallel::mclapply(
            mc.cores = ceiling(parallel::detectCores() / 4),
            seq.int(1, ncol(d1), 1),
            function(y) {
              d2 <- fun.test.in(y)
              return(d2) # nolint
            }
          ),
          names(d1)
        )
      )
    }
    if (Sys.info()[["sysname"]] == "Windows") {
      d_mult <- dplyr::bind_rows(
        setNames(
          lapply(
            seq.int(1, ncol(d1), 1),
            function(y) {
              d2 <- fun.test.in(y)
              return(d2) # nolint
            }
          ),
          names(d1)
        )
      )
    }
    d_mult <- d_mult[d_mult[["Comparison"]] != "1", ]
    return(d_mult) # nolint
  }
  fun.stat.ttest <- function() { # nolint
    fun.test.in <- function(y) { # nolint
      d2 <- tryCatch(
        {
          g1 <- d1[grepl(grpname[[1]], md1[[md_var1]]), ][[y]]
          g2 <- d1[grepl(grpname[[2]], md1[[md_var1]]), ][[y]]
          d2 <- setNames(
            as.data.frame(
              t.test(g1, g2)[["p.value"]]
            ),
            "p.value"
          )
        },
        error = function(e) {
          print("Wilcox Test failed for selected comparison...")
          d2 <- data.frame("p.value" = NA)
          names(d2) <- c("p.value")
          return(d2) # nolint
        }
      )
      d2 <- dplyr::mutate(
        d2,
        "Name" = names(d1)[[y]],
        "Comparison" = paste(
          grpname[[1]],
          grpname[[2]], sep = "-"
        )
      )
      d2 <- d2[, c("Name", "Comparison", "p.value")]
      d2[["FDR"]] <- p.adjust(d2[["p.value"]], method = "fdr")
      return(d2) # nolint
    }
    if (Sys.info()[["sysname"]] != "Windows") {
      d_mult <- dplyr::bind_rows(
        setNames(
          parallel::mclapply(
            mc.cores = ceiling(parallel::detectCores() / 4),
            seq.int(1, ncol(d1), 1),
            function(y) {
              d2 <- fun.test.in(y)
              return(d2) # nolint
            }
          ),
          names(d1)
        )
      )
    }
    if (Sys.info()[["sysname"]] == "Windows") {
      d_mult <- dplyr::bind_rows(
        setNames(
          lapply(
            seq.int(1, ncol(d1), 1),
            function(y) {
              d2 <- fun.test.in(y)
              return(d2) # nolint
            }
          ),
          names(d1)
        )
      )
    }
    d_mult <- d_mult[d_mult[["Comparison"]] != "1", ]
    return(d_mult) # nolint
  }
  # Combine results with fold change
  fun.format.res <- function() { # nolint
    # Fix anova comparisons to match fold change comparisons
    dm_fix <- data.frame(
      "Comparison" = unique(d_mult[["Comparison"]]),
      "Comparison.fc" = unlist(
        lapply(
          seq.int(1, length(unique(d_mult[["Comparison"]])), 1),
          function(x) {
            fx1 <- unique(d_mult[["Comparison"]])[[x]]
            fx2 <- ifelse(
              fx1 %in% unique(fold_all[["Comparison.fc"]]) == FALSE,
              paste(
                unlist(strsplit(fx1, "-"))[[2]],
                unlist(strsplit(fx1, "-"))[[1]],
                sep = "-"
              ),
              fx1
            )
            return(fx2) # nolint
          }
        )
      )
    )
    d_out <- dplyr::left_join(
      dplyr::left_join(
        d_mult,
        dm_fix,
        by = "Comparison"
      ),
      fold_all,
      by = c("Comparison.fc", "Name")
    )
    return(d_out) # nolint
  }
  # Define fold change function (for HiVE comparisons)
  if (mtype == "norm") {
    fun.comb <- function(x, y) {x / y} # nolint
  }
  if (mtype == "scaled") {
    fun.comb <- function(x, y) {x - y} # nolint
  }
  #---- Data precheck ----
  print("Normality Check:")
  d1_check <- fun.check.norm()
  #---- One-way ANOVA: Standard ----
  if (fc_method == "standard") {
    print("Running standard statistical analysis...")
    d1 <- setNames(d1, an1[[cn1]])
  }
  #---- One-way ANOVA: Class ----
  if (fc_method == "class") {
    print("Running class-based statistical analysis...")
    d1 <- aggregate(
      t(as.matrix(d1)),
      list(an1[[cc1]]),
      function(x) mean(x)
    )
    d1 <- setNames(as.data.frame(t(d1[, 2:ncol(d1)])), d1[[1]])
  }
  #---- One-way ANOVA: HiVE----
  if (fc_method == "HiVE") {
    print("Running HiVE-based statistical analysis...")
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
  #---- Run test and combine with fold change ----
  if (is.null(name_group)) {
    d_mult <- fun.stat.anova()
  }
  if (!is.null(name_group)) {
    print("Group name specified; run statistical test with two groups...")
    if (
      (length(d1_check[d1_check[[1]] > (0.05 / nrow(d1_check)), ])) /
        nrow(d1_check) < 0.5
    ) {
      print("Most compound intensities are not normally distributed; Running wilcox rank sum test...") # nolint
      d_mult <- fun.stat.wilcox()
    }
    if (
      (length(d1_check[d1_check[[1]] > (0.05 / nrow(d1_check)), ])) /
        nrow(d1_check) > 0.5
    ) {
      print("Most compound intensities are normally distributed; Running T-test...") # nolint
      d_mult <- fun.stat.ttest()
    }
  }
  fold_all <- ms_stat_fc( # nolint
    exp1 = exp1,
    asy = asy,
    mat_type = mtype,
    col_grp = md_var1,
    col_class = cc1,
    fc_method = fc_method
  )
  print("Combining results...")
  d_out <- fun.format.res()
  print("Statistical analysis complete!")
  return(d_out)
}
