#' Data Input
#'
#' Creates a data frame of processing parameters
#' used for data processing and analysis.
#'
#' @param f Name of a data file in the "data/" folder.
#' Accepts .xlsx, .csv, and .txt files.
#' @param md_num Number of metadata columns in input data file.
#' @param qc_rep Should a QC report be generated for the selected dataset?
#' (Either TRUE or FALSE)
#' @param norm1 Logical indicating if the data have been normalized.
#' (Either TRUE or FALSE)
#' @param ref1 Reference annotation list for assigning compound details and
#' class information (provided as a .txt file).
#' @param ref_alt If ref1 is NA, provide a text file containing
#' a custom reference annotation.
#' @return List containing formatted input data and metadata
#' for downstream analysis.
#' @examples
#'
#' # ms_input(
#' #   "example.txt",
#' #   4,
#' #   FALSE,
#' #   TRUE,
#' #   "reference.txt"
#' # )
#'
#' @export
ms_input <- function( # nolint
  f = NULL,
  md_num,
  qc_rep,
  norm1,
  ref1 = NULL,
  ref_alt = NULL
) {
  if(is.null(f)) { # nolint
    print("No data file was selected; using example dataset...")
    data("example1", package = "Stevens.MSAnalyze")
    d <- example1 # nolint
    a1 <- data.frame(
      "ID" = seq.int(
        1,
        length(names(d[, (md_num + 1):ncol(d)])),
        1
      ),
      "name" = names(d[, (md_num + 1):ncol(d)]),
      "type" = ifelse(
        grepl(
          "iSTD|1_CE 22:1|1_Sphingosine d17:1",
          names(d[, (md_num + 1):ncol(d)])
        ),
        "iSTD",
        "an.comp"
      )
    )
    r1 <- DBI::dbGetQuery(
      db1, # nolint
      'select * from master where "Source" = "Stev_ozaw"'
    )
    lpar <- list(
      "data" = d[, -c(1:md_num, a1[a1[["type"]] == "iSTD", "ID"] + md_num)],
      "meta" = d[, 1:md_num],
      "anno" = a1,
      "istd" = d[, a1[a1[["type"]] == "iSTD", "ID"] + md_num],
      "qc" = qc_rep,
      "nrm" = norm1,
      "anno.ref" = r1
    )
  }
  if(is.na(ref1)) { # nolint
    r1 <- read.table(
      ref_alt,
      sep = "\t",
      header = TRUE
    )
  }
  if(!is.na(ref1)) { # nolint
    r1 <- DBI::dbGetQuery(
      db1, # nolint
      'select * from master' # nolint
    )
    r1 <- r1[r1[["Source"]] == ref1, ]
  }
  if(is.null(ref1) == TRUE) { # nolint
    r1 <- DBI::dbGetQuery(
      db1, # nolint
      'select * from master' # nolint
    )
  }
  if(is.null(f) == FALSE && tools::file_ext(f) == "xlsx") { # nolint
    d <- readxl::read_excel(
      paste(
        "data/",
        f,
        sep = ""
      ),
      sheet = 1
    )
    a1 <- data.frame(
      "ID" = seq.int(
        1,
        length(names(d[, (md_num + 1):ncol(d)])),
        1
      ),
      "name" = names(d[, (md_num + 1):ncol(d)]),
      "type" = ifelse(
        grepl(
          "iSTD|1_CE 22:1|1_Sphingosine d17:1",
          names(d[, (md_num + 1):ncol(d)])
        ),
        "iSTD",
        "an.comp"
      )
    )
    lpar <- list(
      "data" = d[, -c(1:md_num, a1[a1[["type"]] == "iSTD", "ID"] + md_num)],
      "meta" = d[, 1:md_num],
      "anno" = a1,
      "istd" = d[, a1[a1[["type"]] == "iSTD", "ID"] + md_num],
      "qc" = qc_rep,
      "nrm" = norm1,
      "anno.ref" = r1
    )
  }
  if(is.null(f) == FALSE && tools::file_ext(f) == "csv") { # nolint
    d <- read.csv(
      paste(
        "data/",
        f,
        sep = ""
      ),
      check.names = TRUE,
      header = TRUE
    )
    a1 <- data.frame(
      "ID" = seq.int(
        1,
        length(names(d[, (md_num + 1):ncol(d)])),
        1
      ),
      "name" = names(d[, (md_num + 1):ncol(d)]),
      "type" = ifelse(
        grepl(
          "iSTD|1_CE.22.1|1_Sphingosine.d17.1",
          names(d[, (md_num + 1):ncol(d)])
        ),
        "iSTD",
        "an.comp"
      )
    )
    lpar <- list(
      "data" = d[, -c(1:md_num, a1[a1[["type"]] == "iSTD", "ID"] + md_num)],
      "meta" = d[, 1:md_num],
      "anno" = a1,
      "istd" = d[, a1[a1[["type"]] == "iSTD", "ID"] + md_num],
      "qc" = qc_rep,
      "nrm" = norm1,
      "anno.ref" = r1
    )
  }
  if(is.null(f) == FALSE && tools::file_ext(f) == "txt") { # nolint
    d <- read.table(
      paste(
        "data/",
        f,
        sep = ""
      ),
      check.names = TRUE,
      header = TRUE,
      sep = "\t"
    )
    a1 <- data.frame(
      "ID" = seq.int(
        1,
        length(names(d[, (md_num + 1):ncol(d)])),
        1
      ),
      "name" = names(d[, (md_num + 1):ncol(d)]),
      "type" = ifelse(
        grepl(
          "iSTD|1_CE.22.1|1_Sphingosine.d17.1",
          names(d[, (md_num + 1):ncol(d)])
        ),
        "iSTD",
        "an.comp"
      )
    )
    lpar <- list(
      "data" = d[, -c(1:md_num, a1[a1[["type"]] == "iSTD", "ID"] + md_num)],
      "meta" = d[, 1:md_num],
      "anno" = a1,
      "istd" = d[, a1[a1[["type"]] == "iSTD", "ID"] + md_num],
      "qc" = qc_rep,
      "nrm" = norm1,
      "anno.ref" = r1
    )
  }
  return(lpar)
}

#' Add annotations to lipid database
#'
#' Appends annotations from a dataset not already present in the
#' lipid database.
#'
#' @param f Name of a data file in the "data/" folder.
#' Accepts .xlsx, .csv, and .txt files.
#' @param md_num Number of metadata columns in input data file.
#' @param qc_rep Should a QC report be generated for the selected dataset?
#' (Either TRUE or FALSE)
#' @param norm1 Logical indicating if the data have been normalized.
#' (Either TRUE or FALSE)
#' @param ref1 Reference annotation list for assigning compound details and
#' class information (provided as a .txt file).
#' @return List containing formatted input data and metadata
#' for downstream analysis.
#' @examples
#'
#' # ms_input(
#' #   "example.txt",
#' #   4,
#' #   FALSE,
#' #   TRUE,
#' #   "reference.txt"
#' # )
#'
#' @export
ms_add_anno <- function(f) {
  # load dataset
  d_prev3 <- read.table(
    "data/data.pnnl.balf.txt",
    sep = "\t",
    header = TRUE
  )
  # load lipid reference
  amast <- read.table( # nolint
    "ref/0.masterlist.txt",
    sep = "\t",
    header = TRUE
  )
  amast[["input.name"]] <- gsub(
    "\\-|\\(|\\)|\\:|\\/|\\;|\\ |\\|",
    ".",
    amast[["Name"]]
  )
  an1 <- data.frame(
    "Source" = rep("Pnnl_balf", length(names(d_prev3[, 3:ncol(d_prev3)]))),
    "Species" = rep("Human", length(names(d_prev3[, 3:ncol(d_prev3)]))),
    "Organ" = rep("Lung", length(names(d_prev3[, 3:ncol(d_prev3)]))),
    "Matrix" = rep("BALF", length(names(d_prev3[, 3:ncol(d_prev3)]))),
    "Tx" = rep("Normal", length(names(d_prev3[, 3:ncol(d_prev3)]))),
    "input.name" = names(d_prev3[, 3:ncol(d_prev3)])
  )
  # Assign classes to annotations already present in df
  an_pres <- an1[an1[["input.name"]] %in% sort(unique(amast[["input.name"]])), ]
  an_pres <- dplyr::left_join(
    an_pres,
    amast[, -c(1:5)],
    by = "input.name"
  )
  # Join with annotations not found in df
  an_miss <- an1[
    an1[["input.name"]] %in% sort(unique(amast[["input.name"]])) == FALSE,
  ]
  an_out <- dplyr::bind_rows(
    an_pres,
    an_miss
  )
  # Export for assigning remaining classes
  write.table(
    an_out,
    "ref/label.data.txt",
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
  )

}