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
ms_add_anno <- function(
  f
) {
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

#' Data Input (Targeted)
#'
#' Reads input data from exported Skyline raw files.
#'
#' @param dir1 Load all files from the specified directory if not NULL.
#' @return List containing formatted input data and metadata
#' for downstream analysis.
#' @examples
#'
#' # msT_input(
#' #   dir1 = "data/"
#' # )
#'
#' @export
tms_input <- function( # nolint
  dir1 = NULL
) {
  # Load data
  f1 <- list.files(dir1, full.names = TRUE)
  f2 <- f1[grepl("sample", f1) & grepl(".txt", f1)]
  f3 <- f1[grepl("order", f1) & grepl(".txt", f1)]
  f1 <- f1[f1 %in% f2 == FALSE]
  f1 <- f1[f1 %in% f3 == FALSE]
  d1 <- setNames(
    lapply(
      f1,
      function(x) {
        if(grepl(".tsv", x)) { # nolint
          d1a <- read.table(x, sep = "\t", header = TRUE)
        }
        if(grepl(".csv", x)) { # nolint
          d1a <- read.csv(x)
        }
        dref <- gsub(paste(dir1, "/", sep = ""), "", x)
        return(list("data" = d1a, "file" = dref)) # nolint
      }
    ),
    gsub("\\.tsv|\\.csv", "", gsub(paste(dir1, "/", sep = ""), "", f1))
  )
  # Load sample list (if present)
  if(length(f2) == 1) { # nolint
    print("Sample list found in directory; loading sample list along with data")
    d1b <- read.table(f2, sep = "\t", header = TRUE)
  }
  if(length(f3) == 1) { # nolint
    print("Injection order found in directory; loading order along with data")
    d1c <- read.table(f3, sep = "\t", header = TRUE)
  }
  if(exists("d1b") && !exists("d1c")) { # nolint
    ld <- list("data" = d1, "sample_list" = d1b)
  }
  if(!exists("d1b") && exists("d1c")) { # nolint
    ld <- list("data" = d1, "sample_order" = d1c)
  }
  if(exists("d1b") && exists("d1c")) { # nolint
    ld <- list("data" = d1, "sample_list" = d1b, "sample_order" = d1c)
  }
  if(!exists("d1b") && !exists("d1c")) { # nolint
    ld <- list("data" = d1)
  }
  return(ld) # nolint
}

#' Data Format (Targeted)
#'
#' Formats input data from exported Skyline raw files. Must have
#' imported a valid sample list containing group information and from
#' msT_input and the exported sample_order if processing Skyline-based
#' targeted data.
#'
#' @param dat Data object generated by msT_input function.
#' @param intype Exported data type (currently supports "Skyline").
#' @return List containing formatted input data and metadata
#' for downstream analysis.
#' @examples
#'
#' # msT_format(
#' #   dat = d1
#' # )
#'
#' @import dplyr
#' @export
tms_format <- function( # nolint
  dat,
  intype = "Skyline"
) {
  # Load data and consolidate chromatogram and peak intensities
  ld1 <- dat
  ## If input type is Skyline (default)
  if(intype == "Skyline") { # nolint
    # Format chromatogram data
    ld1[["data"]][["chromatograms"]][[1]] <- setNames(ld1[["data"]][["chromatograms"]][[1]][ # nolint
      , c(
        "PeptideModifiedSequence", "PrecursorCharge", "ProductMz",
        "TotalArea", "Times", "Intensities"
      )
    ], c(
      "Name", "Adduct", "Product(mz)",
      "Peak_area_total", "RT_all", "Intensity"
    ))
    ## Important: Skyline exports intensities and chromatograms differently!
    ## Need column indicating sample order in Skyline and count of unique
    ## compounds exported from software.
    # Format intensities
    ld1[["data"]][["intensities"]][[1]] <- setNames(ld1[["data"]][["intensities"]][[1]][ # nolint
      , c(
        "Molecule", "Molecule.List.Name", "Replicate.Name",
        "Precursor.Mz", "Retention.Time", "Area", "Background"
      )
    ], c(
      "Name", "synthesis_pathway", "SampleID",
      "Precursor(mz)", "RT", "Peak_area", "Background"
    ))
    # Format transition list
    ld1[["data"]][["transition_list"]][[1]] <- setNames(ld1[["data"]][["transition_list"]][[1]][ # nolint
      , c(
        "Molecule.List.Name", "Molecule.Name",
        "Precursor.Mz", "Product.Mz"
      )
    ], c(
      "synthesis_pathway", "Name",
      "Precursor(mz)", "Product(mz)"
    ))
    # Merge chromatograms with intensities
    ## Set order of chromatograms to match intensities
    ## Assign chromatogram sample identities
    ld1[["data"]][["chromatograms"]][[1]] <- dplyr::select(
      dplyr::mutate(
        ld1[["data"]][["chromatograms"]][[1]],
        "SampleID" = unlist(
          lapply(
            as.numeric(ld1[["sample_order"]][["Order"]]),
            function(x) {
              rep(
                ld1[["sample_order"]][["SampleID"]][[x]],
                length(
                  unique(
                    ld1[["data"]][["chromatograms"]][[1]][["Name"]]
                  )
                )
              )
            }
          )
        )
      ),
      "SampleID", everything() # nolint
    )
    ld1[["data"]][["merged"]] <- dplyr::left_join(
      ld1[["data"]][["intensities"]][[1]],
      ld1[["data"]][["chromatograms"]][[1]],
      by = c("SampleID", "Name")
    )
    # Create index and assign groups
    ld1[["data"]][["merged"]][["ID"]] <- seq.int(
      1,
      nrow(ld1[["data"]][["merged"]]),
      1
    )
    ld1[["data"]][["merged"]][["Sample"]] <- gsub(
      "_.*", "", ld1[["data"]][["merged"]][["SampleID"]]
    )
    ld1[["sample_list"]][["SampleID"]] <-
      ld1[["sample_list"]][["Sample_ID"]]
    ld1[["sample_list"]][["Group"]] <- ld1[["sample_list"]][["sub_treatment"]]
    ld1[["sample_list"]][["Type"]] <- ld1[["sample_list"]][["treatment"]]
    ld1[["sample_list"]][["Sample"]] <- gsub(
      "_.*", "", ld1[["sample_list"]][["SampleID"]]
    )
    ld1[["data"]][["merged"]] <- dplyr::select(
      dplyr::left_join(
        ld1[["data"]][["merged"]],
        ld1[["sample_list"]][
          ,
          c("Sample", "Code", "Replicate", "Type", "Group")
        ],
        by = "Sample"
      ),
      "ID", "Sample", "Code", "Replicate","Type", "Group", everything() # nolint
    )
    ## Assign quality control types
    ld1[["data"]][["merged"]][["Replicate"]] <- ifelse(
      is.na(ld1[["data"]][["merged"]][["Replicate"]]),
      1,
      ld1[["data"]][["merged"]][["Replicate"]]
    )
    ld1[["data"]][["merged"]][["Code"]] <- ifelse(
      is.na(ld1[["data"]][["merged"]][["Code"]]),
      1,
      ld1[["data"]][["merged"]][["Code"]]
    )
    ld1[["data"]][["merged"]][["Group"]] <- ifelse(
      is.na(ld1[["data"]][["merged"]][["Group"]]) &
        grepl("QC|qc", ld1[["data"]][["merged"]][["SampleID"]]),
      "QC",
      ifelse(
        is.na(ld1[["data"]][["merged"]][["Group"]]) &
          grepl("MB|mb", ld1[["data"]][["merged"]][["SampleID"]]),
        "Blank",
        ifelse(
          is.na(ld1[["data"]][["merged"]][["Group"]]) &
            grepl("cal", ld1[["data"]][["merged"]][["SampleID"]]),
          "Std_curve",
          ifelse(
            is.na(ld1[["data"]][["merged"]][["Group"]]) &
              grepl("Pool|pool", ld1[["data"]][["merged"]][["SampleID"]]),
            "Pool",
            ifelse(
              is.na(ld1[["data"]][["merged"]][["Group"]]) &
                grepl("preinj", ld1[["data"]][["merged"]][["SampleID"]]),
              "Preinjection",
              ld1[["data"]][["merged"]][["Group"]]
            )
          )
        )
      )
    )
    ## Assign sample type
    ld1[["data"]][["merged"]][["Type"]] <- ifelse(
      grepl(
        "QC|Blank|Std_curve|Pool|Preinjection",
        ld1[["data"]][["merged"]][["Group"]]
      ),
      "QC",
      ld1[["data"]][["merged"]][["Type"]]
    )
    ## Rename experimental samples to match code, replicate, and group names
    ## EDIT to allow flexible column input
    ld1[["data"]][["merged"]][["label"]] <- ifelse(
      grepl(
        "QC|Blank|Std_curve|Pool|Preinjection",
        ld1[["data"]][["merged"]][["Group"]]
      ),
      ld1[["data"]][["merged"]][["SampleID"]],
      paste(
        ld1[["data"]][["merged"]][["Code"]],
        ld1[["data"]][["merged"]][["Replicate"]],
        sep = "."
      )
    )
    ## Assign injection order (if present)
    ld1[["data"]][["merged"]] <- dplyr::select(
      dplyr::left_join(
        ld1[["data"]][["merged"]],
        ld1[["sample_order"]],
        by = "SampleID"
      ),
      "Order", "label", everything() # nolint
    )
  }
  return(ld1[["data"]][["merged"]])
}
