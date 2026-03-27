#' Create SummarizedExperiment from Excel
#'
#' Creates a SummarizedExperiment for use with COMb functions from an
#' Excel file. This function produces an output equivalent to ms_postproc,
#' which is useful in cases where processed data is already available
#' for analysis. The input file should contain 3 essential components:
#' the data matrix, compound/feature metadata, and sample metadata.
#' Optionally, the user may provide an additional file containing study
#' metadata for use with specific COMb functions, including ms_qc and
#' ms_qc_report.
#'
#' @param file_data Path to input Excel file.
#' @param file_md (Optional) Path to input Excel file containing study metadata.
#' @param sheetno Numerical value indicating the sheet number containing
#' input data.
#' @param md_no_feat Number of compound/feature metadata columns.
#' @param md_no_samp Number of sample metadata columns.
#' @param format_data If TRUE, format input data matrix to remove NA values
#' and correct data type.
#' @return A SummarizedExperiment containing the data matrix
#' with the final feature and sample metadata
#' for the specified ionization mode(s).
#' @examples
#'
#' # data1 <- ms_create_exp(
#' #   "path/to/excelfile.xlsx",
#' #   md_no_feat = 10,
#' #   md_no_samp = 5
#' # )
#'
#' @export
ms_create_exp <- function(
  file_data,
  file_md = NULL,
  sheetno = 1,
  md_no_feat,
  md_no_samp,
  format_data = TRUE
) {
  # Set params
  file1 <- file_data
  file3 <- file_md
  nosheet <- sheetno
  mdfeat <- md_no_feat
  mdsamp <- md_no_samp
  # Load data
  ## data matrix
  d1 <- readxl::read_excel(
    path = file1,
    sheet = nosheet,
    skip = mdsamp,
    col_names = FALSE
  )
  d1 <- as.matrix(d1[, (mdfeat + 1):ncol(d1)])
  ## feature metadata
  d2 <- readxl::read_excel(
    path = file1,
    sheet = nosheet,
    skip = mdsamp - 1,
    col_names = TRUE
  )[1:mdfeat]
  ## sample metadata
  d3 <- readxl::read_excel(
    path = file1,
    sheet = nosheet,
    col_names = FALSE
  )
  d3 <- d3[1:mdsamp, (mdfeat):ncol(d3)]
  d3 <- setNames(as.data.frame(t(d3)), d3[[1]])[-1, ]
  ## study metadata (if provided)
  if (is.null(file3)) {
    d4 <- list("Study name" = NULL)
  }
  if (!is.null(file3)) {
    d4 <- readxl::read_excel(
      path = file3,
      col_names = TRUE
    )
    d4 <- setNames(
      lapply(
        seq.int(1, nrow(d4), 1),
        function(i) {
          d4[i, ][[2]]
        }
      ), d4[[1]]
    )
  }
  # Create SummarizedExperiment
  d_out <- SummarizedExperiment::SummarizedExperiment(
    assays = list("input" = d1),
    rowData = d2,
    colData = d3,
    metadata = d4
  )
  if (format_data == TRUE) {
    assay(d_out, "format") <- assay(d_out, "input")
    assay(d_out, "format")[grepl("n/d|na", assay(d_out, "format"))] <- NA
    assay(d_out, "format") <- apply(
      assay(d_out, "format"), 2, function(x) as.numeric(x)
    )
    assay(d_out, "format")[is.na(assay(d_out, "format"))] <- 0
  }
  return(d_out) # nolint
}
