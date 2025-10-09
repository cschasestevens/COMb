# COMb v1.04 (20251009)
Comprehensive Observation/analysis of Metabolomics data (COMb)

## Description

Utilizes various R packages to perform processing and analysis of untargeted and targeted metabolomics/lipidomics datasets.    The methods included in this package provide a standardized workflow for data processing, statistical methods, and visualization of metabolomics data.    Most analyses can be run in parallel to expedite computationally heavy analyses.    The package is compatible with all operating systems. However, parallel processing is only available on Mac and Linux OS.

## Getting Started

### Dependencies
* Windows 10-11, WSL Ubuntu 22.04 or higher, Linux Ubuntu 22.04 or higher, or macOS 12.7.1 or higher
* R version 4.5.0 or higher (https://cran.r-project.org/)
* (Optional) RStudio version 2023.06.2 or higher (https://posit.co/download/rstudio-desktop/)
* (Optional) Conda installation of Python 3.9 and umap-learn (Only used for implementing UMAP via Python)
* R-packages (downloaded from CRAN or Bioconductor):
    * Suggests: 
        * knitr,
        * rmarkdown,
        * reticulate,
        * BiocManager
    * Imports: 
        * dplyr,
        * ggplot2,
        * ggsci,
        * viridis,
        * readxl,
        * ggpubr,
        * shadowtext,
        * parallel,
        * reshape2,
        * ggrepel,
        * mvnormtest,
        * grid,
        * EnhancedVolcano,
        * circlize,
        * ComplexHeatmap,
        * ggpmisc,
        * DescTools

### Installation
* Run the following in a new R session on the command line or within R-Studio:

```
devtools::install_github(
  "cschasestevens/COMb", 
  ref = "main", 
  build_vignettes = TRUE
)
```

## Help
* Browse vignettes by running the following:

```
browseVignettes("COMb")
```

* Access function documentation by running the following:

```
# Type function name after package name
?COMb::function_name()
```

## Authors

* Nathanial Chase Stevens, PhD, University of North Carolina at Chapel Hill
* Email: Nathanial_Stevens@med.unc.edu
* Alternate email: cschasestevens@gmail.com
* LinkedIn: https://www.linkedin.com/in/nathanial-chase-stevens-phd-08775180/

## Version History
* 1.04
    * Added internal standard normalization to
    ms_data_norm().
    * Updated ms_stat_crich() function to accept additional parameters for group and column selection.
    * Fixed error in ms_dim_rd() related to duplicate row names.
    * Added ms_data_check_pres() function to calculate the number of samples where a compound is present.
    * Revised ms_data_check() function.
* 1.03
    * Revised LOESS normalization to include MA-based correction
* 1.02
    * Added postprocessing functions for curated MS-DIAL alignment results
    * Added LOESS normalization as a method available in ms_data_norm()
* 1.01
    * Added function for preprocessing alignment results exported from MS-DIAL
* 1.00
    * Initial Release

## License

This project is licensed under the MIT license - see the LICENSE.md file for details.

## Acknowledgments

* ComplexHeatmap package: Gu Z, Eils R, Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics. <doi:10.1093/bioinformatics/btw313>
* Circlize package: Gu, Z. circlize implements and enhances circular visualization in R. Bioinformatics 2014.
