# COMb v2.04 (20260424)
Comprehensive Observation/analysis of Metabolomics data (COMb)

## Description

Utilizes various R packages to perform processing and analysis of untargeted metabolomics/lipidomics datasets.    The methods included in this package provide a standardized workflow for data processing, statistical methods, and visualization of metabolomics data.    Most analyses can be run in parallel to expedite lengthy analyses.    The package is compatible with all operating systems. However, parallel processing is only available on Mac and Linux operating systems.

The HiVE package (https://github.com/cschasestevens/HiVE) extends the plotting capabilities of COMb by enabling generation of detailed lipid network plots for multimodal lipid analysis, including lipidomics and transcriptomics data. See the HiVE package documentation for more information.

## Getting Started

### Dependencies
* Windows 10-11, WSL Ubuntu 22.04 or higher, Linux Ubuntu 22.04 or higher, or macOS 12.7.1 or higher
* R version 4.5.0 or higher (https://cran.r-project.org/)
* (Optional) RStudio version 2023.06.2 or higher (https://posit.co/download/rstudio-desktop/)
* (Optional) Conda installation of Python 3.9 and umap-learn (Only used for implementing UMAP via Python)
* R-packages (downloaded from CRAN or Bioconductor):
    * Suggests: 
        * reticulate,
        * BiocManager,
        * rmarkdown,
        * parallel,
        * DescTools,
        * ggpmisc,
        * plotly,
        * HiVE,
        * ggnewscale,
        * ggrepel,
        * knitr,
        * kableExtra
    * Imports: 
        * ggsci,
        * viridis,
        * RColorBrewer,
        * ggplot2,
        * dplyr,
        * ggpubr,
        * gtools,
        * circlize,
        * ComplexHeatmap,
        * grid,
        * readxl,
        * SummarizedExperiment,
        * magrittr,
        * reshape2,
        * umap,
        * mvnormtest,
        * EnhancedVolcano

### Installation
* Run the following in a new R session on the command line or within R-Studio:

```
pak::pak("cschasestevens/COMb")
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
* 2.04
    * Added additional customization options to ms_plot_vio and ms_plot_bar.
    * Added option to include heatmap annotations in ms_check_pres.
* 2.03
    * Fixed error in ms_qc_report .rmd file related to missing data input tables.
    * Modified ms_create_exp to match output of ms_postproc for preventing errors in downstream QC assessment.
    * Added option to perform missing value imputation with no normalization in ms_data_norm.
    * Added pareto scale as a default option for ms_dim_rd to resolve data input errors when uploading normalized data that still contain missing values.
    * Fixed error in ms_qc related to selecting the compound class column name.
    * Updated installation instructions using pak (install_github is deprecated in newer versions of devtools).
* 2.02
    * Added package tutorial.
    * Updated package documentation.
    * Fixed error in ms_plot_crich to reflect new plot theme.
    * Added ms_create_exp function, which converts Excel files containing feature data, sample data, and a data matrix into a SummarizedExperiment for use with COMb functions.
    * Updated function documentation.
* 2.01
    * Updated documentation for all functions.
    * ms_plot_hm has replaced ms_plot_heat, which is now deprecated.
    * ms_stat_crich updated to allow more flexibility for reference dataset inputs.
    * Integrated wilcox/t-test into ms_stat_uni (ms_stat_univ is now deprecated).
* 2.00b
    * ms_plot_lipidnet has been deprecated; lipid network plots are now created using the HiVE package (https://github.com/cschasestevens/HiVE).
    * Updated ms_stat_uni for compatibility with updated ms_stat_fc function.
    * Updated ms_stat_fc for simplified input and usage with ms_stat_uni.
    * ms_stat_anova has now been superceded by ms_stat_uni.
* 2.00a
    * Processed data are now formatted as SummarizedExperiments for reproducibility.
    * Most functions now also accept SummarizedExperiment objects as inputs.
    * Added ms_qc_report function for generating summary reports from quality control assessment conducted by ms_qc.
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
