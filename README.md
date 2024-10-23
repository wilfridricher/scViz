# scViz

scViz is a RShiny application for exploratory analysis and visualization of scRNAseq. Currently, scViz imports Seurat objects (https://satijalab.org/seurat/) which must be obtained independently.

To use scViz, you need:
- R environment (v4.3.0 or greater)
- Rdata/RDS file including Seurat object containing clustering results and tSNE and/or UMAP dimensionality reduction information.
- R packages: Seurat (v4.3.0 or greater), clustree (v0.5.0 or greater), ggplot2 (v3.4.2 or greater), shiny (v1.7.4 or greater), shinyjs (v2.1.0 or greater), shinyWidgets (v0.7.6 or greater), colourpicker (v1.2.0 or greater).

