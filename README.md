# scViz

scViz is a RShiny application for exploratory analysis and visualization of scRNAseq. Currently, scViz imports Seurat objects (https://satijalab.org/seurat/) which must be obtained independently.

To use scViz, you need:
- R environment (v3.5 or greater)
- Rdata/RDS file including Seurat object containing clustering results and tSNE and/or UMAP dimensionality reduction information.
- R packages: Seurat (v3.0.1 or greater), clustree (v0.2.2 or greater), ggplot2 (v3.1.0 or greater), shiny (v1.2.0 or greater), shinyjs (v1.0 or greater), shinyWidgets (v0.4.8 or greater), colourpicker (v1.0 or greater).

