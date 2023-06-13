# ST data analysis workflow
To understand the spatial cell type distribution in the tumor microenvironment, we first normalized all of the data and used the scRNA-seq data corresponding to the cancer type as a reference, as shown in the figure below. Dimensions, debatch effects, grouping, and annotation were applied before deconvoluting the ST data with STRIDE. We first performed cell type deconvolution at the major lineage level, using Seurat's FindAllMarkers function on scRNA-seq data to compute the top 120 differentially expressed genes for each cell type, which we then utilized as signature genes for running STRIDE on. We selected the top 120 differentially expressed genes for each subtype and utilized them to run [STRIDE](https://github.com/wanglabtongji/STRIDE) in each major lineage cell type independently.
 
![ST data analysis workflow]([图片链接](https://github.com/JialiYue/Pan-Cancer-TLS/blob/main/ST_data_analysis_workflow.png))

## System requirements
* Python (>= 3.7) 
* numpy, pandas, matplotlib, scikit-image (== 0.19.2)


