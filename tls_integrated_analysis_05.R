rm(list = ls())
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
#library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(ggpubr)

RNA.res <- readRDS("/fs/home/yuejiali/ST/BRCA/Analysis/TLS_Finder/BRCA_tls_20230213_res.rds")
RNA.res@meta.data$study <- unlist(lapply(strsplit(rownames(RNA.res@meta.data),"@_"),function(x) x[1]))
RNA.res@meta.data$tls.label <- unlist(lapply(strsplit(rownames(RNA.res@meta.data),"@_"),function(x) x[2]))

RNA.res[["percent.mt"]] <- PercentageFeatureSet(object = RNA.res, pattern = "^MT-")
RNA.res <- subset(RNA.res, nFeature_Spatial > 200 & nFeature_Spatial < 10000 & percent.mt < 10)

RNA.res <- SCTransform(RNA.res,assay = "Spatial",vars.to.regress = 'percent.mt',verbose=TRUE)
RNA.res <- RunPCA(RNA.res)
# UMAP and clustering with top PCs
RNA.res <- RunUMAP(RNA.res, reduction='pca', dims = 1:20)
RNA.res <- FindNeighbors(RNA.res, reduction='pca',dims = 1:12)
p <- DimPlot(RNA.res, reduction = "pca",group.by = "study",label =FALSE,pt.size = 2) #+ NoLegend()
p

#Batch Remove
{
  library(harmony)
  RNA.res <- RunHarmony(RNA.res,'orig.ident', assay.use="SCT")
  # UMAP and clustering with harmonized PCs
  RNA.res <- RunUMAP(RNA.res, reduction='harmony', dims = 1:12)
  RNA.res <- FindNeighbors(RNA.res, reduction='harmony')
  RNA.res <- FindClusters(RNA.res, resolution = 1.4)
  p <- DimPlot(RNA.res, reduction = "harmony",group.by = "seurat_clusters",cols = brewer.pal(n = 4, name = "RdBu"),
               label =FALSE,pt.size = 3) #+ NoLegend()
  p
  #ggsave(file.path(paste0("BRCA_no_batch.pdf")),p,height = 5,width = 6.5)
  
}

#marker gene expression
{
  genes = c("RORC", "IL7R", "LTB","CXCL13")
  p1<- FeaturePlot(RNA.res, features = c("LTB","CXCL13","CD79B","MS4A1"),reduction = 'harmony',pt.size = 0.8,combine = FALSE)
  p2 <- lapply(p1, function (x) x  + scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = 'Spectral'))))
  p = CombinePlots(p2)
  p             
  #ggsave(file.path(paste0("/fs/home/yuejiali/ST/BRCA/Analysis/TLS_Finder/BRCA_tls_pos_4_genes.pdf")),p,height = 6,width = 7)  
  
}

{
  library(ComplexHeatmap)
  library(circlize) ## colo
  
  features=c("CD8A","CD4",#T cells
             'CD79A','CD79B','CD19','MS4A1',#B cells
             'RORC','IL7R',#initiating markers
             'SPIB','LTB','MS4A1','VPREB3','CD79A','CD22','CD37','CXCR5',
             'IRF8','CD79B','LIMD2','BANK1','CD19','ICAM3','TLR10',
             'PPP1R16B','WDFY4','POU2F2','GRAP','FDCSP','RALGPS2','SELL',
             'RAC2','FCMR','RUBCNL','MMP9'#activated    'FCRLA','CR2','FCRL3',  
             
  )
  
  
  row_split = c(rep("T cells",2),rep("B cells",4),rep("Initiating",2),rep("Activated",26))
  row_split = factor(row_split,levels = c("T cells","B cells","Initiating","Activated"))
  length(features)
  length(row_split)
  
  DefaultAssay(RNA.res)<- "Spatial"
  data.features <- FetchData(object = RNA.res, vars = features)
  data.features$id <- RNA.res$SCT_snn_res.1.4
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  ## prepare data.plot
  data.exp <- NULL
  data.plot <- sapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    if(!is.null(data.exp)){
      data.exp <- cbind(data.exp,avg.exp)
    }else{data.exp <- avg.exp}
    return(data.exp)
  })
  col<-colnames(data.plot)
  data.plot = t(apply(data.plot, 1, scale))
  colnames(data.plot)<-col
  # set color gradient
  col_fun <- colorRamp2(c(-1, 0, 5), c("#118ab2", "#fdffb6", "#e63946"))
  # split heatmap
  col_split = c(rep("0",1),rep("1",1),rep("2",1),rep("3",1))
  col_split =factor(col_split,levels = c('0','1','2','3'))
  
  annot = c("T cells","B cells","Initiating","Activated" )
  row_color=c("#e76f51","#f4a261","#0077b6","#ddbea9","#00b4d8","#dc2f02","#fca311","#00CDAC","#8DD7BF")
  ha = HeatmapAnnotation(df = data.frame(Marker=row_split),which = "row",
                         col = list(Marker = c("T cells" = "#e76f51","B cells"="#00CDAC","Initiating" = "#f4a261","Activated"="#0077b6")))
  a<-c('0','1','2','3')
  pdf(paste0("/fs/home/yuejiali/ST/BRCA/Analysis/TLS_Finder/BRCA_marker_heat.pdf"),width = 5,height = 12)
  Heatmap(data.plot[,a], col = col_fun,cluster_columns = FALSE,cluster_rows = FALSE,
          show_column_names = FALSE,show_row_names = TRUE,
          column_names_side = "top",row_names_side = "right",
          row_split = row_split,column_split = col_split,
          row_gap = unit(3, "mm"),column_gap =unit(3, "mm"), 
          left_annotation = ha,
          heatmap_legend_param=list(title = "Z-score",legend_height=unit(3, "cm")))
  dev.off()
}
