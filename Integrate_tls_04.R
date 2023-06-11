{
rm(list = ls())
library(monocle)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Seurat)
library(ggplot2)
library(enrichplot)
library(ggsci)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(CellChat)
library(patchwork)
library(ggpubr)
library(ggbreak)
}

# read obj and integrate tls 
{
  setwd("/fs/home/yuejiali/ST/BRCA/TLS_Predict/BRCA_32572199/TLS_Finder/")
  samples <- c('BC23268_C1','BC23268_C2','BC23268_D1','BC23269_C1',
               'BC23269_D1','BC23270_E2','BC23272_E1','BC23272_E2','BC23287_C1',
               'BC23288_D2','BC23288_E1','BC23288_E2','BC23508_E2','BC23567_E2','BC23803_E1','BC23803_E2',
               'BC23895_C2','BC23903_C2','BC23903_D1','BC24105_C1',
               'BC24223_D2','BC24223_E1','BC24223_E2')
  #c("A1","A2","A5","B2","B3","B4","B5","C1","C2","C3","C4","C5","D1","D2","D3",'D5','D6','E1','E2','E3','F1','F2','F3','G1','G2','G3','H1','H2','H3')
  #c('CID4290','CID4465','CID44971','CID4535')
  
  seurat_tls_list = list()
  for (i in samples){
    print(i)
    coloc <- read.table(paste0("/fs/home/yuejiali/ST_data/BRCA_32572199/ST/decov/",i,"_t_and_B_colocal.tsv"),header = TRUE,sep='\t')
    Spatial.obj = readRDS(file.path( paste0(i, "_min_Spatial_Object.rds")))
    Spatial.obj@meta.data <- cbind(Spatial.obj@meta.data,coloc[rownames(Spatial.obj@meta.data),])
    #rownames(Spatial.res@meta.data) <- Spatial.res@meta.data$barcode
    Spatial.obj@meta.data$isTLS <- as.factor(Spatial.obj@meta.data$isTLS)
    Spatial.obj.tls = subset(Spatial.obj, isTLS != "Stroma") 
    Spatial.obj.tls$TLS.Label = paste(i, Spatial.obj.tls$Label,sep = "_")
    seurat_tls_list[[i]] = Spatial.obj.tls
  }
}

#获取tls的meta.data和平均表达
{
  samples = c('BC23268_C1','BC23268_C2','BC23268_D1','BC23269_C1',
              'BC23269_D1','BC23270_E2','BC23272_E1','BC23272_E2','BC23287_C1',
              'BC23288_D2','BC23288_E1','BC23288_E2','BC23508_E2','BC23567_E2','BC23803_E1','BC23803_E2',
              'BC23895_C2','BC23903_C2','BC23903_D1','BC24105_C1',
              'BC24223_D2','BC24223_E1','BC24223_E2')
  #c("B3","B4","B5","C1","C2","C3","C5","D1","D2",'D6','E1','E2','E3','F1','F2','F3','G1','G2,'G3','H1','H2','H3')
  #c('CID4290','CID4465','CID4535','CID44971')
  genes_intersect <- NULL
  for (i in samples){
    print(i)
    Spatial.obj.tls = seurat_tls_list[[i]]
    DefaultAssay(Spatial.obj.tls) = "Spatial"
    genes_intersect = intersect(genes_intersect, rownames(Spatial.obj.tls))
  }
  
  
  pca_df_list = list()
  tls_expr_list = list()
  tls_frac_list = list()
  
  celltype_level = c("B", "Plasma","CD4T","CD8T","Tprolif", 
                     "Mono.Macro", "DC", "SMC","Endothelial","Fibroblasts", "Malignant","T_B","B_T")
  
  for (sample in c('BC23268_C2','BC23269_C1','BC23268_D1',
                   'BC23269_D1','BC23270_E2','BC23272_E1','BC23272_E2','BC23287_C1',
                   'BC23288_D2','BC23288_E1','BC23288_E2','BC23508_E2','BC23567_E2','BC23803_E1','BC23803_E2',
                   'BC23895_C2','BC23903_C2','BC23903_D1','BC24105_C1',
                   'BC24223_D2','BC24223_E1','BC24223_E2')
  ) {
    # c("B3","B4","B5","C1","C3","C5","D1","D2",'D6','E1','E2','E3','F1','F2','F3','G1','G3','H1','H2','H3')
    #c('CID4290','CID4465','CID4535')
    print(sample)
    Spatial.obj.tls = seurat_tls_list[[sample]]
    DefaultAssay(Spatial.obj.tls) =    "Spatial" 
    expr.tls = GetAssayData(Spatial.obj.tls, slot = "counts")[genes_intersect,]
    TLS_spot_list = split(colnames(Spatial.obj.tls), Spatial.obj.tls$TLS.Label)
    TLS_ct_frac = sapply(TLS_spot_list, function(x){
      tls_deconv_df = Spatial.obj.tls@meta.data[x, celltype_level]#
      tls_deconv_df_top = apply(tls_deconv_df, 1, function(x){
        x_top = x
        x_top[x < sort(x, decreasing = TRUE)[5]] = 0
        return(x_top)
      })
      tls_deconv_df_top = t(tls_deconv_df_top)
      tls_deconv_df_top = as.data.frame(tls_deconv_df_top)
      tls_deconv_df_top = tls_deconv_df_top/rowSums(tls_deconv_df_top)
      return(colMeans(tls_deconv_df_top))
    })
    
    TLS_expr = sapply(TLS_spot_list, function(x){
      return(rowMeans(expr.tls[,x]))
    })
    
    tls_expr_list[[sample]] = TLS_expr
    tls_frac_list[[sample]] = TLS_ct_frac
  }
  tls_expr_all = do.call("cbind", tls_expr_list)
  RNA.res <- CreateSeuratObject(tls_expr_all,project = "BRCA_32572199",assay = "Spatial") 
  #BRCA_32572199  <- RNA.res
  #fac_34650042 <- tls_frac_list
  fac_34650042_all = do.call("cbind",fac_34650042)
  colnames(fac_34650042_all) = paste0("BRCA34650042@_",colnames(fac_34650042_all))
}

#Integrate BRCA all tls object
tls_frac_all = cbind(fac_32572199_all,fac_34650042_all,fac_34493872_all)
tls_frac_all = t(tls_frac_all)
tls_frac_all=as.data.frame(tls_frac_all)
head(tls_frac_all)
RNA.res <- merge(x=BRCA_32572199 ,y=c(BRCA_34650042,BRCA_34493872),
                 add.cell.ids = c("BRCA32572199@","BRCA34650042@","BRCA34493872@"))
str(RNA.res)
RNA.res@meta.data <- cbind(RNA.res@meta.data,tls_frac_all[rownames(RNA.res@meta.data),])
saveRDS(RNA.res,"/fs/home/yuejiali/ST/BRCA/Analysis/TLS_Finder/BRCA_tls_20230213_res.rds")


