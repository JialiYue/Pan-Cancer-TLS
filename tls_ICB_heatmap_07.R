rm(list = ls())
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(GSVA))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(DOSE))
suppressMessages(library(GSEABase))
suppressMessages(library(stringr))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))
suppressMessages(library(rlang))
#suppressMessages(library(survminer))
#suppressMessages(library(survival))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(GeoTcgaData))

setwd("/fs/home/yuejiali/ST/BRCA/Analysis/TLS_Finder/")
RNA.res <- readRDS('/fs/home/yuejiali/ST/BRCA/Analysis/TLS_Finder/BRCA_tls_state_res.rds')
Idents(RNA.res) <- RNA.res@meta.data$state.cluster
marker_gene<- FindAllMarkers(object = RNA.res,assay = "Spatial",  only.pos = TRUE, min.pct = 0.1, thresh.use = 0.25)

ICB_files=list.files("/fs/home/yuejiali/ICB/Immunotherapy_Dataset/Tiger_ICB_Data/")
study <- unlist(lapply(strsplit(ICB_files[grep('.Rds',ICB_files)],'\\.'),function(x) x[1]))
study

df <- list()
df_1 <- list()
for(i in c(1:8,10:13,15:24,26:length(study))){
  print(i)
  clin<-read.table(paste0('/fs/home/yuejiali/ICB/Immunotherapy_Dataset/Tiger_ICB_Data/',study[i],'.Response.tsv'),header = TRUE,sep="\t",fill = TRUE)
  clin<-clin[,c("sample_id",'response_NR')]
  colnames(clin)<-c("Patient","Response")
  R<-clin$Patient[which(clin$Response=="R")]
  NR<-clin$Patient[which(clin$Response=="N")]
  ICB<- as.data.frame(readRDS(paste0("/fs/home/yuejiali/ICB/Immunotherapy_Dataset/Tiger_ICB_Data/",study[i],".Response.Rds")))
  rownames(ICB) <- ICB$GENE_SYMBOL
  ICB <- ICB[,2:dim(ICB)[2]]
  ICB[is.na(ICB)] <- 0
  clus <- c('cluster.1','cluster.2','cluster.3','cluster.4')
  for(j in c(4)){
    marker_gene <-read.table('/fs/home/yuejiali/ST/BRCA/Analysis/TLS_Finder/BRCA_tls_act_4_20230412_markers.tsv',header = TRUE,sep="\t")
    top <-  marker_gene 
    gene_used <- top$gene[which(top$cluster==clus[j])]
    cluster_sig <- gene_used #c("VWF","FLI1","ACKR1","PRCP","HLA-DRA","HLA-DQA1","HLA-DQA2")#
    cluster_gmt <- list(cluster_sig=cluster_sig)
    ICB<-as.matrix(ICB)
    
    res <- gsva(ICB,cluster_gmt,method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
    res<-data.frame(t(res))
    res$Patinet<-unlist(lapply(rownames(res),function(x) gsub("\\.","-",x)))
    res$Response<-"R"
    res$Response[res$Patinet %in% R]<-"R"
    res$Response[res$Patinet %in% NR]<-"NR"
    p <- boxplot(cluster_sig ~ Response, data=res)
    df[[study[i]]] <- data.frame("Response"=p$names,"Median"=p$stats[3,])
    rela.gsva <- df[[study[i]]][which(df[[study[i]]]$Response=="R"),"Median"] - df[[study[i]]][which(df[[study[i]]]$Response!="R"),"Median"]                           
    df_1[[study[i]]] <-  data.frame("Cor_R"=rela.gsva,"Cancer_type"=study[i],"Celltype"=clus[j])
    stat.test=NULL                       
    if(length(unique(res$Response)) >=2){
      stat.test <- res %>%
        wilcox_test(cluster_sig ~ Response)#,alternative = "greater"，
      result <- stat.test$p %>% 
        signif(digits = 3)
      df_1[[study[i]]]$P_value <-  result
    } else 
    {
      df_1[[study[i]]]$P_value <-  1 
    }
  }
}

#需要手动把每个cluster的ICB结果合并在一起
ICB_all = do.call("rbind", df_1)
clus4 <- ICB_all
ICB_all

tls_ICB_all <- rbind(clus1[rownames(clus1),],clus2[rownames(clus1),],clus3[rownames(clus1),],clus4[rownames(clus1),])
#tls_ICB_all <- rbind(clus1,clus2,clus3)
#tls_ICB_all <- tls_ICB_all[-which(is.na(tls_ICB_all$Celltype)),]
dim(tls_ICB_all)
head(tls_ICB_all)

P_value_heatmap=matrix(tls_ICB_all$P_value,ncol=4,byrow=FALSE)
pmt=P_value_heatmap
pmt

Cor_value_heatmap=matrix(tls_ICB_all$Cor_R,ncol=4,byrow=FALSE)
Cor_value_heatmap

colnames(Cor_value_heatmap)=unique(as.vector(unlist(tls_ICB_all$Celltype)))
rownames(Cor_value_heatmap)=unique(tls_ICB_all$Cancer_type)

if (!is.null(pmt)){
  ssmt <- pmt < 0.01
  pmt[ssmt] <-'**'
  smt <- pmt < 0.05 & pmt > 0.01
  pmt[smt] <- '*'
  pmt[!ssmt&!smt]<- ''
} else {
  pmt <- FALSE
}

colnames(pmt)=colnames(Cor_value_heatmap)
rownames(pmt)=rownames(Cor_value_heatmap)

#这里的cell type order 可能得随着pancancer 中的结果改变。 

library(RColorBrewer)

color.use <- tryCatch({
  RColorBrewer::brewer.pal(n = 10, name = "RdBu")
}, error = function(e) {
  (scales::viridis_pal(option = "Spectral", direction = -1))(10)
})

color.use
setwd("/fs/home/yuejiali/ST/BRCA/Analysis/Plot/")
color.use <-c('#67001F','#B2182B','#D6604D','#F4A582','#FDDBC7','#ffe8d6','#D1E5F0','#4393C3','#2166AC','#053061')
library(pheatmap)
a <-c('cluster.1','cluster.2','cluster.3','cluster.4')
p <- pheatmap(Cor_value_heatmap[,a],scale = "none",cluster_row = TRUE, cluster_col = TRUE, border=NA,
              display_numbers = pmt[,a],fontsize_number = 12, number_color = "white",
              cellwidth = 10, cellheight =10,color=rev(color.use))
ggsave(file.path("BRCA_tls_pearson_median_cor_p_value_4act.clustering_20230412.pdf"),p,height = 8,width = 5)


