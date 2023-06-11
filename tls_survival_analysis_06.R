rm(list = ls())
suppressMessages(library(survminer))
suppressMessages(library(survival))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(qpcR))
suppressMessages(library(magrittr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

cluster.1 <- c('BRCA34493872@_CID4465_1','BRCA34493872@_CID44971_1','BRCA34493872@_CID44971_13','BRCA34493872@_CID44971_14','BRCA34493872@_CID44971_15','BRCA34493872@_CID44971_2','BRCA34493872@_CID44971_3','BRCA34493872@_CID44971_4','BRCA34493872@_CID44971_6','BRCA34493872@_CID44971_8','BRCA34493872@_CID44971_9','BRCA34493872@_CID4535_2','BRCA34493872@_CID4535_8')
cluster.2 <- c('BRCA32572199@_BC23268_C2_2','BRCA32572199@_BC23268_D1_1','BRCA32572199@_BC23269_D1_1','BRCA32572199@_BC23270_E2_1','BRCA32572199@_BC23272_E1_2','BRCA32572199@_BC23272_E2_2','BRCA32572199@_BC23272_E2_4','BRCA32572199@_BC23287_C1_1','BRCA32572199@_BC23288_D2_3','BRCA32572199@_BC23288_E1_3','BRCA32572199@_BC23288_E1_4','BRCA32572199@_BC23288_E2_1','BRCA32572199@_BC23288_E2_2','BRCA32572199@_BC23288_E2_3','BRCA32572199@_BC23288_E2_4','BRCA32572199@_BC23567_E2_1','BRCA32572199@_BC23567_E2_2','BRCA32572199@_BC23803_E1_2','BRCA32572199@_BC23903_C2_1','BRCA32572199@_BC23903_D1_4','BRCA32572199@_BC24105_C1_1','BRCA34650042@_C1_1','BRCA34650042@_C2_3','BRCA34650042@_D1_1','BRCA34650042@_D1_2','BRCA34650042@_D1_3','BRCA34650042@_D2_1','BRCA34650042@_D2_2','BRCA34650042@_D6_1','BRCA34650042@_D6_2','BRCA34650042@_D6_3','BRCA34650042@_E1_1','BRCA34650042@_E1_2','BRCA34650042@_E1_3','BRCA34650042@_E1_4','BRCA34650042@_E2_1','BRCA34650042@_E2_2','BRCA34650042@_E2_3','BRCA34650042@_E2_5','BRCA34650042@_E2_6','BRCA34650042@_E2_8','BRCA34650042@_E3_1','BRCA34650042@_E3_2','BRCA34650042@_E3_3','BRCA34650042@_E3_4','BRCA34650042@_F1_1','BRCA34650042@_F1_2','BRCA34650042@_F2_1','BRCA34650042@_F2_3','BRCA34650042@_F3_1','BRCA34650042@_F3_2','BRCA34650042@_H1_2','BRCA34650042@_H1_5','BRCA34650042@_H2_2','BRCA34493872@_CID44971_16','BRCA34493872@_CID4535_6','BRCA34493872@_CID4535_7')
cluster.3 <- c('BRCA32572199@_BC23268_C1_1','BRCA32572199@_BC23268_C2_1','BRCA32572199@_BC23268_C2_3','BRCA32572199@_BC23268_C2_4','BRCA32572199@_BC23269_C1_1','BRCA32572199@_BC23269_C1_2','BRCA32572199@_BC23269_C1_3','BRCA32572199@_BC23269_D1_2','BRCA32572199@_BC23269_D1_3','BRCA32572199@_BC23269_D1_4','BRCA32572199@_BC23272_E1_1','BRCA32572199@_BC23272_E2_1','BRCA32572199@_BC23272_E2_3','BRCA32572199@_BC23272_E2_5','BRCA32572199@_BC23288_D2_1','BRCA32572199@_BC23288_D2_2','BRCA32572199@_BC23288_E1_1','BRCA32572199@_BC23288_E1_2','BRCA32572199@_BC23508_E2_1','BRCA32572199@_BC23803_E1_1','BRCA32572199@_BC23803_E1_3','BRCA32572199@_BC23803_E2_1','BRCA32572199@_BC23895_C2_1','BRCA32572199@_BC23895_C2_2','BRCA32572199@_BC23903_C2_2','BRCA32572199@_BC23903_C2_3','BRCA32572199@_BC23903_D1_1','BRCA32572199@_BC23903_D1_2','BRCA32572199@_BC23903_D1_3','BRCA32572199@_BC24223_D2_1','BRCA32572199@_BC24223_D2_2','BRCA32572199@_BC24223_D2_3','BRCA32572199@_BC24223_E1_1','BRCA32572199@_BC24223_E1_2','BRCA32572199@_BC24223_E1_3','BRCA32572199@_BC24223_E2_1','BRCA32572199@_BC24223_E2_2','BRCA32572199@_BC24223_E2_3','BRCA32572199@_BC24223_E2_4','BRCA32572199@_BC24223_E2_5','BRCA34650042@_C1_2','BRCA34650042@_C2_1','BRCA34650042@_C2_2','BRCA34650042@_C3_1','BRCA34650042@_C3_2','BRCA34650042@_C5_1','BRCA34650042@_D2_3','BRCA34650042@_E2_4','BRCA34650042@_E2_7','BRCA34650042@_E3_5','BRCA34650042@_E3_6','BRCA34650042@_F1_3','BRCA34650042@_F1_4','BRCA34650042@_F2_2','BRCA34650042@_F2_4','BRCA34650042@_F3_3','BRCA34650042@_G1_2','BRCA34650042@_G2_4','BRCA34650042@_G3_3','BRCA34493872@_CID4290_1','BRCA34493872@_CID4290_10','BRCA34493872@_CID4290_11','BRCA34493872@_CID4290_12','BRCA34493872@_CID4290_13','BRCA34493872@_CID4290_14','BRCA34493872@_CID4290_2','BRCA34493872@_CID4290_3','BRCA34493872@_CID4290_4','BRCA34493872@_CID4290_5','BRCA34493872@_CID4290_6','BRCA34493872@_CID4290_7','BRCA34493872@_CID4290_8','BRCA34493872@_CID4290_9','BRCA34493872@_CID44971_10','BRCA34493872@_CID44971_11','BRCA34493872@_CID4535_1','BRCA34493872@_CID4535_3')
cluster.4 <- c('BRCA34650042@_B3_1','BRCA34650042@_B4_1','BRCA34650042@_B5_1','BRCA34650042@_G1_1','BRCA34650042@_G1_3','BRCA34650042@_G2_1','BRCA34650042@_G2_2','BRCA34650042@_G2_3','BRCA34650042@_G2_5','BRCA34650042@_G2_6','BRCA34650042@_G2_7','BRCA34650042@_G3_1','BRCA34650042@_G3_2','BRCA34650042@_G3_4','BRCA34650042@_G3_5','BRCA34650042@_H1_1','BRCA34650042@_H1_3','BRCA34650042@_H1_4','BRCA34650042@_H2_1','BRCA34650042@_H3_1','BRCA34650042@_H3_2')
RNA.res@meta.data$act.cluster <- 'cluster.1'
RNA.res@meta.data$act.cluster[which(rownames(RNA.res@meta.data) %in% cluster.1)] <- 'cluster.1'
RNA.res@meta.data$act.cluster[which(rownames(RNA.res@meta.data) %in% cluster.2)] <- 'cluster.2'
RNA.res@meta.data$act.cluster[which(rownames(RNA.res@meta.data) %in% cluster.3)] <- 'cluster.3'
RNA.res@meta.data$act.cluster[which(rownames(RNA.res@meta.data) %in% cluster.4)] <- 'cluster.4'
Idents(RNA.res) <- as.factor(RNA.res@meta.data$act.cluster)

expmat<-as.matrix(expmat)

#marker_gene$cluster[which(marker_gene$cluster=='T3-FOS')] <-"T3-IL7R"

num<- c(150)

for(n in c(1:length(num))){   
  #paste0("/fs/home/yuejiali/subtype_anno/Subgroup/","Plasma_B","_DiffGenes.tsv")
  marker_gene<-read.table("/fs/home/yuejiali/ST/BRCA/Analysis/TLS_Finder/BRCA_tls_activated_markers_20230412.tsv",header = TRUE,sep="\t")
  a<-levels(RNA.res)
  data<-list()
  top_marker_gene <- marker_gene %>% group_by(cluster) %>% top_n(100, avg_logFC)
  for (i in c(1:length(a))){
    data1<- as.data.frame( top_marker_gene[which( top_marker_gene$cluster==a[i]),c("gene")])
    names(data1)<-a[i]
    data<-qpcR:::cbind.na(data1,data)
  }
  data<-data[,-4]
  TCGA_files=list.files("/fs/home/hanya/Project/Survival_cohorts/TCGA/")
  TCGA_correlation<<-NULL
  temp_list=list.files(paste0("/fs/home/hanya/Project/Survival_cohorts/TCGA/","BRCA"))
  cancer_expMat=readRDS(paste0("/fs/home/hanya/Project/Survival_cohorts/TCGA/","BRCA","/",temp_list[grep("RPKM",temp_list)]))
  Survival_info=t(read.table(paste0("/fs/home/hanya/Project/Survival_cohorts/TCGA/","BRCA","/",temp_list[grep("clin",temp_list)]),header = TRUE,sep=",",row.names = 1))
  #1 代表删失数据
  Survival_info<<-as.data.frame(Survival_info)
  time_index=match("_OS",colnames(Survival_info))
  statu_index=match("_OS_IND",colnames(Survival_info))
  marker_gene<-data
  clus<-colnames(data)
  survival_gsva<-NULL
  for (i in c(1:length(clus))) {
    print(i)
    survival_infor=data.frame(time=Survival_info[,time_index],statu=Survival_info[,statu_index])
    rownames(survival_infor)=gsub("\\.","-",rownames(Survival_info))
    DEgenes_cluster=marker_gene[[as.character(clus[i])]]
    inter_gene=intersect(DEgenes_cluster,rownames(cancer_expMat))
    mini_index=match(inter_gene,rownames(expmat))
    TCGA_index=match(inter_gene,rownames(cancer_expMat))
    mini_sample_index=match(rownames(RNA.res@meta.data)[which(RNA.res@meta.data$act.cluster == clus[i])],colnames(expmat))
    cluster_average=rowMeans(expmat[mini_index,mini_sample_index])
    cluster_average=log2(cluster_average/mean(expmat) +1)
    cor_result=apply(cancer_expMat[TCGA_index,],2,function(m){unlist(cor.test(m,cluster_average))[4]})
    length(intersect(rownames(survival_infor),names(cor_result)))
    sample=intersect(rownames(survival_infor),names(cor_result))
    index1=match(sample,rownames(survival_infor))
    index2=match(sample,names(cor_result))
    if(length(survival_gsva) == 0){
      survival_gsva<<-cbind(survival_infor[index1,],cor_result[index2])
    }else{
      survival_gsva<<-cbind(survival_gsva,cor_result[index2])
    }
    colnames(survival_gsva)[dim(survival_gsva)[2]]<- clus[i]
  }
}