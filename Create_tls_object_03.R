rm(list = ls())
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(plotly)
  library(STutility)
  library(zeallot)
  library(ggplot2)
  library(future)
  library(ggpubr)
  library(MAESTRO)
})

#利用tls finder result得到Seurat object
samples<- c("1160920F","CID4465","CID44971","CID4290","CID4535")
tissue_paths <- list()
position <- list()
RNA.res <- list()
TLS <- list()
for (i in 1:length(samples)) {
  setwd("/fs/home/yuejiali/ST_data/BRCA_34493872/ST/filtered_count_matrices")
  RNA.res[[i]] <- Read10X_h5(paste0(samples[i],"_gene_count.h5"))
  RNA.res[[i]] <- CreateSeuratObject(counts = RNA.res[[i]],assay = "Spatial", project = samples[i], min.cells = 0, min.features = 200) 
  tissue_paths[[i]] <- paste0("/fs/home/yuejiali/ST_data/BRCA_34493872/ST/",samples[i],"/",samples[i],"_spatial/tissue_positions_list.csv")
  position[[i]] <- read.csv(tissue_paths[[i]], header = TRUE )
  rownames(position[[i]]) <- position[[i]]$barcode
  colnames(position[[i]]) <- c("barcode","tissue","y","x","Y","X")
  RNA.res[[i]]@meta.data <-cbind(RNA.res[[i]]@meta.data,position[[i]][rownames(RNA.res[[i]]@meta.data),])
  #载入tls finder result
  df2 <- read.table(paste0("~/ST_data/BRCA_34493872/Multiple_distribution/TLS_Result/",samples[i],"_tls_loc_label_pixel.txt"),header = TRUE,sep="\t")
  #cut掉单独的spots
  a <- unique(df2[,c('row','col')])
  a$Spot <- rownames(a)
  mat <-matrix(data=0, nrow = max(a$row)+1, ncol = max(a$col)+1, byrow = FALSE, dimnames = list(as.character(c(0:max(a$row))),as.character(c(0:max(a$col)))))
  for(m in 1:dim(a)[1]){
    mat[as.character(a[m,"row"]),as.character(a[m,"col"])] <- 1
  }
  coor <- which(mat == 1, arr.ind=TRUE)   
  coor <- cbind(coor, clus = cutree(hclust(dist(coor,"maximum" ), "single"), h = 2))#"manhattan"
  coor <- as.data.frame(coor)
  coor$row <- coor$row -1
  coor$col <- coor$col -1
  colnames(coor) <- c("row","col","Label")
  index <- NULL
  for (xx in 1:dim(coor)[1]){
    r = coor$row[xx]
    c = coor$col[xx]
    check <- c(which((coor$col == c-2) & (coor$row == r)),
               which((coor$col == c-1) & (coor$row == r-1)),
               which((coor$col == c-1) & (coor$row == r+1)),
               which((coor$col == c+2) & (coor$row == r)),
               which((coor$col == c+1) & (coor$row == r-1)),
               which((coor$col == c+1) & (coor$row == r+1)))
    if(length(check)>0){#CID44971 〉1
      index <- c(index,xx)
    }
    
  }
  coor$row <- as.integer(coor$row)
  coor$col <- as.integer(coor$col) 
  Label <- merge(a,coor[index,])
  df2 <- df2[,c("col","row","Y","X","Label_pixel")]
  # df2 <- df2[,c("x","y","pixel_x","pixel_y","Y","X","pixel")]
  tls <- merge(Label,df2,all.x = TRUE,by=c('col','row'))
  TLS[[i]] <- tls
  TLS[[i]] <- TLS[[i]][,c('col','row','Label')]
  colnames(TLS[[i]]) <- c('x','y','Label')
  #标记出哪些是tls
  RNA.res[[i]]@meta.data<- merge(RNA.res[[i]]@meta.data,TLS[[i]],by= c("x","y"),all.x = TRUE)
  RNA.res[[i]]@meta.data$isTLS <- "Stroma"
  RNA.res[[i]]@meta.data$isTLS[which(!is.na(RNA.res[[i]]@meta.data$Label))] <- "TLS"
  #cell<- RNA.res[[i]]@meta.data$barcode[which(is.na(RNA.res[[i]]@meta.data$barcode))]
  #RNA.res[[i]] <- subset(RNA.res[[i]],cells = cell)
  rownames(RNA.res[[i]]@meta.data) <- RNA.res[[i]]@meta.data$barcode
  saveRDS(RNA.res[[i]],paste0("/fs/home/yuejiali/ST/BRCA/TLS_Predict/BRCA_34493872/TLS_Finder/",samples[i],"_min_Spatial_Object.rds"))
}

#label predicted tls
{
  RNA.res <- list()
  TLS <- list()
  sub <- list()
  min_TLS <- c('MS4A1','CXCR5','SELL','CD19','LTB','CD79B','CD37','CD79A','TCL1A') 
  #samples<- 'CID44971'#c('1142243F','1160920F','CID4290','CID4465','CID44971','CID4535')
  #'1142243F','1160920F','CID4290','CID4465',,'CID4535'
  for (i in 1:length(samples)) {
    
    setwd("/fs/home/yuejiali/ST/BRCA/TLS_Predict/BRCA_34493872/TLS_Finder/")
    RNA.res[[i]] <- readRDS(paste0(samples[i],"_min_Spatial_Object.rds"))
    RNA.res[[i]] <- AddModuleScore(
      object = RNA.res[[i]],
      features = list(intersect(min_TLS, rownames(RNA.res[[i]]@assays$Spatial))),
      nbin=10,name = 'TLS_feature')
    expr.data <- RNA.res[[i]]@meta.data
    if(length(names(table(RNA.res[[i]]@meta.data$isTLS))) > 1){
      sub[[i]] <- subset(RNA.res[[i]], isTLS == "TLS")
      df2<- sub[[i]]@meta.data
      df2$Label <- as.character(df2$Label)
      df <- expr.data
      p <-  ggplot(df,aes(x=X,y=-Y,fill=TLS_feature1)) +
        #geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
        geom_point(shape = 21, colour = "black", size = 1.5, stroke = 0.5)+
        #coord_equal() + 
        geom_point(data=df2, aes(x=X, y=-Y,color=Label), pch=21,  size=1.5, stroke=1)+
        #coord_cartesian(expand=FALSE)+
        scale_fill_gradientn(colours = myPalette(100),limits=c(-1,1.2))+                
        xlab("") +
        ylab("") +
        ggtitle(samples[i])+
        labs(fill = "TLS_feature")+
        theme_set(theme_bw(base_size = 14))+
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"),
              axis.text = element_blank(),
              axis.ticks = element_blank()) 
      setwd("/fs/home/yuejiali/ST/BRCA/TLS_feature_score/BRCA_34493872/")
      ggsave(file=file.path(paste0(samples[i],"_TLS_feature_Spotplot_min_TLS.png")),p, width = 6, height = 5)
    } else
      
    {print(samples[i])}
    
  } 
}