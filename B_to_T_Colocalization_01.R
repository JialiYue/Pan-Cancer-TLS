rm(list = ls())
library(Seurat)
library(ggplot2)
library(ggridges)
library(reshape2)
library(RColorBrewer)
library(ggsci)
library(scatterpie)
library(gcookbook)
library(magrittr)
library(egg)
library(ggpubr)
myPalette <- colorRampPalette(brewer.pal(n = 9, name = "GnBu"))

samples<- c("1142243F","1160920F","CID4465","CID44971")#"1142243F","1160920F","CID4465","CID44971"；"CID4290","CID4535"
t <- list()
b <- list()
s <- list()
{
for(i in samples){
  print(i)
  #载入STRIDE的deconvolution的celltype prop信息
  deconv_df = read.table(paste0("/fs/home/yuejiali/ST_data/BRCA_34493872/BRCA_Alex/TNBC/top120/",i,"_spot_celltype_frac.txt"), 
                         row.names = 1, header = TRUE, sep = "\t")
  #rownames(deconv_df) <- unlist(lapply(strsplit(rownames(deconv_df),"X"),function(x) x[2]))
  colnames(deconv_df) = c('B','CD4T','CD8T','DC','Endothelial','Fibroblasts','Malignant','Mono.Macro','Plasma','SMC','Tprolif')
  #载入barcode的坐标信息
  st_loc_df = read.csv(file.path(paste0("/fs/home/yuejiali/ST_data/BRCA_34493872/ST/",i,"/",i,"_spatial/tissue_positions_list.csv")),  header = TRUE )#,sep = "\t"
  #rownames(st_loc_df) <- paste0(st_loc_df$x,"x",st_loc_df$y)
  rownames(st_loc_df) <- st_loc_df$barcode  
  #把celltype prop跟坐标相对应
  st_loc_deconv = cbind(st_loc_df, deconv_df[rownames(st_loc_df), ])
  st_loc_deconv <- st_loc_deconv[!is.na(st_loc_deconv$B),] 
  st_loc_deconv$T = rowSums(st_loc_deconv[,c("CD4T", "CD8T")])#"Plasma"
  st_loc_deconv$Bcell = rowSums(st_loc_deconv[,c("B",'Plasma')])
  setwd("/fs/home/yuejiali/ST_data/BRCA_34493872/Multiple_distribution")  
  #计算每个spot与临近spot的B to T colocalization, Visium 55um的spot呈六边形阵列，不是Visium的，需适当调整
  st_loc_deconv$T_B <- 0
  st_loc_deconv$B_T <- 0
  for(m in 1:dim(st_loc_deconv)[1]){
    x <- st_loc_deconv[m,"col"]
    y <- st_loc_deconv[m,"row"]
    index <- c(which((st_loc_deconv$col == x-2) & (st_loc_deconv$row == y)),
               which((st_loc_deconv$col == x-1) & (st_loc_deconv$row == y-1)),
               which((st_loc_deconv$col == x-1) & (st_loc_deconv$row == y+1)),
               which((st_loc_deconv$col == x+2) & (st_loc_deconv$row == y)),
               which((st_loc_deconv$col == x+1) & (st_loc_deconv$row == y-1)),
               which((st_loc_deconv$col == x+1) & (st_loc_deconv$row == y+1)))
    T_B <- st_loc_deconv[m,"T"]*st_loc_deconv[m,"B"]
    B_T <- st_loc_deconv[m,"T"]*st_loc_deconv[m,"B"]
    for(n in index){
      T_B <- T_B + 0.5*st_loc_deconv[m,"T"]*st_loc_deconv[n,"B"]
      B_T <- B_T + 0.5*st_loc_deconv[m,"B"]*st_loc_deconv[n,"T"] 
    }    
    
    st_loc_deconv$T_B[m] <- T_B
    st_loc_deconv$B_T[m] <- B_T
  }
  write.table(st_loc_deconv[,c('T_B','B_T')],paste0('/fs/home/yuejiali/ST_data/BRCA_34493872/BRCA_Alex/decov/',i,
                                   '_t_and_B_colocal.tsv'),sep = '\t', row.names = TRUE, col.names = TRUE,quote = FALSE)
  
  df <- st_loc_deconv
  #选择适当的background的sd系数、T和B proportion作为cutoff，筛选出具有适当的B_T的colocalization和适当T和B cell prop的spot
  bmean <- mean(df$B_T)
  bsd <- sd(df$B_T)
  colo_min <- bmean + 0.5*bsd
  colo_min
  T_min = mean(df$T) - 0.5*sd(df$T)
  T_min
  B_min = mean(df$B) - 0.5*sd(df$B)
  B_min
  df2 <- st_loc_deconv[which((st_loc_deconv$B_T >=colo_min) &
                               # (st_loc_deconv$B_T <=colo_max) &
                               (st_loc_deconv$T >= T_min) &
                               # (st_loc_deconv$T <= T_max) &
                               (st_loc_deconv$B >= B_min) 
                             #(st_loc_deconv$B <= B_max)
  ),]
  p <-  ggplot(df,aes(x=col,y=-row,fill=B_T)) +
    #geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
    geom_point(shape = 21, colour = "white", size = 1.5, stroke = 0.5)+
    #coord_equal() + 
    #geom_point(data=df3, aes(x=x,y=-y), pch=21, fill=NA, size=4,colour="red", stroke=0.3)+
    geom_point(data=df2, aes(x=col,y=-row), pch=5, fill=NA, size=1.5,colour="#00afb9", stroke=0.3)+
    #coord_cartesian(expand=FALSE)+
    scale_fill_gradientn(colours = myPalette(100),limits=c(0,0.3))+                
    xlab("") +
    ylab("") +
    ggtitle(i)+
    labs(fill = "B_to_T_Colocalization")+
    theme_set(theme_bw(base_size = 14))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank()) 
  
  ggsave(file=file.path(paste0(i,"_TLS_candidate.png")),p, width = 6, height = 5)
  #a <- unique(rbind(df3[,c('y','x')],df2[,c('y','x')]))
  #通过构造矩阵并进行聚类来cut掉单个独立的spot，同上也需要根据spot的呈现形状进行调整
  a <- unique(df2[,c('row','col')])
  a$Spot <- rownames(a)
  mat <-matrix(data=0, nrow = max(a$row)+1, ncol = max(a$col)+1, byrow = FALSE, dimnames = list(as.character(c(0:max(a$row))),as.character(c(0:max(a$col)))))
  for(m in 1:dim(a)[1]){
    mat[as.character(a[m,"row"]),as.character(a[m,"col"])] <- 1
  }
  coor <- which(mat == 1, arr.ind=TRUE)   
  coor <- cbind(coor, clus = cutree(hclust(dist(coor, "maximum"), "single"), h = 1))
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
    if(length(check)>0){
      index <- c(index,xx)
    }
    
  }
  coor$row <- as.integer(coor$row)
  coor$col <- as.integer(coor$col) 
  Label <- merge(a,coor[index,])
  p <-  ggplot(df,aes(x=col,y=-row,fill=B_T)) +
    #geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
    geom_point(shape = 21, colour = "white", size = 1.5, stroke = 0.5)+
    #coord_equal() + 
    # geom_point(data=df2, aes(x=col,y=-row), pch=21, fill=NA, size=1.5,colour="red", stroke=0.3)+
    geom_point(data=Label, aes(x=col,y=-row), pch=5, fill=NA, size=1.5,colour="#00afb9", stroke=0.3)+
    #coord_cartesian(expand=FALSE)+
    scale_fill_gradientn(colours = myPalette(100),limits=c(0,0.3))+                
    xlab("") +
    ylab("") +
    ggtitle(paste0(i," labeled spots"))+
    labs(fill = "B_to_T_Connectivity")+
    theme_set(theme_bw(base_size = 14))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank()) 
  ggsave(file=file.path(paste0(i,"_TLS_cell_fraction_label_mean_sd_T_0.5_","B_0.5_",
                               "coloc_plus_0.5",".png")),p, width = 6, height = 5)
  write.table(Label,paste0(i,'_TLS_label__mean_sd_T_0.5_','B_0.5_',
                           'coloc_plus_0.5','.tsv'),sep = '\t', row.names = FALSE, col.names = TRUE,quote = FALSE)
}
  
}  

{
  #Visium 100um的spot呈正方矩阵
  for(m in 1:dim(st_loc_deconv)[1]){
    x <- st_loc_deconv[m,"x"]
    y <- st_loc_deconv[m,"y"]
    index <- c(which((st_loc_deconv$x == x) & (st_loc_deconv$y == y-1)),
               which((st_loc_deconv$x == x) & (st_loc_deconv$y == y+1)),
               which((st_loc_deconv$x == x-1) & (st_loc_deconv$y == y)),
               which((st_loc_deconv$x == x-1) & (st_loc_deconv$y == y-1)),
               which((st_loc_deconv$x == x-1) & (st_loc_deconv$y == y+1)),
               which((st_loc_deconv$x == x+1) & (st_loc_deconv$y == y)),
               which((st_loc_deconv$x == x+1) & (st_loc_deconv$y == y-1)),
               which((st_loc_deconv$x == x+1) & (st_loc_deconv$y == y+1)))
    T_B <- st_loc_deconv[m,"T"]*st_loc_deconv[m,"B"]
    B_T <- st_loc_deconv[m,"T"]*st_loc_deconv[m,"B"]
    for(n in index){
      T_B <- T_B + 0.5*st_loc_deconv[m,"T"]*st_loc_deconv[n,"B"]
      B_T <- B_T + 0.5*st_loc_deconv[m,"B"]*st_loc_deconv[n,"T"] 
    }    
    
    st_loc_deconv$T_B[m] <- T_B
    st_loc_deconv$B_T[m] <- B_T
    
  }
}