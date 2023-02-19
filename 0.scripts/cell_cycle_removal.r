library(Seurat)
library(ggplot2)
cell_cycle_removal <- function(seurat, cc_gene, mode = 'part', 
                               threshold, per.to.remove=0.1, nPCs=15, assay){
  # cell cycle removal
  df<-Loadings(seurat)
  # calculate cell cycle PC
  dfcc<-df[which(rownames(df) %in% cc_gene),]
  dfcc<-data.frame(vc=colSums(abs(dfcc)), pc=colnames(dfcc))
  dfcc$pc<-factor(dfcc$pc, levels = dfcc$pc)
  
  # manually check cell cycle PC
  p<-ggplot(dfcc,aes(pc,vc))+geom_bar(stat = 'identity')+
    ggtitle('Cell cycle contribution')+theme_bw()+
    theme(plot.title = element_text(size = 14, hjust=0.5,face = 'bold'),
          axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
  if(mode=='part')
  {
    return(p)
  }
  
  # identify PC related to cell cycle
  PC_cc<-as.character(dfcc[which(dfcc$vc>threshold),'pc'])
  cc_related<-c()
  for(i in PC_cc){
    top<-names(sort(abs(df[,i]), decreasing = T))[1:(per.to.remove*dim(df)[1])]
    cc_related<-union(cc_related, top)
  }
  if ('integrated' %in% Assays(seurat)){
    DefaultAssay(seurat) <- 'integrated'
  }else{
    DefaultAssay(seurat) <- 'RNA'
  }
  
  hvg_removeCC <- setdiff(VariableFeatures(seurat), cc_related)
  print(length(hvg_removeCC))
  # with new hvg set, rerun PCA and every step where PCA involves
  seurat <- RunPCA(seurat, features =hvg_removeCC)
  print(c('using ',nPCs))
  seurat <- FindNeighbors(seurat, dims = 1:nPCs)
  seurat <- RunUMAP(seurat, dims = 1:nPCs)
  return(seurat)
}
