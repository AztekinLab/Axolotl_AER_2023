library(clusterProfiler)
library(ggplot2)

func_enri <- function(genes, db, organism, suffix){

    # convert gene symbol to gene id (enrichment function only takes gene id)
    genes<-AnnotationDbi::select(db,
                   columns = c("SYMBOL","ENTREZID"),
                   keytype = "SYMBOL",
                   keys = genes)
    genes<-genes$ENTREZID
    
    # GO term enrichment analysis  
    ego<-enrichGO(OrgDb = db,
                  gene = genes,
                  ont = "BP",
                  pvalueCutoff = 0.01,
                  readable = TRUE)
    write.table(ego,
                paste0("GO-enrich_",suffix,".txt"),
                sep = "\t",
                row.names = F)
    pathway <- as.data.frame(ego)
    pathway = pathway[1:20,]
    pathway = pathway[order(pathway$pvalue,
                            decreasing = T),]
    tt <- factor(pathway$Description, 
                 levels = unique(pathway$Description))
    pp = ggplot(pathway, aes(-1*log10(p.adjust), tt))
    pbubble = pp + geom_point(aes(size = Count, color = -1*log10(p.adjust)))+
      scale_colour_gradient(low = "blue",high = "red") +
      theme_bw() +scale_size(range = c(2, 10)) +
      labs(title = paste("GO-enrich_",suffix,sep = ""))+
      ylab("") + 
      xlim(min(-1*log10(pathway$p.adjust))*0.9, max(-1*log10(pathway$p.adjust))*1.1)+
      theme(axis.text.y = element_text(size = 20, 
                                       family = "Helvetica", 
                                       color = "black",  
                                       angle = 0),
            axis.text.x = element_text(size = 15, 
                                       family = "Helvetica", 
                                       color = "black", 
                                       angle = 0))+
      theme(legend.title = element_text(color = "black", 
                                        size = 20, 
                                        family = "Helvetica"))+
      theme(legend.text = element_text(color="azure4", 
                                       size = 15,  
                                       family = "Helvetica"))+
      theme(axis.title = element_text(color="black", 
                                      size=16, 
                                      family = "Helvetica"))+
      theme(legend.key.size = unit(1.1,'cm'))
    ggsave(file = paste0("GO-enrich_",suffix,".pdf"), 
           plot = pbubble, 
           width = 20, 
           height = 10)
    
    # KEGG enrichment analysis  
    ekk<-enrichKEGG(gene=genes, organism = organism, pvalueCutoff=0.05)
    write.table(ekk, paste0("KEGG-enrich_",suffix,".txt"), sep="\t", row.names = F)
    pathway<-as.data.frame(ekk)
    pathway = pathway[1:20,]
    pathway = pathway[order(pathway$pvalue,decreasing=T),]
    tt <- factor(pathway$Description, levels=unique(pathway$Description))
    
    pp = ggplot(pathway,aes(-1*log10(p.adjust),tt))
    pbubble = pp +geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
      scale_colour_gradient(low="blue",high="red") +theme_bw() +
      scale_size(range=c(2, 10)) +
      labs(title = paste("KEGG-enrich_",suffix,sep=""))+ylab("")+
      theme(axis.text.y = element_text(size = 20, family = "Helvetica", color = "black", angle = 0),
            axis.text.x = element_text(size = 15, family = "Helvetica", color = "black",  angle = 0))+
      theme(legend.title = element_text(color="black", size=20, family = "Helvetica"))+
      theme(legend.text = element_text(color="azure4", size = 15,  family = "Helvetica"))+
      theme(axis.title = element_text(color="black", size=16, family = "Helvetica"))+
      theme(legend.key.size=unit(1.1,'cm'))
    ggsave(file=paste0("KEGG-enrich_",suffix, ".pdf"), plot=pbubble, width=15, height=9)
  }
  
