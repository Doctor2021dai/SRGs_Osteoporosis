rm(list = ls())
library(ConsensusClusterPlus)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(survival)
library(survminer)
library(Rtsne)
library(ggplot2)
library(ggsci)
library(pheatmap)

feExpr <- read.table("02DEGs/GSE56815_DEGs_SSGs_exprSet.txt",sep = "\t",header = T,row.names = 1,check.names = F)

feExpr <- feExpr[,Tumor]

feExpr_matrix <- as.matrix((feExpr) )
max(feExpr_matrix);min(feExpr_matrix)

# The result will output the classification of k from 2 to 7, the clustering method uses hc, 
# the sampling ratio is 0.8, and finally the png image is output
Cluster  <- ConsensusClusterPlus(feExpr_matrix,
                                maxK = 7, # Maximum number of clusters
                                reps = 1000, # Number of times to resample
                                pItem = 0.8, # The resampling ratio of the samples
                                pFeature = 1,
                                clusterAlg = "km", # The clustering algorithm used, you can choose "hc"(hclust), "pam", "km"(k-means)
                                distance = "euclidean",  # The method of calculating the distance, you can choose pearson, spearman, euclidean, binary, maximum, canberra, minkowski
                                seed = 2020,
                                title = "03ConsensusCluster",
                                plot = "png")

maxK = 7
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
for(i in Kvec){
    M = Cluster[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
}#end for i
# The optimal K
optK = Kvec[which.min(PAC)]
optK
# 
# redraw the heatmap
annCol <- data.frame(Cluster = paste0("Cluster",
                                      Cluster[[optK]][["consensusClass"]]),
                     row.names = colnames(feExpr_matrix))
head(annCol)

mycol <- c("#6b86d7","#c7e2a8")
annColors <- list(Cluster = c("Cluster1" = mycol[1],
                              "Cluster2" = mycol[2]))
                           
heatdata <- Cluster[[optK]][["consensusMatrix"]]
dimnames(heatdata) <- list(colnames(feExpr_matrix), colnames(feExpr_matrix))
heatdata[1:3, 1:3]
# draw heatmap
pdf(file="03ConsensusCluster/ConsensusCluster.pdf", width = 5,height = 4.5)
pheatmap(mat = heatdata,
         color = colorRampPalette((c("white", "steelblue")))(100),
         border_color = NA,
         annotation_col = annCol,
         annotation_colors = annColors,
         show_colnames = F,
         show_rownames = F)
dev.off()

# Get the typing result
Cluster_res <- annCol %>%
    as.data.frame() %>%
    rownames_to_column("id")
table(Cluster_res$Cluster)
write.table(Cluster_res,file = "03ConsensusCluster/GSE56815_Cluster_sampleId.txt",row.names = F,quote = F,sep = "\t")



# Load the Rtsne package
library(Rtsne)
set.seed(6)
# Using the Rtsne function for tSNE dimensionality reduction analysis
feExpr1 <- as.matrix(t(feExpr_matrix))
tsne_out <- Rtsne(feExpr1,pca=FALSE,dims=2,
                  perplexity=3,theta=0.0) # Run TSNE
plot(tsne_out$Y,col=factor( Cluster_res$Cluster),asp=1)
tsne_plot  = data.frame(tSNE1  =tsne_out$Y[,1], tSNE2  = tsne_out$Y[,2],Cluster = Cluster_res$Cluster)
tsne_plot$Cluster <- factor(tsne_plot$Cluster ,labels = c("Cluster1","Cluster2"))
head(tsne_plot)

ggplot(tsne_plot,aes(x=tSNE1,y=tSNE2))+
    geom_point(aes(fill = Cluster),shape = 21,color = "black")+
    scale_fill_manual(values=mycol[1:2]) +
    stat_ellipse(aes(color = Cluster,fill = Cluster),
                 geom = 'polygon',
                 level = 0.95,
                 alpha = 0.3,
                 linetype = 2)+
    scale_color_manual(values =c("#7992db","#d7eac2"))+
    scale_fill_manual(values =c("#7992db","#d7eac2"))+
    theme_classic()+
    theme(legend.position = "top")

ggsave("03ConsensusCluster/GSE56815_Cluster_tSNE.pdf",width = 4.2,height = 4.3)
