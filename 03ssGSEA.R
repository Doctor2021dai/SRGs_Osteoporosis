# BiocManager::install("GSVA")

library(GSVA)
library(dplyr)
Sys.setenv(LANGUAGE = "en") # Display English error message
options(stringsAsFactors = FALSE) # Prohibit chr from converting to factor

# Read in immune cell marker genes
immunity <- read.table("04immune/CellMarker.txt", comment.char = "#",sep = "\t",row.names = 1)


geneSet <- list()

for (i in rownames(immunity)) {
  geneSet[[i]] <- as.character( immunity[i,][which(immunity[i,]!="")]) 
  
}

# Use ssGSEA to quantify the level of infiltration
mycol <- c("#6b86d7","#c7e2a8")

library("limma")   
GSE_id <- "GSE56815"
expFile="01rawData/GSE56815_symbol_exprSet.txt" 
expr=read.table(expFile,sep="\t",header=T,check.names=F,row.names = 1)

expr<- expr[,Tumor]

expr[1:3,1:5]
expr_out <- expr %>% rownames_to_column("id")
write.table(expr_out,file=paste0("04immune/" ,GSE_id,"_lowBMD_exprSet.txt"),row.names=F,quote=F,sep = "\t")

class(expr)
expr_matrix <- as.matrix(expr)

# The following step completes ssGSEA
ssgsea <- as.data.frame(t(gsva(expr_matrix, geneSet, method = "ssgsea",kcdf='Gaussian',abs.ranking=TRUE)))
# Define the ssGSEA score correction function Min-Max standardization 
#refers to linear transformation of the original data, mapping the value to [0, 1]
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
# Correct the ssGSEA score
ssgseaOut=normalize(ssgsea)
# str(ssgseaOut)

ssgseaOut_df= ssgseaOut %>% tibble::rownames_to_column(var = "id")
ssgseaOut_df[1:3,1:6]
write.table(ssgseaOut_df,file="04immune/GSE56815_lowBMD_ssGSEA_score.txt",sep="\t",quote=F,row.names=F)


gsva <- t(ssgseaOut)


library(pheatmap)
Cluster_res <- read.table("03ConsensusCluster/GSE56815_Cluster_sampleId.txt",header = T,sep = "\t",check.names = F)
Cluster1 <- Cluster_res$id[Cluster_res$Cluster =="Cluster1"]
Cluster2 <- Cluster_res$id[Cluster_res$Cluster =="Cluster2"]

colnames(gsva) <- substr(colnames(gsva),1,16)
gsva <- gsva[,c(Cluster1,Cluster2)]
annotation_col = data.frame(
  Type = factor(rep(c("Cluster1","Cluster2"), c(length(c(Cluster1)) ,length(Cluster2))),levels = c("Cluster1","Cluster2"))
)
rownames(annotation_col) = colnames(gsva)
head(annotation_col)

ann_colors = list(
  Type = c(`Cluster1` =mycol[1], `Cluster2` = mycol[2] ))
head(ann_colors)



gsva <- gsva[apply(gsva, 1, function(x) sd(x)!=0),]

loc <- order(annotation_col$Type,colSums(gsva),decreasing = T)


pdf(file="04immune/GSE56815_lowBMD_ssGSEA_score_heatmap.pdf",height=8,width=10)
pheatmap(gsva[,loc],
         scale = "row", 
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = F, 
         cluster_cols =F,
         border=FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

dev.off()


# boxplot
boxdata <-  as.data.frame(t(gsva))
risk <- read.table("03ConsensusCluster/GSE56815_Cluster_sampleId.txt",header = T,sep = "\t",check.names = F)
boxdata$id <- rownames(boxdata)
boxdata <- merge(boxdata,risk,by="id")


write.table(boxdata, "04immune/GSE56815_Cluster_ssGSEA_score_group.txt", quote = F, row.names = F,sep = "\t")


boxdata1 <- boxdata[,-1]
library(tidyr)
library(ggpubr)
boxdata1 <- gather(boxdata1,cellType, value,colnames(boxdata1)[-ncol(boxdata1)])
boxdata1$Cluster <- factor(boxdata1$Cluster,levels = c("Cluster1","Cluster2"))
p2 <- ggboxplot(boxdata1, x = "cellType", y = "value",
               color = "Cluster", palette =  mycol,
               add = "jitter", 
               add.params = list(alpha=0.6),)+ rotate_x_text(45) +xlab("")+ylab("Normalized ssGSEA score")
p2 + stat_compare_means(aes(group = Cluster),method="wilcox.test", label = "p.signif")
ggsave("04immune/GSE56815_Cluster_ssGSEA_score_boxplot.pdf",height=6,width=10)

# The expression of different subtypes of Senescence_related_genes was compared using Wilcoxon rank-sum test
expr <-  read.table("04immune/GSE56815_lowBMD_exprSet.txt",header = T,sep = "\t",check.names = F,row.names = 1)

genes <-  read.table("04immune/Senescence_related_genes.txt",header = F,sep = "\t",check.names = F)

group_res <- read.table("03ConsensusCluster/GSE56815_Cluster_sampleId.txt",header = T,sep = "\t",check.names = F)
group1 <- group_res$id[group_res$Cluster =="Cluster1"]
group2 <- group_res$id[group_res$Cluster =="Cluster2"]
scores <- expr[rownames(expr)%in%genes$V1, ]
scores <-as.matrix(scores[,c(group1,group2)])

library(tibble)
boxdata <- as.data.frame(t( (scores) )) %>% rownames_to_column("id")
boxdata <- merge(boxdata,group_res,by="id")

write.table(boxdata, "04immune/Senescence_related_genes_group_exprSet.txt", quote = F, row.names = F,sep = "\t")


boxdata1 <- boxdata[,-1]

library(tidyr)
library(ggpubr)
mycol <- c("#6b86d7","#c7e2a8")
boxdata1 <- gather(boxdata1,cellType, value,colnames(boxdata1)[-ncol(boxdata1)])
boxdata1$Cluster <- factor(boxdata1$Cluster,levels = c("Cluster1","Cluster2"))
boxdata1$cellType <- factor(boxdata1$cellType ,levels = genes$V1)

p2 <- ggboxplot(boxdata1, x = "cellType", y = "value",
                color = "Cluster", palette = mycol,
                add = "jitter",
                add.params = list(alpha=0.6),)+ rotate_x_text(45) +xlab("")+ylab("Gene expression")
p2 + stat_compare_means(aes(group = Cluster),method="wilcox.test", label = "p.signif")
ggsave("04immune/Senescence_related_genes_boxplot.pdf",height=5,width=11)


# ~~~~~~~~~~~~~~~~~~
# 基因跟免疫细胞的相关性 kruskal.test
# 进行spearman 相关性分析，返回相关性系数和p值
#  "SLC8A3"  "CYP3A7"  "TNFRSF8" "P2RY6"   "PKD1L1"  "HRH1"    "COL5A3"  "AGBL4"   "LRP2"    "ACTG2"   "HKDC1"   "SGO1"   
#  "ASTL"    "IL20RB"
# myoverlap <- c("SLC8A3",  "CYP3A7" , "TNFRSF8" ,"P2RY6",   "PKD1L1" , "HRH1" ,   "COL5A3" , "AGBL4",  "LRP2"   , "ACTG2" ,  "HKDC1",   "SGO1",   
#               "ASTL",    "IL20RB" )

tcga_expr <- expr
tcga_gsva <-as.data.frame(t(gsva))
# genelist <- read.csv("05Lasso/feature_lasso.csv")$x 
genelist <- c("RAP1GAP","GNL2","PSD4","COBLL1","CDV3","NR1H3","G6PC")
immuscore <- function(gene){
  y <- as.numeric(tcga_expr[gene,])
  colnames <- colnames(tcga_gsva)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(tcga_gsva[,x]), y , method="pearson")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
