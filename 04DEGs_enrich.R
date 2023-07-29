rm(list = ls())
library("stringr")
library("clusterProfiler")
library("org.Hs.eg.db")
library("GOSemSim")
library("enrichplot")
library("ggplot2")

expr <-  read.table("04immune/GSE56815_lowBMD_exprSet.txt",header = T,sep = "\t",check.names = F,row.names = 1)

group_res <- read.table("03ConsensusCluster/GSE56815_Cluster_sampleId.txt",header = T,sep = "\t",check.names = F)
Cluster1 <- group_res$id[group_res$Cluster =="Cluster1"]
Cluster2 <- group_res$id[group_res$Cluster =="Cluster2"]

group=c(rep('Cluster1',length(Cluster1)),
        rep('Cluster2',length(Cluster2)))  
group=factor(group)

table(group)
# Subset the expression matrix
exprSet=expr[,c(Cluster1,Cluster2)]
# Differentially expressed gene analysis
design=model.matrix(~group)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
res=topTable(fit,adjust='fdr',coef="groupCluster2",number=Inf)
allDiff <- na.omit(res)

logFCCutoff <- 0.5
pvalueCutoff <- 0.05

outDiff=res[(abs(allDiff$logFC)>logFCCutoff & allDiff$P.Value<pvalueCutoff),]
outDiff <- outDiff %>% rownames_to_column(var = "id")
write.table(outDiff,file=paste0("05subtype/" ,GSE_id,"_limma_diff.txt"),row.names=F,quote=F,sep = "\t")
allDiff <- allDiff  %>% rownames_to_column(var = "id")
write.table(allDiff,file=paste0("05subtype/" ,GSE_id,"_limma_all.txt"),row.names=F,quote=F,sep = "\t")

# Volcano map of differentially expressed genes
library(ggplot2)

Significant=ifelse((res$P.Value< 0.05 & abs(res$logFC)> 0.5), ifelse(res$logFC > 0.5,"Up","Down"), "Not")

p = ggplot(res, aes(logFC, -log10(P.Value)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#004BFB", "#BFBFBF", "#F91F10"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-0.5,0.5), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
  
  xlab(expression("log"["2"]*"FC"))+
  ylab(expression("-log"["10"]*"P Value"))+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()

#Save as picture
pdf("05subtype/GSE56815_subtype_vol.pdf",width=5.5,height=5)
print(p)
dev.off()

library(pheatmap)

hmExp= exprSet[outDiff$id,]

max(hmExp)
min(hmExp)

Type= factor(c(rep("Cluster1",25),rep("Cluster2",15)),levels = c("Cluster1","Cluster2"))
names(Type)=colnames(hmExp)
Type=as.data.frame(Type)

loc <- order(Type,colSums(hmExp),decreasing = F)
pdf(file="05subtype/subtype_DEGs_heatmap.pdf",height=8,width=9)
pheatmap::pheatmap(hmExp,
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8,
         border=FALSE)

dev.off()

# Enrichment analysis of differentially expressed genes
rt=read.table("05subtype/GSE56815_subtype_limma_diff.txt",sep="\t",header=T,check.names=F,row.names = 1)

entrezIDs <- mget(rownames(rt), org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out = cbind(symbol=rownames(rt),entrezID=entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
out = out[,c(1,7)]
out = cbind(symbol=rownames(out),out)
write.table(out,file="05subtype/id.xls",quote=F,row.names = F,sep="\t")

rt <- read.table("05subtype/id.xls",sep="\t",header = T)
entrezID_gene <- na.omit(rt$entrezID)

x <- enrichGO(gene = entrezID_gene,
               OrgDb = org.Hs.eg.db,
               pvalueCutoff =0.05,
               qvalueCutoff = 0.05,
               ont="all",
               readable =T
)
write.table(x,file="05subtype/GO.xls",quote=F,row.names = F,sep = "\t")


pdf(file="05subtype/GO_barplot.pdf",width = 9.5,height = 9)
barplot(x, drop = TRUE, showCategory =10,split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free') +scale_y_discrete(labels=function(x) str_wrap(x, width=60)) 
dev.off()


pdf(file="05subtype/GO_bubble.pdf",width = 9.5,height = 9)
dotplot(x,showCategory = 10,split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')+scale_y_discrete(labels=function(x) str_wrap(x, width=60)) 
dev.off()

# install.packages('R.utils')
R.utils::setOption( "clusterProfiler.download.method",'auto' )
kk <- enrichKEGG(gene = entrezID_gene, organism = "hsa", 
                 pvalueCutoff =0.05, 
                 qvalueCutoff =0.05)

y <- setReadable(kk,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
write.table(y,file="05subtype/KEGG.xls",quote=F,row.names = F,sep = "\t")

pdf(file="05subtype/KEGG_barplot.pdf",width = 7,height = 3.5)
barplot(kk, drop = TRUE, showCategory = 10)+scale_y_discrete(labels=function(x) str_wrap(x, width=60)) 
dev.off()

pdf(file="05subtype/KEGG_bubble.pdf",width = 7,height = 3.5)+scale_y_discrete(labels=function(x) str_wrap(x, width=60))
dotplot(kk, showCategory = 10)+scale_y_discrete(labels=function(x) str_wrap(x, width=60)) 

dev.off()
