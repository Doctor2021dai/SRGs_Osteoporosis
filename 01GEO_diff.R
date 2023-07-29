# load the R package
library(GEOquery)
library(limma)
library(tidyverse)
Sys.setenv("VROOM_CONNECTION_SIZE"=9999999)

# GSE56815 40 high (20 pre- and 20 postmanopausal) and 40 low hip BMD (20 pre- and 20 postmanopausal) subjects
GSE_id <- "GSE56815"
# download data
gset = getGEO(GSE_id, destdir="01rawData/",getGPL = F)
# Get the ExpressionSet object, including the expression matrix and grouping information
gset=gset[[1]]

# Get group information
pdata=pData(gset)

pd1=pdata[,apply(pdata, 2, function(x){
    length(unique(x))>=2})] #narrow down
dim(pd1)
write.table( pd1, paste0("01rawData/" ,GSE_id,"_clinical.txt"),sep = "\t",quote = F,row.names = F)
# high BMD low BMD
Normal=rownames(pd1[grepl("high BMD",as.character(pd1$characteristics_ch1.1)),])
Tumor=rownames(pd1[grepl("low BMD",as.character(pd1$characteristics_ch1.1)),])

group=c(rep('high BMD',length(Normal)),
        rep('low BMD',length(Tumor)))  
group=factor(group)

table(group)
exprSet=exprs(gset)
# Subset the expression matrix
exprSet=exprSet[,c(Normal,Tumor)]

boxplot(exprSet,outline=FALSE, notch=T, las=2)

# Need to correct it, the method used is similar to Quntile Normalization
library(limma) 

exprSet=normalizeBetweenArrays(exprSet)

# Determine whether data conversion is required
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}

# Process GEO data gpl gene name annotation file
gpl <- read.table('01rawData/GPL96-57554.txt',header = TRUE,sep = '\t',quote = '', comment.char = '#',check.names = F,fill = T)
colnames(gpl)

library(dplyr)
probe2id <- gpl %>%
    dplyr::select('ID','Gene Symbol') %>%
    tidyr::separate_rows('Gene Symbol',sep = ' /// ') %>%
    dplyr::rename(probe_id = 'ID',symbol='Gene Symbol') %>%
    dplyr::filter(symbol != '')

exprSet1 <- exprSet   %>% as.data.frame() %>%
    tibble::rownames_to_column(var="probe_id") %>% 
    as_tibble() %>% 
    # Merge probe information
    dplyr::inner_join(probe2id,by="probe_id") %>% 
    # Remove redundant information
    dplyr::select(-probe_id) %>% 
    # rearrange
    dplyr::select(symbol,everything()) %>% 
    # Remove NA in symbol
    dplyr::filter(symbol != "NA") %>% 
    # Group by symbol to get the average
    dplyr::group_by(symbol) %>%
    dplyr::summarise_all(mean) %>%
    # column names become row names
    tibble::column_to_rownames(var = "symbol")


exprSet1[1:3,1:5]
exprSet_out <- exprSet1 %>%  tibble::rownames_to_column(var = "symbol")
write.table(exprSet_out,file = paste0("01rawData/" ,GSE_id,"_symbol_exprSet.txt"),row.names = F,quote = F,sep = "\t")
m6A <- read.table("01rawData/CellAge Senescence Genes.csv",header = T,sep=",")

# Differentially expressed gene analysis
design=model.matrix(~group)
fit=lmFit(exprSet1,design)
fit=eBayes(fit) 
res=topTable(fit,adjust='fdr',coef="grouplow BMD",number=Inf)
allDiff <- na.omit(res)

logFCCutoff <- 0.2  
pvalueCutoff <- 0.05

outDiff=res[(abs(allDiff$logFC)>=logFCCutoff & allDiff$P.Value<pvalueCutoff),]
outDiff <- outDiff %>% rownames_to_column(var = "id")
write.table(outDiff,file=paste0("02DEGs/" ,GSE_id,"_limma_diff.txt"),row.names=F,quote=F,sep = "\t")
allDiff <- allDiff  %>% rownames_to_column(var = "id")
write.table(allDiff,file=paste0("02DEGs/" ,GSE_id,"_limma_all.txt"),row.names=F,quote=F,sep = "\t")

# Volcano map of differentially expressed genes
library(ggplot2)

Significant=ifelse((res$P.Value< 0.05 & abs(res$logFC)>= 0), ifelse(res$logFC > 0,"Up","Down"), "Not")
# draw volcano map
p = ggplot(res, aes(logFC, -log10(P.Value)))+
    geom_point(aes(col=Significant))+
    scale_color_manual(values=c("#004BFB", "#BFBFBF", "#F91F10"))+
    labs(title = " ")+
    #geom_vline(xintercept=c(-1,1), colour="black", linetype="dashed")+
    geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
    xlab(expression("log"["2"]*"FC"))+
    ylab(expression("-log"["10"]*"P Value"))+
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()


#Save as picture
pdf("02DEGs/GSE56815_vol.pdf",width=5.5,height=5)
print(p)
dev.off()



# Use the ggvenn package to draw a Venn diagram
SRGs <- read.table("01rawData/CellAge Senescence Genes.csv",header = T,sep=",",check.names = F)

x <- list(DEGs=outDiff$id,SSGs=SRGs$`Gene Symbol`)
library(ggvenn)
mycol <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#ADD1E5",
           '#CCEBC5', '#FCCDE5', '#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")
ggvenn(x,c("DEGs","SSGs"),
       stroke_size = 0.3,show_percentage = FALSE,
       fill_color =mycol[4:8])
ggsave("02DEGs/venn.pdf",width = 4.2,height = 4.5)

intersectGenes<- intersect(outDiff$id,SRGs$`Gene Symbol`)
write.table(intersectGenes,file = "02DEGs/venn_intersectGenes.txt",row.names = F,quote = F,sep = "\t",col.names = F)


# Draw a heatmap of differentially expressed genes
hmGene <- read.table("02DEGs/venn_intersectGenes.txt",header = F)

library(pheatmap)

hmExp= exprSet1[hmGene$V1,]

max(hmExp)
min(hmExp)

Type= factor(c(rep("high BMD",40),rep("low BMD",40)),levels = c("high BMD","low BMD"))
names(Type)=colnames(hmExp)
Type=as.data.frame(Type)

loc <- order(Type,colSums(hmExp),decreasing = F)
pdf(file="02DEGs/SSGs_DEGs_heatmap.pdf",height=8,width=9)
pheatmap( hmExp[,loc],
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

hmExp_out <- hmExp %>% rownames_to_column("id")
write.table(hmExp_out,file=paste0("02DEGs/" ,GSE_id,"_DEGs_SSGs_exprSet.txt"),row.names=F,quote=F,sep = "\t")
