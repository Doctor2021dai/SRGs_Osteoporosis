library(tidyverse)
source('06randomForest/msvmRFE.R')   #文件夹内自带
library(VennDiagram)
library(e1071)
library(caret)
library(randomForest)

library(GEOquery)

Sys.setenv("VROOM_CONNECTION_SIZE"=999999)

# GSE56815 40 high (20 pre- and 20 postmanopausal) and 40 low hip BMD (20 pre- and 20 postmanopausal) subjects
GSE_id <- "GSE56815"
# download data
gset = getGEO(GSE_id, destdir="01rawData/",getGPL = F)
# Get the ExpressionSet object, including the expression matrix and grouping information
gset=gset[[1]]

# Get group information, click to view information
pdata=pData(gset)

pd1=pdata[,apply(pdata, 2, function(x){
  length(unique(x))>=2})] #narrow down
dim(pd1)

Normal=rownames(pd1[grepl("high BMD",as.character(pd1$characteristics_ch1.1)),])
Tumor=rownames(pd1[grepl("low BMD",as.character(pd1$characteristics_ch1.1)),])

                  
rt=read.table("01rawData/GSE56815_symbol_exprSet.txt",sep="\t",header=T,check.names=F,row.names = 1)
tpm <- read.table("02DEGs/venn_intersectGenes.txt",sep="\t",header=F,check.names=F)
data_set <-  as.data.frame(t(rt[tpm$V1,]))

data_set$Type  <- ifelse(rownames(data_set) %in% Normal, 0,1)

data_set$Type <- factor(data_set$Type)
data_set[1:4,1:4]
table(data_set$Type)

set.seed(321)
index <- createDataPartition(data_set$Type, p = 0.7, list = FALSE)
train <- data_set[index, ]
test  <- data_set[-index, ]
table(train$Type)
table(test$Type)

library("DALEX")
# 5-fold cross-validation
control = trainControl(method="repeatedcv", number = 10, savePredictions=TRUE)
# Random Forest
mod_rf = train(Type ~ .,
               data = train, method='rf', trControl = control)

# Support Vector Machines
mod_svm <- train(Type ~.,
                 data = train, method = "svmRadial", prob.model = TRUE, trControl=control)

# Create a custom predict function
p_fun <- function(object, newdata){
  predict(object, newdata=newdata, type="prob")[,2]
}

# Convert the outcome variable to a numeric binary vector
yTest <- as.numeric(as.character(test$Type))
# Create explainer objects for each machine learning model
explainer_rf  <- explain(mod_rf, label = "RF",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
explainer_svm <- explain(mod_svm, label = "SVM",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
# Model Performance
# Calculate model performance and residuals
mp_rf  <- model_performance(explainer_rf)
mp_svm <- model_performance(explainer_svm)
p1 <- plot(mp_rf, mp_svm)
p1
p2 <- plot(mp_rf, mp_svm, geom = "boxplot")
p2
# Activate the patchwork package to combine plots
library("patchwork")
pdf("06randomForest/distribution_of_residuals.pdf",width = 9,height = 5)
p1 + p2
dev.off()

# There are two very important parameters in the function randomForest() of the random forest algorithm, and these two parameters will affect the model
# accuracy, they are mtry and ntree respectively. Generally, the choice of mtry is to try one by one until an ideal value is found.
# The selection of ntree can roughly judge the value when the error in the model is stable through graphics.
# Dependent Variable Independent Variable Construction Formula
form_cls <- as.formula(
  paste0(
    "Type ~ ",paste(colnames(train)[1:ncol(train)-1], collapse = " + ")
  )
)
form_cls

n <- ncol(train)-1
errRate <- NULL
set.seed(321)
for (i in 1:n) {
  mtry_fit <- randomForest(form_cls , data = train,mtry=i,importance = TRUE)
  err <- mean(mtry_fit$err.rate)
  errRate[i] <- err
  print(err)
}
# Choose the mtry with the smallest average error
which.min(errRate)
# m is 8
# Then select the ntree value, ntree specifies the number of decision trees contained in the random forest, the default is 500;
# When it is around 400, the error in the model is basically stable, so take ntree=400.
set.seed(321)
ntree_fit <- randomForest(form_cls, data = train,mtry=12,ntree=1000,importance = TRUE)
pdf("06randomForest/randomForest_ntree.pdf",width = 5,height = 5)
plot(ntree_fit)
legend("top",legend = colnames(ntree_fit$err.rate),
       lty=1:3,
       col = 1:3,
       horiz = T)
dev.off()
# abline( v = 500, col = "gray60")
set.seed(321)
forest <- randomForest(form_cls, data = train, mtry = 12,
                         importance = TRUE,proximity=TRUE,ntree = 800)
forest
pdf("06randomForest/randomForest_ntree_500.pdf",width = 5,height = 5)
plot(forest)
dev.off()

# Partial correlation plot
partialPlot(x=forest,
            pred.data = train,
            x.var = PDPK1,
            which.class = "1",
            ylab = "1")
#Training set self-test
forest_predict_train <- predict(forest, newdata = train)

obs_p_forest_train = data.frame(Predicted=forest_predict_train,Actual=train$Type)
forest_train <- table(train$Type, forest_predict_train, dnn = c('Actual', 'Predicted'))
forest_train
# test set test
forest_predict_test <- predict(forest, test)
obs_p_forest_test = data.frame(Predicted=forest_predict_test,Actual=test$Type)
forest_test <- table(test$Type, forest_predict_test, dnn = c('Actual', 'Predicted'))
forest_test

#10 repetitions of ten-fold cross-validation
set.seed(321)
train.cv <- replicate(10,rfcv(train[-ncol(train)], train$Type, cv.fold = 10,step = 1.5), simplify = FALSE)
train.cv
#Extract validation result drawing
train.cv <- data.frame(sapply(train.cv, '[[', 'error.cv'))
train.cv$genes <- rownames(train.cv)
train.cv <- reshape2::melt(train.cv, id = 'genes')
train.cv$genes <- as.numeric(as.character(train.cv$genes))
train.cv.mean <- aggregate(train.cv$value, by = list(train.cv$genes), FUN = mean)
head(train.cv.mean, 10)
#拟合线图
#拟合线图
library(ggplot2)

p <- ggplot(train.cv.mean, aes(Group.1, x)) +
  geom_line(colour="#6b86d7") +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of genes', y = 'Cross-validation error')
p

#大约提取前 30 个重要的 OTUs
p + geom_vline(xintercept = 10, linetype="dotted")+annotate('text',x = 9,y=0.2,label='n = 10',color='black')
ggsave("06randomForest/Cross-validation error.pdf",width = 4,height = 4)


# View the score indicating the importance of each variable
summary(forest)
importance_otu <- forest$importance
head(importance_otu)

#Or use the function importance()
importance_otu <- data.frame(importance(forest))
head(importance_otu)
#It can be sorted according to the level of importance, for example, according to the "Mean Decrease Accuracy" indicator
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(importance_otu)
# output table
write.table(importance_otu, '06randomForest/RF_importance_gene.txt', sep = '\t', col.names = NA, quote = T)
#Plot display top10 important variables
pdf("06randomForest/randomForest_Top10_importance.pdf",width = 7,height = 6)
varImpPlot(forest, n.var = min(10, nrow(forest$importance)), main = 'Top 10 - genes importance')
dev.off()

importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(importance_otu)

# SVM-RFE ~~~~~~~
#Using 10-fold cross-validation 
train[1:4,1:4]
input <- train[-51]
input$group   <- ifelse(rownames(input) %in% Normal, "high","low")

input$group <- factor(input$group)
input <- input[c(ncol(input),1:ncol(input)-1)]
#Split the data and assign random numbers
svmRFE(input, k = 10, halve.above = 20)

nfold = 10
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
#feature selection
results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=20) 
# View the main variables
top.features = WriteFeatures(results, input, save=F) 
head(top.features)
#Save the features found by SVM-REF to a file
write.table(top.features, '06randomForest/feature_svm_gene.txt', sep = '\t', col.names = NA, quote = T)

# Enable parallel as the backend for foreach parallel computing
library(doParallel)
detectCores()
cl <- makeCluster(10)
registerDoParallel(cl)
featsweep <- foreach(x=1:50, .packages = c("e1071")) %dopar% FeatSweep.wrap(x,results, input)
# end parallel
stopCluster(cl)
save(featsweep,file = "06randomForest/featsweep1020.RData")
load(file = "06randomForest/featsweep1020.RData")
# drawing
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#dev.new(width=4, height=4, bg='white')
pdf("06randomForest/B_svm-error.pdf",width = 5,height = 5)
PlotErrors(errors, no.info=no.info) # View error rate

dev.off()

#dev.new(width=4, height=4, bg='white')
pdf("06randomForest/B_svm-accuracy.pdf",width = 5,height = 5)
Plotaccuracy(1-errors,no.info=no.info) # View accuracy
dev.off()

# The position of the red circle in the figure is the lowest point of error rate
which.min(errors) 
top.features <- top.features[1:which.min(errors), "FeatureName"]
write.table(top.features, '06randomForest/svm_errorlow_gene.txt', sep = '\t', col.names = NA, quote = T)

#ggvenn package to draw Venn diagrams
LYMPHANGIOGENESIS <- read.table("06randomForest/RF_importance_gene.txt",header = T,sep="\t",check.names = F)

x <- list(RF=LYMPHANGIOGENESIS[1:10,1],SVM=top.features)
library(ggvenn)
mycol <- c("#009E73","#F0E442","#0072B2","#D55E00")
ggvenn(x,c("RF","SVM"),
       stroke_size = 0.3,show_percentage = FALSE,
       fill_color =mycol)
ggsave("06randomForest/venn.pdf",width = 4.2,height = 4.5)

intersectGenes<- intersect(top.features,LYMPHANGIOGENESIS[1:10,1])
write.table(intersectGenes,file = "06randomForest/venn_intersectGenes.txt",row.names = F,quote = F,sep = "\t",col.names = F)

svm_model <- svm(Type ~TRIM28+PDPK1+WWP1,data=train)
svm_model
# Check the forecast of the sample
svm_predict_train<- predict(svm_model,train)
obs_p_svm_train = data.frame(Predicted=svm_predict_train,Actual=train$Type)
svm_train <- table(train$Type,svm_predict_train,dnn = c('Actual', 'Predicted'))
svm_train
# test-level prediction output confusion matrix
svm_predict_test<- predict(svm_model,test)
obs_p_svm_test = data.frame(Predicted=svm_predict_test,Actual=test$Type)
svm_test <- table(test$Type,svm_predict_test,dnn = c('Actual', 'Predicted'))
svm_test


set.seed(321)
forest <- randomForest(Type ~ TRIM28+PDPK1+WWP1, data = train,
                       importance = TRUE,proximity=TRUE)
forest

# Training set self-test
forest_predict_train <- predict(forest, train)
obs_p_forest_train = data.frame(Predicted=forest_predict_train,Actual=train$Type)
forest_train <- table(train$Type, forest_predict_train, dnn = c('Actual', 'Predicted'))
forest_train
# test set test
forest_predict_test <- predict(forest, test)
obs_p_forest_test = data.frame(Predicted=forest_predict_test,Actual=test$Type)
forest_test <- table(test$Type, forest_predict_test, dnn = c('Actual', 'Predicted'))
forest_test
# Draw the ROC curve
library(pROC) 
rf_roc <- roc(train$Type,as.numeric(forest_predict_train))
svm_roc <- roc(train$Type,as.numeric(svm_predict_train))
pdf("06randomForest/RF_SVM_auc_train.pdf",width = 5,height = 5)
plot(rf_roc,  # The previously constructed ROC object
     print.auc=TRUE, # Output the AUC value on the image
     col="#377eb8",
     main="trainSet",
     print.auc.x=0.5, print.auc.y=0.5, # AUC value coordinates are (x, y)
     grid= FALSE, # Don't show grid background lines
     legacy.axes=TRUE)  

plot.roc(svm_roc, add=TRUE,  
         col = "#4daf4a", 
         print.auc=TRUE, 
         print.auc.col = "#4daf4a", 
         print.auc.x=0.5, print.auc.y=0.4)  
# Set the legend label
legend("bottomright", 
       legend=c("RF", "SVM"), 
       col=c("#377eb8", "#4daf4a"), 
       lwd=4)
dev.off()
#Draw the ROC curve
library(pROC) 
rf_roc <- roc(test$Type,as.numeric(forest_predict_test))
svm_roc <- roc(test$Type,as.numeric(svm_predict_test))
pdf("06randomForest/RF_SVM_auc_test.pdf",width = 5,height = 5)
plot(rf_roc,  # The previously constructed ROC object
     print.auc=TRUE, # Output the AUC value on the image
     col="#377eb8",
     main="testSet",
     print.auc.x=0.5, print.auc.y=0.5, # AUC value coordinates are (x, y)
     grid= FALSE, 
     legacy.axes=TRUE)  

plot.roc(svm_roc, add=TRUE, # add curve
         col = "#4daf4a", # set curve color
         print.auc=TRUE, # Output AUC on the image
         print.auc.col = "#4daf4a", # Set the color of the AUC text
         print.auc.x=0.5, print.auc.y=0.4)  # The coordinates of AUC are (x, y)

# Set the legend label
legend("bottomright", 
       legend=c("RF", "SVM"), 
       col=c("#377eb8", "#4daf4a"), 
       lwd=4)
dev.off()

#train
roc_train1 <- roc(train$Type,train$PDPK1)
roc_train2 <- roc(train$Type,train$TRIM28)
roc_train3 <- roc(train$Type,train$WWP1)

pdf("06randomForest/trainSet_RF_SVM_auc.pdf",width = 5,height = 5)
plot(roc_train1,  
     print.auc=TRUE, 
     col="#377eb8",
     main="trainSet",
     print.auc.x=0.5, print.auc.y=0.5, 
     grid= FALSE, 
     legacy.axes=TRUE)  

plot.roc(roc_train2, add=TRUE,  
         col = "#4daf4a", 
         print.auc=TRUE, 
         print.auc.col = "#4daf4a", 
         print.auc.x=0.5, print.auc.y=0.4)  
plot.roc(roc_train3, add=TRUE,  
         col = "#ffb549", 
         print.auc=TRUE, 
         print.auc.col = "#ffb549", 
         print.auc.x=0.5, print.auc.y=0.3) 
# Set the legend label
legend("bottomright",
       legend=c("PDPK1", "TRIM28","WWP1"),
       col=c("#377eb8", "#4daf4a","#ffb549"),
       lwd=4)
dev.off()


#test
roc_test1 <- roc(test$Type,test$PDPK1)
roc_test2 <- roc(test$Type,test$TRIM28)
roc_test3 <- roc(test$Type,test$WWP1)

pdf("06randomForest/testSet_RF_SVM_auc.pdf",width = 5,height = 5)
plot(roc_test1,  
     print.auc=TRUE, 
     col="#377eb8",
     main="testSet",
     print.auc.x=0.5, print.auc.y=0.5, 
      grid= FALSE, 
     legacy.axes=TRUE)  

plot.roc(roc_test2, add=TRUE,  
         col = "#4daf4a", 
         print.auc=TRUE, 
         print.auc.col = "#4daf4a", 
         print.auc.x=0.5, print.auc.y=0.4)  
plot.roc(roc_test3, add=TRUE,  
         col = "#ffb549", 
         print.auc=TRUE, 
         print.auc.col = "#ffb549", 
         print.auc.x=0.5, print.auc.y=0.3)  
# Set the legend label
legend("bottomright",
       legend=c("PDPK1", "TRIM28","WWP1"),
       col=c("#377eb8", "#4daf4a","#ffb549"),
       lwd=4)
dev.off()

# TRIM28+PDPK1+WWP1
validset <- read.table("09GEO_valid/GSE62402_symbol_featuregenes_exprSet.txt",header = T,sep = "\t")
library(pROC)

roc_valid1 <- roc(validset$Type,validset$PDPK1)
roc_valid2 <- roc(validset$Type,validset$TRIM28)
roc_valid3 <- roc(validset$Type,validset$WWP1)

pdf("06randomForest/GSE62402_RF_SVM_auc_valid123.pdf",width = 5,height = 5)
plot(roc_valid1,  
     print.auc=TRUE, 
     col="#377eb8",
     main="validSet",
     print.auc.x=0.5, print.auc.y=0.5, 
     grid= FALSE, 
     legacy.axes=TRUE)  

plot.roc(roc_valid2, add=TRUE, 
         col = "#4daf4a", 
         print.auc=TRUE, 
         print.auc.col = "#4daf4a", 
         print.auc.x=0.5, print.auc.y=0.4)  
plot.roc(roc_valid3, add=TRUE,  
         col = "#ffb549", 
         print.auc=TRUE, 
         print.auc.col = "#ffb549", 
         print.auc.x=0.5, print.auc.y=0.3)  
# Set the legend label
legend("bottomright",
       legend=c("PDPK1", "TRIM28","WWP1"),
       col=c("#377eb8", "#4daf4a","#ffb549"),
       lwd=4)
dev.off()
#Validation set self-test
forest_predict_valid <- predict(forest, validset,type = "class")
obs_p_forest_valid  = data.frame(Predicted=forest_predict_valid ,Actual=validset$Type)
forest_valid <- table(validset$Type, forest_predict_valid, dnn = c('Actual', 'Predicted'))
forest_valid


svm_predict_valid<- predict(svm_model,validset)
obs_p_svm_valid = data.frame(Predicted=svm_predict_valid,Actual=validset$Type)
svm_valid <- table(validset$Type,svm_predict_valid,dnn = c('Actual', 'Predicted'))
svm_valid

# Draw random forest importance variable lollipop map
imp <- read.table("06randomForest/RF_importance_gene.txt",header = T,sep = "\t",row.names = 1)
imp_sub=imp[1:10,]
imp_sub$Genes<-rownames(imp_sub)
# add factors to keep order
imp_sub$Genes=factor(imp_sub$Genes,order=T,levels = rev(imp_sub$Genes))
p=ggplot(data = imp_sub, mapping = aes(x=Genes,y=MeanDecreaseAccuracy)) + 
  geom_segment(aes(x = Genes, y = 0, xend = Genes, yend = MeanDecreaseAccuracy),color="grey",size=1.5) +
  geom_point(aes(size = MeanDecreaseAccuracy),color = "#009E73")+
  coord_flip()+theme_bw()+
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(colour="black"),
        axis.title.y = element_blank(), 
        axis.text.y=element_text(colour="black"),
        panel.border=element_rect(colour = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))
p
ggsave("06randomForest/varImpPlot_genes.pdf",width = 3.8,height =4.5 )
