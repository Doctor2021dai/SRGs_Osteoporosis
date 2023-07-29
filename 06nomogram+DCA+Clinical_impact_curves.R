### nomogram
# Use the lrm function to draw the nomogram of the predicted probability
library(regplot)
library(RColorBrewer)
finalhub <- hub[-6]
logicdata <- rtt[,c(finalhub,'Type')]
logicdata$Type <- ifelse(logicdata$Type=='Normal',0,1)

LR1 <- lrm(Type ~PDPK1+TRIM28+WWP1,
           data = data_set,x=T,y=T)
# window save picture
X <- regplot(LR1,observation=data_set["Type"],
             interval = 'confidence',title='Nomogram',
             clickabel=F,showP = FALSE,points=TRUE,center = TRUE)
dev.off()

# Calibration curve
library(rms)
cal1 <- calibrate(LR1, method="boot")
pdf("07nomogram/calibrate.pdf",width = 5,height = 5)
plot(cal1,lwd = 2,lty = 1,
     bty = "l",
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced (%)",ylab = "Observed (%)",
     col = c("#FB8072"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
mtext("")
dev.off()
#2. DCA Decision Curve
library(rmda)


Nomogram <- decision_curve(Type~PDPK1+TRIM28+WWP1,data=data_set,family = binomial(link = 'logit'),
                       thresholds= seq(0,1, by = 0.01),confidence.intervals =0.95,
                       population.prevalence = 0.3)

pdf('07nomogram/DCA.pdf',height = 5.5,width  = 5.5)
plot_decision_curve(#List,
                    Nomogram,
                    curve.names="Nomogram",
                    cost.benefit.axis = T,
                    cost.benefit.axis = FALSE,confidence.intervals = FALSE,
                    standardize = TRUE,col = brewer.pal(7,'Set1'))
dev.off()

# clinical impact curve
# Use the Nomogram model to predict the risk stratification of 1000 people, display the "cost: benefit" axis, assign 8 scales, and display the confidence interval
pdf('07nomogram/CIC.pdf',height = 5.5,width  = 5)
plot_clinical_impact(Nomogram,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits = 8,col = c('red','blue'),confidence.intervals = T)
dev.off()
