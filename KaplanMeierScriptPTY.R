######Caplan Meier Nonsense

#IL-6 All Caplan-Meiyer 20200815 
if (!require('survival')) install.packages('survival'); library(survival)  #????????????????????????????????????????????????????????????????????????????????????!
if (!require('survBootOutliers')) install.packages('survBootOutliers'); library(survBootOutliers)  #????????????????????????????????????????????????????????????????????????????????????!

#working directory ??????????????????
#setwd("G:/Data/EILSA")

#Variable
#moderate vs nonDSM

#Variable=read.csv("PhenotypeData286_20200815.csv",stringsAsFactors = FALSE, row.names = 1) #phenotype


?Surv

#B
Surv(data1$days,data1$DSM)
plot(survfit(Surv(data1$days,data1$DSM)~data1$IL6BHRwhole),lty=c(1,2), xlim=c(0,5))
legend("topright",c("IL-6 Low","IL-6 High"), lty=c(1,2))

plot(survfit(Surv(data1$days,data1$DSM)~data1$IL6BHRwhole),lty=c(1,2), xlim=c(0,5),col=c("black","red"))
lines(survfit(Surv(data1$days,data1$DSM)~data1$IL6BHRwhole,conf.type="log"),mark.T=F,conf.int=TRUE,col=c("blue","red"))
legend("topright",c("IL-6 Low","IL-6 High"),lty=1,col=c("blue","red"))

###
pdf("IL6_B_1.pdf")
plot(survfit(Surv(data1$days,data1$DSM)~data1$IL6BHRwhole),lty=c(1,2), xlim=c(0,5),col=c("black","red"),ylab="Population without delirium", xlab="Postoperative Days", main="IL-6 arrival at ICU")
lines(survfit(Surv(data1$days,data1$DSM)~data1$IL6BHRwhole,conf.type="log"),mark.T=F,conf.int=TRUE,col=c("blue","red"))
legend("bottomleft",c("IL-6 Low","IL-6 High","Cutoff = 1175 pg/ml","P = 3.0E-04"),lty=c(1,1,0,0),col=c("blue","red","black","black"))
dev.off()

text1 <- paste("Cutoff", x)
x <- 123

survdiff(Surv(data1$days,data1$DSM)~data1$IL6BHRwhole)
#Call:
survtest <- survdiff(formula = Surv(data1$days, data1$DSM) ~ data1$IL6BHRwhole)
survdiff(formula = Surv(data1$days, data1$DSM) ~ data1$IL6CHRwhole)
survtest$chisq

N Observed Expected (O-E)^2/E (O-E)^2/V
data1$IL6BHRwhole=0 210       55     69.2      2.92      13.2
data1$IL6BHRwhole=1  76       36     21.8      9.30      13.2

Chisq= 13.2  on 1 degrees of freedom, p= 3e-04 



######lets try and make a function for this shit 
#need to run a plot function w the stuff
#need the varaibles that go into the legend 
###cutoff
###High / Low
###p value 

pchisq(13.2, 1, lower.tail = F)
###still very much working on this
createKaplan <- function(x, y) {
  
  ROCtable <- createNAtable(x) #create dataset w/out NAs specific to each analyte
  x_NA <-x[!is.na(x)]          #create analyte vector w/out NAs
  
  #create ROC and KSplot variables
  newROC <- rocit(score = x_NA, class = ROCtable$DSM, negref = "0", method = "non")
  newKSplot <- ksplot(newROC, legend = F, xlim = 6000)
  
  #get index of cutoff and therefore sens and spec
  cutoffindex <- which(abs(newROC$Cutoff - newKSplot$`KS Cutoff`) 
                       == min(abs(newROC$Cutoff- newKSplot$`KS Cutoff`)))
  
  #create each variable
  newCutoff <- newKSplot$`KS Cutoff`
  
  ########start of Kaplan Meier bit 
  #get the p value 
  survfit1 <- survdiff(formula = Surv(data1$days, data1$DSM) ~ data1$IL6BHRwhole)
  pvalue <- pchisq(survfit1$chisq, 1, lower.tail = F)
  
  #create and return vector w each variable
  ROCvector <- (c(newCutoff, newAUC, newSens, newSpec)) 
  return(ROCvector)
}


###but lets try getting the pvalues out of it 

getKaplanPvalue <- function(x) {
  survfit1 <- survdiff(formula = Surv(data1$days, data1$DSM) ~ x)
  pvalue <- pchisq(survfit1$chisq, 1, lower.tail = F)
  pvalue
}

getKaplanPvalue(data1$TGFb1_timeCsurvival)
survfit1 <- survdiff(formula = Surv(data1$days, data1$DSM) ~ data1$B2M_Asurvival)
pvalue <- pchisq(survfit1$chisq, 1, lower.tail = F)


colnames(data1)
df = NULL
for (i in 260:273) { 
  rocdf <- getKaplanPvalue(data1[,i])
  df <- rbind(df, rocdf)
}
write.csv(df,"~/PTY/Data\\firstROCs4.csv", row.names = FALSE)
#

getKaplanPvalue(data1$IL.1RA_Asurvival)
data1[,194]


plot(survfit(Surv(data1$days,data1$DSM)~data1$B2M_Bsurvival),lty=c(1,2), xlim=c(0,5),col=c("black","red"),ylab="Population without delirium", xlab="Postoperative Days", main="IL-6 arrival at ICU")
lines(survfit(Surv(data1$days,data1$DSM)~data1$B2M_Bsurvival,conf.type="log"),mark.T=F,conf.int=TRUE,col=c("blue","red"))
legend("bottomleft",c("IL-6 Low","IL-6 High","Cutoff = 1175 pg/ml","P = 3.0E-04"),lty=c(1,1,0,0),col=c("blue","red","black","black"))

getKaplanPvalue(data1$IL.1RA_Asurvival)


#######################OK TIME TO MAKE THESE FIGURES


####IL6 B

survdiff(formula = Surv(data1$days, data1$DSM) ~ data1$IL6BHRwhole)
#getROCcutoff(data1$IL6_timeB)

pdf("IL6_B_1.pdf")
plot(survfit(Surv(data1$days,data1$DSM)~data1$IL6BHRwhole),lty=c(1,2), xlim=c(0,5),col=c("black","red"),ylab="Population without delirium", xlab="Postoperative Days", main="IL-6 at Arrival at ICU")
lines(survfit(Surv(data1$days,data1$DSM)~data1$IL6BHRwhole,conf.type="log"),mark.T=F,conf.int=TRUE,col=c("blue","red"))
legend("bottomleft",c("IL-6 Low","IL-6 High","Cutoff = 897.6 pg/ml","P = 0.0003"),lty=c(1,1,0,0),col=c("blue","red","black","black"))
dev.off()


####IL6 C

survdiff(formula = Surv(data1$days, data1$DSM) ~ data1$IL6_timeCsurvival)
#getROCcutoff(data1$IL6_timeC)

pdf("IL6_C_1.pdf")
plot(survfit(Surv(data1$days,data1$DSM)~data1$IL6_timeCsurvival),lty=c(1,2), xlim=c(0,5),col=c("black","red"),ylab="Population without delirium", xlab="Postoperative Days", main="IL-6 at POD1")
lines(survfit(Surv(data1$days,data1$DSM)~data1$IL6_timeCsurvival,conf.type="log"),mark.T=F,conf.int=TRUE,col=c("blue","red"))
legend("bottomleft",c("IL-6 Low","IL-6 High","Cutoff = 201.9 pg/ml","P = 0.001"),lty=c(1,1,0,0),col=c("blue","red","black","black"))
dev.off()


####IL6 D

survdiff(formula = Surv(data1$days, data1$DSM) ~ data1$IL6_timeDsurvival)
#getROCcutoff(data1$IL6_timeD)

pdf("IL6_D_1.pdf")
plot(survfit(Surv(data1$days,data1$DSM)~data1$IL6_timeDsurvival),lty=c(1,2), xlim=c(0,5),col=c("black","red"),ylab="Population without delirium", xlab="Postoperative Days", main="IL-6 at POD4-7")
lines(survfit(Surv(data1$days,data1$DSM)~data1$IL6_timeDsurvival,conf.type="log"),mark.T=F,conf.int=TRUE,col=c("blue","red"))
legend("bottomleft",c("IL-6 Low","IL-6 High","Cutoff = 29.5 pg/ml","P = 0.001"),lty=c(1,1,0,0),col=c("blue","red","black","black"))
dev.off()


####ST2 C

survdiff(formula = Surv(data1$days, data1$DSM) ~ data1$ST2_Csurvival)
#getROCcutoff(data1$ST2_C)

pdf("ST2_C_1.pdf")
plot(survfit(Surv(data1$days,data1$DSM)~data1$ST2_Csurvival),lty=c(1,2), xlim=c(0,5),col=c("black","red"),ylab="Population without delirium", xlab="Postoperative Days", main="ST2 at POD1")
lines(survfit(Surv(data1$days,data1$DSM)~data1$ST2_Csurvival,conf.type="log"),mark.T=F,conf.int=TRUE,col=c("blue","red"))
legend("bottomleft",c("ST2 Low","ST2 High","Cutoff = 481.8 ng/ml","P = 0.05"),lty=c(1,1,0,0),col=c("blue","red","black","black"))
dev.off()


####ST2 D

survdiff(formula = Surv(data1$days, data1$DSM) ~ data1$ST2_Dsurvival)
#getROCcutoff(data1$ST2_D)

pdf("ST2_D_1.pdf")
plot(survfit(Surv(data1$days,data1$DSM)~data1$ST2_Dsurvival),lty=c(1,2), xlim=c(0,5),col=c("black","red"),ylab="Population without delirium", xlab="Postoperative Days", main="ST2 at POD4-7")
lines(survfit(Surv(data1$days,data1$DSM)~data1$ST2_Dsurvival,conf.type="log"),mark.T=F,conf.int=TRUE,col=c("blue","red"))
legend("bottomleft",c("ST2 Low","ST2 High","Cutoff = 61.5 ng/ml","P = 0.01"),lty=c(1,1,0,0),col=c("blue","red","black","black"))
dev.off()


####IL-1RA A

survdiff(formula = Surv(data1$days, data1$DSM) ~ data1$IL.1RA_Asurvival)
#getROCcutoff(data1$IL.1RA_A)

pdf("IL-1RA_A_1.pdf")
plot(survfit(Surv(data1$days,data1$DSM)~data1$IL.1RA_Asurvival),lty=c(1,2), xlim=c(0,5),col=c("black","red"),ylab="Population without delirium", xlab="Postoperative Days", main="IL-1RA at Preoperative Baseline")
lines(survfit(Surv(data1$days,data1$DSM)~data1$IL.1RA_Asurvival,conf.type="log"),mark.T=F,conf.int=TRUE,col=c("blue","red"))
legend("bottomleft",c("ST2 Low","ST2 High","Cutoff = 79.3 pg/ml","P = 0.04"),lty=c(1,1,0,0),col=c("blue","red","black","black"))
dev.off()


####IL-10 C

survdiff(formula = Surv(data1$days, data1$DSM) ~ data1$IL.10_Csurvival)
#getROCcutoff(data1$IL.10_C)

pdf("IL-10_C_1.pdf")
plot(survfit(Surv(data1$days,data1$DSM)~data1$IL.10_Csurvival),lty=c(1,2), xlim=c(0,5),col=c("black","red"),ylab="Population without delirium", xlab="Postoperative Days", main="IL-10 at POD1")
lines(survfit(Surv(data1$days,data1$DSM)~data1$IL.10_Csurvival,conf.type="log"),mark.T=F,conf.int=TRUE,col=c("blue","red"))
legend("bottomleft",c("ST2 Low","ST2 High","Cutoff = 1.24 pg/ml","P = 0.008"),lty=c(1,1,0,0),col=c("blue","red","black","black"))
dev.off()


####IL-17F D

survdiff(formula = Surv(data1$days, data1$DSM) ~ data1$IL.17F_Dsurvival)
#getROCcutoff(data1$IL.17F_D)

pdf("IL-17F_D_1.pdf")
plot(survfit(Surv(data1$days,data1$DSM)~data1$IL.17F_Dsurvival),lty=c(1,2), xlim=c(0,5),col=c("black","red"),ylab="Population without delirium", xlab="Postoperative Days", main="IL-17F at POD4-7")
lines(survfit(Surv(data1$days,data1$DSM)~data1$IL.17F_Dsurvival,conf.type="log"),mark.T=F,conf.int=TRUE,col=c("blue","red"))
legend("bottomleft",c("ST2 Low","ST2 High","Cutoff = 20.0 pg/ml","P = 0.03"),lty=c(1,1,0,0),col=c("blue","red","black","black"))
dev.off()










####Creating a 2x2 block of Kaplan Figure for IL6
pdf("All4_IL6_1.pdf")

par(mfrow=c(2,2))

plot(survfit(Surv(data1$days,data1$DSM)~data1$IL6_timeAsurvival),lty=c(1,2), xlim=c(0,5),col=c("black","red"),ylab="Population without delirium", xlab="Postoperative Days", main="IL-6 at Preoperative Baseline")
lines(survfit(Surv(data1$days,data1$DSM)~data1$IL6_timeAsurvival,conf.type="log"),mark.T=F,conf.int=TRUE,col=c("blue","red"))
legend("bottomleft",c("IL-6 Low","IL-6 High","Cutoff = 7.06 pg/ml","P = 0.1"),lty=c(1,1,0,0),col=c("blue","red","black","black"))

####IL6B
plot(survfit(Surv(data1$days,data1$DSM)~data1$IL6BHRwhole),lty=c(1,2), xlim=c(0,5),col=c("black","red"),ylab="Population without delirium", xlab="Postoperative Days", main="IL-6 at Arrival at ICU")
lines(survfit(Surv(data1$days,data1$DSM)~data1$IL6BHRwhole,conf.type="log"),mark.T=F,conf.int=TRUE,col=c("blue","red"))
legend("bottomleft",c("IL-6 Low","IL-6 High","Cutoff = 897.6 pg/ml","P = 0.0003"),lty=c(1,1,0,0),col=c("blue","red","black","black"))


####IL6 C


plot(survfit(Surv(data1$days,data1$DSM)~data1$IL6_timeCsurvival),lty=c(1,2), xlim=c(0,5),col=c("black","red"),ylab="Population without delirium", xlab="Postoperative Days", main="IL-6 at POD1")
lines(survfit(Surv(data1$days,data1$DSM)~data1$IL6_timeCsurvival,conf.type="log"),mark.T=F,conf.int=TRUE,col=c("blue","red"))
legend("bottomleft",c("IL-6 Low","IL-6 High","Cutoff = 201.9 pg/ml","P = 0.001"),lty=c(1,1,0,0),col=c("blue","red","black","black"))


####IL6 D

plot(survfit(Surv(data1$days,data1$DSM)~data1$IL6_timeDsurvival),lty=c(1,2), xlim=c(0,5),col=c("black","red"),ylab="Population without delirium", xlab="Postoperative Days", main="IL-6 at POD4-7")
lines(survfit(Surv(data1$days,data1$DSM)~data1$IL6_timeDsurvival,conf.type="log"),mark.T=F,conf.int=TRUE,col=c("blue","red"))
legend("bottomleft",c("IL-6 Low","IL-6 High","Cutoff = 29.5 pg/ml","P = 0.001"),lty=c(1,1,0,0),col=c("blue","red","black","black"))

dev.off()



