## new ROC function to do the following:
#1- ROC plot using ROCit and save in folder
#2- call AUC, Sensitivity, and Specificty from cutpointR,
#3- save data in new empty df
#4- return df

#game plan
#create fucntion first that creates roc plot
#create function that uses cutpointR to save data to 3 variables 
#create for loop that can call both functions for i in colnames
#have that for loop save to a dataframe somehow

######ROC Function
#create new datatable w/out NAs
#plot the ROC curve of it 




createNAtable <- function(x) {

  analyte_NA <- data1 %>% filter(!is.na(x))
  
}
createROC_single <- function(x) {
  
  ROCtable <- createNAtable(x) #create dataset w/out NAs
  x_NA <-x[!is.na(x)]          #create analyte vector w/out NAs
  
  #create ROC variable
  newROC <- ?rocit(score = x_NA, class = ROCtable$DSM, negref = "0", method = "non")
  
  #create pdf
  pdf(file = "~/PTY/Data/ROC Curves/ROC%03d.pdf", onefile = F, width = 7, height = 7)
                           
  #create ROC plot
  plot(newROC, legend = F)

  dev.off()
  }
createROC_multiple <- function(x) {
  
  ROCtable <- createNAtable(x) #create dataset w/out NAs
  x_NA <-x[!is.na(x)]          #create analyte vector w/out NAs
  
  #create ROC variable
  newROC <- rocit(score = x_NA, class = ROCtable$DSM, negref = "0", method = "non")
  
  #create ROC plot
  plot(newROC, legend = F)
  
}
createROC_multiple(data1$age)
#createCutpoint <- function(x) {
  #I want this to take a vector argument data1$IL6_timeA
  #I want it to calculate cutoff, AUC, Sens, Spec from cutpointR
  #return these values in a big spreadsheet with the others  
  
  direction = "<="
  #if (median()) {
  #  statement
 # }
  cutpointdata <- cutpointr(data1, x, DSM, na.rm = T,
                            method = oc_youden_normal, use_midpoints = T)
  return(c(cutpointdata$AUC, cutpointdata$optimal_cutpoint,
           cutpointdata$sensitivity, cutpointdata$specificity,
           cutpointdata$acc))
}


createROC_multiple(data1[,143])


###creating ROC curves by column #
pdf(file = "~/PTY/Data/ROC Curves/ROC%03d.pdf", onefile = F, width = 7, height = 7)

for (i in 143:186) {
  createROC_multiple(data1[, i])
}
dev.off()
###

############Ok this time, create a df but out of the ROCit package if i can manage

createROC_df <- function(x) {

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
  newAUC <- newROC$AUC
  newSens <- newROC$TPR[cutoffindex]
  newSpec <- 1 - newROC$FPR[cutoffindex]
 
  #create and return vector w each variable
  ROCvector <- (c(newCutoff, newAUC, newSens, newSpec)) 
  return(ROCvector)
}

getROCcutoff <- function(x) {
  
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
  
  #create and return vector w each variable
  ROCvector <- (c(newCutoff)) 
  return(ROCvector)
}

getROCcutoff(data1$B2M_A)

df = NULL
write.csv(df,"~/PTY/Data\\final ones.csv", row.names = FALSE)
for (i in 184:186) { 
  
  rocdf <- createROC_df(data1[,i])
  df <- rbind(df, rocdf)
}
createROC_df(data1$IL.17F_D)

colnames(data1)



df = NULL
for (i in 89:96) { 
  
  rocdf <- getROCcutoff(data1[,i])
  df <- rbind(df, rocdf)
}
write.csv(df,"~/PTY/Data\\firstROCs.csv", row.names = FALSE)





