####OR Script!!!
#OK for an OR script what do I need to do
#for each variable:
  #create NA vector and NA df
  #create variables for: 
    #number of high il6, delirium
    #number of low il6, delirium
    #number of high il6, no delirium
    #number of low il6, no delirium
  #return them in a 2x2 matrix
  #run the OR package function on it 
  #save the OR and CI as a vector and bind it to a df

createORs <- function(x) {
  #x <- data1$IL6_timeB
  ROCtable <- createNAtable(x) #create dataset w/out NAs specific to each analyte
  x_NA <-x[!is.na(x)]          #create analyte vector w/out NAs
  
      #create ROC and KSplot variables
      newROC <- rocit(score = x_NA, class = ROCtable$DSM, negref = "0", method = "non")
      newKSplot <- ksplot(newROC, legend = F, xlim = 6000)
      #get index of cutoff
      cutoffindex <- which(abs(newROC$Cutoff - newKSplot$`KS Cutoff`) == min(abs(newROC$Cutoff- newKSplot$`KS Cutoff`)))
      newCutoff <- newKSplot$`KS Cutoff`

  #create each variable
  highcases <- subset(x, !is.na(x) & data1$DSM == "1" & x > newCutoff)
  lowcases <- subset(x, !is.na(x) & data1$DSM == "1" & x < newCutoff)
  highnoncases <- subset(x, !is.na(x) & data1$DSM == "0" & x > newCutoff)
  lownoncases <- subset(x, !is.na(x) & data1$DSM == "0" & x < newCutoff)
  
  #create and return 2x2 matrix of sens/spec table 
  ORtable <- matrix(c(length(lownoncases), length(highnoncases),
                      length(lowcases), length(highcases)),nrow = 2, ncol = 2)
  return(ORtable)
}

length(subset(data1$IL6_timeB, data1$DSM == "1", !is.na(data1$IL6_timeB)))

createORs(data1$IL6_timeB)

IL6ORtest <- createORs(data1$MCP.2_C)


IL6ORtest2 <- oddsratio.wald(IL6ORtest)
IL6ORtest2$measure[2,1]

df1 = NULL
write.csv(df1,"~/PTY/Data\\finalonessss.csv", row.names = FALSE)
for (i in 170:170) { 
  
  ORmatrix <- createORs(data1[,i])
  ORdf <- oddsratio.wald(ORmatrix)
  ORvector <- c(ORdf$measure[2,1], ORdf$measure[2,2], ORdf$measure[2,3])
  df1 <- rbind(df1, ORvector)
  return(df1)
}

sqrt()

colnames(data1)
x <- c(1, 2, 3, 4, 5)
y <- c(2, 4, 6, 8, 10)
plot(x, y)

colnames(data1)
createORs(data1$IL.17F_D)

#median test higher or lower 
mediantest <- subset(data1$IL.17F_A, data1$DSM == "1")
mediantest1 <- subset(data1$IL.17F_A, data1$DSM == "0")
mean(mediantest1, na.rm= T)




createSEs <- function(x) {
  #x <- data1$IL6_timeB
  ROCtable <- createNAtable(x) #create dataset w/out NAs specific to each analyte
  x_NA <-x[!is.na(x)]          #create analyte vector w/out NAs
  
  #create ROC and KSplot variables
  newROC <- rocit(score = x_NA, class = ROCtable$DSM, negref = "0", method = "non")
  newKSplot <- ksplot(newROC, legend = F, xlim = 6000)
  #get index of cutoff
  cutoffindex <- which(abs(newROC$Cutoff - newKSplot$`KS Cutoff`) == min(abs(newROC$Cutoff- newKSplot$`KS Cutoff`)))
  newCutoff <- newKSplot$`KS Cutoff`
  
  #create each variable
  highcases <- subset(x, !is.na(x) & data1$DSM == "1" & x > newCutoff)
  lowcases <- subset(x, !is.na(x) & data1$DSM == "1" & x < newCutoff)
  highnoncases <- subset(x, !is.na(x) & data1$DSM == "0" & x > newCutoff)
  lownoncases <- subset(x, !is.na(x) & data1$DSM == "0" & x < newCutoff)
  
  #calculate SE of OR now 
  SE <- sqrt((1/length(highcases)) + (1/length(lowcases)) +
               (1/length(highnoncases)) +(1/length(lownoncases)))

  return(SE)
  
  }
createSEs(data1$IL6_timeA)

df2 <- c()

for (i in 143:186) { 
  
  SEvector <- createSEs(data1[,i])
  
  #SEvector <- 
  df2 <- rbind(df2, SEvector)
  print(SEvector)
}

createSEs(data1$ST2_D)
log(4.72)






