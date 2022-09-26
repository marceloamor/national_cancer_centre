###Analysis for PTY Paper 

setwd("~/PTY/Data")

data1 <- read.csv(file = "Full Cohort ELISA Data (IL6+Validation+Uplex).csv", header = TRUE)

#install.packages("plyr")
#install.packages("dplyr")
#install.packages("cutpointr")
#install.packages("reshape2")
#install.packages("epitools")
#install.packages("ggpubr")
#install.packages('forestplot')
#install.packages("survival")
#install.packages("survBootOutliers")


#library(plyr)
library(dplyr)
library(ROCit)
library(cutpointr)
library(ggplot2)
library(reshape2)
library(epitools)
library(ggpubr)
library(forestplot)
library(survival)
library(survBootOutliers)
#detach("package:plyr")
detach("package:dplyr")
installed.packages()
#increases threshold for scientific notation
options("scipen"=100, "digits"=5)

###summary stats

dplyr::group_by(data1, DSM) %>%
  dplyr::summarise(
    average = mean(age, na.rm = TRUE),
    SDage = sd(age, na.rm = TRUE),
    BMI = mean(BMI, na.rm = TRUE),
    SDBMI = sd(BMI, na.rm = TRUE),
    education = mean(Education, na.rm = TRUE),
    SDedu = sd(Education, na.rm = TRUE)
  )

#first function, returns values grouped by DSM
DSMValues <- function(variable) {
  DSMgrouped <- dplyr::group_by(data1, DSM) %>%
  dplyr::summarise(
    average = mean(variable, na.rm = TRUE),
    SD = sd(variable, na.rm = TRUE))

  CaseAvg <- DSMgrouped$average[2]
  ControlAvg <- DSMgrouped$average[1]
  CaseSD <- DSMgrouped$SD[2]
  ControlSD <- DSMgrouped$SD[1]
  
  return(c(CaseAvg, CaseSD, ControlAvg, ControlSD))
  return(DSMgrouped)
}
DSMValues(data1$Education)

###second function, returns all values, grouped and nongrouped
##OLD DPLYR VERSION OF FUNCTION
#getPvalue1.0 <- function(variable) {
  #Collective values
#  mean1 <- mean(variable, na.rm = TRUE)
#  sdMean <- sd(variable, na.rm = TRUE)
#  MannU <- wilcox.test(variable ~ DSM, data = data1, exact = FALSE)

#  #GroupedBy values 
#  DSMgrouped <- dplyr::group_by(data1, DSM) %>%
#    dplyr::summarise(
#      average = mean(variable, na.rm = TRUE),
#      SD = sd(variable, na.rm = TRUE),
#      .preserve = T)
  
#  CaseAvg <- DSMgrouped$average[2]
#  ControlAvg <- DSMgrouped$average[1]
#  CaseSD <- DSMgrouped$SD[2]
#  ControlSD <- DSMgrouped$SD[1]
  
  #create and return dataframe
#  df1 <- data.frame()
#  labels1 <- c("Mean", "SD", "CaseMean", "CaseSD", "ControlMean", "ControlSD", "MannU")
#  values1 <- c(mean1, sdMean, CaseAvg, CaseSD, ControlAvg, ControlSD, MannU$p.value)
#  df1 <- rbind(df1, data.frame(labels1, values1))
#  return(df1)
#}
 
###NEW version of getPvalue function!
getPvalue <- function(variable) {
  #Collective values
  mean1 <- mean(variable, na.rm = TRUE)
  sdMean <- sd(variable, na.rm = TRUE)
  MannU <- wilcox.test(variable ~ DSM, data = data1, exact = FALSE, alternative = "two.sided")
  
  #GroupedBy values 
  DSMmean <- tapply(variable, data1$DSM, mean, na.rm = T)
  DSMsd <- tapply(variable, data1$DSM, sd, na.rm = T)
  
  CaseAvg <- DSMmean[2]
  ControlAvg <- DSMmean[1]
  CaseSD <- DSMsd[2]
  ControlSD <- DSMsd[1]
  
  #create and return dataframe
  df1 <- data.frame()
  labels1 <- c("Mean", "SD", "CaseMean", "CaseSD", "ControlMean", "ControlSD", "MannU")
  values1 <- c(mean1, sdMean, CaseAvg, CaseSD, ControlAvg, ControlSD, MannU$p.value)
  df1 <- rbind(df1, data.frame(labels1, values1))
  return(df1)
}

getPvalue(data1$Bleeding)

#Mann Whitney U
MannU <- wilcox.test(DSM ~ IL.10_A, data = data1,
                   exact = F, na.rm = T, alternative = "two.sided")
t.test(data1$age, data1$Bleeding, alternative = "two.sided", var.equal = F)
?t.test
t.test(data1$IL6_timeA, data1$age, paired = F)
MannU$p.value
lmIL6Age<-lm(data1$age ~ data1$IL6_timeA)
plot(data1$IL6_timeB, data1$age)

abline(a= lmIL6Age$coefficients[1], b = lmIL6Age$coefficients[2])

#start mass producing stats now 
ageStat <- getPvalue(data1$age)
sexStat <- getPvalue(isFemaleStat$sex)
EducationStat <- getPvalue(data1$Education)
WorkStat <- getPvalue(data1$Work)
MaritalStat <- getPvalue(data1$MaritalStatus)
SportsStat <- getPvalue(data1$Sports)
SittingStat <- getPvalue(data1$Sitting)
MMSEStat <- getPvalue(data1$MMSE_baseline)
AnxietyStat <- getPvalue(data1$HADS)
BleedingStat <- getPvalue(data1$Bleeding)
AlcoholStat <- getPvalue(data1$CAGE)

ttest <- t.test(data1$Sitting, data1$IL6_timeB)
getPvalue(data1$IL.10_A)
ttest$p.value
MannU <- t.test(IL.10_C ~ IL6_timeB, data = data1, exact = FALSE, alternative = "two.sided")

kruskal.test(IL6_timeB ~ DSM, data = data1)


#sexStat fix - isFemaleStat
as.numeric(as.factor(data1$sex, na.rm = T))
sexdata <- data1$sex
sexdata <- revalue(sexdata, c("F"=1))
sexdata <- revalue(sexdata, c("M"=0))
length(data1$sex)
isFemaleStat <- data1 %>%
  mutate(sex = ifelse(MaritalStatus=='2', 1,
                      ifelse(sex == '1', 0, 0)))

summary(sexdata)
sum(isFemaleStat$sex, na.rm = T)

print(chisq.test(table(data1$MaritalStatus, data1$DSM)))

table(sexdata)
###Biomarkers now
#IL6
IL6_A <- getPvalue(data1$IL6_timeA)
IL6_B <- getPvalue(data1$IL6_timeB)
IL6_C <- getPvalue(data1$IL6_timeC)
IL6_D <- getPvalue(data1$IL6_timeD)
IL6_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
       list(IL6_A, IL6_B, IL6_C, IL6_D))
#TGFb1
TGF_A <- getPvalue(data1$TGFb1_timeA)
TGF_B <- getPvalue(data1$TGFb1_timeB)
TGF_C <- getPvalue(data1$TGFb1_timeC)
TGF_D <- getPvalue(data1$TGFb1_timeD)
TGF_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                  list(TGF_A, TGF_B, TGF_C, TGF_D))
#B2M
B2M_A <- getPvalue(data1$B2M_A)
B2M_B <- getPvalue(data1$B2M_B)
B2M_C <- getPvalue(data1$B2M_C)
B2M_D <- getPvalue(data1$B2M_D)
B2M_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                  list(TGF_A, TGF_B, TGF_C, TGF_D))
#Prolactin
Prolactin_A <- getPvalue(data1$Prolactin_A)
Prolactin_B <- getPvalue(data1$Prolactin_B)
Prolactin_C <- getPvalue(data1$Prolactin_C)
Prolactin_D <- getPvalue(data1$Prolactin_D)
Prolactin_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                  list(Prolactin_A, Prolactin_B, Prolactin_C, Prolactin_D))

#ST2
ST2_A <- getPvalue(data1$ST2_A)
ST2_B <- getPvalue(data1$ST2_B)
ST2_C <- getPvalue(data1$ST2_C)
ST2_D <- getPvalue(data1$ST2_D)
ST2_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                  list(ST2_A, ST2_B, ST2_C, ST2_D))

#MIG
MIG_A <- getPvalue(data1$MIG_A)
MIG_B <- getPvalue(data1$MIG_B)
MIG_C <- getPvalue(data1$MIG_C)
MIG_D <- getPvalue(data1$MIG_D)
MIG_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                  list(MIG_A, MIG_B, MIG_C, MIG_D))

#Gas6
Gas6_A <- getPvalue(data1$Gas6_A)
Gas6_B <- getPvalue(data1$Gas6_B)
Gas6_C <- getPvalue(data1$Gas6_C)
Gas6_D <- getPvalue(data1$Gas6_D)
Gas6_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                  list(Gas6_A, Gas6_B, Gas6_C, Gas6_D))

#TSH
TSH_A <- getPvalue(data1$TSH_A)
TSH_B <- getPvalue(data1$TSH_B)
TSH_C <- getPvalue(data1$TSH_C)
TSH_D <- getPvalue(data1$TSH_D)
TSH_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                  list(TSH_A, TSH_B, TSH_C, TSH_D))

#ALK6
ALK6_A <- getPvalue(data1$ALK6_A)
ALK6_B <- getPvalue(data1$ALK6_B)
ALK6_C <- getPvalue(data1$ALK6_C)
ALK6_D <- getPvalue(data1$ALK6_D)
ALK6_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                  list(ALK6_A, ALK6_B, ALK6_C, ALK6_D))

#CD6
CD6_A <- getPvalue(data1$CD6_A)
CD6_B <- getPvalue(data1$CD6_B)
CD6_C <- getPvalue(data1$CD6_C)
CD6_D <- getPvalue(data1$CD6_D)
CD6_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                  list(CD6_A, CD6_B, CD6_C, CD6_D))

#Thyroglobulin
Thyro_A <- getPvalue(data1$Thyro_A)
Thyro_B <- getPvalue(data1$Thyro_B)
Thyro_C <- getPvalue(data1$Thyro_C)
Thyro_D <- getPvalue(data1$Thyro_D)
Thyro_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                  list(Thyro_A, Thyro_B, Thyro_C, Thyro_D))

#G.CSF
G.CSF_A <- getPvalue(data1$G.CSF_A)
G.CSF_B <- getPvalue(data1$G.CSF_B)
G.CSF_C <- getPvalue(data1$G.CSF_C)
G.CSF_D <- getPvalue(data1$G.CSF_D)
G.CSF_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                  list(G.CSF_A, G.CSF_B, G.CSF_C, G.CSF_D))

#GM.CSF
GM.CSF_A <- getPvalue(data1$GM.CSF_A)
GM.CSF_B <- getPvalue(data1$GM.CSF_B)
GM.CSF_C <- getPvalue(data1$GM.CSF_C)
GM.CSF_D <- getPvalue(data1$GM.CSF_D)
GM.CSF_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                    list(GM.CSF_A, GM.CSF_B, GM.CSF_C, GM.CSF_D))

#IL.1RA
IL.1RA_A <- getPvalue(data1$IL.1RA_A)
IL.1RA_B <- getPvalue(data1$IL.1RA_B)
IL.1RA_C <- getPvalue(data1$IL.1RA_C)
IL.1RA_D <- getPvalue(data1$IL.1RA_D)
IL.1RA_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                    list(IL.1RA_A, IL.1RA_B, IL.1RA_C, IL.1RA_D))

#IL.10
IL.10_A <- getPvalue(data1$IL.10_A)
IL.10_B <- getPvalue(data1$IL.10_B)
IL.10_C <- getPvalue(data1$IL.10_C)
IL.10_D <- getPvalue(data1$IL.10_D)
IL.10_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                     list(IL.10_A, IL.10_B, IL.10_C, IL.10_D))

#IL.17F
IL.17F_A <- getPvalue(data1$IL.17F_A)
IL.17F_B <- getPvalue(data1$IL.17F_B)
IL.17F_C <- getPvalue(data1$IL.17F_C)
IL.17F_D <- getPvalue(data1$IL.17F_D)
IL.17F_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                     list(IL.17F_A, IL.17F_B, IL.17F_C, IL.17F_D))

#IL.21
IL.21_A <- getPvalue(data1$IL.21_A)
IL.21_B <- getPvalue(data1$IL.21_B)
IL.21_C <- getPvalue(data1$IL.21_C)
IL.21_D <- getPvalue(data1$IL.21_D)
IL.21_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                     list(IL.21_A, IL.21_B, IL.21_C, IL.21_D))

#IL.33
IL.33_A <- getPvalue(data1$IL.33_A)
IL.33_B <- getPvalue(data1$IL.33_B)
IL.33_C <- getPvalue(data1$IL.33_C)
IL.33_D <- getPvalue(data1$IL.33_D)
IL.33_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                    list(IL.33_A, IL.33_B, IL.33_C, IL.33_D))

#MCP.2
MCP.2_A <- getPvalue(data1$MCP.2_A)
MCP.2_B <- getPvalue(data1$MCP.2_B)
MCP.2_C <- getPvalue(data1$MCP.2_C)
MCP.2_D <- getPvalue(data1$MCP.2_D)
MCP.2_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                     list(MCP.2_A, MCP.2_B, MCP.2_C, MCP.2_D))

#MIP.beta
MIP3.beta_A <- getPvalue(data1$MIP3.beta_A)
MIP3.beta_B <- getPvalue(data1$MIP3.beta_B)
MIP3.beta_C <- getPvalue(data1$MIP3.beta_C)
MIP3.beta_D <- getPvalue(data1$MIP3.beta_D)
MIP3.beta_All <- Reduce(function(x,y) merge(x = x, y = y, by = "labels1"), 
                     list(MIP3.beta_A, MIP3.beta_B, MIP3.beta_C, MIP3.beta_D))


####Quick Call List
ageStat <- getPvalue(data1$age)
sexStat <- getPvalue(data1$sex)
EducationStat <- getPvalue(data1$Education)
WorkStat <- getPvalue(data1$Work)
MaritalStat <- getPvalue(data1$MaritalStatus)
SportsStat <- getPvalue(data1$Sports)
SittingStat <- getPvalue(data1$Sitting)
MMSEStat <- getPvalue(data1$MMSE_baseline)

IL6_All
#write.csv(IL6_All,"~/PTY/Data\\IL6_All.csv", row.names = FALSE)

TGF_All
#write.csv(TGF_All,"~/PTY/Data\\TGF_All.csv", row.names = FALSE)

B2M_All
#write.csv(B2M_All,"~/PTY/Data\\B2M_All.csv", row.names = FALSE)

Prolactin_All
#write.csv(Prolactin_All,"~/PTY/Data\\Prolactin_All.csv", row.names = FALSE)

ST2_All
#write.csv(ST2_All,"~/PTY/Data\\ST2_All.csv", row.names = FALSE)

MIG_All
#write.csv(MIG_All,"~/PTY/Data\\MIG_All.csv", row.names = FALSE)

Gas6_All
#write.csv(Gas6_All,"~/PTY/Data\\Gas6_All.csv", row.names = FALSE)

TSH_All
#write.csv(TGF_All,"~/PTY/Data\\TSH_All.csv", row.names = FALSE)

ALK6_All
#write.csv(ALK6_All,"~/PTY/Data\\ALK6_All.csv", row.names = FALSE)

CD6_All
#write.csv(CD6_All,"~/PTY/Data\\CD6_All.csv", row.names = FALSE)

Thyro_All
#write.csv(Thyro_All,"~/PTY/Data\\Thyro_All.csv", row.names = FALSE)

G.CSF_All
#write.csv(G.CSF_All,"~/PTY/Data\\G.CSF_All.csv", row.names = FALSE)

GM.CSF_All
#write.csv(GM.CSF_All,"~/PTY/Data\\GM.CSF_All.csv", row.names = FALSE)

IL.1RA_All
#write.csv(IL.1RA_All,"~/PTY/Data\\IL.1RA_All.csv", row.names = FALSE)

IL.10_All
#write.csv(IL.10_All,"~/PTY/Data\\IL.10_All.csv", row.names = FALSE)

IL.17F_All
#write.csv(IL.17F_All,"~/PTY/Data\\IL.17F_All.csv", row.names = FALSE)

IL.21_All
#write.csv(IL.21_All,"~/PTY/Data\\IL.21_All.csv", row.names = FALSE)

IL.33_All
#write.csv(IL.33_All,"~/PTY/Data\\IL.33_All.csv", row.names = FALSE)

MCP.2_All
#write.csv(MCP.2_All,"~/PTY/Data\\MCP.2_All.csv", row.names = FALSE)

MIP3.beta_All
#write.csv(MIP3.beta_All,"~/PTY/Data\\MIP3.beta_All.csv", row.names = FALSE)


#################ROC TIMEEEE

rocIL6 <- data1 %>% filter(!is.na(data1$IL6_timeA))
rocIL6A <- rocit(score = rocIL6$IL6_timeA, class = rocIL6$DSM, negref = "0", method = "non")
rocIL6B <- rocit(score = rocIL6$IL6_timeB, class = rocIL6$DSM, negref = "0", method = "non")
rocIL6C <- rocit(score = rocIL6$IL6_timeC, class = rocIL6$DSM, negref = "0", method = "non")
rocIL6D <- rocit(score = rocIL6$IL6_timeD, class = rocIL6$DSM, negref = "0", method = "non")


rocIL6Bnon <- rocit(score = data1$IL6_timeB, class = dataset$DSM,
                    negref = "0", method = "non")
rocIL6Bplot <- plot(rocIL6B, legend = F, add = T)

##ok this syntax works to create and plot each ROC
plot(rocIL6D <- rocit(score = rocIL6$IL6_timeD, class = rocIL6$DSM, negref = "0", method = "non")
, legend = F, add = T)

IL6test1 <- data1 %>% filter(!is.na(data1$B2M_B))
IL6test1$IL6_timeB

B2Mtest1 <- data1 %>% filter(!is.na(data1$B2M_B))

rocB2Mtest1 <- rocit(score = B2Mtest1$B2M_B, class = B2Mtest1$DSM,
                    negref = "0", method = "non")
rocB2Mtest2 <- rocit(score = B2Mtest1$B2M_C, class = B2Mtest1$DSM,
                     negref = "0", method = "non")
rocB2Mtest3 <- rocit(score = B2Mtest1$B2M_D, class = B2Mtest1$DSM,
                     negref = "0", method = "non")
length(rocB2Mtest1$B2M_B)





####


cptest <- cutpointr(data1, age, DSM, na.rm = T,
      method = maximize_metric, metric = sum_sens_spec)

cptestIL6 <- cutpointr(data1, IL6_timeB, DSM, na.rm = T,
                   method = maximize_metric, metric = sum_sens_spec)



