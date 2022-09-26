####lets do some MC4R analysis

#set wd and import data 
setwd("~/Uni/4th year/Dissertation")
mc4r <- read.csv("mc4r_phenotypes.csv", header = T, na.strings = ".")
mc4r2 <- read.csv("mc4r_phenotypes2.csv", header = T, na.strings = ".")
mc4r_noNA <- read.csv("mc4r_noNA.csv", header = T, na.strings = ".")
#lets generate some useful variables

# based on variable values
nonsense <- mc4r[mc4r$nonsense == '1',]
lof <- mc4r[mc4r$function_overall== 'LoF',]
WT <- mc4r[mc4r$function_overall== 'WT',]
gof <- mc4r[mc4r$function_overall== 'GoF',]
barplot()
getmean <- function(x) {mean<- mean(x, na.rm = T)
            return(mean)}
getmean(lof$bmi_mean)
mean(lof$bmi_mean, na.rm = T)

bardata_phenotypes <- cbind(c(getmean(nonsense$obesity_mean), getmean(lof$obesity_mean), getmean(WT$obesity_mean), getmean(gof$obesity_mean)), # Obesity Values
                   c(getmean(nonsense$severe_obesity), getmean(lof$severe_obesity), getmean(WT$severe_obesity), getmean(gof$severe_obesity)), #Severe obesity values
                   c(getmean(nonsense$diabetes_mean), getmean(lof$diabetes_mean), getmean(WT$diabetes_mean), getmean(gof$diabetes_mean))) #diabetes values
colnames(bardata_phenotypes) <- c("Obesity", "Sev. Obes.", "Diabetes")
rownames(bardata_phenotypes) <- c("Nonsense", "LoF", "WT", "GoF")
rownames <- rownames(bardata_phenotypes)

cols <- c("red", "purple", "blue", "green")
barplot(height = bardata_phenotypes*100,
        beside = TRUE,                        # Put the bars next to each other
        col = cols , main = "prevalence of phenotypes between groups", 
        ylab = "%", xlab = "phenotypes")
legend(x="topright", legend=rownames, cex=.8,bty="n", pch=15,
       col = cols, y.intersp = .4)


#bmi stuff
#c(getmean(nonsense$bmi_mean), getmean(lof$bmi_mean), getmean(WT$bmi_mean), getmean(gof$bmi_mean)),# BMI values
bardata_bmi <- c(getmean(nonsense$bmi_mean), getmean(lof$bmi_mean), getmean(WT$bmi_mean), getmean(gof$bmi_mean)) # Obesity Values
colnames(bardata_bmi) <- c("Nonsense", "LoF", "WT", "GoF")

barplot(height = bardata_bmi-min(bardata_bmi),
        beside = TRUE,      # Put the bars next to each other
        col = cols , main = "bmi between groups", 
        ylab = "difference in bmi from minimum", xlab = "groups", ylim = )
legend(x="topright", legend=rownames, cex=.8,bty="n", pch=15,
       col = cols, y.intersp = .4)

#####Histograms
hist(log(na.omit(mc4r2$betaarrestin_mean)))



#forest plots!!
library(forestplot) 

#####BETA ARRESTIN
beta_forest <- 
        structure(list(
                beta  = c(NA, -0.42781, -0.6067,0.59952), 
                lower = c(NA, -0.52068, -1.81, 0.0151656),
                upper = c(NA, -0.33494, 0.599, 1.047388)),
                .Names = c("beta", "lower", "upper"), 
                row.names = c(NA, -4L), 
                class = "data.frame")

beta_tabletext <- cbind(
        c("Class", "GoF","WT", "LoF"),
        c("??", "-0.43", "-0.61", "0.60"),
        c("95% CI", "(-0.52, -0.34)","(-1.81, 0.60)", "(0.015, 1.05)"),
        c("P-value", "1.7e-22", "0.324","8.7e-03"))

forestplot(beta_tabletext, beta_forest$beta, beta_forest$lower, beta_forest$upper,
           new_page = TRUE, xlab = "Beta (95% CI) of BMI",
           is.summary = c(TRUE,F, F,F), 
           boxsize=.04, graph.pos = 4, txt_gp = fpTxtGp(cex=0.75, xlab = gpar(cex=.55), ticks = gpar(cex=.6)),
           xlog = F, hrzl_lines = list("2" = gpar(lty = 1, col = "black")),
           col = fpColors(box = "royalblue", line = "black", summary = "royalblue"))


#####cAMP!!
camp_forest <- 
        structure(list(
                beta  = c(NA, -0.58886, -0.1398,0.997366), 
                lower = c(NA, -0.70855, -.27,0.133173),
                upper = c(NA, -0.46916, -0.0011, 1.861559)),
                .Names = c("beta", "lower", "upper"), 
                row.names = c(NA, -4L), 
                class = "data.frame")

camp_tabletext <- cbind(
        c("Class", "GoF", "WT","LoF"),
        c("??", "-0.59", "-0.14","0.99"),
        c("95% CI", "(-0.71, -0.47)","-0.28, -0.001", "(0.13, 1.86)"),
        c("P-value", "5.3e-22", "4.8e-02","2.4e-02"))

forestplot(camp_tabletext, camp_forest$beta, camp_forest$lower, camp_forest$upper,
           new_page = TRUE, xlab = "Beta (95% CI) of BMI",
           is.summary = c(TRUE,F, F,F), 
           boxsize=.08, graph.pos = 4, txt_gp = fpTxtGp(cex=0.75, xlab = gpar(cex=.55), ticks = gpar(cex=.7)),
           xlog = F, hrzl_lines = list("2" = gpar(lty = 1, col = "black")),
           col = fpColors(box = "royalblue", line = "black", summary = "royalblue"))


#####cell surface expression!! no WT included
cellsurf_forest <- 
        structure(list(
                beta  = c(NA, -0.4151, 0.627), 
                lower = c(NA, -0.99, -0.226),
                upper = c(NA, 0.168, 1.48)),
                .Names = c("beta", "lower", "upper"), 
                row.names = c(NA, -3L), 
                class = "data.frame")

cellsurf_tabletext <- cbind(
        c("Class", "GoF","LoF"),
        c("??", "-0.42", "0.63"),
        c("95% CI", "(-0.99, 0.17)", "(-0.22, 1.48)"),
        c("P-value", "0.16", "0.15"))

forestplot(cellsurf_tabletext, cellsurf_forest$beta, cellsurf_forest$lower, cellsurf_forest$upper,
           new_page = TRUE, xlab = "Beta (95% CI) of BMI",
           is.summary = c(TRUE,F, F), 
           boxsize=.08, graph.pos = 4, txt_gp = fpTxtGp(cex=0.75, xlab = gpar(cex=.55), ticks = gpar(cex=.7)),
           xlog = F, hrzl_lines = list("2" = gpar(lty = 1, col = "black")),
           col = fpColors(box = "royalblue", line = "black", summary = "royalblue"))


#####overall!!
all_forest <- 
        structure(list(
                beta  = c(NA, -0.42841, -0.523, 0.641096), 
                lower = c(NA, -0.52124, -1.831, 0.177572),
                upper = c(NA, -0.33557, 0.78, 1.10462)),
                .Names = c("beta", "lower", "upper"), 
                row.names = c(NA, -4L), 
                class = "data.frame")

all_tabletext <- cbind(
        c("Class", "GoF", "WT","LoF"),
        c("??", "-0.43", "-0.52","0.64"),
        c("95% CI", "(-0.52, -0.34)", "(-1.8, 0.79)","(0.18, 1.10)"),
        c("P-value", "1.5e-19", "0.43","6.7e-03"))

forestplot(all_tabletext, all_forest$beta, all_forest$lower, all_forest$upper,
           new_page = TRUE, xlab = "Beta (95% CI) of BMI",
           is.summary = c(TRUE,F, F,F), 
           boxsize=.08, graph.pos = 4, txt_gp = fpTxtGp(cex=0.75, xlab = gpar(cex=.55), ticks = gpar(cex=.7)),
           xlog = F, hrzl_lines = list("2" = gpar(lty = 1, col = "black")),
           col = fpColors(box = "royalblue", line = "black", summary = "royalblue"))

####Time to make some scatters!! 
###BETA-ARRESTIN
logBeta <- log(mc4r_noNA$betaarrestin_mean)
logBetaLower <- log(mc4r_noNA$betaarrestin_95lower)
logBetaUpper <- log(mc4r_noNA$betaarrestin_95upper)

plot(logBeta,mc4r_noNA$bmi_mean, pch=19, col = "blue", ylim = range(c(15,45)))
arrows(logBeta, mc4r_noNA$bmi_lowerCI, logBeta, mc4r_noNA$bmi_upperCI, length=0.05, angle=90, code=3)
arrows(logBetaLower, mc4r_noNA$bmi_mean, logBetaUpper, mc4r_noNA$bmi_mean, length=0.05, angle=90, code=3)
#beta_lm <- lm(mc4r_noNA$bmi_mean ~ logBeta)
beta_lm_sum <- summary(beta_lm)
beta_r2 <- round((beta_lm_sum$adj.r.squared),digits = 4)
beta_p <- round((beta_lm_sum$coefficients[2,4]), digits =7)
legend("bottomleft", legend= c("Pvalue= ",beta_p, "R2= ", beta_r2), bty="n",cex=.7)
abline(beta_lm, col = "red")
beta_lm_sum$coefficients


### UNLOGGED cAMP!!
Camp <- mc4r_noNA$camp_mean
CampLower <- mc4r_noNA$camp_95lower
CampUpper <- mc4r_noNA$camp_95upper

plot(Camp,mc4r_noNA$bmi_mean, pch=19, col = "blue", ylim = range(c(15,45)))
arrows(Camp, mc4r_noNA$bmi_lowerCI, Camp, mc4r_noNA$bmi_upperCI, length=0.05, angle=90, code=3)
arrows(CampLower, mc4r_noNA$bmi_mean, CampUpper, mc4r_noNA$bmi_mean, length=0.05, angle=90, code=3)
camp_lm <- lm(mc4r_noNA$bmi_mean ~ Camp)
camp_lm_sum <- summary(camp_lm)
camp_r2 <- round((camp_lm_sum$adj.r.squared),digits = 4)
camp_p <- round((camp_lm_sum$coefficients[2,4]), digits =4)
legend("bottomleft", legend= c("Pvalue= ",camp_p, "R2= ", camp_r2), bty="n",cex=.7)
abline(camp_lm, col = "red")


###cell surface expression!!!!
logCellSurf <- log(mc4r_noNA$cellsurface_mean)
logCellSurfLower <- log(mc4r_noNA$cellsurface_95lower)
logCellSurfUpper <- log(mc4r_noNA$cellsurface_95upper)

plot(logCellSurf,mc4r_noNA$bmi_mean, pch=19, col = "blue", ylim = range(c(15,45)))
arrows(logCellSurf, mc4r_noNA$bmi_lowerCI, logCellSurf, mc4r_noNA$bmi_upperCI, length=0.05, angle=90, code=3)
arrows(logCellSurfLower, mc4r_noNA$bmi_mean, logCellSurfUpper, mc4r_noNA$bmi_mean, length=0.05, angle=90, code=3)
#cellsurface_lm <- lm(mc4r_noNA$bmi_mean ~ logCellSurf)
#cellsurface_lm_sum <- summary(cellsurface_lm)
#cellsurface_r2 <- round((cellsurface_lm_sum$adj.r.squared),digits = 4)
#cellsurface_p <- round((cellsurface_lm_sum$coefficients[2,4]), digits =4)
legend("bottomleft", legend= c("Pvalue= ",cellsurface_p, "R2= ", cellsurface_r2), bty="n",cex=.7)
abline(cellsurface_lm, col = "red")


#####################################LOGGED BMI TIME
logBMI <- log(mc4r_noNA$bmi_mean)
logBMIlower <- log(mc4r_noNA$bmi_lowerCI)
logBMIupper <- log(mc4r_noNA$bmi_upperCI)

###BETA-ARRESTIN
logBeta <- log(mc4r_noNA$betaarrestin_mean)
logBetaLower <- log(mc4r_noNA$betaarrestin_95lower)
logBetaUpper <- log(mc4r_noNA$betaarrestin_95upper)

plot(logBeta,logBMI, pch=19, col = "blue", ylim = range(c(3,4)))
arrows(logBeta, logBMIlower, logBeta, logBMIupper, length=0.05, angle=90, code=3)
arrows(logBetaLower, logBMI, logBetaUpper, logBMI, length=0.05, angle=90, code=3)
beta_lm <- lm(logBMI ~ logBeta)
beta_lm_sum <- summary(beta_lm)
beta_r2 <- round((beta_lm_sum$adj.r.squared),digits = 4)
beta_p <- round((beta_lm_sum$coefficients[2,4]), digits =7)
legend("bottomleft", legend= c("Pvalue= ",beta_p, "R2= ", beta_r2), bty="n",cex=.7)
abline(beta_lm, col = "red")
beta_lm_sum$coefficients


### UNLOGGED cAMP!!
Camp <- mc4r_noNA$camp_mean
CampLower <- mc4r_noNA$camp_95lower
CampUpper <- mc4r_noNA$camp_95upper

######YOU GOT TO HERE!!!!!

plot(Camp,logBMI, pch=19, col = "blue", ylim = range(c(2,4)))
arrows(Camp, logBetaLower, Camp, logBetaUpper, length=0.05, angle=90, code=3)
arrows(CampLower, logBMI, CampUpper, logBMI, length=0.05, angle=90, code=3)
camp_lm <- lm(logBMI ~ Camp)
camp_lm_sum <- summary(camp_lm)
camp_r2 <- round((camp_lm_sum$adj.r.squared),digits = 4)
camp_p <- round((camp_lm_sum$coefficients[2,4]), digits =4)
legend("bottomleft", legend= c("Pvalue= ",camp_p, "R2= ", camp_r2), bty="n",cex=.7)
abline(camp_lm, col = "red")


###cell surface expression!!!!
logCellSurf <- log(mc4r_noNA$cellsurface_mean)
logCellSurfLower <- log(mc4r_noNA$cellsurface_95lower)
logCellSurfUpper <- log(mc4r_noNA$cellsurface_95upper)

plot(logCellSurf,mc4r_noNA$bmi_mean, pch=19, col = "blue", ylim = range(c(15,45)))
arrows(logCellSurf, mc4r_noNA$bmi_lowerCI, logCellSurf, mc4r_noNA$bmi_upperCI, length=0.05, angle=90, code=3)
arrows(logCellSurfLower, mc4r_noNA$bmi_mean, logCellSurfUpper, mc4r_noNA$bmi_mean, length=0.05, angle=90, code=3)
#cellsurface_lm <- lm(mc4r_noNA$bmi_mean ~ logCellSurf)
#cellsurface_lm_sum <- summary(cellsurface_lm)
#cellsurface_r2 <- round((cellsurface_lm_sum$adj.r.squared),digits = 4)
#cellsurface_p <- round((cellsurface_lm_sum$coefficients[2,4]), digits =4)
legend("bottomleft", legend= c("Pvalue= ",cellsurface_p, "R2= ", cellsurface_r2), bty="n",cex=.7)
abline(cellsurface_lm, col = "red")


#####time to make some bar plots!
install.packages("ggpubr")
library(ggplot2)
library(ggpubr)

###cAMP!
camp_binned <- c(0,25,50,75,100,125,150)
bmi_cbin <- c(32.44,27.16,28.74,27.71,27.21,27.6,26.79)
SE_cbin <- c(2.119792443,0.444158899,0.872788634,0.716000254,0.069282198,0.332740246,0.057580786)
camp_binned_df<- data.frame(camp_binned, bmi_cbin, 
                      SE_cbin)

camp_binned_df %>%
        ggplot(aes(x = camp_binned, y = bmi_cbin)) +
        geom_col(width=20, fill = "steelblue")+ coord_cartesian(ylim=c(20,35))+
        geom_hline(yintercept=27.4, linetype="dashed")  + 
        geom_errorbar(aes(ymin=bmi_cbin-SE_cbin, ymax=bmi_cbin+SE_cbin),
                      width=5, position=position_dodge(.9), color="black")

###betaarrestin!!
betalog_binned <- c(25, 50, 75, 100, 125, 150)
bmi_betalogbin <- c(28.28,27.68,27.64,26.67,27.2,26.79)
SE_betalogbin <- c(0.950464002,0.526782688,0.290297089,0.649020855,0.069760993,0.05772223)
betalog_bin_df<- data.frame(betalog_binned, bmi_betalogbin, 
                            SE_betalogbin)

betalog_bin_df %>%
        ggplot(aes(x = betalog_binned, y = bmi_betalogbin)) +
        geom_col(width=20, fill = "steelblue")+ coord_cartesian(ylim=c(20,30))+
        geom_hline(yintercept=27.4, linetype="dashed")  + 
        scale_x_continuous(labels = betalog_binned, breaks = c(1:6)) +
        geom_errorbar(aes(ymin=bmi_betalogbin-SE_betalogbin, ymax=bmi_betalogbin+SE_betalogbin),
                      width=5, position=position_dodge(.9), color="black")


###cell surface expression!!
cs_binned <- c(25, 50, 75, 100, 125, 150)
bmi_csbin <- c(28.04,27.84,27.2,26.83,27.39,26.27)
SE_csbin <- c(0.447399264,0.506934269,0.069851362,0.057031396,1.036179987,0.921595838)
cs_bin_df<- data.frame(cs_binned, bmi_csbin, 
                       SE_csbin)

cs_bin_df %>%
        ggplot(aes(x = cs_binned, y = bmi_csbin)) +
        geom_col(width=20, fill = "steelblue")+ coord_cartesian(ylim=c(20,30))+
        geom_hline(yintercept=27.4, linetype="dashed")  + 
        scale_x_continuous(labels = cs_binned, breaks = c(1:6)) +
        geom_errorbar(aes(ymin=bmi_csbin-SE_csbin, ymax=bmi_csbin+SE_csbin),
                      width=5, position=position_dodge(.9), color="black")



######New plots!!!!!!!!
beta_binned <- c("0", "25", "50", "100", "200", "250")
bmi_betabin <- c(28.86,27.23,27.7,26.58,27.2,26.8)
SE_betabin <- c(0.543219771,0.440426509,0.337726821,0.583298178,0.069697843,0.057594808)
beta_bin_df<- data.frame(beta_binned, bmi_betabin, SE_betabin)
beta_bin_df$beta_binned <- factor(beta_bin_df$beta_binned, levels = c("0", "25", "50", "100", "200", "250"))

camp_binned <- c("0","25","50","75","100","125","150")
bmi_cbin <- c(32.44,27.16,28.74,27.71,27.21,27.6,26.79)
SE_cbin <- c(2.119792443,0.444158899,0.872788634,0.716000254,0.069282198,0.332740246,0.057580786)
camp_bin_df<- data.frame(camp_binned, bmi_cbin, SE_cbin)
camp_bin_df$camp_binned <- factor(camp_bin_df$camp_binned, levels = c("0","25","50","75","100","125","150"))


install.packages("gridExtra")
library(gridExtra)

BMI_Yaxis <- expression(paste("BMI (kg/m"^" 2", ")"))
beta_P <- expression(paste("P-value= 1.3x10"^" -7"))
camp_P <- expression(paste("P-value= 1.7x10"^" -8"))
###Beta arrestin NEW plot
BetaPlot <- ggplot(beta_bin_df, aes(beta_binned, bmi_betabin))+
                   geom_col(fill = "steelblue")+ coord_cartesian(ylim=c(20,30))+
                   geom_hline(yintercept=27.4, linetype="dashed")+
                   geom_errorbar(aes(ymin=bmi_betabin-SE_betabin, ymax=bmi_betabin+SE_betabin),
                      width=.2, position=position_dodge(.9), color="black")+
                   labs(y=BMI_Yaxis,x="Beta-Arrestin Recruitment", title = "Beta Arrestin")+
                   theme(plot.title = element_text(hjust = 0.5)) +
                   annotate("text", x= Inf, y= Inf, label = "P-value= 1.3E-7", hjust = 1.3, vjust =2, fontface=2)
CampPlot <- ggplot(camp_bin_df, aes(camp_binned, bmi_cbin))+
                   geom_col(fill = "steelblue")+ coord_cartesian(ylim=c(20,35))+
                   geom_hline(yintercept=27.4, linetype="dashed")+
                   geom_errorbar(aes(ymin=bmi_cbin-SE_cbin, ymax=bmi_cbin+SE_cbin),
                      width=.2, position=position_dodge(.9), color="black")+
                   labs(y=BMI_Yaxis,x="cAMP Production", title = "cAMP")+
                   theme(plot.title = element_text(hjust = 0.5)) +
                   annotate("text", x= Inf, y= Inf, label = "P-value= 1.7E-8", hjust = 1.3, vjust =2, fontface=2)
            
par(mfcol=c(1,2))
BetaPlot
CampPlot
grid.arrange(BetaPlot, CampPlot, ncol=2)
?grid.arrange