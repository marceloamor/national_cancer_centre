############making boxplots
##i want a function then loop that creates my boxplots
#paired boxplots taking a 4 continuous vectors as input
#timepoints ABCD
#and a categorical vector to distinguish cases and non
#and return them side by side AaBbCcDd in a boxplot

#do it in ggplot2???


mmIL6 = melt(data = data1, id.vars = c("DSM"), measure.vars = c("IL6_timeA", 'IL6_timeB', 'IL6_timeC', 'IL6_timeD'), na.rm = T)
mmST2 = melt(data = data1, id.vars = c("DSM"), measure.vars = c("ST2_A", 'ST2_B', 'ST2_C', 'ST2_D'), na.rm = T)

mmIL10 = melt(data = data1, id.vars = c("DSM"), measure.vars = c("IL.10_A", 'IL.10_B', 'IL.10_C', 'IL.10_D'), na.rm = T)
mmIL17F = melt(data = data1, id.vars = c("DSM"), measure.vars = c("IL.17F_A", 'IL.17F_B', 'IL.17F_C', 'IL.17F_D'), na.rm = T)

mmTGFb1 = melt(data = data1, id.vars = c("DSM"), measure.vars = c("TGFb1_timeA", 'TGFb1_timeB', 'TGFb1_timeC', 'TGFb1_timeD'), na.rm = T)

#il6
pdf("IL6_Boxplot_1.pdf")
ggplot(mmIL6, aes(x=variable, y=value, fill=factor(DSM))) +
  geom_boxplot() + scale_y_continuous(limits = c(0, 3500)) + 
  theme_bw() + scale_fill_discrete(name = "Group", labels = c("Non-delirium", "Delirium")) +
  labs(x = "Postoperative Timepoint", y = "Biomarker Concentration (pg/ml)") +
  scale_x_discrete(labels = c("IL6_timeA" = "Baseline", "IL6_timeB" = "aICU", "IL6_timeC" = "POD1", "IL6_timeD" = "POD4-7")) +
  ggtitle("Postoperative IL-6 Concentrations between Groups") + stat_compare_means(aes(group = DSM), label = "p.signif", hide.ns = T) 
dev.off()


#st2
dev.new()
ggplot(mmST2, aes(x=variable, y=value, fill=factor(DSM))) +
  geom_boxplot() + scale_y_continuous(limits = c(0, 3500)) + 
  theme_bw() + scale_fill_discrete(name = "Group", labels = c("Non-delirium", "Delirium")) +
  labs(x = "Postoperative Timepoint", y = "Biomarker Concentration (pg/ml)") +
  scale_x_discrete(labels = c("IL6_timeA" = "Baseline", "IL6_timeB" = "aICU", "IL6_timeC" = "POD1", "IL6_timeD" = "POD4-7")) +
  ggtitle("Postoperative ST2 Concentrations between Groups") + stat_compare_means(aes(group = DSM), label = "p.signif", hide.ns = T) 


#IL10
dev.new()
ggplot(mmIL10, aes(x=variable, y=value, fill=factor(DSM))) +
  geom_boxplot() + scale_y_continuous(limits = c(0, 35)) + 
  theme_bw() + scale_fill_discrete(name = "Group", labels = c("Non-delirium", "Delirium")) +
  labs(x = "Postoperative Timepoint", y = "Biomarker Concentration (pg/ml)") +
  scale_x_discrete(labels = c("IL6_timeA" = "Baseline", "IL6_timeB" = "aICU", "IL6_timeC" = "POD1", "IL6_timeD" = "POD4-7")) +
  ggtitle("Postoperative ST2 Concentrations between Groups") + stat_compare_means(aes(group = DSM), label = "p.signif", hide.ns = T) 

#IL17F
dev.new()
ggplot(mmIL17F, aes(x=variable, y=value, fill=factor(DSM))) +
  geom_boxplot() + scale_y_continuous(limits = c(0, 700)) + 
  theme_bw() + scale_fill_discrete(name = "Group", labels = c("Non-delirium", "Delirium")) +
  labs(x = "Postoperative Timepoint", y = "Biomarker Concentration (pg/ml)") +
  scale_x_discrete(labels = c("IL.10_A" = "Baseline", "IL6_timeB" = "aICU", "IL6_timeC" = "POD1", "IL6_timeD" = "POD4-7")) +
  ggtitle("Postoperative ST2 Concentrations between Groups") + stat_compare_means(aes(group = DSM), label = "p.signif", hide.ns = T) 

il17ftest <- na.omit(data1$IL.17F_D)


?ggtitle
?stat_compare_means()
+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

?element_line()
#tgfb1
tgfb1_plot <- ggplot(mmTGFb1, aes(x=variable, y=value, fill=factor(DSM))) +
  geom_boxplot() + scale_y_continuous(limits = c(0, 20000)) + stat_compare_means(aes(group = DSM), label = "p.signif")

compare_means(value ~ variable, data = mmTGFb1, 
              group.by = "DSM")





ToothGrowth
ToothGrowth$dose <- as.factor(ToothGrowth$dose)



