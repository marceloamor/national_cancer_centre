##########Forest Plot!!!

#install.packages("devtools")
#install.packages("ps")
#devtools::install_github("NightingaleHealth/ggforestplot")
#options(download.file.method = "libcurl")

ORvalues <- read.csv("ORvalues.csv",header=T,row.names=1)
TPAOR <- read.csv("OR_A.csv",header=T,row.names=1) #visuospatial was applied as 1 for OR/CI
TPBOR <- read.csv("OR_B.csv",header=T,row.names=1)
TPCOR <- read.csv("OR_C.csv",header=T,row.names=1)
TPDOR <- read.csv("OR_D.csv",header=T,row.names=1)

?forestplot

forestplot(rownames(TPAOR), 
           txt_gp = fpTxtGp(label = list(gpar(fontfamily = "", cex=0.7)),
                            ticks = gpar(fontfamily = "", cex=0.7),
                            xlab  = gpar(fontfamily = "", cex = 0.7)),
           legend = c("Baseline", "aICU", "POD1"),
           fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI, fpDrawDiamondCI),
           boxsize = .25, # We set the box size to better visualize the type
           line.margin = .02, # We need to add this to avoid crowding
           mean = cbind(TPAOR[,1], TPBOR[,1], TPCOR[,1]),
           lower = cbind(TPAOR[,2], TPBOR[,2], TPCOR[,2]),
           upper = cbind(TPAOR[,3], TPBOR[,3], TPCOR[,3]),
           clip =c(0, max(c(TPAOR[,3], TPBOR[,3], TPCOR[,3], na.rm = TRUE))),
           col=fpColors(box=c("yellow", "red", "orange")),
           grid = structure(c(1), 
                            gp = gpar(lty = 2, col = "#CCCCFF")), 
           xlab="Odds Ratio and 95% CI",)











colnames(df_logodds)

df_logodds <-
  ggforestplot::df_logodds_associations %>%
  dplyr::arrange(name) %>%
  dplyr::left_join(ggforestplot::df_NG_biomarker_metadata, by = "name") %>% 
  dplyr::filter(group == "Amino acids") %>%
  # Set the study variable to a factor to preserve order of appearance
  # Set class to factor to set order of display.
  dplyr::mutate(
    study = factor(
      study,
      levels = c("Meta-analysis", "NFBC-1997", "DILGOM", "FINRISK-1997", "YFS")
    )
  )

ggforestplot::forestplot(
  df = df_logodds,
  estimate = beta,
  logodds = T,
  colour = study,
  shape = study,
  title = "Associations to type 2 diabetes",
  xlab = "Odds ratio for incident type 2 diabetes (95% CI)
  per 1???SD increment in metabolite concentration"
) +
  # You may also want to add a manual shape scale to mark meta-analysis with a
  # diamond shape
  ggplot2::scale_shape_manual(
    values = c(23L, 21L, 21L, 21L, 21L),
    labels = c("Meta-analysis", "NFBC-1997", "DILGOM", "FINRISK-1997", "YFS")
  )

ggforestplot::forestplot(
  df = ORvalues,
  name = biomarker,
  estimate = OR,
  se = log_SE,
  logodds = T,
  colour = timepoint,
  shape = timepoint,
  title = "forestplot that shows marce is smart",
  xlab = "Odds ratio for incident type 2 diabetes (95% CI)
  per 1???SD increment in metabolite concentration"
) +
  # You may also want to add a manual shape scale to mark meta-analysis with a
  # diamond shape
  ggplot2::scale_shape_manual(
    values = c(23L, 21L, 21L, 21L, 21L),
    labels = c("Meta-analysis", "NFBC-1997", "DILGOM", "FINRISK-1997", "YFS")
  )



std.error(data1$IL6_timeA,na.rm = T)

sd(data1$IL6_timeB, na.rm = T)/sqrt(sum(!is.na(data1$IL6_timeB)))



