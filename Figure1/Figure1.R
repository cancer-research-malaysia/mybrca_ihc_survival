library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(dplyr)
library("RColorBrewer")

df <- read.csv("../IntratumoralDataset.csv", header=T, sep=',')

df$time_int_yrs <- df$Dtime_diag/365.25


##--------------Figure 1a. Clinical subtype

clinicalPalette1 <- c("#3C4CA8","#B2DFE8","#ED1F24","#C651AF")

a.df <- df
a.df$clinicalsubtype2 <- factor(a.df$clinicalsubtype2, levels = c("HR+/HER2-","HR+/HER2+","HR-/HER2-","HR-/HER2+"))
a.df$clinicalsubtype2 = relevel(a.df$clinicalsubtype2, ref = "HR+/HER2-")

#cox proportional hazard model
cfit4a <- coxph(Surv(time = a.df$time_int_yrs, event = a.df$Death1) ~ clinicalsubtype2
                +Grade1+Stage1,
                data=a.df)
clinical_os <- ggadjustedcurves(cfit4a, data=a.df, method = "conditional",
                                variable = "clinicalsubtype2",palette=clinicalPalette1,
                                legend=c(0.3,0.3), ylim=c(0.2,1.0),ylab="Overall survival",
                                xlab="Time(years)")
sample_sizes <- table(a.df$clinicalsubtype2)

clinical_os <- clinical_os + guides(color = guide_legend(title = "Clinical Subtype"))
clinical_os <- clinical_os + scale_color_manual(labels = paste0(names(sample_sizes)," (", sample_sizes, ")"),
                                      values = clinicalPalette1)
time_of_interest <- 5
clinical_os <- clinical_os +
  geom_vline(xintercept = time_of_interest, linetype = "dashed",color = "black")+
  annotate("text", x = 5.5, y = 0.2, label = "5 Years")
print(clinical_os)

clinical_fp <- ggforest(cfit4a, data = a.df)
print(clinical_fp)

ggsave("figure-1a.pdf",clinical_os, width=5.71, height=3.90, units="in")
ggsave("suppfigure-1a.pdf",clinical_fp, width=5.71, height=3.90, units="in")



##--------------Figure 1b. Tumour stage

fig1b.df <- df

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")

fit_s <- coxph(Surv(time = fig1b.df$time_int_yrs, event = fig1b.df$Death1) ~ Stage1
               +Grade1+clinicalsubtype2,
               data=fig1b.df)

stage_os <- ggadjustedcurves(fit_s, data=fig1b.df, method = "conditional",variable = "Stage1",palette=cbbPalette, legend=c(0.2,0.3), ylim=c(0.2,1.0),ylab="Overall survival", xlab="Time(years)")
sample_sizes_b <- table(fig1b.df$Stage1)
stage_os <- stage_os + guides(color = guide_legend(title = "Tumor Stage"))
stage_os <- stage_os + scale_color_manual(labels = paste0(names(sample_sizes_b)," (", sample_sizes_b, ")"),
                                      values = cbbPalette)
stage_os <- stage_os +
  geom_vline(xintercept = time_of_interest, linetype = "dashed",color = "black")+
  annotate("text", x = 5.5, y = 0.2, label = "5 Years")
print(stage_os)

stage_fp <- ggforest(fit_s, fig1b.df)
print(stage_fp)

ggsave("figure-1b.pdf",stage_os, width=5.71, height=3.90, units="in")
# ggsave("suppfigure-1b.pdf",stage_fp, width=5.71, height=3.90, units="in")
