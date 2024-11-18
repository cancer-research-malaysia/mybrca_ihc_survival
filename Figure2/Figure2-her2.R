library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

df <- read.csv("../IntratumoralDataset.csv", header = T, sep = ",")

df$time_int_yrs<-df$Dtime_diag/365.25

cbbPalette <- c("#000000", "#E69F00")

her2.df <- filter(df, HER2_Status == "P")


##Cox proportional hazad ratio, adjusted by Grade1, Stage1 and Hormone therapy
##--------------------------------------------------------------------------------

##-------------------HER2 patients, CD3

cd3her2_fit <- coxph(Surv(time = her2.df$time_int_yrs, event = her2.df$Death1) ~ CD3_Status
                     +Grade1+Stage1+HormoneTherapy,
                     data=her2.df)
cd3her2_fit_os <- ggadjustedcurves(cd3her2_fit, data=her2.df, method = "conditional", variable = "CD3_Status", palette=cbbPalette, ylim=c(0.2,1.0), legend = c(0.8, 0.3),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_cd3 <- table(her2.df$CD3_Status)
cd3her2_fit_os <- cd3her2_fit_os + guides(color = guide_legend(title = "CD3 Status"))
cd3her2_fit_os <- cd3her2_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_cd3)," (", sample_sizes_cd3, ")"),
                                                      values = cbbPalette)
cd3her2_fit_fp <- ggforest(cd3her2_fit, data=her2.df)

print(cd3her2_fit_os)
print(cd3her2_fit_fp)




##-------------------HER2 patients, CD4

cd4her2_fit <- coxph(Surv(time = her2.df$time_int_yrs, event = her2.df$Death1) ~ CD4_Status
                     +Grade1+Stage1+HormoneTherapy,
                     data=her2.df)
cd4her2_fit_os <- ggadjustedcurves(cd4her2_fit, data=her2.df, method = "conditional", variable = "CD4_Status", palette=cbbPalette, ylim=c(0.2,1.0), legend = c(0.8, 0.3),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_cd4 <- table(her2.df$CD4_Status)
cd4her2_fit_os <- cd4her2_fit_os + guides(color = guide_legend(title = "CD4 Status"))
cd4her2_fit_os <- cd4her2_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_cd4)," (", sample_sizes_cd4, ")"),
                                                      values = cbbPalette)
cd4her2_fit_fp <- ggforest(cd4her2_fit, data=her2.df)

print(cd4her2_fit_os)
print(cd4her2_fit_fp)




##-------------------HER2 patients, CD8

cd8her2_fit <- coxph(Surv(time = her2.df$time_int_yrs, event = her2.df$Death1) ~ CD8_Status
                     +Grade1+Stage1+HormoneTherapy,
                     data=her2.df)
cd8her2_fit_os <- ggadjustedcurves(cd8her2_fit, data=her2.df, method = "conditional", variable = "CD8_Status", palette=cbbPalette, ylim=c(0.2,1.0), legend = c(0.8, 0.3),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_cd8 <- table(her2.df$CD8_Status)
cd8her2_fit_os <- cd8her2_fit_os + guides(color = guide_legend(title = "CD8 Status"))
cd8her2_fit_os <- cd8her2_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_cd8)," (", sample_sizes_cd8, ")"),
                                                      values = cbbPalette)
cd8her2_fit_fp <- ggforest(cd8her2_fit, data=her2.df)

print(cd8her2_fit_os)
print(cd8her2_fit_fp)




##-------------------HER2 patients, PD-L1

cbbPalette1 <- c("#E69F00","#000000")
pdl1her2_fit <- coxph(Surv(time = her2.df$time_int_yrs, event = her2.df$Death1) ~ PDL1
                      +Grade1+Stage1+HormoneTherapy,
                      data=her2.df)
pdl1her2_fit_os <- ggadjustedcurves(pdl1her2_fit, data=her2.df, method = "conditional", variable = "PDL1", palette=cbbPalette1, ylim=c(0.2,1.0), legend = c(0.8, 0.4),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_pdl1 <- table(her2.df$PDL1)
pdl1her2_fit_os <- pdl1her2_fit_os + guides(color = guide_legend(title = "PDL1 Expression"))
pdl1her2_fit_os <- pdl1her2_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_pdl1)," (", sample_sizes_pdl1, ")"),
                                                        values = cbbPalette1)
pdl1her2_fp <- ggforest(pdl1her2_fit, data=her2.df)

print(pdl1her2_fit_os)
print(pdl1her2_fp)


