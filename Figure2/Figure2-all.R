library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

df <- read.csv("../IntratumoralDataset.csv", header = T, sep = ",")

df$time_int_yrs<-df$Dtime_diag/365.25

cbbPalette <- c("#000000", "#E69F00")


##Cox proportional hazad ratio, adjusted by HR status, HER2 status, Grade1, Stage1
##--------------------------------------------------------------------------------

##-------------------All patients, CD3


cd3status_fit <- coxph(Surv(time = df$time_int_yrs, event = df$Death1) ~ CD3_Status
                       +HR_Status+HER2_Status+Grade1+Stage1,
                       data=df)
cd3status_fit_os <- ggadjustedcurves(cd3status_fit, data=df, method = "conditional", variable = "CD3_Status", palette=cbbPalette, ylim=c(0.2,1.0), legend = c(0.8, 0.3),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_cd3 <- table(df$CD3_Status)
cd3status_fit_os <- cd3status_fit_os + guides(color = guide_legend(title = "CD3 Status"))
cd3status_fit_os <- cd3status_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_cd3)," (", sample_sizes_cd3, ")"),
                                                                    values = cbbPalette)
cd3status_fit_fp <- ggforest(cd3status_fit, data=df)

print(cd3status_fit_os)
print(cd3status_fit_fp)




##-------------------All patients, CD4


cd4status_fit <- coxph(Surv(time = df$time_int_yrs, event = df$Death1) ~ CD4_Status
                       +HR_Status+HER2_Status+Grade1+Stage1,
                       data=df)
cd4status_fit_os <- ggadjustedcurves(cd4status_fit, data=df, method = "conditional", variable = "CD4_Status", palette=cbbPalette, ylim=c(0.2,1.0), legend = c(0.8, 0.3),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_cd4 <- table(df$CD4_Status)
cd4status_fit_os <- cd4status_fit_os + guides(color = guide_legend(title = "CD4 Status"))
cd4status_fit_os <- cd4status_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_cd4)," (", sample_sizes_cd4, ")"),
                                                          values = cbbPalette)
cd4status_fit_fp <- ggforest(cd4status_fit, data=df)

print(cd4status_fit_os)
print(cd4status_fit_fp)




##-------------------All patients, CD8


cd8status_fit <- coxph(Surv(time = df$time_int_yrs, event = df$Death1) ~ CD8_Status
                       +HR_Status+HER2_Status+Grade1+Stage1,
                       data=df)
cd8status_fit_os <- ggadjustedcurves(cd8status_fit, data=df, method = "conditional", variable = "CD8_Status", palette=cbbPalette, ylim=c(0.2,1.0), legend = c(0.8, 0.3),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_cd8 <- table(df$CD8_Status)
cd8status_fit_os <- cd8status_fit_os + guides(color = guide_legend(title = "CD8 Status"))
cd8status_fit_os <- cd8status_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_cd8)," (", sample_sizes_cd8, ")"),
                                                          values = cbbPalette)
cd8status_fit_fp <- ggforest(cd8status_fit, data = df)

print(cd8status_fit_os)
print(cd8status_fit_fp)




##-------------------All patients, PD-L1

cbbPalette1 <- c("#E69F00","#000000")
pdl1_fit <- coxph(Surv(time = df$time_int_yrs, event = df$Death1) ~ PDL1
                  +HR_Status+HER2_Status+Grade1+Stage1,
                  data=df)
pdl1_fit_os <- ggadjustedcurves(pdl1_fit, data=df, method = "conditional", variable = "PDL1", palette=cbbPalette1, ylim=c(0.2,1.0), legend = c(0.8, 0.4),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_pdl1 <- table(df$PDL1)
pdl1_fit_os <- pdl1_fit_os + guides(color = guide_legend(title = "PDL1 Expression"))
pdl1_fit_os <- pdl1_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_pdl1)," (", sample_sizes_pdl1, ")"),
                                                values = cbbPalette1)
PDL1_fp <- ggforest(pdl1_fit, data = df)

print(pdl1_fit_os)
print(PDL1_fp)


