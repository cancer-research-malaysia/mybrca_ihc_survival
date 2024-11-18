library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

df <- read.csv("../IntratumoralDataset.csv", header = T, sep = ",")

df$time_int_yrs<-df$Dtime_diag/365.25

cbbPalette <- c("#000000", "#E69F00")

hr.df <- filter(df, HR_Status == "P")


##Cox proportional hazad ratio, adjusted by Grade1, Stage1 and Hormone therapy
##--------------------------------------------------------------------------------

##-------------------HR+ patients, CD3

cd3hr_fit <- coxph(Surv(time = hr.df$time_int_yrs, event = hr.df$Death1) ~ CD3_Status
                     +Grade1+Stage1+HormoneTherapy,
                     data=hr.df)
cd3hr_fit_os <- ggadjustedcurves(cd3hr_fit, data=hr.df, method = "conditional", variable = "CD3_Status", palette=cbbPalette, ylim=c(0.2,1.0), legend = c(0.8, 0.3),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_cd3 <- table(hr.df$CD3_Status)
cd3hr_fit_os <- cd3hr_fit_os + guides(color = guide_legend(title = "CD3 Status"))
cd3hr_fit_os <- cd3hr_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_cd3)," (", sample_sizes_cd3, ")"),
                                                      values = cbbPalette)
cd3hr_fit_fp <- ggforest(cd3hr_fit, data=hr.df)

# print(cd3hr_fit_os)
# print(cd3hr_fit_fp)



##-------------------HR+ patients, CD4

cd4hr_fit <- coxph(Surv(time = hr.df$time_int_yrs, event = hr.df$Death1) ~ CD4_Status
                     +Grade1+Stage1+HormoneTherapy,
                     data=hr.df)
cd4hr_fit_os <- ggadjustedcurves(cd4hr_fit, data=hr.df, method = "conditional", variable = "CD4_Status", palette=cbbPalette, ylim=c(0.2,1.0), legend = c(0.8, 0.3),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_cd4 <- table(hr.df$CD4_Status)
cd4hr_fit_os <- cd4hr_fit_os + guides(color = guide_legend(title = "CD4 Status"))
cd4hr_fit_os <- cd4hr_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_cd4)," (", sample_sizes_cd4, ")"),
                                                      values = cbbPalette)
cd4hr_fit_fp <- ggforest(cd4hr_fit, data=hr.df)

# print(cd4hr_fit_os)
# print(cd4hr_fit_fp)




##-------------------HR+ patients, CD8

cd8hr_fit <- coxph(Surv(time = hr.df$time_int_yrs, event = hr.df$Death1) ~ CD8_Status
                     +Grade1+Stage1+HormoneTherapy,
                     data=hr.df)
cd8hr_fit_os <- ggadjustedcurves(cd8hr_fit, data=hr.df, method = "conditional", variable = "CD8_Status", palette=cbbPalette, ylim=c(0.2,1.0), legend = c(0.8, 0.3),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_cd8 <- table(hr.df$CD8_Status)
cd8hr_fit_os <- cd8hr_fit_os + guides(color = guide_legend(title = "CD8 Status"))
cd8hr_fit_os <- cd8hr_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_cd8)," (", sample_sizes_cd8, ")"),
                                                      values = cbbPalette)
cd8hr_fit_fp <- ggforest(cd8hr_fit, data=hr.df)

# print(cd8hr_fit_os)
# print(cd8hr_fit_fp)



##-------------------HR+ patients, PD-L1


cbbPalette1 <- c("#E69F00","#000000")
pdl1hr_fit <- coxph(Surv(time = hr.df$time_int_yrs, event = hr.df$Death1) ~ PDL1
                      +Grade1+Stage1+HormoneTherapy,
                      data=hr.df)
pdl1hr_fit_os <- ggadjustedcurves(pdl1hr_fit, data=hr.df, method = "conditional", variable = "PDL1", palette=cbbPalette1, ylim=c(0.2,1.0), legend = c(0.8, 0.4),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_pdl1 <- table(hr.df$PDL1)
pdl1hr_fit_os <- pdl1hr_fit_os + guides(color = guide_legend(title = "PDL1 Expression"))
pdl1hr_fit_os <- pdl1hr_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_pdl1)," (", sample_sizes_pdl1, ")"),
                                                        values = cbbPalette1)
pdl1hr_fp <- ggforest(pdl1hr_fit, data=hr.df)

# print(pdl1hr_fit_os)
# print(pdl1hr_fp)
