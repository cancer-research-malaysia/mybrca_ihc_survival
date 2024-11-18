library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

df <- read.csv("../IntratumoralDataset.csv", header = T, sep = ",")

df$time_int_yrs<-df$Dtime_diag/365.25

cbbPalette <- c("#000000", "#E69F00")

tnbc.df <- filter(df, clinicalsubtype2 == "HR-/HER2-")

##Cox proportional hazad ratio, adjusted by Grade1, Stage1 and Chemotherapy
##--------------------------------------------------------------------------------

##-------------------TNBC patients, CD3


cd3tnbc_fit <- coxph(Surv(time = tnbc.df$time_int_yrs, event = tnbc.df$Death1) ~ CD3_Status
                       +Grade1+Stage1+Chemotherapy,
                       data=tnbc.df)
cd3tnbc_fit_os <- ggadjustedcurves(cd3tnbc_fit, data=tnbc.df, method = "conditional", variable = "CD3_Status", palette=cbbPalette, ylim=c(0.2,1.0), legend = c(0.8, 0.3),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_cd3 <- table(tnbc.df$CD3_Status)
cd3tnbc_fit_os <- cd3tnbc_fit_os + guides(color = guide_legend(title = "CD3 Status"))
cd3tnbc_fit_os <- cd3tnbc_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_cd3)," (", sample_sizes_cd3, ")"),
                                                          values = cbbPalette)
cd3tnbc_fit_fp <- ggforest(cd3tnbc_fit, data=tnbc.df)

# print(cd3tnbc_fit_os)
# print(cd3tnbc_fit_fp)



##-------------------TNBC patients, CD4


cd4tnbc_fit <- coxph(Surv(time = tnbc.df$time_int_yrs, event = tnbc.df$Death1) ~ CD4_Status
                       +Grade1+Stage1+Chemotherapy,
                       data=tnbc.df)
cd4tnbc_fit_os <- ggadjustedcurves(cd4tnbc_fit, data=tnbc.df, method = "conditional", variable = "CD4_Status", palette=cbbPalette, ylim=c(0.2,1.0), legend = c(0.8, 0.3),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_cd4 <- table(tnbc.df$CD4_Status)
cd4tnbc_fit_os <- cd4tnbc_fit_os + guides(color = guide_legend(title = "CD4 Status"))
cd4tnbc_fit_os <- cd4tnbc_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_cd4)," (", sample_sizes_cd4, ")"),
                                                          values = cbbPalette)
cd4tnbc_fit_fp <- ggforest(cd4tnbc_fit, data=tnbc.df)

# print(cd4tnbc_fit_os)
# print(cd4tnbc_fit_fp)



##-------------------TNBC patients, CD8


cd8tnbc_fit <- coxph(Surv(time = tnbc.df$time_int_yrs, event = tnbc.df$Death1) ~ CD8_Status
                       +Grade1+Stage1+Chemotherapy,
                       data=tnbc.df)
cd8tnbc_fit_os <- ggadjustedcurves(cd8tnbc_fit, data=tnbc.df, method = "conditional", variable = "CD8_Status", palette=cbbPalette, ylim=c(0.2,1.0), legend = c(0.8, 0.3),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_cd8 <- table(tnbc.df$CD8_Status)
cd8tnbc_fit_os <- cd8tnbc_fit_os + guides(color = guide_legend(title = "CD8 Status"))
cd8tnbc_fit_os <- cd8tnbc_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_cd8)," (", sample_sizes_cd8, ")"),
                                                          values = cbbPalette)
cd8tnbc_fit_fp <- ggforest(cd8tnbc_fit, data=tnbc.df)

# print(cd8tnbc_fit_os)
# print(cd8tnbc_fit_fp)



##-------------------TNBC patients, PD-L1

cbbPalette1 <- c("#E69F00","#000000")
pdl1tnbc_fit <- coxph(Surv(time = tnbc.df$time_int_yrs, event = tnbc.df$Death1) ~ PDL1
                  +Grade1+Stage1+Chemotherapy,
                  data=tnbc.df)
pdl1tnbc_fit_os <- ggadjustedcurves(pdl1tnbc_fit, data=tnbc.df, method = "conditional", variable = "PDL1", palette=cbbPalette1, ylim=c(0.2,1.0), legend = c(0.8, 0.4),ylab="Overall Survival", xlab="Time(years)")
sample_sizes_pdl1 <- table(tnbc.df$PDL1)
pdl1tnbc_fit_os <- pdl1tnbc_fit_os + guides(color = guide_legend(title = "PDL1 Expression"))
pdl1tnbc_fit_os <- pdl1tnbc_fit_os + scale_color_manual(labels = paste0(names(sample_sizes_pdl1)," (", sample_sizes_pdl1, ")"),
                                                values = cbbPalette1)
pdl1tnbc_fp <- ggforest(pdl1tnbc_fit, data=tnbc.df)

# print(pdl1tnbc_fit_os)
# print(pdl1tnbc_fp)

