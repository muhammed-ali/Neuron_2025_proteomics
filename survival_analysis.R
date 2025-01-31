rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(2022)

library(dplyr)
library(data.table)
library(survival)
library(survminer)

load("/Input_Data/Survival_Analysis_COvsAD_Lasso10_IterativeModel.Rdata")

sfit <- survfit(Surv(Years, censored_status)~Biomarker, data=survival_df)
sfit
# Call: survfit(formula = Surv(Years, censored_status) ~ Biomarker, data = survival_df)
#
#               n events median 0.95LCL 0.95UCL
# Biomarker=0 511     73     12      10      NA
# Biomarker=1 285    241      4       4       5

# Kaplan-Meier plot for Biomarker
sfit_biomarker <- survfit(Surv(Years, censored_status)~Biomarker, data=survival_df)
plot(sfit_biomarker) # plain
ggsurvplot(sfit_biomarker) # colored
surv_pvalue(sfit_biomarker)
#    variable         pval   method   pval.txt
# 1 Biomarker 4.492006e-53 Log-rank p < 0.0001
png("/Figures/Survival_Analysis_COvsAD_Lasso10_IterativeModel_Prediction.png", units="mm", width=250, height=200, res=1000)
ggsurvplot(sfit_biomarker, conf.int=TRUE, pval=FALSE, risk.table=FALSE, combine = TRUE,
           fontsize = 6, xlab = " Time (Years)", ylab = c("Probability of not developing AD"),
           font.x = c(25), font.y = c(25), font.legend = c(25), font.tickslab = c(25),
           legend.labs=c("Negative", "Positive"), legend.title="Proteomic signature",  
           palette=c("green4", "red"))
dev.off()

# Kaplan-Meier plot for AT Status
survival_df$ATbin <- ifelse(survival_df$ATbin == "A-T-", 0, 1)
sfit_AT <- survfit(Surv(Years, censored_status)~ATbin, data=survival_df)
plot(sfit_AT) # plain
ggsurvplot(sfit_AT) # colored
surv_pvalue(sfit_AT)
#   variable         pval   method   pval.txt
# 1    ATbin 3.991115e-52 Log-rank p < 0.0001
png("/Figures/Survival_Analysis_COvsAD_Lasso10_IterativeModel_ATnPred.png", units="mm", width=250, height=200, res=1000)
ggsurvplot(sfit_AT, conf.int=TRUE, pval=TRUE, risk.table=FALSE, combine = TRUE,  
           fontsize = 6, xlab = " Time (Years)", ylab = c("Probability of not developing AD"),
           font.x = c(25), font.y = c(25), font.legend = c(25), font.tickslab = c(25),
           legend.labs=c("Negative", "Positive"), legend.title="AT status",  
           palette=c("green4", "red"))
dev.off()


# Effect of Sex
sex_df <- with(survival_df, data.frame(gender = c("F", "M"), age_at_csf_draw = rep(mean(age_at_csf_draw, na.rm = TRUE), 2), Biomarker = c(0, 0)))
fit <- survfit(res.cox, newdata = sex_df)
png("/Figures/Survival_Analysis_COvsAD_Lasso10_IterativeModel_Prediction_Sex.png", units="mm", width=250, height=200, res=1000)
ggsurvplot(fit, data = survival_df, conf.int = TRUE, legend.title="Sex", legend.labs=c("Female", "Male"), palette=c("green4", "red"), risk.table=FALSE,
  ggtheme = theme_minimal(), fontsize = 6, xlab = " Time (Years)", ylab = c("Probability of not developing AD"), font.x = c(25), font.y = c(25), font.legend = c(25), font.tickslab = c(25))
dev.off()

# Effect of Age
age_df <- with(survival_df, data.frame(gender = rep("F", 5), age_at_csf_draw = quantile(age_at_csf_draw), Biomarker = rep(0, 5)))
row.names(age_df) <- age_df$age_at_csf_draw
fit <- survfit(res.cox, newdata = age_df)
png("/Figures/Survival_Analysis_COvsAD_Lasso10_IterativeModel_Prediction_Age.png", units="mm", width=250, height=200, res=1000)
ggsurvplot(fit, data = survival_df, conf.int = TRUE, legend.title="Age", risk.table=FALSE, ggtheme = theme_minimal(),
  fontsize = 6, xlab = " Time (Years)", ylab = c("Probability of not developing AD"), font.x = c(25), font.y = c(25), font.legend = c(25), font.tickslab = c(25))
dev.off()