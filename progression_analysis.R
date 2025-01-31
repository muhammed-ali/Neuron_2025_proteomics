rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(2022)

library(sjPlot)
library(nlme)

df <- read.table("/Input_Data/ADNI_MAP_Progression_Table_10_IterativeModel.txt", sep="\t", header=T, stringsAsFactors=F)
dim(df) # 1364   19
df$sex <- as.factor(df$sex)
df$Proteomic_signature <- as.factor(df$Proteomic_signature)
fit <- lme(cdrsum ~ years + age_at_csf_draw + sex + Proteomic_signature + cdrsum_first + years*Proteomic_signature, data=df, random=~1|ID)
sum_fit <- summary(fit)
sum_fit_coef <- coef(sum_fit)
sum_fit_coef
#                                         Value  Std.Error   DF   t-value      p-value
# (Intercept)                       -1.49863143 0.92939392 1029 -1.612483 1.071636e-01
# years                              0.52309467 0.04579327 1029 11.422961 1.549549e-28
# age_at_csf_draw                    0.02231688 0.01264724  328  1.764565 7.856756e-02
# sexMale                           -0.21877521 0.18979031  328 -1.152721 2.498644e-01
# Proteomic_signaturePositive        0.27996944 0.24435376  328  1.145755 2.527319e-01
# cdrsum_first                       1.01824571 0.05889290  328 17.289786 4.588661e-48
# years:Proteomic_signaturePositive  0.41277429 0.05650667 1029  7.304877 5.554228e-13

# plotting progression slopes
png("/Figures/Progression_Analysis_Slopes.png", units="mm", width=120, height=120, res=1000)
plot_model(fit, type = "pred", terms = c("years", "Proteomic_signature"), colors=c('red', 'green4')) + theme(text = element_text(size=20),
  panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour = "black"), 
  legend.title=element_text(size=10), legend.text=element_text(size=10), legend.position="top") + ggtitle("")
dev.off()