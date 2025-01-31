rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(2022)

library(caret)
library(pROC)
library(dplyr)
library(tibble)
library(tidyr)

##
## Feature filtering (Remove corr > 0.8 features and those not present in Stanford ADRC) 
##

# read proteomic data for Discovery (Stage 1)
load("/Input_Data/Discovery_MAP_Ruiz_ADvsCO_log10_Proteomic_Data_7029_DAA.RData")
dim(discovery_data_all) # 1170 7035

# filter for analytes that passed bonferroni correction in meta-analysis
ADCO_significant <- read.table("/Input_Data/Discovery_ReplicationOfDiscoveryFDR_ATN_ageSexPlate_lm_sv_Metaanalysis_Bonf_Filtered.txt", sep="\t", header=T, stringsAsFactors=F, quote="", fill=T)
dim(ADCO_significant) # 2173   20
selected_cols <- c("age_at_csf_draw", "gender", "PlateId", "level", "sv1", "sv2", ADCO_significant$Analyte)
discovery_data_all <- discovery_data_all[,selected_cols]
dim(discovery_data_all) # 1170 2179

# filter for analytes not present in the Stanford ADRC
load("/Input_Data/proteomics_stanford_107sample_4734analytes.RData")
dim(proteomics_stanford_final_sv) # 107 4741
analytes_discovery <- unique(colnames(discovery_data_all[7:2179]))
length(analytes_discovery) # 2173
analytes_stanford <- unique(colnames(proteomics_stanford_final_sv[8:4741]))
length(analytes_stanford) # 4734
table(analytes_discovery %in% analytes_stanford)
# FALSE  TRUE
#   721  1452
overlapping_analytes <- intersect(analytes_discovery, analytes_stanford)
length(overlapping_analytes) # 1452
selected_cols <- c("age_at_csf_draw", "gender", "PlateId", "level", "sv1", "sv2", overlapping_analytes)
discovery_stanford_filt <- discovery_data_all[,selected_cols]
dim(discovery_stanford_filt) # 1170 1458

# remove highly correlated analytes (cor > 0.80)
d2 <- discovery_stanford_filt
d2$sex <- ifelse(d2$gender == "F", 0, 1)
d2$gender <- NULL
d2$PlateId <- NULL
d2$sv1 <- NULL
d2$sv2 <- NULL
dim(d2) # 1170 1455
d2_cor <- d2 %>% as.matrix %>% cor %>% as.data.frame %>% rownames_to_column(var = 'var1') %>% gather(var2, value, -var1)
dim(d2_cor) # 2117025       3
d2_cor <- na.omit(d2_cor)
dim(d2_cor) # 2114117        3
d2_cor <- d2_cor[order(d2_cor$value),]
d2_cor <- d2_cor[!(d2_cor$var1 == d2_cor$var2),] # remove corr b/w same variables
tail(d2_cor)
dim(d2_cor) # 2112662        3
nrow(d2_cor[d2_cor$value > 0.80,]) # 9826 different pairs of analytes for which cor > 0.80.

# remove one of the correlated analytes (cor > 0.80) from the analyte list
correlated_Analytes <- d2_cor[d2_cor$value > 0.80,]
dim(correlated_Analytes) # 9826    3
length(unique(correlated_Analytes$var1)) # 424
length(unique(correlated_Analytes$var2)) # 424
table(correlated_Analytes$var1 %in% correlated_Analytes$var2)
#  TRUE
#  9826
# remove corr b/w mirror combinations, keep one.
correlated_Analytes$temp <- apply(correlated_Analytes, 1, function(x) paste(sort(x), collapse=""))
dim(correlated_Analytes) # 9826     4
correlated_Analytes <- correlated_Analytes[!duplicated(correlated_Analytes$temp), 1:3]
dim(correlated_Analytes) # 4913    3
length(unique(correlated_Analytes$var1)) # 361
length(unique(correlated_Analytes$var2)) # 338
table(correlated_Analytes$var1 %in% correlated_Analytes$var2)
# FALSE  TRUE
#   282  4631
analytes_noncor <- unique(overlapping_analytes[!(overlapping_analytes %in% correlated_Analytes$var1)])
length(analytes_noncor) # 1094
1452 - 1094 # 358 correlated (corr > 0.8) features are removed

# create expression matrix for non-corr analytes that are also present in Stanford ADRC
selected_cols <- c("age_at_csf_draw", "gender", "PlateId", "level", "sv1", "sv2", analytes_noncor)
discovery_stanford_filt_noncor <- discovery_stanford_filt[,selected_cols]
dim(discovery_stanford_filt_noncor) # 1170 1100
analyte_mapping <- ADCO_significant[,c("Analyte", "EntrezGeneSymbol.x")]
analyte_mapping$Analyte_Symbol <- paste0(analyte_mapping$Analyte, "_", analyte_mapping$EntrezGeneSymbol.x, sep="")
head(analyte_mapping, 2)
dim(analyte_mapping) # 2173    3
table(names(discovery_stanford_filt_noncor)[7:1100] %in% analyte_mapping$Analyte)
# TRUE
# 1094
analyte_mapping <- analyte_mapping[analyte_mapping$Analyte %in% names(discovery_stanford_filt_noncor)[7:1100],]
dim(analyte_mapping) # 1094    3
names(discovery_stanford_filt_noncor)[7:1100] <- analyte_mapping$Analyte_Symbol[match(names(discovery_stanford_filt_noncor)[7:1100],analyte_mapping$Analyte)]


##
## LASSO Regression (50 iterations) to create AD prediction model
##

# define a "not in" opperator
`%!in%` <- Negate(`%in%`)
  
# expression matrix
quantMatrix <- as.matrix(discovery_stanford_filt_noncor[,7:1100])
dim(quantMatrix) # 1170 1094

# phenotype
pheno <- discovery_stanford_filt_noncor[,1:6]
pheno$sample_ID <- row.names(pheno)
dim(pheno) # 1170    7

# list the proteins
proteins <- colnames(quantMatrix)
length(proteins) # 1094
prot <- as.data.frame(proteins)
names(prot) <- "Protein"
df <- as.data.frame(prot)
colnames(df) <- c("Proteins")
dim(df) # 1094    1
aucdf <- data.frame()

for (i in 1:50)
  {
    set.seed(i)
    # set up training (70%) and testing (30%) datasets
    trainSamples <- c(sample(pheno$sample_ID[pheno$level=='1'], round(length(pheno$level[pheno$level=='1'])*0.7)),
                    sample(pheno$sample_ID[pheno$level=='0'], round(length(pheno$level[pheno$level=='0'])*0.7)))
    trainX <- quantMatrix[rownames(quantMatrix)%in%trainSamples,]
    trainX <- as.matrix(trainX)
    trainY <- as.numeric(pheno$level[pheno$sample_ID %in% trainSamples])
    # build a lasso model
    set.seed(567)
    cvOut <- cv.glmnet(trainX,trainY,family="binomial",alpha=1)
    lamMin <- cvOut$lambda.min
    # rebuild the model using best lambda
    lassoBest <- glmnet(trainX,trainY,family="binomial",alpha=1, lambda = lamMin)
    coefProt <- coef(lassoBest) %>% as.matrix() %>% as_tibble(rownames = "Proteins")
    sigProt <- coefProt[coefProt$s0 != 0,]
    df <- merge(df,sigProt, by="Proteins", all.x = T)
    testX <- quantMatrix[rownames(quantMatrix)%!in%trainSamples,]
    testY <- as.numeric(pheno$level[pheno$sample_ID %!in% trainSamples])
    preds <- as.data.frame(predict(lassoBest, newx = testX, type = "response"))
    t <- ci.auc(testY, as.numeric(preds$s0),conf.level = 0.9)
    t <- data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Proteins = (nrow(sigProt)-1))
    aucdf <- rbind(aucdf,t)
    print(i)
}

# save predictions and performance obtained from lasso model training
dim(df) # 1094   51
dim(aucdf) # 50  4
write.csv(aucdf, "/Input_Data/ModelPerformance_50runs.csv", row.names = F)
col_names <- c("Proteins", seq(from=1,to=50))
colnames(df) <- col_names
df$count <- rowSums(!is.na(select(df, 2:51)))
dim(df) # 1094   52
df <- df[order(df$count, decreasing=T),]
write.csv(df,"/Input_Data/ProteinWeights_Sorted_50runs.csv", row.names = F)
df$Analyte <- sapply(strsplit(as.character(df$Proteins), "_"), '[', 1)

# save the final AD prediction model containing top 10 proteins appearing in the most Lasso models
iterativeModel_10analytes <- df[1:10,]$Analyte
save(iterativeModel_10analytes, file="/Input_Data/iterativeModel_10analytes.RData")


##
## A+T+ vs. A-T-
##

# Discovery (Knight-ADRC)

load("/Input_Data/Discovery_MAP_Ruiz_ADvsCO_log10_Proteomic_Data_7029_DAA.RData")
selected_cols <- c("age_at_csf_draw", "gender", "level", iterativeModel_10analytes)
disc_auc_data <- as.data.frame(discovery_data_all)[, selected_cols]
dim(disc_auc_data) # 1170   13
table(disc_auc_data$level) # 0 = A-T-; 1 = A+T+
#   0   1
# 680 490

# Creating a 70/30 split for train/test in the training dataset
index = sample(1:nrow(disc_auc_data), size = 0.70 * nrow(disc_auc_data))
train_discovery = disc_auc_data[index, ] # 819  13
test_discovery = disc_auc_data[-index, ] # 351  13

# for Baseline (without any analyte)
bl_disc <- glm(as.formula(level ~ age_at_csf_draw + gender), data = train_discovery, family = 'binomial')
coefficients_baseline <- as.data.frame(summary(bl_disc)[13])
out.prob <- predict(bl_disc, newdata = test_discovery, type='response')
auc_discovery_washU_25_bl <- roc(as.formula(test_discovery$level ~ out.prob), plot = FALSE, print.auc = FALSE)
npv_ppv_best <- coords(auc_discovery_washU_25_bl, x = "best", best.method="youden", input = "threshold", ret = c("threshold","acc","npv","ppv","sensitivity","specificity"))

# for all analytes (N=10)
all_disc <- glm(level ~ ., data = train_discovery, family = 'binomial')
coefficients_all_disc <- as.data.frame(summary(all_disc)[13])
out_prob_stage1 <- predict(all_disc, newdata = test_discovery, type='response')
auc_discovery_washU_25_analytes <- roc(as.formula(test_discovery$level ~ out_prob_stage1), plot = FALSE, print.auc = FALSE)
npv_ppv_best <- coords(auc_discovery_washU_25_analytes, x = "best", best.method="youden", input = "threshold", ret = c("threshold","acc","npv","ppv","sensitivity","specificity"))
fixed_threshold <- npv_ppv_best$threshold # 0.7060383
save(all_disc, bl_disc, train_discovery, test_discovery, iterativeModel_10analytes, fixed_threshold, file="/Input_Data/ADvsCO_10Analytes_Iterative_Model.RData")

# Replication (Knight-ADRC)

load("/Input_Data/Replication_ADNI_Pau_ADvsCO_log10_Proteomic_Data_7029_DAA.RData")
length(selected_cols) # 13
repl_auc_data <- as.data.frame(replication_data_all)[, selected_cols]
dim(repl_auc_data) # 593  13
table(repl_auc_data$level)
#   0   1
# 235 358

# for Baseline (without any analyte)
out.prob <- predict(bl_disc, newdata = repl_auc_data, type='response')
auc_replication_washU_25_bl <- roc(as.formula(repl_auc_data$level ~ out.prob), plot = FALSE, print.auc = FALSE)
npv_ppv_best <- coords(auc_replication_washU_25_bl, x = "best", best.method="youden", input = "threshold", ret = c("threshold","acc","npv","ppv","sensitivity","specificity"))
round(auc_replication_washU_25_bl$auc, 2) # 0.59

# for all analytes (N=10)
out_prob_stage2 <- predict(all_disc, newdata = repl_auc_data, type='response')
auc_replication_washU_25_analytes <- roc(as.formula(repl_auc_data$level ~ out_prob_stage2), plot = FALSE, print.auc = FALSE)
npv_ppv_best <- coords(auc_replication_washU_25_analytes, x = "all", best.method="youden", input = "threshold", ret = c("threshold","acc","npv","ppv","sensitivity","specificity"))
round(auc_replication_washU_25_analytes$auc, 2) # 0.98

# Replication (Stanford-ADRC)

load("/Input_Data/proteomics_stanford_138sample_4734analytes.RData")
dim(proteomics_stanford_final) # 138 4742
table(proteomics_stanford_final$at_status)
# 00 01 10 11
# 80  6 25 27
proteomics_stanford_final <- proteomics_stanford_final[proteomics_stanford_final$at_status %in% c("11", "00"),]
table(proteomics_stanford_final$at_status)
# 00 11
# 80 27
proteomics_stanford_final$level <- ifelse(proteomics_stanford_final$at_status == "00", 0, 1)
replication_data_all <- proteomics_stanford_final
length(selected_cols) # 13
repl_auc_data <- as.data.frame(replication_data_all)[, selected_cols]
dim(repl_auc_data) # 107  13
table(repl_auc_data$level)
#  0  1
# 80 27

# For Baseline (without any analyte)
out.prob <- predict(bl_disc, newdata = repl_auc_data, type='response')
auc_replication_stanford_25_bl <- roc(as.formula(repl_auc_data$level ~ out.prob), plot = FALSE, print.auc = FALSE)
npv_ppv_best <- coords(auc_replication_stanford_25_bl, x = "best", best.method="youden", input = "threshold", ret = c("threshold","acc","npv","ppv","sensitivity","specificity"))
round(auc_replication_stanford_25_bl$auc, 2) # 0.57

# For all analytes
out_prob_stage3 <- predict(all_disc, newdata = repl_auc_data, type='response')
auc_replication_stanford_25_analytes <- roc(as.formula(repl_auc_data$level ~ out_prob_stage3), plot = FALSE, print.auc = FALSE)
npv_ppv_best <- coords(auc_replication_stanford_25_analytes, x = "all", best.method="youden", input = "threshold", ret = c("threshold","acc","npv","ppv","sensitivity","specificity"))
round(auc_replication_stanford_25_analytes$auc, 2) # 0.99

# ROC plot
png("/Figures/AUC_DiscRepl_COvsAD_Lasso_10_IterativeModel.png", units="mm", width=120, height=120, res=1000)
plot_curve <- plot.roc(auc_discovery_washU_25_bl, xlim=c(1,0), ylim=c(0, 1), col=alpha("blue", 0.5), lty=2, main = c("A-T- vs A+T+"))
text(0.6, 0.20, paste("Stage 1 Baseline: ", round(auc_discovery_washU_25_bl$auc, 2), sep=""), col=alpha("blue", 0.5), cex=0.75, pos=4)
plot_curve <- plot.roc(auc_discovery_washU_25_analytes, xlim=c(1,0), ylim=c(0, 1), lty=1, lwd=6, col=alpha("blue", 0.5), add=TRUE)
text(0.6, 0.25, paste("Stage 1 (30%; Testing): ",  round(auc_discovery_washU_25_analytes$auc, 2), sep=""), col=("blue"), cex=0.85, pos=4)
plot_curve <- plot.roc(auc_replication_washU_25_bl, xlim=c(1,0), ylim=c(0, 1), col=alpha("red", 0.5), lty=2, add=TRUE)
text(0.6, 0.10, paste("Stage 2 Baseline: ", round(auc_replication_washU_25_bl$auc, 2), sep=""), col=alpha("red", 0.5), cex=0.75, pos=4)
plot_curve <- plot.roc(auc_replication_washU_25_analytes, xlim=c(1,0), ylim=c(0, 1), lty=1, lwd=6, col=alpha("red", 0.5), add=TRUE)
text(0.6, 0.15, paste("Stage 2 Testing: ",  round(auc_replication_washU_25_analytes$auc, 2), sep=""), col=("red"), cex=0.85, pos=4)
plot_curve <- plot.roc(auc_replication_stanford_25_bl, xlim=c(1,0), ylim=c(0, 1), col=alpha("green4", 0.5),lty=2, add=TRUE)
text(0.6, 0.0, paste("Stanford Baseline: ", round(auc_replication_stanford_25_bl$auc, 2), sep=""), col=alpha("green4",0.5), cex=0.75, pos=4)
plot_curve <- plot.roc(auc_replication_stanford_25_analytes, xlim=c(1,0), ylim=c(0, 1), lty=1, lwd=6, col=alpha("green4", 0.5), add=TRUE)
text(0.6, 0.05, paste("Stanford Validation: ",  round(auc_replication_stanford_25_analytes$auc, 2), sep=""), col=("green4"), cex=0.80, pos=4)
dev.off()

save(auc_discovery_washU_25_bl, auc_discovery_washU_25_analytes, auc_replication_washU_25_bl, auc_replication_washU_25_analytes, auc_replication_stanford_25_bl, 
	auc_replication_stanford_25_analytes, file="/Input_Data/AUC_ATpos_ATneg/AUC_DiscRepl_COvsAD_Lasso_10_IterativeModel.RData")