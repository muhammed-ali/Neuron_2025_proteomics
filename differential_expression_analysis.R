rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(2022)

library(a4Base)
library(sva)
library(ggrepel)
library(data.table)
library(ggplot2)
library(Biobase)
library(dplyr)
library(EnhancedVolcano)
library(readxl)
library(scran)


##
## Discovery Linear Regression (Stage 1): A-T- vs. A+T+
##

load("/Input_Data/Discovery_MAP_Ruiz_edata_pdata_withSV.RData")
dim(edata_pdata_disc) # 1545 7039

# replace the level according to ATN and keep the group we are interested in (AD and CO)
edata_pdata_disc$ATN <- ifelse(edata_pdata_disc$ATbin == "A+T+", "AD",  
  ifelse(edata_pdata_disc$ATbin == "A-T-", "CO",
  ifelse(edata_pdata_disc$ATbin == "A+T-", "Asymp", "Exclude")))

# reorder columns
edata_pdata_disc <- edata_pdata_disc[,-c("ATbin")]
dim(edata_pdata_disc) #  1545 7039
edata_pdata_disc <- edata_pdata_disc %>% select(ExtIdentifier, UniquePhenoID, age_at_csf_draw, gender, PlateId, status_for_analysis, ATN, level, sv1, sv2, everything())
dim(edata_pdata_disc) #  1545 7039

# Subset based on CC status from actual pheno so that DLB, PD, and ADAD are not included
edata_pdata_disc <- edata_pdata_disc[edata_pdata_disc$status_for_analysis %in% c("AD", "CO", "OT"),]
dim(edata_pdata_disc) #  1517 7039

# Subset based on ATN status - only keeping AD (A+T+) and CO (A-T-)
edata_pdata_disc <- edata_pdata_disc[edata_pdata_disc$ATN %in% c("AD", "CO"),]
dim(edata_pdata_disc) #  1170 7039

edata_pdata_disc$level <- ifelse(edata_pdata_disc$ATN == 'CO',0,1)
dim(edata_pdata_disc) #  1170 7039
table(edata_pdata_disc$level)
#   0   1
# 680 490
ADCO_disc <- edata_pdata_disc[,c("ExtIdentifier", "UniquePhenoID", "ATN")]
dim(ADCO_disc) # 1170    3

# get analytes
analyte <- as.data.frame(edata_pdata_disc[,11:ncol(edata_pdata_disc)])
dim(analyte) #  1170 7029
analyte_name <- as.data.frame(colnames(analyte))
head(analyte_name)
# replace the NaN and Inf values in our data frame
edata_pdata_disc[is.na(edata_pdata_disc) | edata_pdata_disc == "Inf"] <- NA

# perform linear regression with age, sex, plate, sv1, sv2 as covariates
sumstats_disc <- data.frame()
# run model (for AD vs CO)
for (i in 1:ncol(analyte)) {
  analyte_ID <- analyte_name[i,1]
  protein <- analyte[,i]
  model <- lm(protein ~ as.factor(level) + as.numeric(age_at_csf_draw) + gender + as.factor(PlateId) + as.numeric(sv1) + as.numeric(sv2), data = edata_pdata_disc)
  output <- as.data.frame(summary(model)[4])
  result <- cbind(as.character(analyte_ID), output[2,1], output[2,2], output[2,4])
  sumstats_disc <- rbind(sumstats_disc, result) 
}
colnames(sumstats_disc) <- c("Analyte","Estimate", "Standard_error","Pvalue")
dim(sumstats_disc) # 7029    4
any(is.na(sumstats_disc)) # FALSE

# multiple test correction
cols = c(2,3,4)    
sumstats_disc[,cols] = apply(sumstats_disc[,cols], 2, function(x) as.numeric(as.character(x)))
sumstats_disc$FDR <- p.adjust(sumstats_disc$Pvalue, method = "fdr", n = nrow(sumstats_disc))
sumstats_disc$bonferroni <- p.adjust(sumstats_disc$Pvalue, method = "bonferroni", n = nrow(sumstats_disc))

# analyte annotation
sumstats_disc$newID2 <- sumstats_disc$Analyte
sumstats_disc <- merge(sumstats_disc, protein_names,by= 'newID2')
sumstats_disc$newID2 <- NULL 
sumstats_disc <- sumstats_disc[c("Analyte", "Target", "EntrezGeneSymbol", "Estimate", "Standard_error", "Pvalue", "FDR", "bonferroni")]
dim(sumstats_disc) # 7029    8
cols = c(4,5,6,7,8)    
sumstats_disc[,cols] = apply(sumstats_disc[,cols], 2, function(x) as.numeric(as.character(x)))
sumstats_disc$diffexpressed <- "NO"
sumstats_disc$diffexpressed[sumstats_disc$Pvalue < 0.05] <- "Significant"
sumstats_disc$diffexpressed[sumstats_disc$Pvalue >= 0.05] <- "Not-Significant"
sumstats_disc <- sumstats_disc[order(sumstats_disc$Pvalue),]
write.table(sumstats_disc, file="/Input_Data/Discovery_MAP_Ruiz_ATN_ADvsCO_ageSexPlate_lm_sv_Results.txt", sep="\t", row.names=F, quote=F)

# volcano plot of Discovery
top_10_labels <- unique(c(sumstats_disc[sumstats_disc$Estimate > 0,]$EntrezGeneSymbol[1:7], sumstats_disc[sumstats_disc$Estimate < 0,]$EntrezGeneSymbol[1:5]))
png("/Figures/Volcano_Plot_Discovery.png", units="mm", width=190, height=142, res=1000)
EnhancedVolcano(sumstats_disc,
    lab = sumstats_disc$EntrezGeneSymbol,
    selectLab = top_10_labels,
    x = 'Estimate',
    y = 'FDR',
    title= NULL,
    subtitle = NULL,
    pCutoff = 0.05,
    FCcutoff = 0,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    pointSize = 2,
    labSize = 5,
    colAlpha = 1/5,
    labFace = 'bold',
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = "",
    ylab = "FDR", xlab = "Estimate",
    legendLabels = c("NS", "NS", "FDR", "FDR"),
    xlim = c(min(sumstats_disc$Estimate, na.rm = TRUE), max(sumstats_disc$Estimate, na.rm = TRUE)))
dev.off()



##
## Replication Linear Regression (Stage 2): A-T- vs. A+T+
##

load("/Input_Data/Replication_ADNI_Pau_edata_pdata_withSV.RData")
dim(edata_pdata_repl) # 930 7039
edata_pdata_repl[1:4,c(1:8,7037:7039)]
table(edata_pdata_repl$ATbin)
# A-T- A-T+ A+T- A+T+
#  240   39  268  383

# replace the level according to ATN and keep the group we are interested in (AD and CO)
edata_pdata_repl$ATN <- ifelse(edata_pdata_repl$ATbin == "A+T+", "AD",  
  ifelse(edata_pdata_repl$ATbin == "A-T-", "CO",
  ifelse(edata_pdata_repl$ATbin == "A+T-", "Asymp", "Exclude")))

# reorder columns
edata_pdata_repl <- edata_pdata_repl[,-c("ATbin")]
dim(edata_pdata_repl) #  930 7039
edata_pdata_repl <- edata_pdata_repl %>% select(ExtIdentifier, UniquePhenoID, age_at_csf_draw, gender, PlateId, status_for_analysis, ATN, level, sv1, sv2, everything())
dim(edata_pdata_repl) #  930 7039

# Subset based on CC status from actual pheno so that DLB, PD, and ADAD are not included
edata_pdata_repl <- edata_pdata_repl[edata_pdata_repl$status_for_analysis %in% c("AD", "CO", "OT"),]
dim(edata_pdata_repl) #  864 7039

# Subset based on ATN status - only keeping AD (A+T+) and CO (A-T-)
edata_pdata_repl <- edata_pdata_repl[edata_pdata_repl$ATN %in% c("AD", "CO"),]
dim(edata_pdata_repl) #  593 7039

edata_pdata_repl$level <- ifelse(edata_pdata_repl$ATN == 'CO',0,1)
dim(edata_pdata_repl) #  593 7039
table(edata_pdata_repl$level)
#   0   1
# 235 358
ADCO_repl <- edata_pdata_repl[,c("ExtIdentifier", "UniquePhenoID", "ATN")]
dim(ADCO_repl) # 593   3
ADCO <- rbind(ADCO_disc, ADCO_repl)
dim(ADCO) # 1763    3
write.table(ADCO, file="/Input_Data/ADCO_samples_disc_repl.txt", sep="\t", row.names=F, quote=F)

# filter for analytes that passed FDR-correction in discovery
discovery_FDR <- read.table("/Input_Data/Discovery_MAP_Ruiz_ATN_ADvsCO_ageSexPlate_lm_sv_Results.txt", sep="\t", header=T, stringsAsFactors=F, quote="")
discovery_FDR <- discovery_FDR[discovery_FDR$FDR < 0.05,]
dim(discovery_FDR) # 3565    9

# get analytes
analyte <- as.data.frame(edata_pdata_repl[,11:ncol(edata_pdata_repl)])
dim(analyte) #  593 7029
analyte_FDR <- analyte[,colnames(analyte) %in% discovery_FDR$Analyte]
dim(analyte_FDR) # 593 3565
analyte_name_FDR <- as.data.frame(colnames(analyte_FDR))
table(colnames(analyte_FDR) %in% discovery_FDR$Analyte)
# TRUE
# 3565
# replace the NaN and Inf values in our data frame
edata_pdata_repl[is.na(edata_pdata_repl) | edata_pdata_repl == "Inf"] <- NA

sumstats_repl <- data.frame()
# run model (for AD vs CO)
for (i in 1:ncol(analyte_FDR)) {
  analyte_ID <- analyte_name[i,1]
  protein <- analyte_FDR[,i]
  model <- lm(protein ~ as.factor(level) + as.numeric(age_at_csf_draw) + gender + as.factor(PlateId) + as.numeric(sv1) + as.numeric(sv2), data = edata_pdata_repl)
  output <- as.data.frame(summary(model)[4])
  result <- cbind(as.character(analyte_ID), output[2,1], output[2,2], output[2,4])
  sumstats_repl <- rbind(sumstats_repl, result) 
}
colnames(sumstats_repl) <- c("Analyte","Estimate", "Standard_error","Pvalue")
dim(sumstats_repl) # 3565    4
any(is.na(sumstats_repl)) # FALSE

# multiple test correction
cols = c(2,3,4)    
sumstats_repl[,cols] = apply(sumstats_repl[,cols], 2, function(x) as.numeric(as.character(x)))
sumstats_repl$FDR <- p.adjust(sumstats_repl$Pvalue, method = "fdr", n = nrow(sumstats_repl))
sumstats_repl$bonferroni <- p.adjust(sumstats_repl$Pvalue, method = "bonferroni", n = nrow(sumstats_repl))

# analyte annotation
sumstats_repl$newID2 <- sumstats_repl$Analyte
sumstats_repl <- merge(sumstats_repl, protein_names,by= 'newID2')
sumstats_repl$newID2 <- NULL 
sumstats_repl <- sumstats_repl[c("Analyte", "Target", "EntrezGeneSymbol", "Estimate", "Standard_error", "Pvalue", "FDR", "bonferroni")]
dim(sumstats_repl) # 3565    8
cols = c(4,5,6,7,8)    
sumstats_repl[,cols] = apply(sumstats_repl[,cols], 2, function(x) as.numeric(as.character(x)))
sumstats_repl$diffexpressed <- "NO"
sumstats_repl$diffexpressed[sumstats_repl$Pvalue < 0.05] <- "Significant"
sumstats_repl$diffexpressed[sumstats_repl$Pvalue >= 0.05] <- "Not-Significant"
sumstats_repl <- sumstats_repl[order(sumstats_repl$Pvalue),]
write.table(sumstats_repl, file="/Input_Data/Replication_ADNI_Pau_ATN_ADvsCO_ageSexPlate_lm_sv_Results.txt", sep="\t", row.names=F, quote=F)

# volcano plot of Replication
top_10_labels <- unique(c(sumstats_repl[sumstats_repl$Estimate > 0,]$EntrezGeneSymbol[1:7], sumstats_repl[sumstats_repl$Estimate < 0,]$EntrezGeneSymbol[1:5]))
png("/Figures/Volcano_Plot_Replication.png", units="mm", width=190, height=142, res=1000)
EnhancedVolcano(sumstats_repl,
    lab = sumstats_repl$EntrezGeneSymbol,
    selectLab = top_10_labels,
    x = 'Estimate',
    y = 'FDR', 
    title= NULL, 
    subtitle = NULL, 
    pCutoff = 0.05,
    FCcutoff = 0,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    pointSize = 2,
    labSize = 5,
    colAlpha = 1/5,
    labFace = 'bold',
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = "",
    ylab = "FDR", xlab = "Estimate",
    legendLabels = c("NS", "NS", "FDR", "FDR"),
    xlim = c(min(sumstats_repl$Estimate, na.rm = TRUE), max(sumstats_repl$Estimate, na.rm = TRUE)))
dev.off()



##
## Meta-analysis (Stage 3): A-T- vs. A+T+
##

discovery_DAA_FDR <- sumstats_disc[sumstats_disc$FDR < 0.05,]
dim(discovery_DAA_FDR) # 3565    9 (FDR in Discovery)
replication_DAA_FDR <- replication_DAA[replication_DAA$FDR < 0.05,]
dim(replication_DAA_FDR) # 2608    9 (FDR in replication)

# concordance check
disc_repl_fdr <- inner_join(discovery_DAA_FDR, replication_DAA_FDR, by="Analyte")
dim(disc_repl_fdr) # 2608   17
disc_repl_fdr_up <- disc_repl_fdr[disc_repl_fdr$Estimate.x > 0 & disc_repl_fdr$Estimate.y > 0, ] # 1659
disc_repl_fdr_down <- disc_repl_fdr[disc_repl_fdr$Estimate.x < 0 & disc_repl_fdr$Estimate.y < 0, ] # 884
discRepl <- rbind(disc_repl_fdr_up, disc_repl_fdr_down)
dim(discRepl) # 2543   17 (FDR in Discovery & Replication + Same direction)

# combine p-values based on weighted-z method for meta-analysis
discRepl$weightedZ <- combinePValues(discRepl$Pvalue.x, discRepl$Pvalue.y, method="z", weights=1:2)
dim(discRepl) # 2543   18
discRepl[1:3,c(1,4,6,7,8,12,14,15,16,18)]
discRepl$weightedZ_FDR <- p.adjust(discRepl$weightedZ, method = "fdr", n = nrow(discRepl))
discRepl$weightedZ_bonf <- p.adjust(discRepl$weightedZ, method = "bonferroni", n = nrow(discRepl))
discRepl_meta_bonf <- discRepl[discRepl$weightedZ < 3.4e-05,]
dim(discRepl_meta_bonf) # 2173   20 (Bonf_meta + FDR in Discovery & Replication + Same direction)
write.table(discRepl_meta_bonf, file="/Input_Data/Discovery_ReplicationOfDiscoveryFDR_ATN_ageSexPlate_lm_sv_Metaanalysis_Bonf_Filtered.txt", sep="\t", row.names=F, quote=F)

# volcano plot of meta-analysis
discRepl_meta_bonf <- discRepl_meta_bonf[order(discRepl_meta_bonf$weightedZ_bonf),]
top_10_labels <- unique(c(discRepl_meta_bonf[discRepl_meta_bonf$Estimate.x > 0,]$EntrezGeneSymbol.x[1:7], discRepl_meta_bonf[discRepl_meta_bonf$Estimate.x < 0,]$EntrezGeneSymbol.x[1:5]))
png("/Figures/Volcano_Plot_Metaanalysis.png", units="mm", width=190, height=142, res=1000)
EnhancedVolcano(discRepl_meta_bonf,
    lab = discRepl_meta_bonf$EntrezGeneSymbol.x,
    selectLab = top_10_labels,
    x = 'Estimate.x',
    y = 'weightedZ_bonf',
    title= NULL,
    subtitle = NULL,
    pCutoff = 0.05,
    FCcutoff = 0,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    pointSize = 2,
    labSize = 5,
    colAlpha = 1/5,
    labFace = 'bold',
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = "",
    ylab = "P-Bonferroni", xlab = "Estimate",
    legendLabels = c("NS", "NS", "P-Bonferroni", "P-Bonferroni"),
    xlim = c(min(discRepl_meta_bonf$Estimate.x, na.rm = TRUE), max(discRepl_meta_bonf$Estimate.x, na.rm = TRUE)))
dev.off()