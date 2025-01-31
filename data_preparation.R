rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(2022)

library(sva)
library(ggrepel)
library(data.table)
library(ggplot2)
library(Biobase)
library(dplyr)
library(EnhancedVolcano)
library(readxl)

##
## Discovery dataset preparation
##

# read proteomics data from Discovery cohorts (Knight-ADRC {MAP} and FACE {Ruiz})
MAP_expression <- readRDS("/Input_Data/MAP_931sample_7029analytes.rds")
Ruiz_expression <- readRDS("/Input_Data/Ruiz_621sample_7029analytes.rds")

# CSF AT dichotomization data
AT <- read.table("/Input_Data/ATdichotomization.allsubjects.20220322.txt", header=T, sep="\t", stringsAsFactors=F)
dim(AT) # 3040    9

AT <- AT[,c("UniquePhenoID", "ATbin")]
AT <- na.omit(AT)
dim(AT) # 3021    2

# read phenotype file
pheno_original <- data.table(read_excel("/Input_Data/Cohort5_CSF_pheno_data_summary_DataFreezeNov2021.xlsx", sheet = "Pheno_data_Carlos"))
dim(pheno_original) # 3065   82
length(unique(pheno_original$UniquePhenoID)) # 3064

pheno_original <- pheno_original[,c("UniquePhenoID", "CC_at_LP_draw...25")]
colnames(pheno_original) <- c("UniquePhenoID", "CC_at_LP_draw")
dim(pheno_original) # 3065    2

# convert CC_Status to AD, CO, OT, and Exclude.
pheno_original$status_for_analysis <- ifelse(pheno_original$CC_at_LP_draw == "AD", "AD",  
  ifelse(pheno_original$CC_at_LP_draw == "CO", "CO",
  ifelse(pheno_original$CC_at_LP_draw == "OT", "OT",
  ifelse(pheno_original$CC_at_LP_draw == "Neuro_OT", "OT",
  ifelse(pheno_original$CC_at_LP_draw == "OT(CO)", "CO", "Exclude")))))
table(pheno_original$status_for_analysis)
#     AD      CO Exclude      OT
#   1198    1017     598     252

# drop original CC_Status column and append ATbin.
pheno_original <- pheno_original[,-2]
pheno_original <- inner_join(pheno_original, AT, by="UniquePhenoID")
dim(pheno_original) # 3022    3
table(pheno_original$status_for_analysis)
#     AD      CO Exclude      OT
#   1188    1007     577     250

# combine expression data with phenotype
discovery_expr <- combineTwoExpressionSet(MAP_expression, Ruiz_expression)
dim(discovery_expr)
# Features  Samples
#     7029     1552

pheno <- data.table(pData(discovery_expr))
dim(pheno) # 1552   46
write.table(pheno, file="/Input_Data/Discovery_MAP_Ruiz_Pheno_from_expression.txt", sep="\t", row.names=F, quote=F)

# only keep required columns from the pheno
pheno_subset <- pheno[,list(PlateId,ExtIdentifier,UniquePhenoID,gender,age_at_csf_draw)]
pheno_subset <- inner_join(pheno_subset, pheno_original, by="UniquePhenoID")
# pheno_subset <- pheno_subset[,-3]
dim(pheno_subset) # 1545    7

# extract expression matrix
edata <- exprs(discovery_expr)
edata.log <- log10(edata)

# transpose matrix for analysis and append pheno columns at the end
transposed_edata <- (t(edata.log))
transposed_edata <- cbind(ExtIdentifier = colnames(edata), transposed_edata)
edata_pdata_disc <- merge(transposed_edata, pheno_subset, by ='ExtIdentifier')
dim(edata_pdata_disc) # 1545 7036

edata_pdata_disc <- edata_pdata_disc %>% select(ExtIdentifier, UniquePhenoID, status_for_analysis, ATbin, age_at_csf_draw, gender, PlateId, everything())
edata_pdata_disc <- data.table(edata_pdata_disc)
dim(edata_pdata_disc) # 1545 7036

# assign factor level to CA/CO, CO <- 0, CA <- 1
edata_pdata_disc <- edata_pdata_disc[status_for_analysis %in% c('AD','CO','OT','Exclude')]
dim(edata_pdata_disc) #  1545 7036

edata_pdata_disc$level<-ifelse(edata_pdata_disc$status_for_analysis=='CO',0,1)
dim(edata_pdata_disc) #  1545 7037


# reorder columns
edata_pdata_disc <- edata_pdata_disc %>% select(ExtIdentifier, UniquePhenoID, age_at_csf_draw, gender, PlateId, status_for_analysis, ATbin, level, everything())
dim(edata_pdata_disc) #  1545 7037

analyte <- as.data.frame(edata_pdata_disc[,9:ncol(edata_pdata_disc)])
dim(analyte) #  1545 7029
analyte_name <- as.data.frame(colnames(analyte))


## SVA calculation for using them as covariates
edata_sva <- analyte
row.names(edata_sva) <- edata_pdata_disc$ExtIdentifier
edata_sva <- data.matrix(edata_sva, rownames.force = TRUE)
edata_sva <- t(edata_sva)
edata_sva[is.na(edata_sva)] <- 0
mod = model.matrix(~as.factor(status_for_analysis), data=edata_pdata_disc)
row.names(mod) <- edata_pdata_disc$ExtIdentifier
mod[,2] <- ifelse(mod[,2] > 0, 0, 1)
mod0 <- model.matrix(~1, data=edata_pdata_disc)
row.names(mod0) <- edata_pdata_disc$ExtIdentifier
n.sv = num.sv(edata_sva, mod, method="be", B = 20, seed = 2022) 
svobj <- sva(edata_sva, mod, mod0, n.sv = 2)
dim(svobj$sv) # 1545    2
SVs <- as.data.frame(svobj$sv)
colnames(SVs) <- c("sv1", "sv2")
SVs$ExtIdentifier <- edata_pdata_disc$ExtIdentifier
edata_pdata_disc <- inner_join(edata_pdata_disc, SVs, by="ExtIdentifier")
dim(edata_pdata_disc) # 1545 7039
## SVA

feature.dt <- data.table(fData(discovery_expr))
dim(feature.dt) # 7029   93
## modify analyte name
feature.dt[, newID := gsub('-', '.', SeqId)]
feature.dt[, newID2 := paste0('X', newID)]
protein_names <- feature.dt[,list(newID2,Target,EntrezGeneSymbol, UniProt)]
dim(protein_names) # 7029    4
save(edata_pdata_disc, protein_names, file="/Input_Data/Discovery_MAP_Ruiz_edata_pdata_withSV.RData")


##
## Replication dataset preparation
##

# read proteomics data from replication cohorts (ADNI, Barcelona-1 {Pau})
ADNI_expression <- readRDS("/Input_Data/Proteomics/ADNI_737sample_7029analytes.rds")
Pau_expression <- readRDS("/Input_Data/Pau_207sample_7029analytes.rds")

replication_expr <- combineTwoExpressionSet(ADNI_expression, Pau_expression)
dim(replication_expr)
# Features  Samples
#     7029      944

pheno <- data.table(pData(replication_expr))
dim(pheno) # 944   46
head(pheno, 2)
write.table(pheno, file="/Input_Data/Replication_ADNI_Pau_Pheno_from_expression.txt", sep="\t", row.names=F, quote=F)

# only keep required columns from the pheno
pheno_subset <- pheno[,list(PlateId,ExtIdentifier,UniquePhenoID,gender,age_at_csf_draw)]
pheno_subset <- inner_join(pheno_subset, pheno_original, by="UniquePhenoID")
# pheno_subset <- pheno_subset[,-3]
dim(pheno_subset) # 930   7

# extract expression matrix
edata <- exprs(replication_expr)
edata.log <- log10(edata)

# transpose matrix for analysis and append pheno columns at the end
transposed_edata <- (t(edata.log))
transposed_edata <- cbind(ExtIdentifier = colnames(edata), transposed_edata)
edata_pdata_repl <- merge(transposed_edata, pheno_subset, by ='ExtIdentifier')
dim(edata_pdata_repl) # 930 7036
edata_pdata_repl <- edata_pdata_repl %>% select(ExtIdentifier, UniquePhenoID, status_for_analysis, ATbin, age_at_csf_draw, gender, PlateId, everything())
edata_pdata_repl <- data.table(edata_pdata_repl)
dim(edata_pdata_repl) # 930 7036

# assign factor level to CA/CO, CO <- 0, CA <- 1
edata_pdata_repl <- edata_pdata_repl[status_for_analysis %in% c('AD','CO','OT','Exclude')]
dim(edata_pdata_repl) #  930 7036
edata_pdata_repl$level<-ifelse(edata_pdata_repl$status_for_analysis=='CO',0,1)
dim(edata_pdata_repl) #  930 7037

# reorder columns
edata_pdata_repl <- edata_pdata_repl %>% select(ExtIdentifier, UniquePhenoID, age_at_csf_draw, gender, PlateId, status_for_analysis, ATbin, level, everything())
dim(edata_pdata_repl) #  930 7037

analyte <- as.data.frame(edata_pdata_repl[,9:ncol(edata_pdata_repl)])
dim(analyte) #  930 7029
analyte_name <- as.data.frame(colnames(analyte))

## SVA calculation for using them as covariates
edata_sva <- analyte
row.names(edata_sva) <- edata_pdata_repl$ExtIdentifier
edata_sva <- data.matrix(edata_sva, rownames.force = TRUE)
edata_sva <- t(edata_sva)
edata_sva[is.na(edata_sva)] <- 0
mod = model.matrix(~as.factor(status_for_analysis), data=edata_pdata_repl)
row.names(mod) <- edata_pdata_repl$ExtIdentifier
mod[,2] <- ifelse(mod[,2] > 0, 0, 1) # while converting "status_for_analysis" to factor it convert AD to 0 and CO to 1, so reverting that
mod0 <- model.matrix(~1, data=edata_pdata_repl)
row.names(mod0) <- edata_pdata_repl$ExtIdentifier
n.sv = num.sv(edata_sva, mod, method="be", B = 20, seed = 2022) 
svobj <- sva(edata_sva, mod, mod0, n.sv = 2)
dim(svobj$sv) # 930   2
SVs <- as.data.frame(svobj$sv)
colnames(SVs) <- c("sv1", "sv2")
head(SVs, 3)
SVs$ExtIdentifier <- edata_pdata_repl$ExtIdentifier
edata_pdata_repl <- inner_join(edata_pdata_repl, SVs, by="ExtIdentifier")
dim(edata_pdata_repl) # 930 7039
## SVA

feature.dt <- data.table(fData(replication_expr))
dim(feature.dt) # 7029   93
## modify analyte name
feature.dt[, newID := gsub('-', '.', SeqId)]
feature.dt[, newID2 := paste0('X', newID)]
protein_names <- feature.dt[,list(newID2, Target, EntrezGeneSymbol, UniProt)]
dim(protein_names) # 7029    4
head(protein_names, 2)

save(edata_pdata_repl, protein_names, file="/Input_Data/Replication_ADNI_Pau_edata_pdata_withSV.RData")