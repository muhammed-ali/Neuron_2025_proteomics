rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(2022)

library(dplyr)

# read sumstats file from meta-analysis
ADCO <- read.table("/Input_Data/Discovery_ReplicationOfDiscoveryFDR_ATN_ageSexPlate_lm_sv_Metaanalysis_Bonf_Filtered.txt", sep="\t", header=T, stringsAsFactors=F, fill=T)
dim(ADCO) # 2173   20
length(unique(ADCO$Analyte)) # 2173
length(unique(ADCO$EntrezGeneSymbol.x)) # 2029
gene_mapping <- ADCO[,c("Analyte", "EntrezGeneSymbol.x", "Target.x")]


# read sumstats file containing "estimate" and "p_value", "average_expression" of each analyte from:
# 1) ADCO (A-T- vs. A+T+)
# 2) sympCO (A-T- vs. A+T-)
# 3) sympAD (A+T- vs. A+T+)

newdf <- read.table("/Input_Data/Estimate_Pval_of_ADCO_sympAD_sympCO_2759_MeanAbundance.txt", sep="\t", header=T, stringsAsFactors=F)
dim(newdf) # 2759   10
newdf <- newdf[newdf$Analyte %in% ADCO$Analyte,]
dim(newdf) # 2173   10

threshold <- 0.05

# proteins displaying significant increase from A-T- to A+T- but no significant difference from A+T- to A+T+
c1 <- newdf[newdf$Pvalue.sympCO < threshold & newdf$Estimate.sympCO > 0 & newdf$Pvalue.sympAD > threshold,]
dim(c1) # 16 10

# proteins displaying significant increase from A-T- to A+T- and significant decrease from A+T- to A+T+
c23 <- newdf[newdf$Pvalue.sympCO < threshold & newdf$Estimate.sympCO > 0 & newdf$Pvalue.sympAD < threshold & newdf$Estimate.sympAD < 0,]
dim(c23) # 486  10
c2 <- c23[c23$Estimate.ADCO > 0,]
dim(c2) # 4 10
#head(c2)
c3 <- c23[c23$Estimate.ADCO < 0,]
dim(c3) # 482  10

# proteins displaying linear significant increase from A-T- to A+T- and from A+T- to A+T+ 
c4 <- newdf[newdf$Pvalue.sympCO < threshold & newdf$Estimate.sympCO > 0 & newdf$Pvalue.sympAD < threshold & newdf$Estimate.sympAD > 0,]
dim(c4) # 46 10

# proteins displaying linear decrease from A-T- to A+T- and from A+T- to A+T+ 
c5 <- newdf[newdf$Pvalue.sympCO < threshold & newdf$Estimate.sympCO < 0 & newdf$Pvalue.sympAD > threshold,]
dim(c5) # 16 10
c67 <- newdf[newdf$Pvalue.sympCO < threshold & newdf$Estimate.sympCO < 0 & newdf$Pvalue.sympAD < threshold & newdf$Estimate.sympAD > 0,]
dim(c67) # 1074  10
c6 <- c67[c67$Estimate.ADCO < 0,]
dim(c6) # 38 10

# proteins displaying significant decrease from A-T- to A+T- and significant increase from A+T- to A+T+
c7 <- c67[c67$Estimate.ADCO > 0,]
dim(c7) # 1036   10

# proteins displaying significant linear decrease from A-T- to A+T- and from A+T- to A+T+ 
c8 <- newdf[newdf$Pvalue.sympCO < threshold & newdf$Estimate.sympCO < 0 & newdf$Pvalue.sympAD < threshold & newdf$Estimate.sympAD < 0,]
dim(c8) # 2 10

# proteins displaying non-significant differences between A-T- to A+T- but significant increase from A+T- to A+T+ 
c9 <- newdf[newdf$Pvalue.sympCO > threshold & newdf$Pvalue.sympAD < threshold & newdf$Estimate.sympAD > 0,]
dim(c9) # 395  10

# proteins displaying non-significant differences between A-T- to A+T- but significant decrease from A+T- to A+T+ 
c10 <- newdf[newdf$Pvalue.sympCO > threshold & newdf$Pvalue.sympAD < threshold & newdf$Estimate.sympAD < 0,]
dim(c10) # 115  10

c1$cluster <- 1
c2$cluster <- 2
c3$cluster <- 3
c4$cluster <- 4
c5$cluster <- 5
c6$cluster <- 6
c7$cluster <- 7
c8$cluster <- 8
c9$cluster <- 9
c10$cluster <- 10

all_clusters <- rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)
dim(all_clusters) # 2150   11
head(all_clusters, 3)
table(all_clusters$cluster)
#   1    2    3    4    5    6    7    8    9   10
#  16    4  482   46   16   38 1036    2  395  115

missing <- newdf[!(newdf$Analyte %in% all_clusters$Analyte),]
dim(missing) # 23  10
head(missing)

# none is significant beside CO_AD
nrow(missing[missing$Pvalue.sympCO < threshold,]) # 0
nrow(missing[missing$Pvalue.sympAD < threshold,]) # 0

nrow(missing[missing$Estimate.sympCO > 0,]) # 10
nrow(missing[missing$Estimate.sympAD > 0,]) # 11
nrow(missing[missing$Estimate.sympCO > 0 & missing$Estimate.sympAD > 0,]) # 9

nrow(missing[missing$Estimate.sympCO < 0,]) # 13
nrow(missing[missing$Estimate.sympAD < 0,]) # 12
nrow(missing[missing$Estimate.sympCO < 0 & missing$Estimate.sympAD < 0,]) # 11

# proteins displaying non-significant overall increase from A-T- to A+T- and from A+T- to A+T+
c11 <- newdf[newdf$Pvalue.sympCO > threshold & newdf$Pvalue.sympAD > threshold & newdf$Estimate.sympCO > 0,]
dim(c11) # 10  10

# proteins displaying non-significant overall decrease from A-T- to A+T- and from A+T- to A+T+
c12 <- newdf[newdf$Pvalue.sympCO > threshold & newdf$Pvalue.sympAD > threshold & newdf$Estimate.sympCO < 0,]
dim(c12) # 13  10

c11$cluster <- 11
c12$cluster <- 12

all_clusters <- rbind(all_clusters,c11,c12)
dim(all_clusters) # 2173   11
head(all_clusters, 3)
table(all_clusters$cluster)
#    1    2    3    4    5    6    7    8    9   10   11   12
#   16    4  482   46   16   38 1036    2  395  115   10   13
length(unique(all_clusters$Analyte)) # 2173

all_clusters$signif_status <- ifelse(all_clusters$Pvalue.sympCO < threshold & all_clusters$Pvalue.sympAD < threshold, "11", 0)
all_clusters$signif_status <- ifelse(all_clusters$Pvalue.sympCO > threshold & all_clusters$Pvalue.sympAD < threshold, "01", all_clusters$signif_status)
all_clusters$signif_status <- ifelse(all_clusters$Pvalue.sympCO < threshold & all_clusters$Pvalue.sympAD > threshold, "10", all_clusters$signif_status)
all_clusters$signif_status <- ifelse(all_clusters$Pvalue.sympCO > threshold & all_clusters$Pvalue.sympAD > threshold, "00", all_clusters$signif_status)
table(all_clusters$signif_status)
#  00   01   10   11
#  23  510   32 1608
23+510+32+1608 # 2173

all_clusters$cluster <- paste0("c", all_clusters$cluster)
table(all_clusters$cluster, all_clusters$signif_status)
#        00   01   10   11
#  c1     0    0   16    0
#  c10    0  115    0    0
#  c11   10    0    0    0
#  c12   13    0    0    0
#  c2     0    0    0    4
#  c3     0    0    0  482
#  c4     0    0    0   46
#  c5     0    0   16    0
#  c6     0    0    0   38
#  c7     0    0    0 1036
#  c8     0    0    0    2
#  c9     0  395    0    0

all_clusters <- inner_join(all_clusters, gene_mapping, by="Analyte")
dim(all_clusters) # 2173   14
colnames(all_clusters)[13] <- "Gene"
head(all_clusters, 3)
all_clusters[198:200,]
# split rows when multiple symbols are given, e.g. ITGA2B|ITGB3
library(tidyr)
all_clusters_unique <- separate_rows(all_clusters, Gene, sep="\\|", convert = TRUE)
all_clusters_unique <- as.data.frame(all_clusters_unique)
dim(all_clusters_unique) # 2194   14


# keep the first gene associated with a protein (Benoit Lehallier et al. Nature Medicine, 2019)
# because several individual proteins (33 of 2,925) were mapped to multiple gene symbols, we kept only the first gene symbol provided by SomaLogic to prevent false-positive enrichment.
all_clusters_unique <- all_clusters_unique[!duplicated(all_clusters_unique$Analyte),]
dim(all_clusters_unique) # 2173   14

# replace gene name for rows which have missing ("None") Gene column
nrow(all_clusters_unique[all_clusters_unique$Gene == "None",]) # 4
all_clusters_unique$Gene <- ifelse(all_clusters_unique$Gene == "None", all_clusters_unique$Target.x, all_clusters_unique$Gene)
all_clusters_unique$Target.x <- NULL # remove this colum, don't need it anymore

# convert Symbol to Entrez ID
library(org.Hs.eg.db)
all_clusters_unique$Entrez <- mapIds(org.Hs.eg.db, all_clusters_unique$Gene, 'ENTREZID', 'SYMBOL')
all_clusters_unique <- all_clusters_unique %>% mutate(Entrez = sapply(Entrez, toString))
dim(all_clusters_unique) # 2173   14
length(unique(all_clusters_unique$Entrez)) # 2026 (including one unique "NA")

# 4 genes do not have Entrez ID. Let's search them manually.
all_clusters_unique[all_clusters_unique$Entrez == "NA",]

# TSTA3 = GFUS = 7264 (https://www.ncbi.nlm.nih.gov/gene/7264)
# MPP6 = PALS2 = 51678 (https://www.ncbi.nlm.nih.gov/gene/51678)
# MAGEA5 = MAGEA5P = 4104 (https://www.ncbi.nlm.nih.gov/gene/4104)
# GGT2 = GGT2P = 728441 (https://www.ncbi.nlm.nih.gov/gene/728441)
# MENT = C1orf56 = 54964 (https://www.ncbi.nlm.nih.gov/gene/54964)

all_clusters_unique$Entrez <- ifelse(all_clusters_unique$Gene == "TSTA3", "7264", all_clusters_unique$Entrez)
all_clusters_unique$Entrez <- ifelse(all_clusters_unique$Gene == "MPP6", "51678", all_clusters_unique$Entrez)
all_clusters_unique$Entrez <- ifelse(all_clusters_unique$Gene == "MAGEA5", "4104", all_clusters_unique$Entrez)
all_clusters_unique$Entrez <- ifelse(all_clusters_unique$Gene == "GGT2", "728441", all_clusters_unique$Entrez)
all_clusters_unique$Entrez <- ifelse(all_clusters_unique$Gene == "MENT", "54964", all_clusters_unique$Entrez)
any(is.na(all_clusters_unique$Entrez)) # FALSE

# get groups based on expression change patterns
all_clusters_unique$group <- ifelse(all_clusters_unique$cluster %in% c("c1", "c2", "c4", "c9","c11"), "g1_up", "unk")
all_clusters_unique$group <- ifelse(all_clusters_unique$cluster %in% c("c3"), "g2_upDown", all_clusters_unique$group)
all_clusters_unique$group <- ifelse(all_clusters_unique$cluster %in% c("c5","c6","c8","c10","c12"), "g3_down", all_clusters_unique$group)
all_clusters_unique$group <- ifelse(all_clusters_unique$cluster %in% c("c7"), "g4_downUp", all_clusters_unique$group)
table(all_clusters_unique$group)
#    g1_up g2_upDown   g3_down g4_downUp
#      471       482       184      1036

dim(all_clusters_unique) # 2173   15
length(unique(all_clusters_unique$Analyte)) # 2173
length(unique(all_clusters_unique$Gene)) # 2030
length(unique(all_clusters_unique$Entrez)) # 2030
max(all_clusters_unique$Pvalue.ADCO) # 0.0250375

write.table(all_clusters_unique, file="/Input_Data/Metaanalysis_Bonf_PsuedoTrajctories.txt", sep="\t", row.names=F, quote=F)