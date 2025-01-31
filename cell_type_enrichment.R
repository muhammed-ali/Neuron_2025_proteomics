rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(2022)

library(data.table)
library(dplyr)
library(pheatmap)

# read the cell type composition file (https://doi.org/https://doi.org/10.1016/j.neuron.2015.11.013.)
ct_composition <- data.table(read_excel("/Input_Data/brain_cell_type_specificity.xlsx", sheet = "cell_types_with_percentiles"))
dim(ct_composition) # 23225    52
ct_composition <- ct_composition[,c(1,52)]
ct_composition <- ct_composition[3:nrow(ct_composition),]
colnames(ct_composition) <- c("Gene", "Cell_Type")
ct_composition[is.na(ct_composition)] <- "Unknown" # Replace missing gene annotation with "Unknown"
ct_composition <- unique(ct_composition)
dim(ct_composition) # 21709     2
head(ct_composition, 3)
#      Gene                  Cell_Type
# 1:   APOE    Human mature astrocytes
# 2:    APP                    Unknown
# 3: MS4A4A Human Microglia/Macrophage
ct_composition <- unique(ct_composition[ct_composition$Cell_Type != "Unknown",])
dim(ct_composition) # 4993    2
head(ct_composition, 3)
#      Gene                  Cell_Type
# 1:   APOE    Human mature astrocytes
# 2: MS4A4A Human Microglia/Macrophage
# 3: MS4A6A Human Microglia/Macrophage

# read the complete proteomics gene lists
ADCO_complete <- read.table("/Input_Data/Discovery_MAP_Ruiz_ATN_ADvsCO_ageSexPlate_lm_sv_Results.txt", sep="\t", header=T, stringsAsFactors=F, fill=T, quote="")
dim(ADCO_complete) # 7029   9
colnames(ADCO_complete)[3] <- "Gene"
# Replace gene name for rows which have missing ("None") Gene column
ADCO_complete$Gene <- ifelse(ADCO_complete$Gene == "None", ADCO_complete$Target, ADCO_complete$Gene)
# Split rows when multiple symbols are given, e.g. PPP3CA|PPP3R1
ADCO_complete <- as.data.frame(separate_rows(ADCO_complete, Gene, sep="\\|", convert = TRUE))
ADCO_complete <- unique(ADCO_complete[!duplicated(ADCO_complete$Gene),])
dim(ADCO_complete) # 6162    9

# Read the significant proteins gene lists
ADCO_significant <- read.table("/Input_Data/Metaanalysis_Bonf_PsuedoTrajctories.txt", sep="\t", header=T, stringsAsFactors=F, quote="", fill=T)
ADCO_significant <- as.data.frame(separate_rows(ADCO_significant, Gene, sep="\\|", convert = TRUE))
ADCO_significant <- unique(ADCO_significant[!duplicated(ADCO_significant$Gene),])
dim(ADCO_significant) # 2030   15
ADCO_significant <- left_join(ADCO_significant, ct_composition, by=c("Gene"="Gene"))
ADCO_significant$Cell_Type[is.na(ADCO_significant$Cell_Type)] <- "Unknown"
dim(ADCO_significant) # 2030   16

# get enrichment p-value
background <- length(unique(ADCO_complete$Gene)) # 6162
CT_df <- as.data.frame(table(ADCO_significant$Cell_Type))
colnames(CT_df) <- c("cell_type", "CT_labels")
# CT_df
#                    cell_type CT_labels
# 1          Human Endothelial        64
# 2    Human mature astrocytes       141
# 3 Human Microglia/Macrophage       166
# 4              Human Neurons       236
# 5     Human Oligodendrocytes        48
# 6                    Unknown      1375

# function for getting CT enrichment
CT_enrich_stats <- function(cluster, ADCO_significant, background){
 genes_in_group <- unique(ADCO_significant[ADCO_significant$group %in% cluster,]$Gene)
 pathway_CT_groups <- as.data.frame(table(ADCO_significant[ADCO_significant$group %in% cluster,]$Cell_Type)) # Entrez
 colnames(pathway_CT_groups) <- c("cell_type", "freq")
 pathway_CT_groups <- pathway_CT_groups[order(pathway_CT_groups$cell_type),]
 pathway_CT_groups <<- inner_join(pathway_CT_groups, CT_df, by="cell_type")
 pathway_CT_groups <- pathway_CT_groups[pathway_CT_groups$cell_type != "Unknown",]
 pathway_CT_groups$enrichment_P <- phyper(pathway_CT_groups$freq, pathway_CT_groups$CT_labels, background-pathway_CT_groups$CT_labels, length(genes_in_group), lower.tail = FALSE)
 print(pathway_CT_groups)
}

table(ADCO_significant$group)
#    g1_up g2_upDown   g3_down g4_downUp
#      426       475       176       953
      
CT_enrich_stats("g1_up", ADCO_significant, background)
#                   cell_type freq CT_labels enrichment_P
#1          Human Endothelial   19        64 4.580748e-09
#2    Human mature astrocytes   32       141 2.414670e-10
#3 Human Microglia/Macrophage   39       166 9.759877e-13
#4              Human Neurons   46       236 1.493975e-11
#5     Human Oligodendrocytes    9        48 1.354573e-03
CT_enrich_stats("g2_upDown", ADCO_significant, background)
#                   cell_type freq CT_labels enrichment_P
#1          Human Endothelial   12        64 9.848551e-04
#2    Human mature astrocytes   19       141 5.439404e-03
#3 Human Microglia/Macrophage   46       166 6.922606e-16
#4              Human Neurons   24       236 6.295029e-02
#5     Human Oligodendrocytes    6        48 7.287790e-02
CT_enrich_stats("g3_down", ADCO_significant, background)
#                   cell_type freq CT_labels enrichment_P
#1          Human Endothelial    5        64 0.0094872132
#2    Human mature astrocytes    9       141 0.0068441937
#3 Human Microglia/Macrophage   12       166 0.0008627106
#4              Human Neurons   14       236 0.0029148839
#5     Human Oligodendrocytes    3        48 0.0472433539
CT_enrich_stats("g4_downUp", ADCO_significant, background)
#                   cell_type freq CT_labels enrichment_P
#1          Human Endothelial   28        64 1.219625e-08
#2    Human mature astrocytes   81       141 5.461365e-32
#3 Human Microglia/Macrophage   69       166 5.109273e-17
#4              Human Neurons  152       236 1.755202e-69
#5     Human Oligodendrocytes   30        48 1.427749e-14
 
# get enrichment p-value
ct_enrich_calc <- function(foreground, background){
	foreground <- foreground[foreground$Gene %in% background$Gene,]
	merge_data <- inner_join(foreground, background, by="Gene")
	ct_count_table <- as.data.frame(table(merge_data$Cell_Type))
	colnames(ct_count_table) <- c("CT", "Count")
	ct_count_table$Percent <- round((ct_count_table$Count/nrow(foreground))*100, 1)
	return(ct_count_table)
}
ct_enrich_all <- ct_enrich_calc(ADCO_complete, ct_composition)
colnames(ct_enrich_all) <- c("CT", "All_n", "All_percent")
ct_enrich_all
#                           CT All_n All_percent
# 1          Human Endothelial   242        14.1
# 2    Human mature astrocytes   317        18.4
# 3 Human Microglia/Macrophage   472        27.4
# 4              Human Neurons   554        32.2
# 5     Human Oligodendrocytes   137         8.0

ct_enrich_signif <- ct_enrich_calc(ADCO_significant, ct_composition)
colnames(ct_enrich_signif) <- c("CT", "Signif_n", "Signif_percent")
ct_enrich_signif
#                           CT Signif_n Signif_percent
# 1          Human Endothelial       64            9.8
# 2    Human mature astrocytes      141           21.5
# 3 Human Microglia/Macrophage      166           25.3
# 4              Human Neurons      236           36.0
# 5     Human Oligodendrocytes       48            7.3

ct_enrich_g1 <- ct_enrich_calc(ADCO_significant[ADCO_significant$group == "g1_up",], ct_composition)
colnames(ct_enrich_g1) <- c("CT", "G1_n", "G1_percent")
ct_enrich_g1
#                           CT G1_n G1_percent
# 1          Human Endothelial   19       13.1
# 2    Human mature astrocytes   32       22.1
# 3 Human Microglia/Macrophage   39       26.9
# 4              Human Neurons   46       31.7
# 5     Human Oligodendrocytes    9        6.2

ct_enrich_g2 <- ct_enrich_calc(ADCO_significant[ADCO_significant$group == "g2_upDown",], ct_composition)
colnames(ct_enrich_g2) <- c("CT", "G2_n", "G2_percent")
ct_enrich_g2
#                           CT G2_n G2_percent
# 1          Human Endothelial   12       11.2
# 2    Human mature astrocytes   19       17.8
# 3 Human Microglia/Macrophage   46       43.0
# 4              Human Neurons   24       22.4
# 5     Human Oligodendrocytes    6        5.6


ct_enrich_g3 <- ct_enrich_calc(ADCO_significant[ADCO_significant$group == "g3_down",], ct_composition)
colnames(ct_enrich_g3) <- c("CT", "G3_n", "G3_percent")
ct_enrich_g3
#                           CT G3_n G3_percent
# 1          Human Endothelial    5       11.6
# 2    Human mature astrocytes    9       20.9
# 3 Human Microglia/Macrophage   12       27.9
# 4              Human Neurons   14       32.6
# 5     Human Oligodendrocytes    3        7.0


ct_enrich_g4 <- ct_enrich_calc(ADCO_significant[ADCO_significant$group == "g4_downUp",], ct_composition)
colnames(ct_enrich_g4) <- c("CT", "G4_n", "G4_percent")
ct_enrich_g4
#                           CT G4_n G4_percent
# 1          Human Endothelial   28        7.8
# 2    Human mature astrocytes   81       22.5
# 3 Human Microglia/Macrophage   69       19.2
# 4              Human Neurons  152       42.2
# 5     Human Oligodendrocytes   30        8.3


# merge all data frames in list
df_list <- list(ct_enrich_all, ct_enrich_signif, ct_enrich_g1, ct_enrich_g2, ct_enrich_g3, ct_enrich_g4)
ct_enrich_df <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
dim(ct_enrich_df) #  5 13
ct_enrich_df <- ct_enrich_df[,c("CT", "All_percent", "Signif_percent", "G1_percent", "G2_percent", "G3_percent", "G4_percent")]
ct_enrich_df$Ratio_Signif <- round(ct_enrich_df$Signif_percent / ct_enrich_df$All_percent, 1)
ct_enrich_df$Ratio_G1 <- round(ct_enrich_df$G1_percent / ct_enrich_df$All_percent, 1)
ct_enrich_df$Ratio_G2 <- round(ct_enrich_df$G2_percent / ct_enrich_df$All_percent, 1)
ct_enrich_df$Ratio_G3 <- round(ct_enrich_df$G3_percent / ct_enrich_df$All_percent, 1)
ct_enrich_df$Ratio_G4 <- round(ct_enrich_df$G4_percent / ct_enrich_df$All_percent, 1)
ct_enrich_df$FC_Signif <- round(foldchange(ct_enrich_df$Signif_percent, ct_enrich_df$All_percent), 1)
ct_enrich_df$FC_G1 <- round(foldchange(ct_enrich_df$G1_percent, ct_enrich_df$All_percent), 1)
ct_enrich_df$FC_G2 <- round(foldchange(ct_enrich_df$G2_percent, ct_enrich_df$All_percent), 1)
ct_enrich_df$FC_G3 <- round(foldchange(ct_enrich_df$G3_percent, ct_enrich_df$All_percent), 1)
ct_enrich_df$FC_G4 <- round(foldchange(ct_enrich_df$G4_percent, ct_enrich_df$All_percent), 1)
ct_enrich_df <- ct_enrich_df[,c(1,8:17)]
ct_enrich_df_subset <- ct_enrich_df[,c("CT", "FC_G1", "FC_G2", "FC_G3", "FC_G4")]
colnames(ct_enrich_df_subset) <- c("CT", "G1_up", "G2_upDown", "G3_down", "G4_downUp")

# prepare FC and annotation data for heatmap
n <- ct_enrich_df_subset$CT
# transpose all but the first column (name)
t_df <- as.data.frame(t(ct_enrich_df_subset[,-1]))
colnames(t_df) <- n
t_df$Group <- factor(row.names(t_df))
Annotation_col <- t_df$Group
Annotation_col <- as.data.frame(Annotation_col)
colnames(Annotation_col) <- "Group"
row.names(Annotation_col) <- Annotation_col$Group

# plotting
png("/Figures/AD_vs_CO_cell_type_composition_Heatmap_FC.png", units="mm", width=160, height=80, res=1000)
pheatmap(
  t(t_df[, -6]),
  annotation_col = Annotation_col, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE, 
  show_colnames = FALSE, 
  annotation_colors = list(Group = c("G1_up"="darkturquoise", "G2_upDown"="green4", "G3_down"="darkgoldenrod", "G4_downUp"="darkmagenta")),
  annotation_legend = TRUE,
  fontsize = 10, 
  fontsize_row = 10,
  border_color = NA,
  na_col = "grey88")
dev.off()