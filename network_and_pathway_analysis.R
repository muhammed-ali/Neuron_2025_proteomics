rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(2022)

library(igraph)
library(enrichplot)
library(clusterProfiler)


##
## Generate input files for network analysis and visualization
##

DEA <- read.table("/Input_Data/Metaanalysis_Bonf_PsuedoTrajctories.txt", sep="\t", header=T, stringsAsFactors=F)
dim(DEA) # 2173   15

table(DEA$group)
#    g1_up g2_upDown   g3_down g4_downUp
#      471       482       184      1036

#
# g1_up
#

g1 <- DEA[DEA$group == "g1_up",]
dim(g1) # 471  15
length(unique(g1$Analyte)) # 471
length(unique(g1$Gene)) # 454
g1$Expression <- ifelse(g1$Estimate.ADCO > 0, 1, 0)
g1_expr <- g1[,c("Gene", "Expression")]
g1_expr <- unique(na.omit(g1_expr))
dim(g1_expr) # 454   2

dir.create("g1_up")
# write GRN input files
write.table(g1_expr, file="./g1_up/expression.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(g1_expr$Gene, file="./g1_up/geneList.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(g1_expr, file="./g1_up/N1_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
g1_expr_N2 <- g1_expr
g1_expr_N2$Expression <- ifelse(g1_expr_N2$Expression < 1, 1, 0)
write.table(g1_expr_N2, file="./g1_up/N2_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(g1_expr_N2, file="./g1_up/expression_N2.txt", sep="\t", row.names=F, quote=F, col.names=F)



#
# g2_upDown
#

g2 <- DEA[DEA$group == "g2_upDown",]
dim(g2) # 482  15
length(unique(g2$Analyte)) # 482
length(unique(g2$Gene)) # 475
g2$Expression <- ifelse(g2$Estimate.ADCO > 0, 1, 0)
g2_expr <- g2[,c("Gene", "Expression")]
g2_expr <- unique(na.omit(g2_expr))
dim(g2_expr) # 475   2
table(g2_expr$Expression)
#   0
# 475

dir.create("g2_upDown")
# write GRN input files
write.table(g2_expr, file="./g2_upDown/expression.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(g2_expr$Gene, file="./g2_upDown/geneList.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(g2_expr, file="./g2_upDown/N1_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
g2_expr_N2 <- g2_expr
g2_expr_N2$Expression <- ifelse(g2_expr_N2$Expression < 1, 1, 0)
write.table(g2_expr_N2, file="./g2_upDown/N2_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(g2_expr_N2, file="./g2_upDown/expression_N2.txt", sep="\t", row.names=F, quote=F, col.names=F)


#
# g3_down
#

g3 <- DEA[DEA$group == "g3_down",]
dim(g3) # 184  15
length(unique(g3$Analyte)) # 184
length(unique(g3$Gene)) # 182
g3$Expression <- ifelse(g3$Estimate.ADCO > 0, 1, 0)
g3_expr <- g3[,c("Gene", "Expression")]
g3_expr <- unique(na.omit(g3_expr))
dim(g3_expr) # 182   2


dir.create("g3_down")
# write GRN input files
write.table(g3_expr, file="./g3_down/expression.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(g3_expr$Gene, file="./g3_down/geneList.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(g3_expr, file="./g3_down/N1_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
g3_expr_N2 <- g3_expr
head(g3_expr)
g3_expr_N2$Expression <- ifelse(g3_expr_N2$Expression < 1, 1, 0)
head(g3_expr_N2)
write.table(g3_expr_N2, file="./g3_down/N2_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(g3_expr_N2, file="./g3_down/expression_N2.txt", sep="\t", row.names=F, quote=F, col.names=F)


#
# g4_downUp
#

g4 <- DEA[DEA$group == "g4_downUp",]
dim(g4) # 1036   15
length(unique(g4$Analyte)) # 1036
length(unique(g4$Gene)) # 962
g4$Expression <- ifelse(g4$Estimate.ADCO > 0, 1, 0)
g4_expr <- g4[,c("Gene", "Expression")]
g4_expr <- unique(na.omit(g4_expr))
dim(g4_expr) # 962   2
table(g4_expr$Expression)
#   1
# 962

dir.create("g4_downUp")
# write GRN input files
write.table(g4_expr, file="./g4_downUp/expression.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(g4_expr$Gene, file="./g4_downUp/geneList.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(g4_expr, file="./g4_downUp/N1_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
g4_expr_N2 <- g4_expr
head(g4_expr)
g4_expr_N2$Expression <- ifelse(g4_expr_N2$Expression < 1, 1, 0)
head(g4_expr_N2)
write.table(g4_expr_N2, file="./g4_downUp/N2_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(g4_expr_N2, file="./g4_downUp/expression_N2.txt", sep="\t", row.names=F, quote=F, col.names=F)


# Metacore Filter applied are:
# Species = Homo Sapiens
# Interaction type = Binding, CrT, Regulation, influence on expression, transcriptional regulation
# Additional Filter: Functional + Binding interactions (also used for network building).


##
## Network reconstruction and perturbation analysis
##

cat('java -jar /JARs/Preprocessor.jar . nodemap.txt interactions.txt geneList.txt
java -jar /JARs/DifferentialNetworkAnalysis.jar expression.txt adjacency.txt GAResult.txt 0 true 1000 50 .
java -jar /JARs/CommonNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt CommonNetworkGenerator_Output.txt
java -jar /JARs/DifferentialNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt DifferentialNetworkGenerator_Output.txt
java -jar /JARs/ComputeCycles.jar CommonNetworkGenerator_Output.txt expression.txt pos.txt neg.txt
java -jar /JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype1.txt 1 SteadyStateCalculatorN1.txt
java -jar /JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype2.txt 2 SteadyStateCalculatorN2.txt
java -jar /JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN1.txt PerturbagenListGeneratorN1.txt
java -jar /JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN2.txt PerturbagenListGeneratorN2.txt
java -jar /JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 1 500000 BruteForcePerturbationsUpdatedN1_1.txt
java -jar /JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 2 500000 BruteForcePerturbationsUpdatedN1_2.txt
java -jar /JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 3 500000 BruteForcePerturbationsUpdatedN1_3.txt
java -jar /JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 4 500000 BruteForcePerturbationsUpdatedN1_4.txt
sort -rn BruteForcePerturbationsUpdatedN1_1.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_1.txt
sort -rn BruteForcePerturbationsUpdatedN1_2.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_2.txt
sort -rn BruteForcePerturbationsUpdatedN1_3.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_3.txt
sort -rn BruteForcePerturbationsUpdatedN1_4.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_4.txt' > run_network_analysis.sh)

# Strongly Connected Component (SCC) analysis for g1_up
n1 <- read.table("./g1_up/NetworkPhenotype1.txt", sep="\t", header=F)
dim(n1) # 1141    3

df <- n1[, -2]
head(df)	# tab separated file e.g GENE1	GENE2
g <- graph.data.frame(df)
# plot(g)
SCC <- clusters(g, mode="strong")
lcc <- induced.subgraph(g, V(g)[which(SCC$membership == which.max(SCC$csize))])
lcc
a <- V(g)[which(SCC$membership == which.max(SCC$csize))]
b <- as.matrix(a)
c <- data.frame(b)
g1_scc_prot <- row.names(c)

N1_filtered <- n1[n1$V1 %in% g1_scc_prot,]
dim(N1_filtered) # 863   3
N1_filtered <- N1_filtered[N1_filtered$V3 %in% g1_scc_prot,]
dim(N1_filtered) # 636   3
write.table(N1_filtered, file="./g1_up/NetworkPhenotype1_filtered.txt", sep="\t", row.names=F, quote=F, col.names=F)


##
## KEGG Enrichment of g1_Up_SCC
##

ADCO_significant <- read.table("/Input_Data/Metaanalysis_Bonf_PsuedoTrajctories_CT_info.txt", sep="\t", header=T, stringsAsFactors=F, quote="", fill=T)
dim(ADCO_significant) # 2173   16

CT_df <- as.data.frame(table(ADCO_significant$Cell_Type))
colnames(CT_df) <- c("cell_type", "CT_labels")
CT_df
#              cell_type CT_labels
# 1           Astrocytes       153
# 2          Endothelial        70
# 3 Microglia/Macrophage       182
# 4              Neurons       267
# 5     Oligodendrocytes        56
# 6              Unknown      1316

background <- length(unique(ADCO_significant$Gene)) # 2029

ADCO_significant_g1_filt <- ADCO_significant[ADCO_significant$Gene %in% g1_scc_prot,]
dim(ADCO_significant_g1_filt) # 168  16

ADCO_enrichKEGG_g1 <- enrichKEGG(
  unique(ADCO_significant_g1_filt$Entrez),
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  universe = as.character(unique(ADCO_complete$Entrez)),
  minGSSize = 3,
  maxGSSize = 1000,
  use_internal_data = FALSE
)
dim(ADCO_enrichKEGG_g1) # 62  9

library(scales) # for breaking pathway names
g1_pathway_fig <- as.data.frame(cbind(ADCO_enrichKEGG_g1@result$Description[1:10], ADCO_enrichKEGG_g1@result$p.adjust[1:10]))
g1_pathway_fig$V2 <- as.numeric(g1_pathway_fig$V2)
g1_pathway_fig <- g1_pathway_fig[order(g1_pathway_fig$V2),]
g1_pathway_fig$V1 <- factor(g1_pathway_fig$V1, levels=c(g1_pathway_fig$V1))
color_code <- c("red", "dodgerblue3", "green3", "purple", "orange", "yellowgreen", "peru", "palevioletred1", "skyblue", "black")
png("/Figures/AD_vs_CO_enrichKEGG_g1_SCC_barPlot.png",  units="mm", width=170, height=140, res=1000)
ggplot(g1_pathway_fig, aes(x=V1, y=V2)) + geom_bar(stat = "identity", fill=color_code) + coord_flip() + 
scale_y_continuous(expand = c(0, 0)) + theme(axis.text.y = element_text(face="bold", color="black", size=14), 
	axis.text.x = element_text(face="bold", color="black", size=10), axis.title.x = element_text(face="bold", color="black", size=14), 
	axis.title.y = element_text(face="bold", color="black", size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
	panel.background = element_blank()) + xlab("Pathways") + ylab("P (FDR)") + scale_x_discrete(labels = wrap_format(30))
dev.off()
write.csv(ADCO_enrichKEGG_g1, file="/Input_Data/AD_vs_CO_enrichKEGG_g1_SCC.csv", row.names=F, quote=F)

# Repeat the same network and pathway enrichment analyses for each protein group.
