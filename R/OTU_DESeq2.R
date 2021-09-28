library("DESeq2")
library(ggplot2)
library("BiocParallel")
register(MulticoreParam(8))


###################
#
# Colors
# 
# Assign colors to the different phyla to make figures consistent.
#
#################################

AcidobacteriaCol<-"#999999"
ActinobacteriaCol<-"#E69F00"
BacteroidetesCol<-"#56B4E9"
CyanobacteriaCol<-"#0072B2"
PlanctomycetesCol<-"#CAB2D6"
ProteobacteriaCol<-"#6A3D9A"

###########

otuTaxonomyPath<-"../Data/cbas_otu_taxonomy.csv"
input_counts<-"../Data/cbas_tempVSctrl.otutab.csv"
input_sample_information<-"../Data/cbas_tempVSctrl.otutab.info"

otuTaxonomy<-read.csv(otuTaxonomyPath, sep="\t")
OTU_DE<-read.csv(input_counts, head=T, sep="\t")
OTU_DE_INFO<-read.csv(input_sample_information, head=T)

core<-c("OTU_1","OTU_11","OTU_13","OTU_18","OTU_2","OTU_21","OTU_3","OTU_388","OTU_4","OTU_47","OTU_48"
  ,"OTU_5","OTU_50","OTU_59","OTU_6","OTU_8","OTU_89") 


coreOTUTaxonomy<-otuTaxonomy[otuTaxonomy$sequence_identifier %in% core,]
OTU_DE<-OTU_DE[OTU_DE$X.OTU.ID %in% core,]


row.names(OTU_DE)<-OTU_DE$X.OTU.ID
OTU_DE$X.OTU.ID<-NULL

#Remove genes that have no counts over all samples
OTU_DE <- floor(OTU_DE)

OTU_DeSeq<-DESeqDataSetFromMatrix(countData=OTU_DE, colData=OTU_DE_INFO, design=~condition)

OTU_DeSeq$condition<-relevel(OTU_DeSeq$condition,ref="Control")

#estimate size factors
#Thus, if all size factors are roughly equal to one, the libraries have been sequenced equally deeply.
OTU_DeSeq<-estimateSizeFactors(OTU_DeSeq)
OTU_DeSeq<-estimateDispersions(OTU_DeSeq, fitType="local")

#wald test
OTU_DeSeq<-nbinomWaldTest(OTU_DeSeq)
OTU_DeSeq_Results<-results(OTU_DeSeq, pAdjustMethod = "BH")

OTU_DeSeq_Results[!is.na(OTU_DeSeq_Results$pvalue<0.05) & OTU_DeSeq_Results$pvalue<0.05,]

lcfDF<-data.frame(OTU=row.names(OTU_DeSeq_Results), BaseMean=OTU_DeSeq_Results$baseMean, LFC=OTU_DeSeq_Results$log2FoldChange, LFC_SE=OTU_DeSeq_Results$lfcSE)

lcfDF<-bind_cols(lcfDF[order(lcfDF$OTU),], coreOTUTaxonomy[order(coreOTUTaxonomy$sequence_identifier),])

#lfc of out core by phylum
svg("../Plots/logFoldChangeBYOTU.svg")
#pdf("../Plots/logFoldChangeBYOTU.pdf")
ggplot(lcfDF, aes(x=reorder(OTU, BaseMean), y=LFC)) + geom_bar(aes(fill=Phylum), stat="identity") + scale_fill_manual(values = c(AcidobacteriaCol, ActinobacteriaCol, BacteroidetesCol, CyanobacteriaCol, PlanctomycetesCol, ProteobacteriaCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))
dev.off()

svg("../Plots/baseMeanBYOTU.svg")
#pdf("../Plots/baseMeanBYOTU.pdf")
ggplot(lcfDF, aes(x=reorder(OTU, BaseMean), y=BaseMean)) + geom_bar(aes(fill=Phylum), stat="identity") + scale_fill_manual(values = c(AcidobacteriaCol, ActinobacteriaCol, BacteroidetesCol, CyanobacteriaCol, PlanctomycetesCol, ProteobacteriaCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1)) + scale_y_log10()
dev.off()

#correlations
cor.test(log10(OTU_DeSeq_Results$baseMean),OTU_DeSeq_Results$log2FoldChange)
cor.test(log10(OTU_DeSeq_Results$baseMean[-16]),OTU_DeSeq_Results$log2FoldChange[-16])

cor.test(log10(OTU_DeSeq_Results$baseMean),OTU_DeSeq_Results$log2FoldChange, method = "spearman")
cor.test(log10(OTU_DeSeq_Results$baseMean[-16]),OTU_DeSeq_Results$log2FoldChange[-16], method = "spearman")

#lm
logBYMean<-lm(OTU_DeSeq_Results$log2FoldChange~log10(OTU_DeSeq_Results$baseMean))
anova(logBYMean)
svg("../Plots/regressionCooksDbyLeverage.svg")
#pdf("../Plots/regressionCooksDbyLeverage.pdf")
plot(logBYMean,6)
dev.off()

logBYMeanNo16<-lm(OTU_DeSeq_Results$log2FoldChange[-16]~log10(OTU_DeSeq_Results$baseMean[-16]))
anova(logBYMeanNo16)

svg("../Plots/logFoldChangeBYbaseMean.svg")
plot(OTU_DeSeq_Results$log2FoldChange~log10(OTU_DeSeq_Results$baseMean), ylim=c(-2,3.5))
abline(h=0, col="plum", lwd=3)
abline(logBYMean$coefficients, col="slateblue", lwd=3)
abline(logBYMeanNo16$coefficients, col="goldenrod1", lwd=3)
dev.off()

svg("../Plots/logFoldChangeBYbaseMeanNo16.svg")
plot(OTU_DeSeq_Results$log2FoldChange[-16]~log10(OTU_DeSeq_Results$baseMean[-16]), ylim=c(-2,3.5))
abline(h=0, col="plum", lwd=3)
abline(logBYMean$coefficients, col="slateblue", lwd=3)
abline(logBYMeanNo16$coefficients, col="goldenrod1", lwd=3)
dev.off()




