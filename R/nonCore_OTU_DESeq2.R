library("DESeq2")
library(ggplot2)
library("BiocParallel")
register(MulticoreParam(8))
library(dplyr)


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
ChlamydiaeCol<-"#009E73"
ChloroflexiCol<-"#F0E442"
CyanobacteriaCol<-"#0072B2"
DadabacteriaCol<-"#D55E00"
DeferribacteresCol<-"#CC79A7"
DeinococcusThermusCol<-"#A6CEE3"
ElusimicrobiaCol<-"#1F78B4"
FirmicutesCol<-"#B2DF8A"
GemmatimonadetesCol<-"#33A02C"
HydrogenedentesCol<-"#FB9A99"
LatescibacteriaCol<-"#E31A1C"
NitrospinaeCol<-"#FDBF6F"
NitrospiraeCol<-"#FF7F00"
PlanctomycetesCol<-"#CAB2D6"
ProteobacteriaCol<-"#6A3D9A"
SpirochaetesCol<-"#8DD3C7"
ThaumarchaeotaCol<-"#FFFFB3"
VerrucomicrobiaCol<-"#FB8072"
WPS2Col<-"#999999"

allBacteriaColorVector<-c(AcidobacteriaCol,ActinobacteriaCol,BacteroidetesCol, ChlamydiaeCol,
                          ChloroflexiCol, CyanobacteriaCol, DadabacteriaCol, DeferribacteresCol,
                          DeinococcusThermusCol, ElusimicrobiaCol, FirmicutesCol, GemmatimonadetesCol,
                          HydrogenedentesCol, LatescibacteriaCol, NitrospinaeCol, NitrospiraeCol,
                          PlanctomycetesCol, ProteobacteriaCol, SpirochaetesCol, ThaumarchaeotaCol,
                          VerrucomicrobiaCol, WPS2Col)


###########

otuTaxonomyPath<-"./Data/cbas_otu_taxonomy.csv"
input_counts<-"./Data/cbas_tempVSctrl.otutab.csv"
input_sample_information<-"./Data/cbas_tempVSctrl.otutab.info"

otuTaxonomy<-read.csv(otuTaxonomyPath, sep="\t")
OTU_DE<-read.csv(input_counts, head=T, sep="\t")
OTU_DE_INFO<-read.csv(input_sample_information, head=T)

core<-c("OTU_1","OTU_11","OTU_13","OTU_18","OTU_2","OTU_21","OTU_3","OTU_388","OTU_4","OTU_47","OTU_48"
        ,"OTU_5","OTU_50","OTU_59","OTU_6","OTU_8","OTU_89") 


coreOTUTaxonomy<-otuTaxonomy[!otuTaxonomy$sequence_identifier %in% core,]
OTU_DE<-OTU_DE[!OTU_DE$X.OTU.ID %in% core,]

withPhylumData<-coreOTUTaxonomy$Phylum != ""
OTU_DE<-OTU_DE[withPhylumData,]
coreOTUTaxonomy<-coreOTUTaxonomy[coreOTUTaxonomy$Phylum != "",]

row.names(OTU_DE)<-OTU_DE$X.OTU.ID
OTU_DE$X.OTU.ID<-NULL
length(OTU_DE)
#Remove genes that have no counts over all samples
OTU_DE <- floor(OTU_DE)

OTU_DeSeq<-DESeqDataSetFromMatrix(countData=OTU_DE, colData=OTU_DE_INFO, design=~condition)
??DESeqDataSetFromMatrix
OTU_DeSeq$condition<-relevel(OTU_DeSeq$condition,ref="Control")

#estimate size factors
#Thus, if all size factors are roughly equal to one, the libraries have been sequenced equally deeply.
??estimateSizeFactors
OTU_DeSeq<-estimateSizeFactors(OTU_DeSeq)
OTU_DeSeq<-estimateDispersions(OTU_DeSeq)

#wald test
#??nbinomWaldTest 
#tests for significance of coefficients in a Negative Binomial GLM, using previously calculated sizeFactors (or normalizationFactors) and dispersion estimates
OTU_DeSeq<-nbinomWaldTest(OTU_DeSeq)
OTU_DeSeq_Results<-results(OTU_DeSeq, pAdjustMethod = "BH") #extracts a result table from a DESeq analysis giving base means across samples
# ?p.adjust #"BH" a less conservative p adjustemnt method by Benjamini & Hochberg (1995)

OTU_DeSeq_Results[!is.na(OTU_DeSeq_Results$pvalue<0.05) & OTU_DeSeq_Results$pvalue<0.05,]#search - ¿ exclude missing vals and combine significant points)
#A fold change describes the ratio of two values
#if there is a two fold increase (fold change=2, Log2FC=1) between A and B, then A is twice as big as B (or A is 200% of B).
#If there is a two fold decrease (fold change = 0.5, Log2FC= -1) between A and B, then A is half as big as B (or B is twice as big as A, or A is 50% of B).

lcfDF<-data.frame(OTU=row.names(OTU_DeSeq_Results), BaseMean=OTU_DeSeq_Results$baseMean, LFC=OTU_DeSeq_Results$log2FoldChange, LFC_SE=OTU_DeSeq_Results$lfcSE)

lcfDF<-bind_cols(lcfDF[order(lcfDF$OTU),], coreOTUTaxonomy[order(coreOTUTaxonomy$sequence_identifier),])

#lfc of out core by phylum
#svg("../Plots/logFoldChangeBYOTU.svg")
#pdf("../Plots/logFoldChangeBYOTU.pdf")

ggplot(lcfDF, aes(x=BaseMean, y=LFC)) + geom_point(aes(colour=Phylum)) + scale_color_manual(values = allBacteriaColorVector) + theme_bw() + scale_x_log10()

#dev.off()

