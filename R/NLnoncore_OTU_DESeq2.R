library(RColorBrewer)
library("DESeq2")
library(ggplot2)
library("BiocParallel")
register(MulticoreParam(8))

uncorrectedCountsPath<-"./Data/all.otutab_raw.csv"
#tbc bactLoadCorrectedCountsPath<-"./Data/all.otutab_Corrected_valuesOnly.csv"
otuTaxonomyPath<-"./Data/all_taxa.csv"
input_sample_information<-"./Data/NL_OTU_DE_INFO.csv"

otuTaxonomy<-read.csv(otuTaxonomyPath, sep=";")
OTU_DE<-read.csv(uncorrectedCountsPath, sep="\t")
OTU_DE_INFO<-read.csv(input_sample_information, head=T)

str(otuTaxonomy$Phylum) #PHylum has 34 levels


coreOTUs<-c("OTU_1", "OTU_10","OTU_102","OTU_105","OTU_1068", "OTU_1074", "OTU_1095", "OTU_11","OTU_110","OTU_12", "OTU_128", "OTU_129","OTU_13",
            "OTU_1379","OTU_14","OTU_151","OTU_16","OTU_17","OTU_170","OTU_179","OTU_19","OTU_2","OTU_20","OTU_2075","OTU_2078","OTU_21","OTU_210",
            "OTU_22", "OTU_235","OTU_24", "OTU_241","OTU_25","OTU_26","OTU_261","OTU_27","OTU_295","OTU_30","OTU_31", "OTU_32", "OTU_325","OTU_326",
            "OTU_33", "OTU_34", "OTU_35", "OTU_36", "OTU_374", "OTU_39", "OTU_41", "OTU_42", "OTU_43", "OTU_455",  "OTU_500", 
            "OTU_504","OTU_525","OTU_53", "OTU_532","OTU_54", "OTU_56","OTU_57","OTU_60","OTU_64","OTU_66", "OTU_7","OTU_716","OTU_718",
            "OTU_73", "OTU_8","OTU_809","OTU_81", "OTU_85", "OTU_86", "OTU_87", "OTU_9","OTU_907","OTU_92", "OTU_931")
noncoreOTUTaxonomy<-otuTaxonomy[!otuTaxonomy$sequence_identifier %in% coreOTUs,]
noncoreOTU_DE<-OTU_DE[!OTU_DE$X.OTU.ID %in% coreOTUs,]

noncorewithPhylumData<-noncoreOTUTaxonomy$Phylum != ""
noncoreOTU_DE<-noncoreOTU_DE[noncorewithPhylumData,]
noncoreOTUTaxonomy<-noncoreOTUTaxonomy[noncoreOTUTaxonomy$Phylum != "",]
#MAke ID name row names rather than column
row.names(noncoreOTU_DE)<-noncoreOTU_DE$X.OTU.ID
noncoreOTU_DE$X.OTU.ID<-NULL
#HOw many phyla
length(unique(noncoreOTUTaxonomy$Phylum))#33
help(RColorBrewer)
cols<-brewer.pal(12, "Set3")
cols<- colorRampPalette(x)(33)
names(cols)<-levels(noncoreOTUTaxonomy$Phylum)
colScale <- scale_colour_manual(name = "Phylum",values = cols)
###
#FOr adding colour scale later
#One plot with all the data
#p <- ggplot(dat,aes(x,y,colour = grp)) + geom_point()
#p1 <- p + colScale

#A second plot with only four of the levels
#p2 <- p %+% droplevels(subset(dat[4:10,])) + colScale
####
#Remove genes that have no counts over all samples
noncoreOTU_DE <- floor(noncoreOTU_DE)
??floor
noncoreOTU_DeSeq<-DESeqDataSetFromMatrix(countData=noncoreOTU_DE, colData=OTU_DE_INFO, design=~species)

#OTU_DeSeq$species<-relevel(OTU_DeSeq$condition,ref="Control")

#estimate size factors
#Thus, if all size factors are roughly equal to one, the libraries have been sequenced equally deeply.
noncoreOTU_DeSeq<-estimateSizeFactors(noncoreOTU_DeSeq)
noncoreOTU_DeSeq<-estimateDispersions(OTU_DeSeq)

#wald test
OTU_DeSeq<-nbinomWaldTest(OTU_DeSeq)
OTU_DeSeq_Results<-results(OTU_DeSeq, pAdjustMethod = "BH")

OTU_DeSeq_Results[!is.na(OTU_DeSeq_Results$pvalue<0.05) & OTU_DeSeq_Results$pvalue<0.05,]

lcfDF<-data.frame(OTU=row.names(OTU_DeSeq_Results), BaseMean=OTU_DeSeq_Results$baseMean, LFC=OTU_DeSeq_Results$log2FoldChange, LFC_SE=OTU_DeSeq_Results$lfcSE)

lcfDF<-bind_cols(lcfDF[order(lcfDF$OTU),], coreOTUTaxonomy[order(coreOTUTaxonomy$sequence_identifier),])

#lfc of out core by phylum
#svg("../Plots/logFoldChangeBYOTU.svg")
#pdf("../Plots/logFoldChangeBYOTU.pdf")

ggplot(lcfDF, aes(x=BaseMean, y=LFC)) + geom_point(aes(colour=Phylum)) + scale_color_manual(values = allBacteriaColorVector) + theme_bw() + scale_x_log10()

#dev.off()

