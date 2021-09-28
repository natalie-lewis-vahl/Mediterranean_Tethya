library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(vegan)
library(reshape2)

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

uncorrectedCountsPath<-"../Data/cbas_tempVSctrl.otutab.csv"
bactLoadCorrectedCountsPath<-"./cbas_tempVSctrl.otutab.bactloadCorrected_valuesOnly.csv"
otuTaxonomyPath<-"../Data/cbas_otu_taxonomy.csv"

otuTaxonomy<-read.csv(otuTaxonomyPath, sep="\t")

countsDF<-read.csv(uncorrectedCountsPath, sep="\t")

OTUPresence<-as.data.frame(countsDF[,-1]>0)
OTUPresence<-cbind(countsDF[,1], OTUPresence)

head(OTUPresence)
controlSamples<-c(2:4,7:9,11:13,17:19)
treatmentSamples<-c(5:6,10,14:16,20:21)

head(OTUPresence[,controlSamples])
head(OTUPresence[,treatmentSamples])

#taxa gained, i.e. not in control but in at least one treated sample
absentInControls<-apply(OTUPresence[,controlSamples], 1, sum) == 0
table(absentInControls)
taxaGained<-otuTaxonomy[otuTaxonomy$sequence_identifier %in% OTUPresence[absentInControls,1],]
dim(taxaGained)
taxaGained<-taxaGained[taxaGained$Phylum != "",]
dim(taxaGained)
gainedOTUCountPerPhylum<-table(taxaGained$Phylum)

#taxa lost, i.e. not in treatment but in at least one control sample
absentInTreatments<-apply(OTUPresence[,treatmentSamples], 1, sum) == 0
table(absentInTreatments)
taxaLost<-otuTaxonomy[otuTaxonomy$sequence_identifier %in% OTUPresence[absentInTreatments,1],]
dim(taxaLost)
taxaLost<-taxaLost[taxaLost$Phylum != "",]
dim(taxaLost)
lostOTUCountPerPhylum<-table(taxaLost$Phylum)

gainedOTUCountPerPhylum
lostOTUCountPerPhylum
otuTotals<-table(otuTaxonomy[otuTaxonomy$Phylum != "",]$Phylum)
otuGain<-data.frame(otuTotalRichness = as.vector(otuTotals), gainedOTUs = as.vector(gainedOTUCountPerPhylum), lostOTUs = as.vector(lostOTUCountPerPhylum), Phylum = names(otuTotals))

otuGain<-otuGain[otuGain$Phylum != "",]
otuGain$relativeGainPerPhylum<-(otuGain$gainedOTUs-otuGain$lostOTUs)/otuGain$otuTotalRichness
otuGain$absoluteGainPerPhylum<-(otuGain$gainedOTUs-otuGain$lostOTUs)

ggplot(otuGain, aes(x=Phylum, y=absoluteGainPerPhylum)) + geom_bar(aes(fill=Phylum), stat="identity") + scale_fill_manual(values = allBacteriaColorVector) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))
ggplot(otuGain, aes(x=Phylum, y=relativeGainPerPhylum)) + geom_bar(aes(fill=Phylum), stat="identity") + scale_fill_manual(values = allBacteriaColorVector) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))


#some check points
#table(absentInTreatments)
#table(absentInControls)
#table(absentInTreatments | absentInControls)
#table(absentInTreatments & absentInControls)


  
  
  