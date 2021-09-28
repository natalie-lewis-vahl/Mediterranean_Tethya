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
CyanobacteriaCol<-"#0072B2"
PlanctomycetesCol<-"#CAB2D6"
ProteobacteriaCol<-"#6A3D9A"

###########

uncorrectedCountsPath<-"../Data/cbas_tempVSctrl.otutab.csv"
bactLoadCorrectedCountsPath<-"./Data/cbas_tempVSctrl.otutab.bactloadCorrected_valuesOnly.csv"
otuTaxonomyPath<-"./Data/cbas_otu_taxonomy.csv"

otuTaxonomy<-read.csv(otuTaxonomyPath, sep="\t")

countsDF<-read.csv(bactLoadCorrectedCountsPath, sep="\t")

OTUPresence<-as.data.frame(countsDF[,-1]>0)#turns into logical value, true or false if present or not
OTUPresence<-cbind(countsDF[,1], OTUPresence)

dim(OTUPresence)
presentIn<-apply(OTUPresence[,-1], 1, sum)#For all columns except 1 count sum of rows/how many samples the bacteria is present in
OTUPresence<-cbind(OTUPresence,presentIn/length(OTUPresence[,-1]))

OTUPresenceWithTaxonomy<-bind_cols(OTUPresence[order(OTUPresence$`countsDF[, 1]`),], otuTaxonomy[order(otuTaxonomy$sequence_identifier),])
coreCommunity<-OTUPresenceWithTaxonomy[OTUPresenceWithTaxonomy$`presentIn/length(OTUPresence[, -1])` == 1.00,]
dim(coreCommunity)

#get rid of unclassified OTUs
coreCommunity<-coreCommunity[coreCommunity$Phylum != "",]

#countSum<-apply(countsDF[-1],1,sum)
#countsDF<-cbind(countsDF, countSum)
countsWithTaxonomy<-bind_cols(countsDF[order(countsDF$X.OTU.ID),], otuTaxonomy[order(otuTaxonomy$sequence_identifier),])

#sum(countsDF[countsDF$X.OTU.ID %in% coreCommunity$`countsDF[, 1]`,]$countSum)/sum(countsDF$countSum)

coreCountsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$X.OTU.ID %in% coreCommunity$`countsDF[, 1]`,]
corePropsWithTaxonomy<-coreCountsWithTaxonomy
corePropsWithTaxonomy[,2:16]<-apply(coreCountsWithTaxonomy[,2:16],2,function(x) x/sum(x))

corePropsWithTaxonomyMelted<-melt(corePropsWithTaxonomy)

levels(corePropsWithTaxonomyMelted$variable)<-c("LauraG1C1S15","LauraG1C2S19","LauraG1C3S23","LauraG2C1S16","LauraG2C2S20","LauraG2C3S24",
                                                "LauraP1C1S13","LauraP1C2S17","LauraP1C3S21","LauraP2C1S14","LauraP2C2S18","LauraP2C3S22",
                                                "LauraG1T2S7","LauraG1T3S11","LauraG2T3S12","LauraP1T1S1","LauraP1T2S5","LauraP1T3S9",
                                                "LauraP2T2S6","LauraP2T3S10")

levels(corePropsWithTaxonomyMelted$variable)<-c("G1C1","G1C2","G1C3","G2C1","G2C2","G2C3","P1C1","P1C2","P1C3","P2C1","P2C2","P2C3",
                                                "G1T2","G1T3","G2T3","P1T1","P1T2","P1T3","P2T2","P2T3")

levels(corePropsWithTaxonomyMelted$variable)<-c("C","C","C","C","C","C","C","C","C","C","C","C",
                                                "T","T","T","T","T","T","T","T")

ggplot(corePropsWithTaxonomyMelted, aes(variable, y=value)) + geom_bar(aes(fill=Phylum), stat="identity", position="stack") + scale_fill_manual(values = c(AcidobacteriaCol, ActinobacteriaCol, BacteroidetesCol, CyanobacteriaCol, PlanctomycetesCol, ProteobacteriaCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))

ggplot(corePropsWithTaxonomyMelted, aes(X.OTU.ID, y=value)) + geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(variable)) + scale_fill_manual(values = c(AcidobacteriaCol, ActinobacteriaCol, BacteroidetesCol, CyanobacteriaCol, PlanctomycetesCol, ProteobacteriaCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))

ggplot(corePropsWithTaxonomyMelted[corePropsWithTaxonomyMelted$X.OTU.ID %in% c("OTU_1","OTU_2","OTU_3","OTU_4","OTU_8"),], aes(X.OTU.ID, y=value)) + geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(variable)) + scale_fill_manual(values = c(BacteroidetesCol, CyanobacteriaCol, ProteobacteriaCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))
ggplot(corePropsWithTaxonomyMelted[corePropsWithTaxonomyMelted$X.OTU.ID %in% c("OTU_1","OTU_2","OTU_3","OTU_4","OTU_8"),], aes(variable, y=value)) + geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(X.OTU.ID)) + scale_fill_manual(values = c(BacteroidetesCol, CyanobacteriaCol, ProteobacteriaCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))

ggplot(corePropsWithTaxonomyMelted[!(corePropsWithTaxonomyMelted$X.OTU.ID %in% c("OTU_1","OTU_2","OTU_3","OTU_4","OTU_8")),], aes(X.OTU.ID, y=value)) + geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(variable)) + scale_fill_manual(values = c(AcidobacteriaCol, ActinobacteriaCol, BacteroidetesCol, CyanobacteriaCol, PlanctomycetesCol, ProteobacteriaCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))
ggplot(corePropsWithTaxonomyMelted[!(corePropsWithTaxonomyMelted$X.OTU.ID %in% c("OTU_1","OTU_2","OTU_3","OTU_4","OTU_8")),], aes(variable, y=value)) + geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(X.OTU.ID)) + scale_fill_manual(values = c(AcidobacteriaCol, ActinobacteriaCol, BacteroidetesCol, CyanobacteriaCol, PlanctomycetesCol, ProteobacteriaCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))


#core
ggplot(corePropsWithTaxonomyMelted, aes(X.OTU.ID, y=value)) + geom_bar(aes(fill=Phylum), stat="identity") + scale_fill_manual(values = c(AcidobacteriaCol, ActinobacteriaCol, BacteroidetesCol, CyanobacteriaCol, PlanctomycetesCol, ProteobacteriaCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))




