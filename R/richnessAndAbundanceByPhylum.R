#!/usr/bin/env Rscript

###########################################
#
# abundanceModels.R
#
# Copyright (C) 2019 Sergio Vargas
#
# Contact: sergio.vargas@lmu.de
#
# This program is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation,
# either version 3 of the License, or any later version.
#
# This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program.
#
# If not, see <https://www.gnu.org/licenses/>
#
################################################################


###########################################
#
# Usage:
#
# This script is a companion for Vargas and Leiva.
# 
# Set the wd to the R folder provided with the package
# once this has been done sourcing should work
#
#
################################################################

################
#Load required libraries
###########################

library(ggplot2)
library(dplyr)
library(reshape2)

######################
#Private functions
#
####################################################


################
#Input/Output paths
###########################

#bactLoadCorrectedCountsPath<-"../Data/cbas_tempVSctrl.otutab.bactloadCorrected_valuesOnly.csv"
uncorrectedCountsPath<-"../Data/cbas_tempVSctrl.otutab.csv"
otuTaxonomyPath<-"../Data/cbas_otu_taxonomy.csv"

plotsPath<-"../Plots/"

################
#User defined variables
###########################

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

plotSVG<-TRUE
plotPDF<-FALSE

################
#Main Source
###########################

countsDF<-read.csv(uncorrectedCountsPath, sep="\t")
otuTaxonomy<-read.csv(otuTaxonomyPath, sep="\t")

#get OTU counts and add them to DF
countSum<-apply(countsDF[-1],1,sum)
countsDF<-cbind(countsDF, countSum)

#bind the taxonomy columns to the otu counts DF
countsWithTaxonomy<-bind_cols(countsDF[order(countsDF$X.OTU.ID),], otuTaxonomy[order(otuTaxonomy$sequence_identifier),])

#remove OTUs with sum zero if any. Uncomment dim lines if check wanted.
#dim(countsWithTaxonomy)
countsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$countSum != 0,]
#dim(countsWithTaxonomy)


#######################
#
# OTU Richness by Phylum
#
#
###########################################

#get a table with the number of OTUs annotated to each phylum
phylumRichness<-table(countsWithTaxonomy$Phylum)

allPhylaRichnessPie<-ggplot(as.data.frame(phylumRichness[-1]), aes(x="", y=Freq, fill=Var1)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = allBacteriaColorVector) + theme_bw()


#richness pie chart

if(plotPDF){
  
  richnessPiePath<-paste(plotsPath,"allPhylaRichnessPie.pdf", sep="")
  pdf(richnessPiePath)
  print(allPhylaRichnessPie)
  dev.off()
  
}

if(plotSVG){
  
  richnessPiePath<-paste(plotsPath,"allPhylaRichnessPie.svg", sep="")
  svg(richnessPiePath)
  print(allPhylaRichnessPie)
  dev.off()
  
}

#split in low and high richness
highRichnessPhyla<-phylumRichness[phylumRichness>10]
highRichnessPhylaPie<-ggplot(as.data.frame(highRichnessPhyla[-1]), aes(x="", y=Freq, fill=Var1)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = c(AcidobacteriaCol, ActinobacteriaCol, BacteroidetesCol, ChloroflexiCol, CyanobacteriaCol, FirmicutesCol, PlanctomycetesCol, ProteobacteriaCol, VerrucomicrobiaCol)) + theme_bw()

lowRichnessPhyla<-phylumRichness[phylumRichness<=10]
lowRichnessPhylaPie<-ggplot(as.data.frame(lowRichnessPhyla), aes(x="", y=Freq, fill=Var1)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = c(ChlamydiaeCol, DadabacteriaCol, DeferribacteresCol, DeinococcusThermusCol, ElusimicrobiaCol, GemmatimonadetesCol, HydrogenedentesCol, LatescibacteriaCol, NitrospinaeCol, NitrospiraeCol, SpirochaetesCol, ThaumarchaeotaCol, WPS2Col)) + theme_bw()

if(plotPDF){
  
  highRichnessPhylaPiePath<-paste(plotsPath,"highRichnessPhylaPie.pdf", sep="")
  pdf(highRichnessPhylaPiePath)
  print(highRichnessPhylaPie)
  dev.off()
  
  lowRichnessPhylaPiePath<-paste(plotsPath,"lowRichnessPhylaPie.pdf", sep="")
  pdf(lowRichnessPhylaPiePath)
  print(lowRichnessPhylaPie)
  dev.off()
  
}

if(plotSVG){
  
  highRichnessPhylaPiePath<-paste(plotsPath,"highRichnessPhylaPie.svg", sep="")
  svg(highRichnessPhylaPiePath)
  print(highRichnessPhylaPie)
  dev.off()
  
  lowRichnessPhylaPiePath<-paste(plotsPath,"lowRichnessPhylaPie.svg", sep="")
  svg(lowRichnessPhylaPiePath)
  print(lowRichnessPhylaPie)
  dev.off()
  
}

#############################################
#
#
# abundance per phylum
#
#
#########################################################################

#get abundance by phylum
phylumCounts<-countsWithTaxonomy %>% group_by(Phylum) %>% summarise(sum=sum(countSum))

allPhylaAbundancePie<-ggplot(phylumCounts[-1,], aes(x="", y=sum, fill=Phylum)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = allBacteriaColorVector) + theme_bw()

if(plotPDF){
  
  allPhylaAbundancePiePath<-paste(plotsPath,"allPhylaAbundancePie.pdf", sep="")
  pdf(allPhylaAbundancePiePath)
  print(allPhylaAbundancePie)
  dev.off()
  
}

if(plotSVG){
  
  allPhylaAbundancePiePath<-paste(plotsPath,"allPhylaAbundancePie.svg", sep="")
  svg(allPhylaAbundancePiePath)
  print(allPhylaAbundancePie)
  dev.off()
  
}

#split in high and low abundance
highCountsPhyla<-phylumCounts[phylumCounts$sum > 1000,]
lowCountsPhyla<-phylumCounts[phylumCounts$sum <= 1000,]

#plot
highAbundancePhylaPie<-ggplot(highCountsPhyla[-1,], aes(x="", y=sum, fill=Phylum)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = c(AcidobacteriaCol, ActinobacteriaCol, BacteroidetesCol, CyanobacteriaCol, FirmicutesCol, PlanctomycetesCol, ProteobacteriaCol, VerrucomicrobiaCol)) + theme_bw()
lowAbundancePhylaPie<-ggplot(lowCountsPhyla, aes(x="", y=sum, fill=Phylum)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = c(ChlamydiaeCol, ChloroflexiCol, DadabacteriaCol, DeferribacteresCol, DeinococcusThermusCol, ElusimicrobiaCol, GemmatimonadetesCol, HydrogenedentesCol, LatescibacteriaCol, NitrospinaeCol, NitrospiraeCol, SpirochaetesCol, ThaumarchaeotaCol, WPS2Col)) + theme_bw()

if(plotPDF){
  
  highAbundancePhylaPiePath<-paste(plotsPath,"highAbundancePhylaPie.pdf", sep="")
  pdf(highAbundancePhylaPiePath)
  print(highAbundancePhylaPie)
  dev.off()
  lowAbundancePhylaPiePath<-paste(plotsPath,"lowAbundancePhylaPie.pdf", sep="")
  pdf(lowAbundancePhylaPiePath)
  print(lowAbundancePhylaPie)
  dev.off()
}

if(plotSVG){
  
  highAbundancePhylaPiePath<-paste(plotsPath,"highAbundancePhylaPie.svg", sep="")
  svg(highAbundancePhylaPiePath)
  print(highAbundancePhylaPie)
  dev.off()
  lowAbundancePhylaPiePath<-paste(plotsPath,"lowAbundancePhylaPie.svg", sep="")
  svg(lowAbundancePhylaPiePath)
  print(lowAbundancePhylaPie)
  dev.off()
  
}


#EOF




