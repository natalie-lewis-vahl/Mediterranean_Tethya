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
library(vegan)
library(reshape2)

######################
#Private functions
#####################################################

#turn counts into presence absence
presentOTU<-function(x){
  
  return((x > 0)*1) 
}

################
#Input/Output paths
###########################

#bactLoadCorrectedCountsPath<-"../Data/cbas_tempVSctrl.otutab.bactloadCorrected_valuesOnly.csv"
uncorrectedCountsPath<-"./Data/cbas_tempVSctrl.otutab.csv"
otuTaxonomyPath<-"./Data/cbas_otu_taxonomy.csv"

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
dim(countsWithTaxonomy)

#################################
#
# Richness per Phylum per sample plots
#
#######################################################################

#turn abundance into richness and count otus by phylum, melt and plot if required (TOOK OUT FILTER FOR UNCLASSIFIED %>% filter(Phylum != ""))
otuCountByPhylumBySample<-countsWithTaxonomy %>% group_by(Phylum)  %>% mutate_at(colnames(countsWithTaxonomy)[2:21], presentOTU) %>% mutate_at(colnames(countsWithTaxonomy)[2:21],sum) %>% distinct_at(vars(colnames(countsWithTaxonomy)[c(2:21,25)])) 
countsWithTaxonomy
otuCountByPhylumBySampleMelted<-melt(otuCountByPhylumBySample)

otuCountByPhylumBySampleMelted<-otuCountByPhylumBySampleMelted%>% 
  filter(Phylum != "")
levels(otuCountByPhylumBySampleMelted$variable)<-c("LauraG1C1S15","LauraG1C2S19","LauraG1C3S23","LauraG2C1S16","LauraG2C2S20","LauraG2C3S24","LauraP1C1S13","LauraP1C2S17","LauraP1C3S21","LauraP2C1S14","LauraP2C2S18","LauraP2C3S22","LauraG1T2S7","LauraG1T3S11","LauraG2T3S12","LauraP1T1S1","LauraP1T2S5","LauraP1T3S9","LauraP2T2S6","LauraP2T3S10")

richnessByPhylumBySampleBarPlot<-ggplot(otuCountByPhylumBySampleMelted, aes(x=variable, y=value)) + geom_bar(position = "stack", stat="identity", aes(fill=Phylum)) + scale_fill_manual(values = allBacteriaColorVector) + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))

if(plotPDF){
  
  richnessByPhylumBySampleBarPlotPath<-paste(plotsPath,"richnessByPhylumBySampleBarPlot.pdf",sep="")
  pdf(richnessByPhylumBySampleBarPlotPath)
  print(richnessByPhylumBySampleBarPlot)
  dev.off()

}

if(plotSVG){
  
  richnessByPhylumBySampleBarPlotPath<-paste(plotsPath,"richnessByPhylumBySampleBarPlot.svg",sep="")
  svg(richnessByPhylumBySampleBarPlotPath)
  print(richnessByPhylumBySampleBarPlot)
  dev.off()
  
}



####################################
#
# Abundance per Phylum per sample plots
#
#####################################################

abundanceByPhylumBySample<-countsWithTaxonomy %>% group_by(Phylum) %>% filter(Phylum != "") %>% mutate_at(colnames(countsWithTaxonomy)[2:21],sum) %>% distinct_at(vars(colnames(countsWithTaxonomy)[c(2:21,25)])) 

abundanceByPhylumBySample<-ungroup(abundanceByPhylumBySample) %>% mutate_at(colnames(abundanceByPhylumBySample)[1:20],function(x)x/sum(x))

abundanceByPhylumBySampleMelted<-melt(abundanceByPhylumBySample)
levels(abundanceByPhylumBySampleMelted$variable)<-c("LauraG1C1S15","LauraG1C2S19","LauraG1C3S23","LauraG2C1S16","LauraG2C2S20","LauraG2C3S24","LauraP1C1S13","LauraP1C2S17","LauraP1C3S21","LauraP2C1S14","LauraP2C2S18","LauraP2C3S22","LauraG1T2S7","LauraG1T3S11","LauraG2T3S12","LauraP1T1S1","LauraP1T2S5","LauraP1T3S9","LauraP2T2S6","LauraP2T3S10")

abundanceByPhylumBySampleBarPlot<-ggplot(abundanceByPhylumBySampleMelted, aes(x=variable, y=value)) + geom_bar(position = "stack", stat="identity", aes(fill=Phylum)) + scale_fill_manual(values = allBacteriaColorVector) + theme_bw()  + theme(axis.text.x=element_text(angle=90, hjust=1))


if(plotPDF){
  
  abundanceByPhylumBySampleBarPlotPath<-paste(plotsPath,"abundanceByPhylumBySampleBarPlot.pdf",sep="")
  pdf(abundanceByPhylumBySampleBarPlotPath)
  print(abundanceByPhylumBySampleBarPlot)
  dev.off()
  
}

if(plotSVG){
  
  abundanceByPhylumBySampleBarPlotPath<-paste(plotsPath,"abundanceByPhylumBySampleBarPlot.svg",sep="")
  svg(abundanceByPhylumBySampleBarPlotPath)
  print(abundanceByPhylumBySampleBarPlot)
  dev.off()
  
}

#EOF