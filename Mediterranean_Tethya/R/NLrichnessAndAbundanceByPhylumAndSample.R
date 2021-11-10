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
# This script is a modification of the companion script for Vargas and Leiva.
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
uncorrectedCountsPath<-"./Data/all.otutab_raw.csv"
otuTaxonomyPath<-"./Data/all_taxa.csv"

plotsPath<-"./Figures/richnessAndAbundanceByPhylumAndSample/"

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

AcidobacteriotaCol<-"yellowgreen"
ActinobacteriotaCol<-"yellow2"
BacteroidotaCol<-"wheat3"
BdellovibrionotaCol<-"navyblue"
CalditrichotaCol<-"ivory1"
CampylobacterotaCol<-"violetred4"
ChloroflexiCol<-"violetred2"
CrenarchaeotaCol<-"violetred"
CyanobacteriaCol<-"violet"
DadabacteriaCol<-"turquoise4"
DeferrisomatotaCol<-"turquoise2"
DeinococcotaCol<-"turquoise"
DependentiaeCol<-"tomato3"                 
DesulfobacterotaCol<-"tomato1"   
EntotheonellaeotaCol<-"thistle4"   
FirmicutesCol<-"thistle4"                    
FusobacteriotaCol<-"tan4"               
GemmatimonadotaCol<-"tan3"
HydrogenedentesCol<-"tan1"
LatescibacterotaCol<-"steelblue4"
MyxococcotaCol<-"steelblue2"
NB1jCol<-"springgreen4"
NitrospinotaCol<-"springgreen2"
NitrospirotaCol<-"slategray2"
PAUC34fCol<-"slategray"
PlanctomycetotaCol<-"slateblue4"
ProteobacteriaCol<-"salmon"
SAR324cladeCol<-"royalblue"
SpirochaetotaCol<-"rosybrown3"
SumerlaeotaCol<-"red2"
UnclassifiedCol<-"purple4"
VerrucomicrobiotaCol<-"purple"
WPS2Col<-"plum3"
WS1Col<-"orange"
######
################

####################
allBacteriaColorVector<-c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,BdellovibrionotaCol,
                          CalditrichotaCol,CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
                          DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DependentiaeCol,DesulfobacterotaCol,
                          EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,GemmatimonadotaCol,HydrogenedentesCol,
                          LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,PAUC34fCol,
                          PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,SpirochaetotaCol,SumerlaeotaCol,
                          UnclassifiedCol,VerrucomicrobiotaCol,WPS2Col,WS1Col)
#################

plotSVG<-TRUE
plotPDF<-FALSE

################
#Main Source
###########################

countsDF<-read.csv(uncorrectedCountsPath, sep="\t")
otuTaxonomy<-read.csv(otuTaxonomyPath, sep=";")

#get OTU counts and add them to DF
countSum<-apply(countsDF[-1],1,sum)
countsDF<-cbind(countsDF, countSum)
#Make cumulative percentage to filter out OTUs which make up less than the 5% cummulative sample count
pct<-countSum/sum(countSum)

countsDF<-cbind(countsDF, pct)
countsDF<-countsDF%>%
  mutate(cumpct=cumsum(pct))
dim(countsDF)

#bind the taxonomy columns to the otu counts DF
countsWithTaxonomy<-bind_cols(countsDF[order(countsDF$X.OTU.ID),], otuTaxonomy[order(otuTaxonomy$sequence_identifier),])
#Filter up to 95 % additive abundance
countsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$cumpct < 0.95,]
dim(countsWithTaxonomy)
#delete column again
countsWithTaxonomy=select(countsWithTaxonomy,-cumpct)
#remove OTUs with sum zero if any. Uncomment dim lines if check wanted.
#dim(countsWithTaxonomy)
countsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$countSum != 0,]
dim(countsWithTaxonomy)
#include NAs as unclassified
countsWithTaxonomy[countsWithTaxonomy==""]<-"Unclassified"
#################################
#
# Richness per Phylum per sample plots
#
#######################################################################

#turn abundance into richness and count otus by phylum, melt and plot if required (TOOK OUT FILTER FOR UNCLASSIFIED %>% filter(Phylum != ""))
otuCountByPhylumBySample<-countsWithTaxonomy %>% group_by(Phylum)  %>% mutate_at(colnames(countsWithTaxonomy)[2:22], presentOTU) %>% mutate_at(colnames(countsWithTaxonomy)[2:22],sum) %>% distinct_at(vars(colnames(countsWithTaxonomy)[c(2:22,25)])) 

otuCountByPhylumBySampleMelted<-melt(otuCountByPhylumBySample)

levels(otuCountByPhylumBySampleMelted$variable)<-c("GW1941","GW1942","GW1944","GW1945","GW1946","GW1947","GW1948","GW1952","GW1953","GW1954","GW1955","GW1957","GW1958","GW1959","GW1964","GW1967","GW1968","GW1969","GW1973","GW1979","GW1982")

richnessByPhylumBySampleBarPlot<-ggplot(otuCountByPhylumBySampleMelted, aes(x=variable, y=value)) + geom_bar(position = "stack", stat="identity", aes(fill=Phylum)) + scale_fill_manual(values = allBacteriaColorVector) + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))+labs(xlab="Species sample",ylab="Number of different OTU's per phylum")

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

########WITHOUT UNCLASSIFIED ####
#turn abundance into richness and count otus by phylum, melt and plot if required (TOOK OUT FILTER FOR UNCLASSIFIED %>% filter(Phylum != ""))
xotuCountByPhylumBySample<-countsWithTaxonomy %>% group_by(Phylum)%>%filter(Phylum !="Unclassified")  %>% mutate_at(colnames(countsWithTaxonomy)[2:22], presentOTU) %>% mutate_at(colnames(countsWithTaxonomy)[2:22],sum) %>% distinct_at(vars(colnames(countsWithTaxonomy)[c(2:22,25)])) 

xotuCountByPhylumBySampleMelted<-melt(xotuCountByPhylumBySample)

levels(xotuCountByPhylumBySampleMelted$variable)<-c("GW1941","GW1942","GW1944","GW1945","GW1946","GW1947","GW1948","GW1952","GW1953","GW1954","GW1955","GW1957","GW1958","GW1959","GW1964","GW1967","GW1968","GW1969","GW1973","GW1979","GW1982")

xrichnessByPhylumBySampleBarPlot<-ggplot(xotuCountByPhylumBySampleMelted, aes(x=variable, y=value)) + geom_bar(position = "stack", stat="identity", aes(fill=Phylum)) + scale_fill_manual(values = allBacteriaColorVector) + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))+labs(xlab="Species sample",ylab="Number of different OTU's per phylum")

if(plotPDF){
  
  xrichnessByPhylumBySampleBarPlotPath<-paste(plotsPath,"xrichnessByPhylumBySampleBarPlot.pdf",sep="")
  pdf(xrichnessByPhylumBySampleBarPlotPath)
  print(xrichnessByPhylumBySampleBarPlot)
  dev.off()
  
}

if(plotSVG){
  
  xrichnessByPhylumBySampleBarPlotPath<-paste(plotsPath,"xrichnessByPhylumBySampleBarPlot.svg",sep="")
  svg(xrichnessByPhylumBySampleBarPlotPath)
  print(xrichnessByPhylumBySampleBarPlot)
  dev.off()
  
}

####################################
#
# Abundance per Phylum per sample plots
#
#####################################################
########WITHOUT UNCLASSIFIED ####
xabundanceByPhylumBySample<-countsWithTaxonomy %>% group_by(Phylum)%>%filter(Phylum !="Unclassified") %>% mutate_at(colnames(countsWithTaxonomy)[2:22],sum) %>% distinct_at(vars(colnames(countsWithTaxonomy)[c(2:22,25)])) 
#calculating the relative abundance of the OTU in each sample
xabundanceByPhylumBySample<-ungroup(xabundanceByPhylumBySample) %>% mutate_at(colnames(xabundanceByPhylumBySample)[1:21],function(x)x/sum(x))
xabundanceByPhylumBySampleMelted<-melt(xabundanceByPhylumBySample)
levels(xabundanceByPhylumBySampleMelted$variable)<-c("GW1941","GW1942","GW1944","GW1945","GW1946","GW1947","GW1948","GW1952","GW1953","GW1954","GW1955","GW1957","GW1958","GW1959","GW1964","GW1967","GW1968","GW1969","GW1973","GW1979","GW1982")

xabundanceByPhylumBySampleBarPlot<-ggplot(xabundanceByPhylumBySampleMelted, aes(x=variable, y=value)) + geom_bar(position = "stack", stat="identity", aes(fill=Phylum)) + scale_fill_manual(values = allBacteriaColorVector) + theme_bw()  + theme(axis.text.x=element_text(angle=90, hjust=1))+labs(xlab="Species sample",ylab="Relative abundance of OTU's across phyla per sample")


if(plotPDF){
  
  xabundanceByPhylumBySampleBarPlotPath<-paste(plotsPath,"xabundanceByPhylumBySampleBarPlot.pdf",sep="")
  pdf(xabundanceByPhylumBySampleBarPlotPath)
  print(xabundanceByPhylumBySampleBarPlot)
  dev.off()
  
}

if(plotSVG){
  
  xabundanceByPhylumBySampleBarPlotPath<-paste(plotsPath,"xabundanceByPhylumBySampleBarPlot.svg",sep="")
  svg(xabundanceByPhylumBySampleBarPlotPath)
  print(xabundanceByPhylumBySampleBarPlot)
  dev.off()
  
}

#EOF