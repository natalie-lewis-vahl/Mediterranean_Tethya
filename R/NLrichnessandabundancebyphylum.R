###########################################
#
# Usage:
#
# This script is a companion for Vargas and Leiva.
library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
######################
#Private functions
#
####################################################


################
#Input/Output paths
###########################

#bactLoadCorrectedCountsPath<-"./Data/cbas_tempVSctrl.otutab.bactloadCorrected_valuesOnly.csv"
uncorrectedCountsPath<-"./Data/all.otutab_raw.csv"
otuTaxonomyPath<-"./Data/all_taxa.csv"

plotsPath<-"./Figures/richnessAndAbundanceByPhylum/"

#Main Source
###########################

countsDF<-read.csv(uncorrectedCountsPath, sep="\t")
otuTaxonomy<-read.csv(otuTaxonomyPath, sep=";")
###########################
otuTaxonomy[otuTaxonomy==""]<-"Unclassified"
x<-arrange(otuTaxonomy,Phylum)
unique(x$Phylum)
#y<-arrange(all_countsWithTaxonomy)
#unique(y$Phylum)
# Colors
#Cols<-(c("yellowgreen","yellow2","wheat3","navyblue","ivory1","violetred4","violetred2","violetred","violet",
#       "turquoise4","turquoise2","turquoise","tomato3","tomato1","thistle4","thistle4","tan4","tan3","tan1",
#       "steelblue4","steelblue2","springgreen4","springgreen2","slategray2","slategray","slateblue4","salmon",
#       "royalblue","rosybrown3","red2","purple4","purple","plum3","orange"))
#names(Cols)<-levels(otuTaxonomy$Phylum)
#colours_set <- scale_fill_manual(name = "Var1",values = Cols)

################
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

#Summary columns for counts of each sp (every 7 columns, start function at column 2 and make the last one form 16)
dim(countsDF)
x<-sapply(seq(2,16,by=7),function(i) rowSums(countsDF[,i:(i+6)]))
#add them to DF
countsDF<-cbind(countsDF, x)
#Suma all columns and add
allSum<-apply(countsDF[,2:22],1,sum)
countsDF<-cbind(countsDF,allSum)
#rename cols
names(countsDF)
countsDF<-rename(countsDF,sumTau="1",sumTme="2",sumTci="3")

#Make cumulative percentage to filter out OTUs which make up less than the 5% cummulative sample count
pct<-allSum/sum(allSum)

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

#remove OTUs with sum zero for all species
dim(countsWithTaxonomy)
all_countsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$allSum!=0,]

#remove OTUs with sum zero for Tau
tau_countsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$sumTau!=0,]
dim(tau_countsWithTaxonomy)#912

#remove OTUs with sum zero for Tme
tme_countsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$sumTme!=0,]
dim(tme_countsWithTaxonomy)#257

#remove OTUs with sum zero for Tci
tci_countsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$sumTci!=0,]
dim(tci_countsWithTaxonomy)#2240

#######################
#
# OTU Richness by Phylum
#FOR ALL SPECIES
#
###########################################

#get a table with the number of OTUs annotated to each phylum for all sp
phylumRichness<-table(all_countsWithTaxonomy$Phylum)
allsp_allPhylaRichnessPie<-ggplot(as.data.frame(phylumRichness), aes(x="", y=Freq,fill=Var1)) + geom_bar(stat = "identity") + coord_polar("y") +scale_fill_manual(values=c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,BdellovibrionotaCol,
                                                                                                                                                                              CalditrichotaCol,CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
                                                                                                                                                                              DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DependentiaeCol,DesulfobacterotaCol,
                                                                                                                                                                              EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,GemmatimonadotaCol,HydrogenedentesCol,
                                                                                                                                                                              LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,PAUC34fCol,
                                                                                                                                                                              PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,SpirochaetotaCol,SumerlaeotaCol,
                                                                                                                                                                              UnclassifiedCol,VerrucomicrobiotaCol,WPS2Col,WS1Col))+ theme_bw()
#Get legend
legend<-get_legend(allsp_allPhylaRichnessPie)
#Remove legend from plot
allsp_allPhylaRichnessPie<-ggplot(as.data.frame(phylumRichness), aes(x="", y=Freq,fill=Var1)) + geom_bar(stat = "identity",show.legend=FALSE) + 
  coord_polar("y") +scale_fill_manual(values=c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,BdellovibrionotaCol, CalditrichotaCol,CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DependentiaeCol,DesulfobacterotaCol,EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,GemmatimonadotaCol,HydrogenedentesCol,LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,PAUC34fCol,PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,SpirochaetotaCol,SumerlaeotaCol,UnclassifiedCol,VerrucomicrobiotaCol,WPS2Col,WS1Col))+ theme_bw()


#richness pie chart

if(plotPDF){
  
  richnessPiePath<-paste(plotsPath,"allsp_allPhylaRichnessPie.pdf", sep="")
  pdf(allsp_richnessPiePath)
  print(allsp_allPhylaRichnessPie)
  dev.off()
  
}

if(plotSVG){
  
  allsp_richnessPiePath<-paste(plotsPath,"allsp_allPhylaRichnessPie.svg", sep="")
  svg(allsp_richnessPiePath)
  print(allsp_allPhylaRichnessPie)
  dev.off()
  
}

#split in low and high richness
highRichnessPhyla<-phylumRichness[phylumRichness>10]
names(highRichnessPhyla)
allsp_highRichnessPhylaPie<-ggplot(as.data.frame(highRichnessPhyla), aes(x="", y=Freq, fill=Var1)) + 
  geom_bar(stat = "identity",show.legend=FALSE) + coord_polar("y") +
  scale_fill_manual(values = c(AcidobacteriotaCol, ActinobacteriotaCol, BacteroidotaCol, 
                               BdellovibrionotaCol, ChloroflexiCol,CyanobacteriaCol, DesulfobacterotaCol, 
                               GemmatimonadotaCol,LatescibacterotaCol,MyxococcotaCol,NB1jCol, PlanctomycetotaCol, 
                               ProteobacteriaCol, UnclassifiedCol, VerrucomicrobiotaCol)) + theme_bw()

lowRichnessPhyla<-phylumRichness[phylumRichness<=10]
names(lowRichnessPhyla)
allsp_lowRichnessPhylaPie<-ggplot(as.data.frame(lowRichnessPhyla), aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat = "identity",show.legend=FALSE) + coord_polar("y") +
  scale_fill_manual(values = c(CalditrichotaCol,CampylobacterotaCol,CrenarchaeotaCol,DadabacteriaCol,DeferrisomatotaCol,
                               DeinococcotaCol,DependentiaeCol,EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,HydrogenedentesCol,
                               NitrospinotaCol,NitrospirotaCol,PAUC34fCol,SAR324cladeCol,SpirochaetotaCol,SumerlaeotaCol,WPS2Col,WS1Col)) +
  theme_bw()

if(plotPDF){
  
  allsp_highRichnessPhylaPiePath<-paste(plotsPath,"allsp_highRichnessPhylaPie.pdf", sep="")
  pdf(allsp_highRichnessPhylaPiePath)
  print(allsp_highRichnessPhylaPie)
  dev.off()
  
  allsp_lowRichnessPhylaPiePath<-paste(plotsPath,"allsp_lowRichnessPhylaPie.pdf", sep="")
  pdf(allsp_lowRichnessPhylaPiePath)
  print(allsp_lowRichnessPhylaPie)
  dev.off()
  
}

if(plotSVG){
  
 allsp_highRichnessPhylaPiePath<-paste(plotsPath,"allsp_highRichnessPhylaPie.svg", sep="")
 svg(allsp_highRichnessPhylaPiePath)
  print(allsp_highRichnessPhylaPie)
  dev.off()
  
  allsp_lowRichnessPhylaPiePath<-paste(plotsPath,"allsp_lowRichnessPhylaPie.svg", sep="")
  svg(allsp_lowRichnessPhylaPiePath)
  print(allsp_lowRichnessPhylaPie)
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
phylumCounts<-all_countsWithTaxonomy %>% group_by(Phylum) %>% summarise(sum=sum(allSum))
head(phylumCounts$Phylum)
allsp_allPhylaAbundancePie<-ggplot(phylumCounts, aes(x="", y=sum, fill=Phylum)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = allBacteriaColorVector) + theme_bw()
allsp_allPhylaAbundancePie
if(plotPDF){
  
  allsp_allPhylaAbundancePiePath<-paste(plotsPath,"allsp_allPhylaAbundancePie.pdf", sep="")
  pdf(allsp_allPhylaAbundancePiePath)
  print(allsp_allPhylaAbundancePie)
  dev.off()
  
}

if(plotSVG){
  
  allsp_allPhylaAbundancePiePath<-paste(plotsPath,"allsp_allPhylaAbundancePie.svg", sep="")
  svg(allsp_allPhylaAbundancePiePath)
  print(allsp_allPhylaAbundancePie)
  dev.off()
  
}

#split in high and low abundance
highCountsPhyla<-phylumCounts[phylumCounts$sum > 1000,]
lowCountsPhyla<-phylumCounts[phylumCounts$sum <= 1000,]
unique(highCountsPhyla$Phylum)
#plot
allsp_highAbundancePhylaPie<-ggplot(highCountsPhyla, aes(x="", y=sum, fill=Phylum)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,CrenarchaeotaCol,CyanobacteriaCol,
                                                                                                                                                                   DesulfobacterotaCol,PlanctomycetotaCol,ProteobacteriaCol,UnclassifiedCol,VerrucomicrobiotaCol)) + theme_bw()
unique(lowCountsPhyla$Phylum)
allsp_lowAbundancePhylaPie<-ggplot(lowCountsPhyla, aes(x="", y=sum, fill=Phylum)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = c(BdellovibrionotaCol,CalditrichotaCol,CampylobacterotaCol,ChloroflexiCol,DadabacteriaCol,DeferrisomatotaCol,
                                                                                                                                                            DeinococcotaCol,DependentiaeCol,EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,GemmatimonadotaCol,
                                                                                                                                                            HydrogenedentesCol,LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,
                                                                                                                                                            PAUC34fCol,SAR324cladeCol,SpirochaetotaCol,SumerlaeotaCol,WPS2Col,WS1Col)) + theme_bw()

if(plotPDF){
  
  allsp_highAbundancePhylaPiePath<-paste(plotsPath,"allsp_highAbundancePhylaPie.pdf", sep="")
  pdf(allsp_highAbundancePhylaPiePath)
  print(allsp_highAbundancePhylaPie)
  dev.off()
  allsp_lowAbundancePhylaPiePath<-paste(plotsPath,"allsp_lowAbundancePhylaPie.pdf", sep="")
  pdf(allsp_lowAbundancePhylaPiePath)
  print(allsp_lowAbundancePhylaPie)
  dev.off()
}

if(plotSVG){
  
  allsp_highAbundancePhylaPiePath<-paste(plotsPath,"allsp_highAbundancePhylaPie.svg", sep="")
  svg(allsp_highAbundancePhylaPiePath)
  print(allsp_highAbundancePhylaPie)
  dev.off()
  allsp_lowAbundancePhylaPiePath<-paste(plotsPath,"allsp_lowAbundancePhylaPie.svg", sep="")
  svg(allsp_lowAbundancePhylaPiePath)
  print(allsp_lowAbundancePhylaPie)
  dev.off()
  
}

#######################
#
# OTU Richness by Phylum
#FOR Tethya aurantium
#
###########################################

#get a table with the number of OTUs annotated to each phylum for all sp
tauphylumRichness<-table(tau_countsWithTaxonomy$Phylum)
dim(tauphylumRichness)
names(tauphylumRichness)
tau_allPhylaRichnessPie<-ggplot(as.data.frame(tauphylumRichness), aes(x="", y=Freq, fill=Var1)) + geom_bar(stat = "identity") + coord_polar("y") +scale_fill_manual(values=c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,BdellovibrionotaCol,
                                                                                                                                                                              CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
                                                                                                                                                                              DadabacteriaCol,DeinococcotaCol,DesulfobacterotaCol,
                                                                                                                                                                              EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,GemmatimonadotaCol,
                                                                                                                                                                              LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,
                                                                                                                                                                              PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,
                                                                                                                                                                              UnclassifiedCol,VerrucomicrobiotaCol))+ theme_bw()
tau_allPhylaRichnessPie

#richness pie chart

if(plotPDF){
  
  tau_richnessPiePath<-paste(plotsPath,"tau_allPhylaRichnessPie.pdf", sep="")
  pdf(tau_richnessPiePath)
  print(tau_allPhylaRichnessPie)
  dev.off()
  
}

##if(plotSVG){

#  richnessPiePath<-paste(plotsPath,"allPhylaRichnessPie.svg", sep="")
#  svg(richnessPiePath)
#  print(allPhylaRichnessPie)
#  dev.off()

#}

#split in low and high richness
tauhighRichnessPhyla<-tauphylumRichness[tauphylumRichness>10]
names(tauhighRichnessPhyla)
tauhighRichnessPhylaPie<-ggplot(as.data.frame(tauhighRichnessPhyla), aes(x="", y=Freq, fill=Var1)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = c(AcidobacteriotaCol, ActinobacteriotaCol, BacteroidotaCol,
                                                                                                                                                                             NB1jCol, PlanctomycetotaCol, ProteobacteriaCol, UnclassifiedCol, VerrucomicrobiotaCol)) + theme_bw()
tauhighRichnessPhylaPie
taulowRichnessPhyla<-tauphylumRichness[tauphylumRichness<=10]
names(taulowRichnessPhyla)
taulowRichnessPhylaPie<-ggplot(as.data.frame(taulowRichnessPhyla), aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat = "identity") + coord_polar("y") +
  scale_fill_manual(values = c(BdellovibrionotaCol,CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,
                               CyanobacteriaCol,DadabacteriaCol,DeinococcotaCol,DesulfobacterotaCol,
                               EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,GemmatimonadotaCol,
                               LatescibacterotaCol,MyxococcotaCol,NitrospinotaCol,NitrospirotaCol,SAR324cladeCol)) +
  theme_bw()
taulowRichnessPhylaPie

if(plotPDF){
  
  tauhighRichnessPhylaPiePath<-paste(plotsPath,"tauhighRichnessPhylaPie.pdf", sep="")
  pdf(tauhighRichnessPhylaPiePath)
  print(tauhighRichnessPhylaPie)
  dev.off()
  
  taulowRichnessPhylaPiePath<-paste(plotsPath,"taulowRichnessPhylaPie.pdf", sep="")
  pdf(taulowRichnessPhylaPiePath)
  print(taulowRichnessPhylaPie)
  dev.off()
  
}

if(plotSVG){
  
  tauhighRichnessPhylaPiePath<-paste(plotsPath,"tauhighRichnessPhylaPie.svg", sep="")
  svg(tauhighRichnessPhylaPiePath)
  print(tauhighRichnessPhylaPie)
  dev.off()
  
  taulowRichnessPhylaPiePath<-paste(plotsPath,"taulowRichnessPhylaPie.svg", sep="")
  svg(taulowRichnessPhylaPiePath)
  print(taulowRichnessPhylaPie)
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
tauphylumCounts<-tau_countsWithTaxonomy %>% group_by(Phylum) %>% summarise(sum=sum(sumTau))
head(tauphylumCounts$Phylum)
tauPhylaAbundancePie<-ggplot(tauphylumCounts, aes(x="", y=sum, fill=Phylum)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,BdellovibrionotaCol,
                                                                                                                                                           CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
                                                                                                                                                           DadabacteriaCol,DeinococcotaCol,DesulfobacterotaCol,
                                                                                                                                                           EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,GemmatimonadotaCol,
                                                                                                                                                           LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,
                                                                                                                                                           PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,
                                                                                                                                                           UnclassifiedCol,VerrucomicrobiotaCol)) + theme_bw()
tauPhylaAbundancePie
if(plotPDF){
  
  tauPhylaAbundancePiePath<-paste(plotsPath,"tauPhylaAbundancePie.pdf", sep="")
  pdf(tauPhylaAbundancePiePath)
  print(tauPhylaAbundancePie)
  dev.off()
  
}

if(plotSVG){
  
  tauPhylaAbundancePiePath<-paste(plotsPath,"tauPhylaAbundancePie.svg", sep="")
  svg(tauPhylaAbundancePiePath)
  print(tauPhylaAbundancePie)
  dev.off()
  
}

#split in high and low abundance
tauhighCountsPhyla<-tauphylumCounts[tauphylumCounts$sum > 1000,]
taulowCountsPhyla<-tauphylumCounts[tauphylumCounts$sum <= 1000,]
unique(tauhighCountsPhyla$Phylum)
#plot
tauhighAbundancePhylaPie<-ggplot(tauhighCountsPhyla, aes(x="", y=sum, fill=Phylum)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = c(ActinobacteriotaCol,BacteroidotaCol,CyanobacteriaCol,
                                                                                                                                                              PlanctomycetotaCol,ProteobacteriaCol,UnclassifiedCol,VerrucomicrobiotaCol)) + theme_bw()
unique(taulowCountsPhyla$Phylum)
taulowAbundancePhylaPie<-ggplot(taulowCountsPhyla, aes(x="", y=sum, fill=Phylum)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = c(AcidobacteriotaCol,BdellovibrionotaCol,
                                                                                                                                                                  CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,
                                                                                                                                                                  DadabacteriaCol,DeinococcotaCol,DesulfobacterotaCol,
                                                                                                                                                                  EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,GemmatimonadotaCol,
                                                                                                                                                                  LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,
                                                                                                                                                                  SAR324cladeCol)) + theme_bw()

if(plotPDF){
  
  tauhighAbundancePhylaPiePath<-paste(plotsPath,"tauhighAbundancePhylaPie.pdf", sep="")
  pdf(tauhighAbundancePhylaPiePath)
  print(tauhighAbundancePhylaPie)
  dev.off()
  taulowAbundancePhylaPiePath<-paste(plotsPath,"taulowAbundancePhylaPie.pdf", sep="")
  pdf(taulowAbundancePhylaPiePath)
  print(taulowAbundancePhylaPie)
  dev.off()
}

if(plotSVG){
  
  tauhighAbundancePhylaPiePath<-paste(plotsPath,"tauhighAbundancePhylaPie.svg", sep="")
  svg(tauhighAbundancePhylaPiePath)
  print(tauhighAbundancePhylaPie)
  dev.off()
  taulowAbundancePhylaPiePath<-paste(plotsPath,"taulowAbundancePhylaPie.svg", sep="")
  svg(taulowAbundancePhylaPiePath)
  print(taulowAbundancePhylaPie)
  dev.off()
  
}

#######################
#
# OTU Richness by Phylum
#FOR Tethya meloni
#
###########################################

#get a table with the number of OTUs annotated to each phylum for all sp
tmephylumRichness<-table(tme_countsWithTaxonomy$Phylum)
dim(tmephylumRichness)
names(tmephylumRichness)
tme_allPhylaRichnessPie<-ggplot(as.data.frame(tmephylumRichness), aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat = "identity") + coord_polar("y") +
  scale_fill_manual(values=c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,
                               CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,DependentiaeCol,
                               DesulfobacterotaCol,EntotheonellaeotaCol,GemmatimonadotaCol,MyxococcotaCol,NB1jCol,
                               NitrospinotaCol,NitrospirotaCol,PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,UnclassifiedCol,VerrucomicrobiotaCol))+ theme_bw()
tme_allPhylaRichnessPie

#richness pie chart

if(plotPDF){
  
  tme_richnessPiePath<-paste(plotsPath,"tme_allPhylaRichnessPie.pdf", sep="")
  pdf(tme_richnessPiePath)
  print(tme_allPhylaRichnessPie)
  dev.off()
  
}

if(plotSVG){

  tme_richnessPiePath<-paste(plotsPath,"tme_allPhylaRichnessPie.svg", sep="")
  svg(tme_richnessPiePath)
  print(tme_allPhylaRichnessPie)
  dev.off()

#}

#split in low and high richness
tmehighRichnessPhyla<-tmephylumRichness[tmephylumRichness>10]
names(tmehighRichnessPhyla)
tmehighRichnessPhylaPie<-ggplot(as.data.frame(tmehighRichnessPhyla), aes(x="", y=Freq, fill=Var1)) + 
  geom_bar(stat = "identity") + coord_polar("y") +
  scale_fill_manual(values = c(ActinobacteriotaCol, BacteroidotaCol,PlanctomycetotaCol, ProteobacteriaCol, UnclassifiedCol)) + theme_bw()
tmehighRichnessPhylaPie
tmelowRichnessPhyla<-tmephylumRichness[tmephylumRichness<=10]
names(tmelowRichnessPhyla)
tmelowRichnessPhylaPie<-ggplot(as.data.frame(tmelowRichnessPhyla), aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat = "identity") + coord_polar("y") +
  scale_fill_manual(values = c(AcidobacteriotaCol,
                               CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,DependentiaeCol,
                               DesulfobacterotaCol,EntotheonellaeotaCol,GemmatimonadotaCol,MyxococcotaCol,NB1jCol,
                               NitrospinotaCol,NitrospirotaCol,SAR324cladeCol,VerrucomicrobiotaCol)) +
  theme_bw()

tmelowRichnessPhylaPie

if(plotPDF){
  
  tmehighRichnessPhylaPiePath<-paste(plotsPath,"tmehighRichnessPhylaPie.pdf", sep="")
  pdf(tmehighRichnessPhylaPiePath)
  print(tmehighRichnessPhylaPie)
  dev.off()
  
  tmelowRichnessPhylaPiePath<-paste(plotsPath,"tmelowRichnessPhylaPie.pdf", sep="")
  pdf(tmelowRichnessPhylaPiePath)
  print(tmelowRichnessPhylaPie)
  dev.off()
  
}

if(plotSVG){
  
  tmehighRichnessPhylaPiePath<-paste(plotsPath,"tmehighRichnessPhylaPie.svg", sep="")
  svg(tmehighRichnessPhylaPiePath)
  print(tmehighRichnessPhylaPie)
  dev.off()
  
  tmelowRichnessPhylaPiePath<-paste(plotsPath,"tmelowRichnessPhylaPie.svg", sep="")
  svg(tmelowRichnessPhylaPiePath)
  print(tmelowRichnessPhylaPie)
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
tmephylumCounts<-tme_countsWithTaxonomy %>% group_by(Phylum) %>% summarise(sum=sum(sumTme))
names(tmephylumCounts$Phylum)
tmePhylaAbundancePie<-ggplot(tmephylumCounts, aes(x="", y=sum, fill=Phylum)) + geom_bar(stat = "identity") + coord_polar("y") +scale_fill_manual(values = c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,
                                                                                                                                                             CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,DependentiaeCol,
                                                                                                                                                             DesulfobacterotaCol,EntotheonellaeotaCol,GemmatimonadotaCol,MyxococcotaCol,NB1jCol,
                                                                                                                                                             NitrospinotaCol,NitrospirotaCol,PlanctomycetotaCol,ProteobacteriaCol,
                                                                                                                                                             SAR324cladeCol,UnclassifiedCol,VerrucomicrobiotaCol)) + theme_bw()
tmePhylaAbundancePie
if(plotPDF){
  
  tmePhylaAbundancePiePath<-paste(plotsPath,"tmePhylaAbundancePie.pdf", sep="")
  pdf(tmePhylaAbundancePiePath)
  print(tmePhylaAbundancePie)
  dev.off()
  
}

if(plotSVG){
  
  tmePhylaAbundancePiePath<-paste(plotsPath,"tmePhylaAbundancePie.svg", sep="")
  svg(tmeallPhylaAbundancePiePath)
  print(tmeallPhylaAbundancePie)
  dev.off()
  
}

#split in high and low abundance
tmehighCountsPhyla<-tmephylumCounts[tmephylumCounts$sum > 1000,]
tmelowCountsPhyla<-tmephylumCounts[tmephylumCounts$sum <= 1000,]
unique(tmehighCountsPhyla$Phylum)
#plot
tmehighAbundancePhylaPie<-ggplot(tmehighCountsPhyla, aes(x="", y=sum, fill=Phylum)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = c(ActinobacteriotaCol,ProteobacteriaCol,UnclassifiedCol)) + theme_bw()
unique(tmelowCountsPhyla$Phylum)
tmelowAbundancePhylaPie<-ggplot(tmelowCountsPhyla, aes(x="", y=sum, fill=Phylum)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = c(AcidobacteriotaCol,BacteroidotaCol,
                                                                                                                                                                  CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,DependentiaeCol,
                                                                                                                                                                  DesulfobacterotaCol,EntotheonellaeotaCol,GemmatimonadotaCol,MyxococcotaCol,NB1jCol,
                                                                                                                                                                  NitrospinotaCol,NitrospirotaCol,PlanctomycetotaCol,
                                                                                                                                                                  SAR324cladeCol,VerrucomicrobiotaCol)) + theme_bw()

if(plotPDF){
  
  tmehighAbundancePhylaPiePath<-paste(plotsPath,"tmehighAbundancePhylaPie.pdf", sep="")
  pdf(tmehighAbundancePhylaPiePath)
  print(tmehighAbundancePhylaPie)
  dev.off()
  tmelowAbundancePhylaPiePath<-paste(plotsPath,"tmelowAbundancePhylaPie.pdf", sep="")
  pdf(tmelowAbundancePhylaPiePath)
  print(tmelowAbundancePhylaPie)
  dev.off()
}

if(plotSVG){
  
  tmehighAbundancePhylaPiePath<-paste(plotsPath,"tmehighAbundancePhylaPie.svg", sep="")
  svg(tmehighAbundancePhylaPiePath)
  print(tmehighAbundancePhylaPie)
  dev.off()
  tmelowAbundancePhylaPiePath<-paste(plotsPath,"tmelowAbundancePhylaPie.svg", sep="")
  svg(tmelowAbundancePhylaPiePath)
  print(tmelowAbundancePhylaPie)
  dev.off()
  
}

#######################
#
# OTU Richness by Phylum
#FOR Tethya citroni
#
###########################################

#get a table with the number of OTUs annotated to each phylum for T citroni
tciphylumRichness<-table(tci_countsWithTaxonomy$Phylum)
dim(tciphylumRichness)
names(tciphylumRichness)
tci_allPhylaRichnessPie<-ggplot(as.data.frame(tciphylumRichness), aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat = "identity") + coord_polar("y") +
  scale_fill_manual(values=c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,BdellovibrionotaCol,CalditrichotaCol,
                               CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,DadabacteriaCol,DeferrisomatotaCol,
                               DeinococcotaCol,DependentiaeCol,DesulfobacterotaCol,EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,
                               GemmatimonadotaCol,HydrogenedentesCol,LatescibacterotaCol,MyxococcotaCol,NB1jCol,
                               NitrospinotaCol,NitrospirotaCol,PAUC34fCol,PlanctomycetotaCol,ProteobacteriaCol,
                               SAR324cladeCol,SpirochaetotaCol,SumerlaeotaCol,UnclassifiedCol,VerrucomicrobiotaCol,WPS2Col,WS1Col))+ theme_bw()
tci_allPhylaRichnessPie

#richness pie chart

if(plotPDF){
  
  tci_richnessPiePath<-paste(plotsPath,"tci_allPhylaRichnessPie.pdf", sep="")
  pdf(tci_richnessPiePath)
  print(tci_allPhylaRichnessPie)
  dev.off()
  
}

if(plotSVG){

  tci_richnessPiePath<-paste(plotsPath,"tci_allPhylaRichnessPie.svg", sep="")
  svg(tci_richnessPiePath)
  print(tci_allPhylaRichnessPie)
  dev.off()

}

#split in low and high richness
tcihighRichnessPhyla<-tciphylumRichness[tciphylumRichness>10]
names(tcihighRichnessPhyla)
tcihighRichnessPhylaPie<-ggplot(as.data.frame(tcihighRichnessPhyla), aes(x="", y=Freq, fill=Var1)) + 
  geom_bar(stat = "identity") + coord_polar("y") +
  scale_fill_manual(values = c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,BdellovibrionotaCol,ChloroflexiCol,
                               DesulfobacterotaCol,GemmatimonadotaCol,LatescibacterotaCol,MyxococcotaCol,NB1jCol,
                              PlanctomycetotaCol,ProteobacteriaCol,UnclassifiedCol,VerrucomicrobiotaCol)) + theme_bw()
tcihighRichnessPhylaPie
tcilowRichnessPhyla<-tciphylumRichness[tciphylumRichness<=10]
names(tcilowRichnessPhyla)
tcilowRichnessPhylaPie<-ggplot(as.data.frame(tcilowRichnessPhyla), aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat = "identity") + coord_polar("y") +
  scale_fill_manual(values = c(CalditrichotaCol,
                               CampylobacterotaCol,CrenarchaeotaCol,CyanobacteriaCol,DadabacteriaCol,DeferrisomatotaCol,
                               DeinococcotaCol,DependentiaeCol,EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,
                               HydrogenedentesCol,NitrospinotaCol,NitrospirotaCol,PAUC34fCol,
                               SAR324cladeCol,SpirochaetotaCol,SumerlaeotaCol,WPS2Col,WS1Col)) +
  theme_bw()

tcilowRichnessPhylaPie

if(plotPDF){
  
  tcihighRichnessPhylaPiePath<-paste(plotsPath,"tcihighRichnessPhylaPie.pdf", sep="")
  pdf(tcihighRichnessPhylaPiePath)
  print(tcihighRichnessPhylaPie)
  dev.off()
  
  tcilowRichnessPhylaPiePath<-paste(plotsPath,"tcilowRichnessPhylaPie.pdf", sep="")
  pdf(tcilowRichnessPhylaPiePath)
  print(tcilowRichnessPhylaPie)
  dev.off()
  
}

if(plotSVG){
  
  tcihighRichnessPhylaPiePath<-paste(plotsPath,"tcihighRichnessPhylaPie.svg", sep="")
  svg(tcihighRichnessPhylaPiePath)
  print(tcihighRichnessPhylaPie)
  dev.off()
  
  tcilowRichnessPhylaPiePath<-paste(plotsPath,"tcilowRichnessPhylaPie.svg", sep="")
  svg(tcilowRichnessPhylaPiePath)
  print(tcilowRichnessPhylaPie)
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
tciphylumCounts<-tci_countsWithTaxonomy %>% group_by(Phylum) %>% summarise(sum=sum(sumTci))
head(tciphylumCounts$Phylum)
tciPhylaAbundancePie<-ggplot(tciphylumCounts, aes(x="", y=sum, fill=Phylum)) + geom_bar(stat = "identity") + coord_polar("y") +
  scale_fill_manual(values = c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,BdellovibrionotaCol,CalditrichotaCol,
                               CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,DadabacteriaCol,DeferrisomatotaCol,
                               DeinococcotaCol,DependentiaeCol,DesulfobacterotaCol,EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,
                               GemmatimonadotaCol,HydrogenedentesCol,LatescibacterotaCol,MyxococcotaCol,NB1jCol,
                               NitrospinotaCol,NitrospirotaCol,PAUC34fCol,PlanctomycetotaCol,ProteobacteriaCol,
                               SAR324cladeCol,SpirochaetotaCol,SumerlaeotaCol,UnclassifiedCol,VerrucomicrobiotaCol,WPS2Col,WS1Col)) + theme_bw()
tciPhylaAbundancePie

if(plotPDF){
  
  tciPhylaAbundancePiePath<-paste(plotsPath,"tciPhylaAbundancePie.pdf", sep="")
  pdf(tciPhylaAbundancePiePath)
  print(tciPhylaAbundancePie)
  dev.off()
  
}

if(plotSVG){
  
  tciPhylaAbundancePiePath<-paste(plotsPath,"tciPhylaAbundancePie.svg", sep="")
  svg(tciallPhylaAbundancePiePath)
  print(tciallPhylaAbundancePie)
  dev.off()
  
}

#split in high and low abundance
tcihighCountsPhyla<-tciphylumCounts[tciphylumCounts$sum > 1000,]
tcilowCountsPhyla<-tciphylumCounts[tciphylumCounts$sum <= 1000,]
unique(tcihighCountsPhyla$Phylum)
#plot
tcihighAbundancePhylaPie<-ggplot(tcihighCountsPhyla, aes(x="", y=sum, fill=Phylum)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,CrenarchaeotaCol,CyanobacteriaCol,
                                                                                                                                                                    DesulfobacterotaCol,PlanctomycetotaCol,ProteobacteriaCol,UnclassifiedCol)) + theme_bw()
unique(tcilowCountsPhyla$Phylum)
tcilowAbundancePhylaPie<-ggplot(tcilowCountsPhyla, aes(x="", y=sum, fill=Phylum)) + geom_bar(stat = "identity") + coord_polar("y") + scale_fill_manual(values = c(BdellovibrionotaCol,CalditrichotaCol,CampylobacterotaCol,ChloroflexiCol,
                                                                                                                                                                  DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DependentiaeCol,
                                                                                                                                                                  EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,GemmatimonadotaCol,HydrogenedentesCol,
                                                                                                                                                                  LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,PAUC34fCol,
                                                                                                                                                                  SAR324cladeCol,SpirochaetotaCol,SumerlaeotaCol,VerrucomicrobiotaCol,WPS2Col,WS1Col)) + theme_bw()

if(plotPDF){
  
  tcihighAbundancePhylaPiePath<-paste(plotsPath,"tcihighAbundancePhylaPie.pdf", sep="")
  pdf(tcihighAbundancePhylaPiePath)
  print(tcihighAbundancePhylaPie)
  dev.off()
  tcilowAbundancePhylaPiePath<-paste(plotsPath,"tcilowAbundancePhylaPie.pdf", sep="")
  pdf(tcilowAbundancePhylaPiePath)
  print(tcilowAbundancePhylaPie)
  dev.off()
}

if(plotSVG){
  
  tcihighAbundancePhylaPiePath<-paste(plotsPath,"tcihighAbundancePhylaPie.svg", sep="")
  svg(tcihighAbundancePhylaPiePath)
  print(tcihighAbundancePhylaPie)
  dev.off()
  tcilowAbundancePhylaPiePath<-paste(plotsPath,"tcilowAbundancePhylaPie.svg", sep="")
  svg(tcilowAbundancePhylaPiePath)
  print(tcilowAbundancePhylaPie)
  dev.off()
  
}


##Combine allsp_ richness plots together
tau_allPhylaRichnessPie<-tau_allPhylaRichnessPie+theme(legend.position="none")
tauhighRichnessPhylaPie<-tauhighRichnessPhylaPie+theme(legend.position="none")
taulowRichnessPhylaPie<-taulowRichnessPhylaPie+theme(legend.position="none")
tme_allPhylaRichnessPie<-tme_allPhylaRichnessPie+theme(legend.position="none")
tmehighRichnessPhylaPie<-tmehighRichnessPhylaPie+theme(legend.position="none")
tmelowRichnessPhylaPie<-tmelowRichnessPhylaPie+theme(legend.position="none")
tci_allPhylaRichnessPie<-tci_allPhylaRichnessPie+theme(legend.position="none")
tcihighRichnessPhylaPie<-tcihighRichnessPhylaPie+theme(legend.position="none")
tcilowRichnessPhylaPie<-tcilowRichnessPhylaPie+theme(legend.position="none")

allspecies_richnessplots<-list(allsp_allPhylaRichnessPie,allsp_highRichnessPhylaPie,allsp_lowRichnessPhylaPie,tau_allPhylaRichnessPie,tauhighRichnessPhylaPie,taulowRichnessPhylaPie,tme_allPhylaRichnessPie,tmehighRichnessPhylaPie,tmelowRichnessPhylaPie,tci_allPhylaRichnessPie,tcihighRichnessPhylaPie,tcilowRichnessPhylaPie,legend)

ggsave("allspeciesrichnessplots.jpg", path=plotsPath, grid.arrange(grobs = allspecies_richnessplots,layout_matrix=rbind(c(1,2,3,13),c(4,5,6,13),c(7,8,9,13),c(10,11,12,13))))

##Combine allsp_ abundance plots together
allsp_allPhylaAbundancePie<-allsp_allPhylaAbundancePie+theme(legend.position="none")
allsp_highAbundancePhylaPie<-allsp_highAbundancePhylaPie+theme(legend.position="none")
allsp_lowAbundancePhylaPie<-allsp_lowAbundancePhylaPie+theme(legend.position="none")
tauPhylaAbundancePie<-tauPhylaAbundancePie+theme(legend.position="none")
tauhighAbundancePhylaPie<-tauhighAbundancePhylaPie+theme(legend.position="none")
taulowAbundancePhylaPie<-taulowAbundancePhylaPie+theme(legend.position="none")
tmePhylaAbundancePie<-tmePhylaAbundancePie+theme(legend.position="none")
tmehighAbundancePhylaPie<-tmehighAbundancePhylaPie+theme(legend.position="none")
tmelowAbundancePhylaPie<-tmelowAbundancePhylaPie+theme(legend.position="none")
tciPhylaAbundancePie<-tciPhylaAbundancePie+theme(legend.position="none")
tcihighAbundancePhylaPie<-tcihighAbundancePhylaPie+theme(legend.position="none")
tcilowAbundancePhylaPie<-tcilowAbundancePhylaPie+theme(legend.position="none")

allspecies_Abundanceplots<-list(allsp_allPhylaAbundancePie,allsp_highAbundancePhylaPie,allsp_lowAbundancePhylaPie,tauPhylaAbundancePie,tauhighAbundancePhylaPie,taulowAbundancePhylaPie,tmePhylaAbundancePie,tmehighAbundancePhylaPie,tmelowAbundancePhylaPie,tciPhylaAbundancePie,tcihighAbundancePhylaPie,tcilowAbundancePhylaPie,legend)

ggsave("allspeciesAbundanceplots.jpg", path=plotsPath, grid.arrange(grobs = allspecies_Abundanceplots,layout_matrix=rbind(c(1,2,3,13),c(4,5,6,13),c(7,8,9,13),c(10,11,12,13))))
