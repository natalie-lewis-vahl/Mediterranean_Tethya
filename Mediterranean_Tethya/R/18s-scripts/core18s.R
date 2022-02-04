library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(vegan)
library(reshape2)
#ENter data and order the data frames
taxa<-read.csv("./Data/18s/taxa_fixed18s.csv",sep=";")
unique(taxa$Phylum)
otus<-read.csv("./Data/18s/all.otutab.csv", sep=";")
otus<-arrange(otus, X.OTU.ID)
taxa<-arrange(taxa,sequence_identifier)

oat<-bind_cols(otus[order(otus$X.OTU.ID),], taxa[order(taxa$sequence_identifier),])
dim(taxa)

countSum<-apply(oat[2:45],1,sum)#sum the rows
oatsum<-cbind(oat,countSum)
#Filter for OTUs with more than 5 reads
oat<-oatsum[oatsum$countSum>=5,]
dim(oat)
#Filtering for 95 percent
#pct<-countSum/sum(countSum)
#oatsum<-cbind(oatsum, pct)
#oatsum[desc(oatsum$pct),]
#oatx<-oatsum%>%
#  arrange(desc(pct))
#oat<-oatx%>%
#  mutate(cumpct=cumsum(pct))
#oat<-oat[oat$cumpct < 0.95,]
#remove columns for cumpct and pct
#oat<-oat[,-c(58,59)]
#remove column cumsum
oat<-oat[,-57]

#FIlter out bacteria
oat<-oat%>%filter(Domain!="Bacteria")
dim(oat)
unique(oat$Phylum)
head(oat)
#remove otu 1 and 2 because of co-amplification
#FInd which rows to delete
which(oat$X.OTU.ID=="OTU_4076")
which(oat$X.OTU.ID=="OTU_4077")
dim(oat)
#Delete row 1 and 2
oat<-oat[-c(1,2),]
oat[oat==""]<-"Unclassified"
#Find last sample column
which(colnames(oat)=="GW1984")
dim(oat)
otus<-oat[,-c(46:56)]
head(otus)
##Presence data for calculating core communities
OTUpresence<-as.data.frame(otus[,-1]>0)#turns into logical value, true or false if present or not
OTUpresence<-cbind(oat[1], OTUpresence,oat[50:56])
head(oat)
head(OTUpresence)
dim(OTUpresence)
#look at phyla distrib
countsperphyla<-OTUpresence%>%
  group_by(Phylum)%>%
  tally()
View(countsperphyla)
#sum of presence for each sp
OTUpresence<-OTUpresence %>%
  rowwise() %>%
  mutate(
    sumall =sum(c_across(2:45))/44,
    sumtau = sum(c_across(c(2:12,45)))/12,
    sumtme = sum(c_across(c(13:23,44,43)))/13,
    sumtci = sum(c_across(24:42))/19
  )
######
#####Inspect number of OTUs per sp and spread of richness
OTUntau<-OTUpresence%>%
  count(OTUpresence$sumtau>0)
OTUntme<-OTUpresence%>%  
  count(OTUpresence$sumtme>0)
OTUntci<-OTUpresence%>%  
  count(OTUpresence$sumtci>0)
#####################1018 OTUs pre filter
otun_all<-OTUpresence%>%
  mutate("which_sp"="all")
nrow(otun_all)#295
otun_tau<-OTUpresence%>%
  filter(sumtau>0)%>%
  mutate("which_sp"="T. au")#227
nrow(otun_tau)
otun_tme<-OTUpresence%>%
  filter(sumtme>0)%>%
  mutate("which_sp"="T. me")
nrow(otun_tme)#153
otun_tci<-OTUpresence%>%
  filter(sumtci>0)%>%
  mutate("which_sp"="T. ci")#136
nrow(otun_tci)
####
library(VennDiagram)
head(otun_tau)
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(otun_tau$`otus[, 1]`, otun_tme$`otus[, 1]`, otun_tci$`otus[, 1]`),
  category.names = c("T. aurantium" , "T. meloni " , "T. citrina"),
  filename = 'Figures/18splots/spreadofrichness_venn_diagramm.png',
  output=TRUE,
  imagetype="png",
  height = 500 , 
  width = 500 , 
  resolution = 300,
  compression = "lzw",
  main="18s OTU richness across and between species",
  main.fontfamily="sans",
  main.cex=0.35,
  main.fontface = "bold",
  
  # Circles
  lwd = 0.5,
  lty = 1,
  fill = myCol,
  
  # Numbers
  cex = 0.4,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.4,
  cat.fontface = "italic",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)
###############
########
#Core communities
core_community<-OTUpresence%>%
  dplyr::filter(sumtau>=0.90 |sumtme>=0.90 | sumtci>=0.90)
nrow(core_community)#13
#OR
cc_all<-OTUpresence%>%
  filter(sumall>= 0.90)%>%
  mutate("which_sp"="all")#0
nrow(cc_all)
cc_tau<-OTUpresence%>%
  filter(sumtau>=0.90)%>%
  mutate("which_sp"="T. au")#11
nrow(cc_tau)
cc_tme<-OTUpresence%>%
  filter(sumtme>= 0.90)%>%
  mutate("which_sp"="T. me")
nrow(cc_tme)#3
cc_tci<-OTUpresence%>%
  filter(sumtci>= 0.90)%>%
  mutate("which_sp"="T. ci")
nrow(cc_tci)#0
head(cc_tau)
#####
#Ven Diagram
library(VennDiagram)

myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(cc_tau$X.OTU.ID, cc_tme$X.OTU.ID, cc_tci$X.OTU.ID),
  category.names = c("T. aurantium" , "T. meloni " , "T. citrina"),
  filename = 'Figures/18splots/core18s_vendiagram.png',
  output=TRUE,
  imagetype="png",
  height = 500 , 
  width = 500 , 
  resolution = 300,
  compression = "lzw",
  main="18s core communities across and between species",
  main.fontfamily="sans",
  main.cex=0.35,
  main.fontface = "bold",
  
  # Circles
  lwd = 0.5,
  lty = 1,
  fill = myCol,
  
  # Numbers
  cex = 0.4,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.4,
  cat.fontface = "italic",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

