library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(vegan)
library(reshape2)
#ENter data and order the data frames
taxa<-read.csv("./Data/18s/taxa_fixed18s.csv",sep=",")
unique(taxa$Phylum)
otus<-read.csv("./Data/18s/all.otutab.csv", sep="\t")
otus<-arrange(otus, X.OTU.ID)
taxa<-arrange(taxa,sequence_identifier)

##Presence data for calculating core communities
OTUpresence<-as.data.frame(otus[,-1]>0)#turns into logical value, true or false if present or not
OTUpresence<-cbind(otus[,1], OTUpresence)
View(OTUpresence)
dim(OTUpresence)
#sum of presence for each sp
OTUpresence<-OTUpresence %>%
  rowwise() %>%
  mutate(
    sumall =sum(c_across(2:45))/44,
    sumtau = sum(c_across(2:12))/11,
    sumtme = sum(c_across(13:23))/11,
    sumtci = sum(c_across(24:45))/22
  )
#Assign taxonomy

otu_and_taxa<-bind_cols(OTUpresence[order(OTUpresence$`otus[, 1]`),], taxa[order(taxa$sequence_identifier),])
unique(otu_and_taxa$Phylum)
xotu_and_taxa<-otu_and_taxa%>%filter(Domain!="Bacteria")
unique(xotu_and_taxa$Phylum)
head(xotu_and_taxa)
#remove otu 1 and 2 because of coamplification
which(xotu_and_taxa$`otus[, 1]`=="OTU_1")
which(xotu_and_taxa$`otus[, 1]`=="OTU_2")
dim(xotu_and_taxa)
xotu_and_taxa<-xotu_and_taxa[-c(1,95),]
unique(xotu_and_taxa$Phylum)
unique(xotu_and_taxa$`otus[, 1]`)
xotu_and_taxa[xotu_and_taxa==""]<-"Unclassified"
View(xotu_and_taxa)
#Core communities
core_community<-xotu_and_taxa%>%
  dplyr::filter(sumtau>=0.90 |sumtme>=0.90 | sumtci>=0.90)
nrow(core_community)#17
#OR
cc_all<-xotu_and_taxa%>%
  filter(sumall>= 0.90)%>%
  mutate("which_sp"="all")#0
nrow(cc_all)
cc_tau<-xotu_and_taxa%>%
  filter(sumtau>=0.90)%>%
  mutate("which_sp"="T. au")#11
nrow(cc_tau)
cc_tme<-xotu_and_taxa%>%
  filter(sumtme>= 0.90)%>%
  mutate("which_sp"="T. me")
nrow(cc_tme)#7
cc_tci<-xotu_and_taxa%>%
  filter(sumtci>= 0.90)%>%
  mutate("which_sp"="T. ci")
nrow(cc_tci)#0
#####
#Ven Diagram
library(VennDiagram)

myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(cc_tau$sequence_identifier, cc_tme$sequence_identifier, cc_tci$sequence_identifier),
  category.names = c("T. aurantium" , "T. meloni " , "T. citroni"),
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
