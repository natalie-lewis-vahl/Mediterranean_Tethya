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

#New colours
AcidobacteriaCol<-"#8dd3c7"
ActinobacteriaCol<-"#ffffb3"
BacteroidetesCol<-"#bebada"
CrenarchaeotaCol<-"#fb8072"
CyanobacteriaCol<-"#80b1d3"
DesulfobacterotaCol<-"#fdb462"
NitrospirotaCol<-"#b3de69"
PlanctomycetesCol<-"#fccde5"
ProteobacteriaCol<-"#ccebc5"
VerrucomicrobiotaCol<-"#bc80bd"
UnclassifiedCol<-"#d9d9d9"

uncorrectedCountsPath<-"./Data/all.otutab_raw.csv"
#tbc bactLoadCorrectedCountsPath<-"./Data/all.otutab_Corrected_valuesOnly.csv"
otuTaxonomyPath<-"./Data/all_taxa.csv"

otuTaxonomy<-read.csv(otuTaxonomyPath, sep=";")

countsDF<-read.csv(uncorrectedCountsPath, sep="\t")

OTUPresence<-as.data.frame(countsDF[,-1]>0)#turns into logical value, true or false if present or not
OTUPresence<-cbind(countsDF[,1], OTUPresence)
dim(OTUPresence)

#sum of presence across samples across and between sp

#ASSUMING each species has an equal number of samples
x<-21 #Number of total samples
y<-7#Number of samples per species
OTUPresence<-OTUPresence %>%
  rowwise() %>%
  mutate(
    sumall =sum(c_across(2:22))/x,
    sumtau = sum(c_across(2:8))/y,
    sumtme = sum(c_across(9:15))/y,
    sumtci = sum(c_across(16:22))/y
  )
#Assign taxonomy

OTUPresence_WithTaxonomy<-bind_cols(OTUPresence[order(OTUPresence$`countsDF[, 1]`),], otuTaxonomy[order(otuTaxonomy$sequence_identifier),])


#calculate core communities

core_community<-OTUPresence_WithTaxonomy%>%
  dplyr::filter(sumtau==1.00 |sumtme==1.00 | sumtci==1.00) #76 in total
list(core_community$sequence_identifier)
#or
cc_all<-OTUPresence_WithTaxonomy%>%
  filter(sumall== 1.00)%>%
  mutate("which_sp"="all")#5 core
cc_tau<-OTUPresence_WithTaxonomy%>%
  filter(sumtau==1.00)%>%
  mutate("which_sp"="T. au")#43 core
cc_tme<-OTUPresence_WithTaxonomy%>%
  filter(sumtme== 1.00)%>%
  mutate("which_sp"="T. me")#25 core
cc_tci<-OTUPresence_WithTaxonomy%>%
  filter(sumtci== 1.00)%>%
  mutate("which_sp"="T. ci")#28 core
#Combinations of shared core bacteria
#Tau and Tme
cc_tau_tme<-OTUPresence_WithTaxonomy%>%
  filter(sumtau==1.00 & sumtme==1.00 )#8 shared
#Tau and Tci
cc_tau_tci<-OTUPresence_WithTaxonomy%>%
  filter(sumtau==1.00 & sumtci==1.00 )#11 shared
#Tme and Tci
cc_tme_tci<-OTUPresence_WithTaxonomy%>%
  filter(sumtme==1.00 & sumtci==1.00 )#6 shared
#Ven Diagram
library(VennDiagram)
??venn.diagram
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(cc_tau$sequence_identifier, cc_tme$sequence_identifier, cc_tci$sequence_identifier),
  category.names = c("T. aurantium" , "T. meloni " , "T. citroni"),
  filename = '#overlapping_cc_venn_diagramm.png',
  output=TRUE,
  imagetype="png",
  height = 500 , 
  width = 500 , 
  resolution = 300,
  compression = "lzw",
  main="Core bacterial communities across and between species",
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

#get rid of unclassified OTUs- but is this potentially important info?
#core_community<-core_community[core_community$Phylum != "",]

#Bind counts with taxonomy to calc the abundance of the core community bacteria 
#(rather than proportion of presence across samples)
countsWithTaxonomy<-bind_cols(countsDF[order(countsDF$X.OTU.ID),], otuTaxonomy[order(otuTaxonomy$sequence_identifier),])

#FIll blank boxes with "unclassified"
names(countsWithTaxonomy)

countsWithTaxonomy[countsWithTaxonomy==""]<-"Unclasssified"
#sum(countsDF[countsDF$X.OTU.ID %in% coreCommunity$`countsDF[, 1]`,]$countSum)/sum(countsDF$countSum)

coreCountsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$X.OTU.ID %in% core_community$`countsDF[, 1]`,]
corePropsWithTaxonomy<-coreCountsWithTaxonomy
corePropsWithTaxonomy[,2:22]<-apply(coreCountsWithTaxonomy[,2:22],2,function(x) x/sum(x))


corePropsWithTaxonomyMelted<-melt(corePropsWithTaxonomy)

levels(corePropsWithTaxonomyMelted$variable)<-c("GW1941","GW1942" ,"GW1944", "GW1945", "GW1946", "GW1947", "GW1948",
                                                "GW1952", "GW1953", "GW1954", "GW1955", "GW1957", "GW1958", "GW1959",
                                                "GW1964", "GW1967", "GW1968", "GW1969", "GW1973", "GW1979", "GW1982")

levels(corePropsWithTaxonomyMelted$variable)<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau",
                                                "Tme","Tme","Tme","Tme","Tme","Tme","Tme",
                                                "Tci","Tci","Tci","Tci","Tci","Tci","Tci")

names(corePropsWithTaxonomyMelted$Phylum)

ggplot(corePropsWithTaxonomyMelted, aes(variable, y=value)) + geom_bar(aes(fill=Phylum), stat="identity", position="stack") + scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol, CrenarchaeotaCol, CyanobacteriaCol, DesulfobacterotaCol, NitrospirotaCol, PlanctomycetesCol, ProteobacteriaCol,UnclassifiedCol,VerrucomicrobiotaCol)) +
  labs(x="Species",y="Relative abundance")+scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citroni)))+theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
#FOr all species
ggplot(corePropsWithTaxonomyMelted, aes(X.OTU.ID, y=value)) + geom_boxplot(aes(fill=Phylum)) + 
  geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(variable)) + 
  scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol, CrenarchaeotaCol, CyanobacteriaCol, DesulfobacterotaCol, NitrospirotaCol, PlanctomycetesCol, ProteobacteriaCol, UnclassifiedCol,VerrucomicrobiotaCol)) +
  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1),legend.position="top",legend.direction="horizontal")
#FOr just one sp?
Tau_corePropsWithTaxonomyMelted<-subset(corePropsWithTaxonomyMelted, variable=="Tau")

ggplot(Tau_corePropsWithTaxonomyMelted, aes(X.OTU.ID, y=value)) + geom_boxplot(aes(fill=Phylum)) + 
  geom_jitter(colour="gray80", alpha=0.65) + 
  scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol, CrenarchaeotaCol, CyanobacteriaCol, DesulfobacterotaCol, NitrospirotaCol, PlanctomycetesCol, ProteobacteriaCol, UnclassifiedCol,VerrucomicrobiotaCol)) +
  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1),legend.position="top",legend.direction="horizontal")

#High abundance cc (above 0.1 relative abundance in any sample)
ggplot(corePropsWithTaxonomyMelted[corePropsWithTaxonomyMelted$X.OTU.ID %in% c("OTU_1","OTU_2","OTU_7","OTU_8","OTU_9","OTU_10","OTU_14","OTU_16","OTU_20","OTU_21","OTU_24","OTU_33","OTU_35"),], aes(X.OTU.ID, y=value)) + 
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(variable)) + scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol, CrenarchaeotaCol, CyanobacteriaCol,PlanctomycetesCol, ProteobacteriaCol, UnclassifiedCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))

ggplot(corePropsWithTaxonomyMelted[corePropsWithTaxonomyMelted$X.OTU.ID %in% c("OTU_1","OTU_2","OTU_7","OTU_8","OTU_9","OTU_10","OTU_14","OTU_16","OTU_20","OTU_21","OTU_24","OTU_33","OTU_35"),], aes(variable, y=value)) + geom_boxplot(aes(fill=Phylum))+ geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(X.OTU.ID)) + scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol, CrenarchaeotaCol, CyanobacteriaCol,PlanctomycetesCol, ProteobacteriaCol, UnclassifiedCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))

#LOw abundance- the remaining cc's
ggplot(corePropsWithTaxonomyMelted[!(corePropsWithTaxonomyMelted$X.OTU.ID %in% c("OTU_1","OTU_2","OTU_7","OTU_8","OTU_9","OTU_10","OTU_14","OTU_16","OTU_20","OTU_21","OTU_24","OTU_33","OTU_35")),], aes(X.OTU.ID, y=value)) +
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(variable)) + scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol,CyanobacteriaCol, DesulfobacterotaCol, NitrospirotaCol, PlanctomycetesCol, ProteobacteriaCol, UnclassifiedCol,VerrucomicrobiotaCol)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1),legend.position="bottom",legend.direction="horizontal")

ggplot(corePropsWithTaxonomyMelted[!(corePropsWithTaxonomyMelted$X.OTU.ID %in% c("OTU_1","OTU_2","OTU_7","OTU_8","OTU_9","OTU_10","OTU_14","OTU_16","OTU_20","OTU_21","OTU_24","OTU_33","OTU_35")),], aes(variable, y=value)) +
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(X.OTU.ID)) + scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol, CyanobacteriaCol, DesulfobacterotaCol, NitrospirotaCol, PlanctomycetesCol, ProteobacteriaCol, UnclassifiedCol,VerrucomicrobiotaCol)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1),legend.position="bottom",legend.direction="horizontal")


#core
ggplot(corePropsWithTaxonomyMelted, aes(X.OTU.ID, y=value)) + geom_bar(aes(fill=Phylum), stat="identity") + scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol, CrenarchaeotaCol, CyanobacteriaCol, DesulfobacterotaCol, NitrospirotaCol, PlanctomycetesCol, ProteobacteriaCol, UnclassifiedCol,VerrucomicrobiotaCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))




