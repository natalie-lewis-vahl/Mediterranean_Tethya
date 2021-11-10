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
head(countsDF)
data<-cbind(countsDF[order(countsDF$X.OTU.ID),],otuTaxonomy[order(otuTaxonomy$sequence_identifier),])
head(data)

#CAlculate pct abundance to filter out the last 5% 
countSum<-apply(countsDF[2:22],1,sum) #sum the rows
pct<-countSum/sum(countSum)

xdata<-cbind(data, pct)
xdata[desc(xdata$pct),]
xdata<-xdata%>%
  mutate(cumpct=cumsum(pct))

xdata<-xdata%>%
  arrange(desc(pct))
xdata<-xdata[xdata$cumpct < 0.95,]
#Delete the cumpct and pct rows again
head(xdata)
xdata=select(xdata,-pct,-cumpct)
#Make into presence absence data
OTUPresence<-as.data.frame(xdata[,2:22]>0)#turns into logical value, true or false if present or not
OTUPresence<-cbind(xdata[,-c(2:22)], OTUPresence)
head(OTUPresence)

#ASSUMING each species has an equal number of samples
x<-21 #Number of total samples
y<-7#Number of samples per species
OTUPresence<-OTUPresence %>%
  rowwise() %>%
  mutate(
    sumall =sum(c_across(9:29))/x,
    sumtau = sum(c_across(9:15))/y,
    sumtme = sum(c_across(16:22))/y,
    sumtci = sum(c_across(23:29))/y
  )

#calculate core communities
OTUPresence[OTUPresence==""]<-"Unclasssified"
core_community<-OTUPresence%>%
  dplyr::filter(sumtau==1.00 |sumtme==1.00 | sumtci==1.00) #72 in total
list(core_community$sequence_identifier)
nrow(core_community)
#or
cc_all<-OTUPresence%>%
  filter(sumall== 1.00)%>%
  mutate("which_sp"="all")#5 core
nrow(cc_all)
cc_tau<-OTUPresence%>%
  filter(sumtau==1.00)%>%
  mutate("which_sp"="T. au")#39core
nrow(cc_tau)
cc_tme<-OTUPresence%>%
  filter(sumtme== 1.00)%>%
  mutate("which_sp"="T. me")#25 core
nrow(cc_tme)
cc_tci<-OTUPresence%>%
  filter(sumtci== 1.00)%>%
  mutate("which_sp"="T. ci")#27 core
head(cc_tci)
#Combinations of shared core bacteria
#Tau and Tme
cc_tau_tme<-OTUPresence%>%
  filter(sumtau==1.00 & sumtme==1.00 )#8 shared
#Tau and Tci
cc_tau_tci<-OTUPresence%>%
  filter(sumtau==1.00 & sumtci==1.00 )#10 shared
#Tme and Tci
cc_tme_tci<-OTUPresence%>%
  filter(sumtme==1.00 & sumtci==1.00 )#6 shared
#Ven Diagram
library(VennDiagram)
??venn.diagram
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(cc_tau$sequence_identifier, cc_tme$sequence_identifier, cc_tci$sequence_identifier),
  category.names = c("T. aurantium" , "T. meloni " , "T. citroni"),
  filename = 'overlapping_cc_venn_diagramm.png',
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

#Calculating the raltive abundance of phyla cross species

#Calculating the relative abundance of phyla cross species

#FIll blank boxes with "unclassified"
names(xdata)

xdata[xdata==""]<-"Unclasssified"
#Melt data set and remove unnecassary rows
datamelt<-melt(xdata[-c(1,23,24,26:29)])

#DAta with just the core OTUs and also with just non core OTUs
##core
coredata<-xdata[xdata$X.OTU.ID %in% core_community$X.OTU.ID,]
head(coredata)
coredatamelt<-melt(coredata[-c(1,23,24,26:29)])
##non-core
noncoredata<-xdata[!(xdata$X.OTU.ID %in% core_community$X.OTU.ID),]#HEREEE
noncoredatamelt<-melt(noncoredata[-c(1,23,24,26:29)])
#
#Set levels

levels(datamelt$variable)<-c("GW1941","GW1942" ,"GW1944", "GW1945", "GW1946", "GW1947", "GW1948",
                             "GW1952", "GW1953", "GW1954", "GW1955", "GW1957", "GW1958", "GW1959",
                             "GW1964", "GW1967", "GW1968", "GW1969", "GW1973", "GW1979", "GW1982")

levels(datamelt$variable)<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau",
                             "Tme","Tme","Tme","Tme","Tme","Tme","Tme",
                             "Tci","Tci","Tci","Tci","Tci","Tci","Tci")

levels(coredatamelt$variable)<-c("GW1941","GW1942" ,"GW1944", "GW1945", "GW1946", "GW1947", "GW1948",
                                 "GW1952", "GW1953", "GW1954", "GW1955", "GW1957", "GW1958", "GW1959",
                                 "GW1964", "GW1967", "GW1968", "GW1969", "GW1973", "GW1979", "GW1982")

levels(coredatamelt$variable)<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau",
                                 "Tme","Tme","Tme","Tme","Tme","Tme","Tme",
                                 "Tci","Tci","Tci","Tci","Tci","Tci","Tci")
levels(noncoredatamelt$variable)<-c("GW1941","GW1942" ,"GW1944", "GW1945", "GW1946", "GW1947", "GW1948",
                                    "GW1952", "GW1953", "GW1954", "GW1955", "GW1957", "GW1958", "GW1959",
                                    "GW1964", "GW1967", "GW1968", "GW1969", "GW1973", "GW1979", "GW1982")

levels(noncoredatamelt$variable)<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau",
                                    "Tme","Tme","Tme","Tme","Tme","Tme","Tme",
                                    "Tci","Tci","Tci","Tci","Tci","Tci","Tci")
#Group by species and calc relative abundance of each OTU for each species
xdatamelt<-datamelt%>%
  group_by(variable)%>%
  summarise(relative_abundance=value/sum(value))
dim(xdatamelt)
head(xdatamelt)
#Bind the "relative_abundance" col (for each sp seperately) to the main data set
dataset<-cbind(datamelt,xdatamelt[2])

##REpeat for core dataset
#Group by species and calc relative abundance of each OTU for each species
xcoredatamelt<-coredatamelt%>%
  group_by(variable)%>%
  summarise(relative_abundance=value/sum(value))
#Bind the "relative_abundance" col (for each sp seperately) to the main data set
coredataset<-cbind(coredatamelt,xcoredatamelt[2])
##REpeat for noncore dataset
#Group by species and calc relative abundance of each OTU for each species
xnoncoredatamelt<-noncoredatamelt%>%
  group_by(variable)%>%
  summarise(relative_abundance=value/sum(value))
#Bind the "relative_abundance" col (for each sp seperately) to the main data set
noncoredataset<-cbind(noncoredatamelt,xnoncoredatamelt[2])
##########

unique(dataset$Phylum)
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
DeferrisomatotaCol<-"turquoise3"
DeinococcotaCol<-"turquoise"
DependentiaeCol<-"tomato3"                 
DesulfobacterotaCol<-"tomato1"   
EntotheonellaeotaCol<-"thistle4"   
FirmicutesCol<-"black"                    
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
head(dataset)

ggplot(dataset, aes(variable, y=relative_abundance)) + geom_bar(aes(fill=Phylum), stat="identity", position="stack") + 
  scale_fill_manual(values = c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,BdellovibrionotaCol,
                               CalditrichotaCol,CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
                               DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DependentiaeCol,DesulfobacterotaCol,
                               EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,GemmatimonadotaCol,HydrogenedentesCol,
                               LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,PAUC34fCol,
                               PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,SpirochaetotaCol,SumerlaeotaCol,
                               UnclassifiedCol,VerrucomicrobiotaCol,WPS2Col,WS1Col)) +
  labs(x="Species",y="Relative abundance")+
  scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citroni)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
#Core OTUs
ggplot(coredataset, aes(variable, y=relative_abundance)) + geom_bar(aes(fill=Phylum), stat="identity", position="stack") + 
  scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol, CrenarchaeotaCol, CyanobacteriaCol, DesulfobacterotaCol, NitrospirotaCol, PlanctomycetesCol, ProteobacteriaCol,UnclassifiedCol,VerrucomicrobiotaCol)) +
  labs(x="Species",y="Relative abundance")+
  scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citroni)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
#non core phyla
ggplot(noncoredataset, aes(variable, y=relative_abundance)) + geom_bar(aes(fill=Phylum), stat="identity", position="stack") + 
  scale_fill_manual(values = c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,BdellovibrionotaCol,
                               CalditrichotaCol,CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
                               DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DependentiaeCol,DesulfobacterotaCol,
                               EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,GemmatimonadotaCol,HydrogenedentesCol,
                               LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,PAUC34fCol,
                               PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,SpirochaetotaCol,SumerlaeotaCol,
                               UnclassifiedCol,VerrucomicrobiotaCol,WPS2Col,WS1Col)) +
  labs(x="Species",y="Relative abundance")+
  scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citroni)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))

#TO look at OTUs across sp and sample in the core OTUs
coredataprop<-coredata
coredataprop[,2:22]<-apply(coredata[,2:22],2,function(x) x/sum(x))

mcore_otu<-melt(coredataprop)

levels(mcore_otu$variable)<-c("GW1941","GW1942" ,"GW1944", "GW1945", "GW1946", "GW1947", "GW1948",
                              "GW1952", "GW1953", "GW1954", "GW1955", "GW1957", "GW1958", "GW1959",
                              "GW1964", "GW1967", "GW1968", "GW1969", "GW1973", "GW1979", "GW1982")

levels(mcore_otu$variable)<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau",
                              "Tme","Tme","Tme","Tme","Tme","Tme","Tme",
                              "Tci","Tci","Tci","Tci","Tci","Tci","Tci")
head(mcore_otu)
#FOr all species
ggplot(mcore_otu, aes(x=X.OTU.ID, y=value)) + geom_boxplot(aes(fill=Phylum)) + 
  geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(variable)) + 
  scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol, CrenarchaeotaCol, CyanobacteriaCol, DesulfobacterotaCol, NitrospirotaCol, PlanctomycetesCol, ProteobacteriaCol, UnclassifiedCol,VerrucomicrobiotaCol)) +
  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1),legend.position="top",legend.direction="horizontal")

#FOr just one sp?
Tau_corePropsWithTaxonomyMelted<-subset(mcore_otu, variable=="Tau")

ggplot(Tau_corePropsWithTaxonomyMelted, aes(X.OTU.ID, y=value)) + geom_boxplot(aes(fill=Phylum)) + 
  geom_jitter(colour="gray80", alpha=0.65) + 
  scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol, CrenarchaeotaCol, CyanobacteriaCol, DesulfobacterotaCol, NitrospirotaCol, PlanctomycetesCol, ProteobacteriaCol, UnclassifiedCol,VerrucomicrobiotaCol)) +
  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1),legend.position="top",legend.direction="horizontal")

#High abundance cc (above 0.1 relative abundance in any sample)
core_high<-mcore_otu%>%
  arrange(desc(value))%>%
  filter(value>0.1)
head(core_high)
ggplot(core_high, aes(X.OTU.ID, y=value)) + 
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(variable)) + scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol, CrenarchaeotaCol, CyanobacteriaCol,PlanctomycetesCol, ProteobacteriaCol, UnclassifiedCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))

ggplot(core_high, aes(variable, y=value)) + geom_boxplot(aes(fill=Phylum))+ geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(X.OTU.ID)) + scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol, CrenarchaeotaCol, CyanobacteriaCol,PlanctomycetesCol, ProteobacteriaCol, UnclassifiedCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))

#LOw abundance- the remaining cc's
core_low<-mcore_otu%>%
  arrange(desc(value))%>%
  filter(value<=0.1)

ggplot(core_low, aes(X.OTU.ID, y=value)) +
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(variable)) + scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol,CrenarchaeotaCol, CyanobacteriaCol, DesulfobacterotaCol, NitrospirotaCol, PlanctomycetesCol, ProteobacteriaCol, UnclassifiedCol,VerrucomicrobiotaCol)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1),legend.position="bottom",legend.direction="horizontal")

ggplot(core_low, aes(variable, y=value)) +
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(X.OTU.ID)) + scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol,CrenarchaeotaCol, CyanobacteriaCol, DesulfobacterotaCol, NitrospirotaCol, PlanctomycetesCol, ProteobacteriaCol, UnclassifiedCol,VerrucomicrobiotaCol)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1),legend.position="bottom",legend.direction="horizontal")


#core
ggplot(mcore_otu, aes(X.OTU.ID, y=value)) + geom_bar(aes(fill=Phylum), stat="identity") + scale_fill_manual(values = c(ActinobacteriaCol, BacteroidetesCol, CrenarchaeotaCol, CyanobacteriaCol, DesulfobacterotaCol, NitrospirotaCol, PlanctomycetesCol, ProteobacteriaCol, UnclassifiedCol,VerrucomicrobiotaCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))




