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

uncorrectedCountsPath<-"./Data/16s_allsamples_otu.csv"
#tbc bactLoadCorrectedCountsPath<-"./Data/all.otutab_Corrected_valuesOnly.csv"
otuTaxonomyPath<-"./Data/16s_allsamples_taxa.csv"

otuTaxonomy<-read.csv(otuTaxonomyPath, sep=";")
countsDF<-read.csv(uncorrectedCountsPath, sep="\t")
head(otuTaxonomy)
data<-cbind(countsDF[order(countsDF$X.OTU.ID),],otuTaxonomy[order(otuTaxonomy$sequence_identifier),])
head(data)
unique(data$Phylum)
#If  you want to filter up to 95% abuncance 
#Alculate pct abundance to filter out the last 5% 
countSum<-apply(countsDF[2:45],1,sum) #sum the rows
pct<-countSum/sum(countSum)

xdata<-cbind(data,countSum, pct)
xdata<-xdata%>%
  arrange(desc(pct))
xdata<-xdata%>%
  mutate(cumpct=cumsum(pct))
#There are more otus with more or equal to 50 reads over the
#95% cumpct abundance threshold so use at least 50 reads as 
#filter point instead

xdata<-xdata[xdata$countSum >=50,]
#Delete the cumpct and pct rows again
head(xdata)
xdata=select(xdata,-countSum,-pct,-cumpct)
#Make into presence absence data
OTUPresence<-as.data.frame(xdata[,2:45]>0)#turns into logical value, true or false if present or not
OTUPresence<-cbind(xdata[,-c(2:45)], OTUPresence)
head(OTUPresence)
which(colnames(OTUPresence)=="GW1984")#column 53
which(colnames(OTUPresence)=="GW1941")#column 10

#To calculate in how what proportion of samples each each otu is 
#separately for the species and for all sp. calc sum of reads
#across the samples and divide by the number of samples
OTUPresence<-OTUPresence %>%
  rowwise() %>%
  mutate(
    sumall =sum(c_across(10:53))/44,
    sumtau = sum(c_across(c(10:20,53)))/12,
    sumtme = sum(c_across(c(21:31,51,52)))/13,
    sumtci = sum(c_across(32:50))/19
  )
#define core community as present in 90% of samples
#calculate core communities
OTUPresence[OTUPresence==""]<-"Unclasssified"
core_community<-OTUPresence%>%
  dplyr::filter(sumtau>=0.90 |sumtme>=0.90 | sumtci>=0.90) #52 in total
list(core_community$sequence_identifier)
nrow(core_community)
#or
cc_all<-OTUPresence%>%
  filter(sumall>= 0.90)%>%
  mutate("which_sp"="all")#6 core
nrow(cc_all)
cc_tau<-OTUPresence%>%
  filter(sumtau>=0.90)%>%
  mutate("which_sp"="T. au")#9core
nrow(cc_tau)
cc_tme<-OTUPresence%>%
  filter(sumtme>=0.90)%>%
  mutate("which_sp"="T. me")#6 core
nrow(cc_tme)
cc_tci<-OTUPresence%>%
  filter(sumtci>=0.90)%>%
  mutate("which_sp"="T. ci")#48   core
nrow(cc_tci) 
#Combinations of shared core bacteria
#Tau and Tme
cc_tau_tme<-OTUPresence%>%
  filter(sumtau>=0.90 & sumtme>=0.90 )#6 shared
#Tau and Tci
cc_tau_tci<-OTUPresence%>%
  filter(sumtau>=0.90 & sumtci>=0.90 )#5 shared
#Tme and Tci
cc_tme_tci<-OTUPresence%>%
  filter(sumtme>=0.90 & sumtci>=0.90 )#4 shared
#Ven Diagram
library(VennDiagram)

myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(cc_tau$sequence_identifier, cc_tme$sequence_identifier, cc_tci$sequence_identifier),
  category.names = c("T. aurantium" , "T. meloni " , "T. citrina"),
  filename = 'Figures/AS_16splots/overlapping_cc_venn_diagramm.png',
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
datamelt<-melt(xdata[-c(1,46,47,49:53)])

#Set levels

levels(datamelt$variable)<-c("GW1941","GW1942" ,"GW1943","GW1944", "GW1945", "GW1946", "GW1947", "GW1948","GW1949","GW1950","GW1951",
                             "GW1952","GW1953", "GW1954", "GW1955", "GW1956", "GW1957", "GW1958", "GW1959","GW1960","GW1961","GW1962",
                             "GW1963","GW1964","GW1965", "GW1966", "GW1967", "GW1968", "GW1969","GW1970","GW1971","GW1972", "GW1973","GW1974", "GW1975", "GW1976", "GW1977", "GW1978", "GW1979","GW1980","GW1981", "GW1982","GW1983", "GW1984")

levels(datamelt$variable)<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau",
                             "Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme",
                             "Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")
#Group by species and calc relative abundance of each OTU for each species
xdatamelt<-datamelt%>%
  group_by(variable)%>%
  summarise(relative_abundance=value/sum(value))
dim(xdatamelt)
head(datamelt)
#Bind the "relative_abundance" col (for each sp seperately) to the main data set
dataset<-cbind(datamelt,xdatamelt[2])
#creates a bunch of 0 values because creates a category for every OTU and sp


#DAta with just the core OTUs and also with just non core OTUs
##core
coredata<-xdata[xdata$X.OTU.ID %in% core_community$X.OTU.ID,]
head(coredata)
coredatamelt<-melt(coredata[-c(1,46,47,49:53)])

levels(coredatamelt$variable)<-c("GW1941","GW1942" ,"GW1943","GW1944", "GW1945", "GW1946", "GW1947", "GW1948","GW1949","GW1950","GW1951",
                             "GW1952","GW1953", "GW1954", "GW1955", "GW1956", "GW1957", "GW1958", "GW1959","GW1960","GW1961","GW1962",
                             "GW1963","GW1964","GW1965", "GW1966", "GW1967", "GW1968", "GW1969","GW1970","GW1971","GW1972", "GW1973","GW1974", "GW1975", "GW1976", "GW1977", "GW1978", "GW1979","GW1980","GW1981", "GW1982","GW1983", "GW1984")

levels(coredatamelt$variable)<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau",
                             "Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme",
                             "Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")
##for core dataset
#Group by species and calc relative abundance of each OTU for each species
xcoredatamelt<-coredatamelt%>%
  group_by(variable)%>%
  summarise(relative_abundance=value/sum(value))
#Bind the "relative_abundance" col (for each sp seperately) to the main data set
coredataset<-cbind(coredatamelt,xcoredatamelt[2])

##non-core
noncoredata<-xdata[!(xdata$X.OTU.ID %in% core_community$X.OTU.ID),]#HEREEE
noncoredatamelt<-melt(noncoredata[-c(1,46,47,49:53)])
#
levels(noncoredatamelt$variable)<-c("GW1941","GW1942" ,"GW1943","GW1944", "GW1945", "GW1946", "GW1947", "GW1948","GW1949","GW1950","GW1951",
                             "GW1952","GW1953", "GW1954", "GW1955", "GW1956", "GW1957", "GW1958", "GW1959","GW1960","GW1961","GW1962",
                             "GW1963","GW1964","GW1965", "GW1966", "GW1967", "GW1968", "GW1969","GW1970","GW1971","GW1972", "GW1973","GW1974", "GW1975", "GW1976", "GW1977", "GW1978", "GW1979","GW1980","GW1981", "GW1982","GW1983", "GW1984")

levels(noncoredatamelt$variable)<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau",
                             "Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme",
                             "Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")


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
#BdellovibrionotaCol<-"navyblue"
CalditrichotaCol<-"tomato3"
#CampylobacterotaCol<-"ivory1"
ChloroflexiCol<-"violetred2"
CrenarchaeotaCol<-"violetred"
CyanobacteriaCol<-"violet"
DadabacteriaCol<-"turquoise4"
DeferrisomatotaCol<-"turquoise3"
DeinococcotaCol<-"turquoise"
#DependentiaeCol<-"tomato3"                 
DesulfobacterotaCol<-"tomato1"   
EntotheonellaeotaCol<-"thistle4"   
#FirmicutesCol<-"black"                    
#FusobacteriotaCol<-"rosybrown3"               
GemmatimonadotaCol<-"tan3"
#HydrogenedentesCol<-"tan1"
LatescibacterotaCol<-"steelblue4"
MyxococcotaCol<-"steelblue2"
NB1jCol<-"springgreen4"
NitrospinotaCol<-"springgreen2"
NitrospirotaCol<-"slategray2"
#PAUC34fCol<-"slateblue"
PlanctomycetotaCol<-"slategray"
ProteobacteriaCol<-"salmon"
SAR324cladeCol<-"royalblue"
SpirochaetotaCol<-"tan4"
#SumerlaeotaCol<-"red2"
UnclassifiedCol<-"purple4"
VerrucomicrobiotaCol<-"purple"
#WPS2Col<-"plum3"
#WS1Col<-"orange"
######
png("./Figures/AS_16splots/relativeabundanceallphyla.png",height=10,width=20,units="cm",res=300)
ggplot(dataset, aes(variable, y=relative_abundance)) + geom_bar(aes(fill=Phylum),stat="identity", position="fill") + 
  scale_fill_manual(values = c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,
                               CalditrichotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
                               DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DesulfobacterotaCol,
                               EntotheonellaeotaCol,GemmatimonadotaCol,
                               LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,
                               PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,SpirochaetotaCol,
                               UnclassifiedCol,VerrucomicrobiotaCol)) +
  labs(x="Species",y="Relative abundance")+
  scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
dev.off()
#Core OTUs
unique(coredataset$Phylum)
png("./Figures/AS_16splots/relativeabundance-corephyla.png",height=10,width=20,units="cm",res=300)
ggplot(coredataset, aes(variable, y=relative_abundance)) + geom_bar(aes(fill=Phylum), stat="identity", position="fill") + 
  scale_fill_manual(values = c(AcidobacteriotaCol, ActinobacteriotaCol, BacteroidotaCol, CrenarchaeotaCol, CyanobacteriaCol,DadabacteriaCol, DesulfobacterotaCol, NitrospirotaCol, PlanctomycetotaCol, ProteobacteriaCol,UnclassifiedCol,VerrucomicrobiotaCol)) +
  labs(x="Species",y="Relative abundance")+
  scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
dev.off()
#non core phyla
unique(noncoredataset$Phylum)
png("./Figures/AS_16splots/relativeabundance-non-corephyla.png",height=10,width=20,units="cm",res=300)
ggplot(noncoredataset, aes(variable, y=relative_abundance)) + geom_bar(aes(fill=Phylum), stat="identity", position="fill") + 
  scale_fill_manual(values = c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,
                               CalditrichotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
                               DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DesulfobacterotaCol,
                               EntotheonellaeotaCol,GemmatimonadotaCol,
                               LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,
                               PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,SpirochaetotaCol,
                               UnclassifiedCol,VerrucomicrobiotaCol)) +
  labs(x="Species",y="Relative abundance")+
  scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
dev.off()
#TO look at OTUs across sp and sample in the core OTUs
coredataprop<-coredata
coredataprop[,2:45]<-apply(coredata[,2:45],2,function(x) x/sum(x))

mcore_otu<-melt(coredataprop)

levels(mcore_otu$variable)<-c("GW1941","GW1942" ,"GW1943","GW1944", "GW1945", "GW1946", "GW1947", "GW1948","GW1949","GW1950","GW1951",
                                    "GW1952","GW1953", "GW1954", "GW1955", "GW1956", "GW1957", "GW1958", "GW1959","GW1960","GW1961","GW1962",
                                    "GW1963","GW1964","GW1965", "GW1966", "GW1967", "GW1968", "GW1969","GW1970","GW1971","GW1972", "GW1973","GW1974", "GW1975", "GW1976", "GW1977", "GW1978", "GW1979","GW1980","GW1981", "GW1982","GW1983", "GW1984")

levels(mcore_otu$variable)<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau",
                                    "Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme",
                                    "Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")


unique(mcore_otu$Phylum)
#FOr all species
png("./Figures/AS_16splots/relativeabundance-otus-betweensp.png",height=10,width=30,units="cm",res=300)
ggplot(mcore_otu, aes(x=X.OTU.ID, y=value)) + geom_boxplot(aes(fill=Phylum)) + 
  geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(variable)) + 
  scale_fill_manual(values = c(AcidobacteriotaCol,ActinobacteriotaCol, BacteroidotaCol, CrenarchaeotaCol, CyanobacteriaCol,DadabacteriaCol,DesulfobacterotaCol, NitrospirotaCol, PlanctomycetotaCol, ProteobacteriaCol, UnclassifiedCol,VerrucomicrobiotaCol)) +
  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1),legend.position="top",legend.direction="horizontal")
dev.off()
#FOr just one sp?
Tau_corePropsWithTaxonomyMelted<-subset(mcore_otu, variable=="Tau")
unique(Tau_corePropsWithTaxonomyMelted$Phylum)
ggplot(Tau_corePropsWithTaxonomyMelted, aes(X.OTU.ID, y=value)) + geom_boxplot(aes(fill=Phylum)) + 
  geom_jitter(colour="gray80", alpha=0.65) + 
  scale_fill_manual(values = c(AcidobacteriotaCol, ActinobacteriotaCol, BacteroidotaCol, CrenarchaeotaCol, CyanobacteriaCol, DadabacteriaCol, DesulfobacterotaCol,NitrospirotaCol, PlanctomycetotaCol, ProteobacteriaCol, UnclassifiedCol,VerrucomicrobiotaCol)) +
  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1),legend.position="top",legend.direction="horizontal")

#High abundance cc (above 0.1 relative abundance in any sample)
core_high<-mcore_otu%>%
  arrange(desc(value))%>%
  filter(value>0.1)
unique(core_high$Phylum)
png("./Figures/AS_16splots/highabundance-otusv1.png",height=10,width=30,units="cm",res=300)
ggplot(core_high, aes(X.OTU.ID, y=value)) + 
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(variable)) + scale_fill_manual(values = c(ActinobacteriotaCol, BacteroidotaCol, CrenarchaeotaCol, CyanobacteriaCol,NitrospirotaCol,PlanctomycetotaCol, ProteobacteriaCol, UnclassifiedCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))
dev.off()
png("./Figures/AS_16splots/highabundance-otusv2.png",height=10,width=30,units="cm",res=300)
ggplot(core_high, aes(variable, y=value)) + geom_boxplot(aes(fill=Phylum))+ geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(X.OTU.ID)) + scale_fill_manual(values = c(ActinobacteriotaCol, BacteroidotaCol, CrenarchaeotaCol, CyanobacteriaCol,NitrospirotaCol,PlanctomycetotaCol, ProteobacteriaCol, UnclassifiedCol)) + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1))
dev.off()
#LOw abundance- the remaining cc's

core_low<-mcore_otu%>%
  arrange(desc(value))%>%
  filter(value<=0.1)
unique(core_low$Phylum)

png("./Figures/AS_16splots/lowabundance-otusv1.png",height=10,width=30,units="cm",res=300)
ggplot(core_low, aes(X.OTU.ID, y=value)) +
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(variable)) + scale_fill_manual(values = c(AcidobacteriotaCol,ActinobacteriotaCol,DadabacteriaCol,DesulfobacterotaCol, BacteroidotaCol,CrenarchaeotaCol, CyanobacteriaCol, NitrospirotaCol, PlanctomycetotaCol, ProteobacteriaCol, UnclassifiedCol,VerrucomicrobiotaCol)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1),legend.position="bottom",legend.direction="horizontal")
dev.off()
png("./Figures/AS_16splots/lowabundance-otusv2.png",height=10,width=30,units="cm",res=300)
ggplot(core_low, aes(variable, y=value)) +
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(X.OTU.ID)) + scale_fill_manual(values = c(AcidobacteriotaCol,ActinobacteriotaCol,DadabacteriaCol,DesulfobacterotaCol, BacteroidotaCol,CrenarchaeotaCol, CyanobacteriaCol, NitrospirotaCol, PlanctomycetotaCol, ProteobacteriaCol, UnclassifiedCol,VerrucomicrobiotaCol)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1),legend.position="bottom",legend.direction="horizontal")
dev.off()

#core
unique(mcore_otu$Phylum)
png("./Figures/AS_16splots/coreotusbargraph.png",height=10,width=20,units="cm",res=300)
ggplot(mcore_otu, aes(X.OTU.ID, y=value)) + geom_bar(aes(fill=Phylum), stat="identity") + scale_fill_manual(values = c(AcidobacteriotaCol,ActinobacteriotaCol, BacteroidotaCol, CrenarchaeotaCol, CyanobacteriaCol,DadabacteriaCol ,DesulfobacterotaCol, NitrospirotaCol, PlanctomycetotaCol, ProteobacteriaCol, UnclassifiedCol,VerrucomicrobiotaCol)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1))
dev.off()

#continues into AS_NLbubbleplots.R

