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
countsDF<-read.csv(uncorrectedCountsPath, sep=";")
head(otuTaxonomy)
data<-cbind(countsDF,otuTaxonomy)
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
dim(xdata)
xdata<-xdata[xdata$countSum >=50,]
dim(xdata)
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

OTUPresence[OTUPresence==""]<-"Unclasssified"
#Inspect number of OTU per phyla overall
countsperphyla<-OTUPresence%>%
  group_by(Phylum)%>%
  tally()
######
#####Inspect number of OTUs per sp and spread of richness
OTUntau<-OTUPresence%>%
  count(OTUPresence$sumtau>0)
OTUntme<-OTUPresence%>%  
  count(OTUPresence$sumtme>0)
OTUntci<-OTUPresence%>%  
  count(OTUPresence$sumtci>0)
####################
otun_all<-OTUPresence%>%
  mutate("which_sp"="all")#411 OTUs
nrow(otun_all)
otun_tau<-OTUPresence%>%
  filter(sumtau>0)%>%
  mutate("which_sp"="T. au")#328
nrow(otun_tau)
otun_tme<-OTUPresence%>%
  filter(sumtme>0)%>%
  mutate("which_sp"="T. me")#279
nrow(otun_tme)
otun_tci<-OTUPresence%>%
  filter(sumtci>0)%>%
  mutate("which_sp"="T. ci")#384
####
library(VennDiagram)

myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(otun_tau$sequence_identifier, otun_tme$sequence_identifier, otun_tci$sequence_identifier),
  category.names = c("T. aurantium" , "T. meloni " , "T. citrina"),
  filename = 'Figures/AS_16splots/spreadofrichness_venn_diagramm.png',
  output=TRUE,
  imagetype="png",
  height = 500 , 
  width = 500 , 
  resolution = 300,
  compression = "lzw",
  main="16s OTU richness across and between species",
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
#define core community as present in 90% of samples
#calculate core communities
core_community<-OTUPresence%>%
  dplyr::filter(sumtau>=0.90 |sumtme>=0.90 | sumtci>=0.90) #52 in total
list(core_community$sequence_identifier)
nrow(core_community)
View(core_community)
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
  main="16s Core communities across and between species",
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

#Calculating the relative abundance of phyla cross species

#Calculating the relative abundance of phyla cross species

#FIll blank boxes with "unclassified"
head(xdata)
xdata[xdata==""]<-"Unclasssified"
unique(xdata$Phylum)
ydata<-xdata%>%
  rowwise() %>%
  mutate(
    sumtau = sum(c_across(c(2:12,45))),
    sumtme = sum(c_across(c(13:23,43,44))),
    sumtci = sum(c_across(24:42))
  )
ydata<-ydata[,c(46,48,54:56)]
ydatamelt<-melt(ydata)
head(ydatamelt)
#Melt data set and remove unnecassary rows
#Set levels

levels(ydatamelt$variable)<-c("sumtau","sumtme","sumtci")
levels(ydatamelt$variable)<-c("Tau","Tme","Tci")
#Group by species and calc relative abundance of each OTU for each species
###THE ISSUE HERE
xdatamelt<-ydatamelt%>%
  group_by(variable)%>%
  summarise(relative_abundance=value/sum(value))
dim(xdatamelt)
head(xdatamelt)
#Bind the "relative_abundance" col (for each sp seperately) to the main data set
dataset<-cbind(ydatamelt,xdatamelt[2])
#creates a bunch of 0 values because creates a category for every OTU and sp


#DAta with just the core OTUs and also with just non core OTUs
##core
coredata<-xdata[xdata$X.OTU.ID %in% core_community$X.OTU.ID,]
#Prev line selects all core OTUs for all species so need to calculate
#percentage of samples present again to filter later
coredatamelt1<-coredata%>%
  rowwise() %>%
  mutate(
    tau = sum(c_across(c(2:12,45))),
    tme = sum(c_across(c(13:23,43,44))),
    tci = sum(c_across(24:42))
  )
coredatamelt2<-melt(coredatamelt1[,c(1,46,48,54:56)])
dim(coredatamelt2)                   
OTUpresence<-as.data.frame(coredata[,2:45]>0)#turns into logical value, true or false if present or not
OTUpresence<-cbind(coredata[,1],OTUpresence)
cclimits<-OTUpresence%>%
  rowwise() %>%
  mutate(
    sumtau = sum(c_across(c(2:12,45))/12),
    sumtme = sum(c_across(c(13:23,43,44))/13),
    sumtci = sum(c_across(24:42)/19))
cclimits2<-cclimits[,c(1,46:48)]
cclimitmelt<-melt(cclimits2)
dim(cclimitmelt)
names(cclimitmelt)[names(cclimitmelt) == 'value'] <- 'limit'
names(cclimitmelt)[names(cclimitmelt) == 'variable'] <- 'sample'

coredatacomb<-cbind(coredatamelt2,cclimitmelt)
coredatamelt<-coredatacomb%>%
  filter(limit>=0.9)
coredatamelt<-coredatamelt[,-c(2,6:8)]
head(coredatamelt)
dim(coredatamelt)
###################

#Set levels
levels(coredatamelt$variable)<-c("Tau","Tme","Tci")

##for core dataset
#Group by species and calc relative abundance of each OTU for each species
xcoredatamelt<-coredatamelt%>%
  group_by(variable)%>%
  summarise(relative_abundance=value/sum(value))
#Bind the "relative_abundance" col (for each sp seperately) to the main data set
coredataset<-cbind(coredatamelt,xcoredatamelt[2])

##non-core
noncoredata<-xdata[!(xdata$X.OTU.ID %in% coredataset$X.OTU.ID),]#HEREEE
ynoncoredata<-noncoredata%>%
  rowwise() %>%
  mutate(
    sumtau = sum(c_across(c(2:12,45))),
    sumtme = sum(c_across(c(13:23,43,44))),
    sumtci = sum(c_across(24:42))
  )

noncoredatamelt2<-melt(ynoncoredata[,c(1,48,54:56)])
OTUpresence<-as.data.frame(noncoredata[,2:45]>0)#turns into logical value, true or false if present or not
OTUpresence<-cbind(noncoredata[,1],OTUpresence)
cclimits<-OTUpresence%>%
  rowwise() %>%
  mutate(
    sumtau = sum(c_across(c(2:12,45))/12),
    sumtme = sum(c_across(c(13:23,43,44))/13),
    sumtci = sum(c_across(24:42)/19))
cclimits2<-cclimits[,c(1,46:48)]
cclimitmelt<-melt(cclimits2)
dim(cclimitmelt)
names(cclimitmelt)[names(cclimitmelt) == 'value'] <- 'limit'
names(cclimitmelt)[names(cclimitmelt) == 'variable'] <- 'sample'

noncoredatacomb<-cbind(noncoredatamelt2,cclimitmelt)
noncoredatamelt<-noncoredatacomb%>%
  filter(!limit>=0.9)
noncoredatamelt<-noncoredatamelt[,-c(5:7)]

#Set levels
head(noncoredatamelt)
levels(noncoredatamelt$variable)<-c("Tau","Tme","Tci")

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
BdellovibrionotaCol<-"darkslategray"
CalditrichotaCol<-"tomato3"
CampylobacterotaCol<-"wheat1"
ChloroflexiCol<-"violetred2"
CrenarchaeotaCol<-"violetred"
CyanobacteriaCol<-"violet"
DadabacteriaCol<-"turquoise4"
DeferrisomatotaCol<-"turquoise3"
DeinococcotaCol<-"turquoise"
#DependentiaeCol<-"tomato3"                 
DesulfobacterotaCol<-"tomato1"   
EntotheonellaeotaCol<-"Violetred4"   
#FirmicutesCol<-"black"                    
#FusobacteriotaCol<-"rosybrown3"               
GemmatimonadotaCol<-"tan3"
#HydrogenedentesCol<-"tan1"
LatescibacterotaCol<-"steelblue4"
MyxococcotaCol<-"steelblue2"
NB1jCol<-"springgreen4"
NitrospinotaCol<-"springgreen2"
NitrospirotaCol<-"slategray2"
PAUC34fCol<-"darkkhaki"
PlanctomycetotaCol<-"purple4"
ProteobacteriaCol<-"salmon"
SAR324cladeCol<-"royalblue"
SpirochaetotaCol<-"tan4"
#SumerlaeotaCol<-"red2"
UnclassifiedCol<-"slategray"
VerrucomicrobiotaCol<-"purple"
WPS2Col<-"plum3"
#WS1Col<-"orange"
######
library
unique(dataset$Phylum)
head(dataset)
#for reordering
#dataset_n<-dataset
#dataset_n$Phylum<-factor(dataset_n$Phylum,levels=c(unique(dataset$Phylum)))
#change order to : values = c(ProteobacteriaCol,UnclassifiedCol,CrenarchaeotaCol,ActinobacteriotaCol,BacteroidotaCol,CyanobacteriaCol,PlanctomycetotaCol,VerrucomicrobiotaCol,DesulfobacterotaCol,AcidobacteriotaCol,NitrospirotaCol,SpirochaetotaCol,NitrospinotaCol,DeferrisomatotaCol,DadabacteriaCol,MyxococcotaCol,SAR324cladeCol,NB1jCol,DeinococcotaCol,ChloroflexiCol,CalditrichotaCol,EntotheonellaeotaCol,LatescibacterotaCol,GemmatimonadotaCol))
  
png("./Figures/AS_16splots/relativeabundanceallphyla.png",height=10,width=20,units="cm",res=300)
ggplot(dataset, aes(variable, y=relative_abundance,fill=Phylum)) + geom_bar(aes(fill=Phylum),stat="identity", position="fill") + 
  scale_fill_manual(values = c(AcidobacteriotaCol, ActinobacteriotaCol, BacteroidotaCol,CalditrichotaCol,ChloroflexiCol, CrenarchaeotaCol,
                               CyanobacteriaCol, DadabacteriaCol,DeferrisomatotaCol, DeinococcotaCol, DesulfobacterotaCol, EntotheonellaeotaCol,
                               GemmatimonadotaCol,LatescibacterotaCol, MyxococcotaCol, NB1jCol ,NitrospinotaCol, NitrospirotaCol,PlanctomycetotaCol,
                               ProteobacteriaCol, SAR324cladeCol, SpirochaetotaCol, UnclassifiedCol,VerrucomicrobiotaCol)) +
  labs(x="Species",y="Relative abundance")+
  scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
dev.off()
#Core OTUs
unique(coredataset$Phylum)
#for reordering
#coredataset_n<-coredataset
#coredataset_n$Phylum<-factor(coredataset_n$Phylum,levels=c(unique(coredataset$Phylum)))
#order for cols: CrenarchaeotaCol,UnclassifiedCol, ActinobacteriotaCol, BacteroidotaCol, CyanobacteriaCol,PlanctomycetotaCol, ProteobacteriaCol, DesulfobacterotaCol, NitrospirotaCol,AcidobacteriotaCol,DadabacteriaCol,VerrucomicrobiotaCol

png("./Figures/AS_16splots/relativeabundance-corephyla.png",height=10,width=20,units="cm",res=300)
ggplot(coredataset, aes(variable, y=relative_abundance)) + geom_bar(aes(fill=Phylum), stat="identity", position="fill") + 
  scale_fill_manual(values = c(AcidobacteriotaCol, ActinobacteriotaCol, BacteroidotaCol,CrenarchaeotaCol, CyanobacteriaCol, DadabacteriaCol, DesulfobacterotaCol,NitrospirotaCol,PlanctomycetotaCol, ProteobacteriaCol,UnclassifiedCol, VerrucomicrobiotaCol)) +
  labs(x="Species",y="Relative abundance")+
  scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
dev.off()
#####
##Grid panel for stacked bars
legend<-get_legend(ggplot(dataset, aes(variable, y=relative_abundance,fill=Phylum)) + geom_bar(aes(fill=Phylum),stat="identity", position="fill") + 
  scale_fill_manual(values = c(AcidobacteriotaCol, ActinobacteriotaCol, BacteroidotaCol,CalditrichotaCol,ChloroflexiCol, CrenarchaeotaCol,
                               CyanobacteriaCol, DadabacteriaCol,DeferrisomatotaCol, DeinococcotaCol, DesulfobacterotaCol, EntotheonellaeotaCol,
                               GemmatimonadotaCol,LatescibacterotaCol, MyxococcotaCol, NB1jCol ,NitrospinotaCol, NitrospirotaCol,PlanctomycetotaCol,
                               ProteobacteriaCol, SAR324cladeCol, SpirochaetotaCol, UnclassifiedCol,VerrucomicrobiotaCol)) +
  labs(x="Species",y="Relative abundance")+
  scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))+
  guides(fill=guide_legend(ncol=1)))

plot1<-ggplot(dataset, aes(variable, y=relative_abundance,fill=Phylum)) + geom_bar(aes(fill=Phylum),stat="identity", position="fill") + 
  scale_fill_manual(values = c(AcidobacteriotaCol, ActinobacteriotaCol, BacteroidotaCol,CalditrichotaCol,ChloroflexiCol, CrenarchaeotaCol,
                               CyanobacteriaCol, DadabacteriaCol,DeferrisomatotaCol, DeinococcotaCol, DesulfobacterotaCol, EntotheonellaeotaCol,
                               GemmatimonadotaCol,LatescibacterotaCol, MyxococcotaCol, NB1jCol ,NitrospinotaCol, NitrospirotaCol,PlanctomycetotaCol,
                               ProteobacteriaCol, SAR324cladeCol, SpirochaetotaCol, UnclassifiedCol,VerrucomicrobiotaCol)) +
  labs(x="Species",y="Relative abundance")+ggtitle("Including all OTUs")+
  scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9),legend.position="none")

plot2<-ggplot(coredataset, aes(variable, y=relative_abundance)) + geom_bar(aes(fill=Phylum), stat="identity", position="fill") + 
  scale_fill_manual(values = c(AcidobacteriotaCol, ActinobacteriotaCol, BacteroidotaCol,CrenarchaeotaCol, CyanobacteriaCol, DadabacteriaCol, DesulfobacterotaCol,NitrospirotaCol,PlanctomycetotaCol, ProteobacteriaCol,UnclassifiedCol, VerrucomicrobiotaCol)) +
  labs(x="Species",y="Relative abundance")+ggtitle("Including only core OTUs")+
  scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9),legend.position="none")

png("Figures/AS_16splots/stackedbar16s.png",width=14,height=7,units="in",res=300)
plot_grid(plot1,plot2,legend,nrow=1,rel_widths=c(2/5,2/5,1/5))
dev.off()
#non core phyla
unique(noncoredataset$Phylum)
#reordering
#noncoredataset_n<-noncoredataset
#noncoredataset_n$Phylum<-factor(noncoredataset_n$Phylum,levels=c(unique(noncoredataset$Phylum)))
#ProteobacteriaCol,UnclassifiedCol,ActinobacteriotaCol,BacteroidotaCol,VerrucomicrobiotaCol,DesulfobacterotaCol,AcidobacteriotaCol,SpirochaetotaCol,CyanobacteriaCol,PlanctomycetotaCol,NitrospinotaCol,CrenarchaeotaCol,DeferrisomatotaCol,MyxococcotaCol,SAR324cladeCol,NB1jCol,DeinococcotaCol,NitrospirotaCol,DadabacteriaCol,ChloroflexiCol,CalditrichotaCol,EntotheonellaeotaCol,LatescibacterotaCol,GemmatimonadotaCol

png("./Figures/AS_16splots/relativeabundance-non-corephyla.png",height=10,width=20,units="cm",res=300)
ggplot(noncoredataset, aes(variable, y=relative_abundance)) + geom_bar(aes(fill=Phylum), stat="identity", position="fill") + 
  scale_fill_manual(values = c(AcidobacteriotaCol, ActinobacteriotaCol, BacteroidotaCol,CalditrichotaCol,ChloroflexiCol, CrenarchaeotaCol,
                               CyanobacteriaCol, DadabacteriaCol,DeferrisomatotaCol, DeinococcotaCol, DesulfobacterotaCol, EntotheonellaeotaCol,
                               GemmatimonadotaCol,LatescibacterotaCol, MyxococcotaCol, NB1jCol ,NitrospinotaCol, NitrospirotaCol,PlanctomycetotaCol,
                               ProteobacteriaCol, SAR324cladeCol, SpirochaetotaCol, UnclassifiedCol,VerrucomicrobiotaCol)) +
  labs(x="Species",y="Relative abundance")+
  scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
dev.off()
#TO look at OTUs across sp and sample in the core OTUs
head(coredata)
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
nmcore_otu<-mcore_otu
nmcore_otu$Phylum<-factor(nmcore_otu$Phylum,levels=c(unique(nmcore_otu$Phylum)))
head(nmcore_otu)
nmcore_otu<-nmcore_otu[!nmcore_otu$value == 0,]

png("./Figures/AS_16splots/relativeabundance-coreotus-betweensp.png",height=20,width=30,units="cm",res=300)
ggplot(nmcore_otu, aes(x=X.OTU.ID, y=value)) + geom_boxplot(aes(fill=Phylum)) + 
  geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(variable)) + 
  scale_fill_manual(values = c(CrenarchaeotaCol,UnclassifiedCol, ActinobacteriotaCol, BacteroidotaCol, CyanobacteriaCol,PlanctomycetotaCol, ProteobacteriaCol, DesulfobacterotaCol, NitrospirotaCol,AcidobacteriotaCol,DadabacteriaCol,VerrucomicrobiotaCol)) +
  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1),legend.position="right",legend.direction="vertical",axis.title.x = element_text(vjust=-2))+
  facet_wrap(~variable,nrow=3)

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
ncore_high<-core_high
ncore_high$Phylum<-factor(ncore_high$Phylum,levels=c(unique(ncore_high$Phylum)))
ncore_high<-ncore_high[!ncore_high$value == 0,]

png("./Figures/AS_16splots/highabundance-coreotusv1.png",height=20,width=20,units="cm",res=300)
ggplot(ncore_high, aes(X.OTU.ID, y=value)) + 
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(rows=vars(variable)) +
  scale_fill_manual(values = c(CrenarchaeotaCol,UnclassifiedCol,PlanctomycetotaCol,CyanobacteriaCol,ActinobacteriotaCol, BacteroidotaCol, ProteobacteriaCol, NitrospirotaCol)) +
  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1),axis.title.x = element_text(vjust=-2))
dev.off()
png("./Figures/AS_16splots/highabundance-coreotusv2.png",height=10,width=30,units="cm",res=300)
ggplot(ncore_high, aes(variable, y=value)) + geom_boxplot(aes(fill=Phylum))+ geom_jitter(colour="gray80", alpha=0.65) +
  facet_grid(cols=vars(X.OTU.ID)) + scale_fill_manual(values = c(CrenarchaeotaCol,UnclassifiedCol,PlanctomycetotaCol,CyanobacteriaCol,ActinobacteriotaCol, BacteroidotaCol, ProteobacteriaCol, NitrospirotaCol)) +
  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1),axis.title.x = element_text(vjust=-2))
dev.off()
#LOw abundance- the remaining cc's

core_low<-mcore_otu%>%
  arrange(desc(value))%>%
  filter(value<=0.1)
unique(core_low$Phylum)
ncore_low<-core_low
ncore_low$Phylum<-factor(ncore_low$Phylum, levels=c(unique(ncore_low$Phylum)))
ncore_low<-ncore_low[!ncore_low$value == 0,]
png("./Figures/AS_16splots/lowabundance-coreotusv1.png",height=20,width=30,units="cm",res=300)
ggplot(ncore_low, aes(X.OTU.ID, y=value)) +
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) +
  facet_grid(cols=vars(variable)) +  facet_wrap(~variable,nrow=3)+
  scale_fill_manual(values = c(UnclassifiedCol,BacteroidotaCol,PlanctomycetotaCol,CyanobacteriaCol,
                               ProteobacteriaCol,DesulfobacterotaCol,ActinobacteriotaCol,CrenarchaeotaCol,
                               NitrospirotaCol,AcidobacteriotaCol,DadabacteriaCol,VerrucomicrobiotaCol)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1),axis.title.x = element_text(vjust=-2),legend.position="right",legend.direction="vertical")
dev.off()
png("./Figures/AS_16splots/lowabundance-coreotusv2.png",height=10,width=80,units="cm",res=300)
ggplot(ncore_low, aes(variable, y=value)) +
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) + facet_grid(cols=vars(X.OTU.ID)) +
  scale_fill_manual(values = c(UnclassifiedCol,BacteroidotaCol,PlanctomycetotaCol,CyanobacteriaCol,
                               ProteobacteriaCol,DesulfobacterotaCol,ActinobacteriotaCol,CrenarchaeotaCol,
                               NitrospirotaCol,AcidobacteriotaCol,DadabacteriaCol,VerrucomicrobiotaCol)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1),axis.title.x = element_text(vjust=-2),legend.position="bottom",legend.direction="horizontal")
dev.off()
??fct_reorder
#core
library(forcats)
unique(mcore_otu$Phylum)
png("./Figures/AS_16splots/coreotusbargraph.png",height=10,width=20,units="cm",res=300)
ggplot(nmcore_otu, aes(x=fct_reorder(as.factor(X.OTU.ID),value,.desc=TRUE), y=value)) + geom_bar(aes(fill=Phylum), stat="identity") +
  scale_fill_manual(values = c(CrenarchaeotaCol,UnclassifiedCol, ActinobacteriotaCol, BacteroidotaCol, CyanobacteriaCol,PlanctomycetotaCol, ProteobacteriaCol, DesulfobacterotaCol, NitrospirotaCol,AcidobacteriotaCol,DadabacteriaCol,VerrucomicrobiotaCol)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1))
dev.off()
################Not just core   OTUs
#high abund for OTUs
dataprop<-xdata
dataprop[,2:45]<-apply(xdata[,2:45],2,function(x) x/sum(x))

m_otu<-melt(dataprop)

levels(m_otu$variable)<-c("GW1941","GW1942" ,"GW1943","GW1944", "GW1945", "GW1946", "GW1947", "GW1948","GW1949","GW1950","GW1951",
                              "GW1952","GW1953", "GW1954", "GW1955", "GW1956", "GW1957", "GW1958", "GW1959","GW1960","GW1961","GW1962",
                              "GW1963","GW1964","GW1965", "GW1966", "GW1967", "GW1968", "GW1969","GW1970","GW1971","GW1972", "GW1973","GW1974", "GW1975", "GW1976", "GW1977", "GW1978", "GW1979","GW1980","GW1981", "GW1982","GW1983", "GW1984")

levels(m_otu$variable)<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau",
                              "Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme",
                              "Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")


unique(m_otu$Phylum)
nm_otu<-m_otu
nm_otu$Phylum<-factor(nm_otu$Phylum,levels=c(unique(nm_otu$Phylum)))
head(nm_otu)
data_high<-m_otu%>%
  arrange(desc(value))%>%
  filter(value>0.05)
unique(data_high$Phylum)
ndata_high<-data_high
ndata_high$Phylum<-factor(ndata_high$Phylum,levels=c(unique(ndata_high$Phylum)))
head(ndata_high)
ndata_high<-ndata_high[!ndata_high$value == 0,]

png("./Figures/AS_16splots/highabundance_otusgraph.png",height=20,width=30,units="cm",res=300)

ggplot(ndata_high, aes(sequence_identifier, y=value)) +
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) +
  facet_grid(cols=vars(variable)) +  facet_wrap(~variable,nrow=3)+
  scale_fill_manual(values = c( ProteobacteriaCol, UnclassifiedCol,CrenarchaeotaCol, 
                                BacteroidotaCol,CyanobacteriaCol,ActinobacteriotaCol)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1),axis.title.x = element_text(vjust=-2),legend.position="right",legend.direction="vertical")+
  labs(x="OTU",y="Relative abundance",title="OTUs > 0.05 Relative abundance")
dev.off()

###Low abundance
data_low<-m_otu%>%
  arrange(desc(value))%>%
  filter(value<=0.05)
unique(data_low$Phylum)
ndata_low<-data_low
ndata_low$Phylum<-factor(ndata_low$Phylum,levels=c(unique(ndata_low$Phylum)))
unique(ndata_low$Phylum)
ndata_low<-ndata_low[!ndata_low$value == 0,]
png("./Figures/AS_16splots/lowabundance_otusgraphlong.png",height=30,width=100,units="cm",res=300)

ggplot(ndata_low, aes(sequence_identifier, y=value)) +
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) +
  facet_grid(cols=vars(variable)) +  facet_wrap(~variable,nrow=3)+
  scale_fill_manual(values = c(ActinobacteriotaCol,ProteobacteriaCol,UnclassifiedCol,VerrucomicrobiotaCol,PlanctomycetotaCol,
                               BacteroidotaCol,AcidobacteriotaCol,CyanobacteriaCol,DesulfobacterotaCol,
                               NitrospirotaCol, CrenarchaeotaCol,NitrospinotaCol,
                               DadabacteriaCol,MyxococcotaCol,DeferrisomatotaCol,CalditrichotaCol,SpirochaetotaCol,
                               ChloroflexiCol,DeinococcotaCol,NB1jCol,EntotheonellaeotaCol,SAR324cladeCol,LatescibacterotaCol,GemmatimonadotaCol)) + theme_bw() +guides(fill=guide_legend(ncol=1))+
  theme(axis.text.x=element_text(angle=90,hjust=1,size=5),axis.title.x = element_text(vjust=-2),legend.position="right",legend.direction="vertical")+
  labs(x="OTU",y="Relative abundance",title="OTUs < 0.05 Relative abundance")
dev.off()

png("./Figures/AS_16splots/lowabundance_otusgraph.png",height=30,width=60,units="cm",res=300)

ggplot(ndata_low, aes(sequence_identifier, y=value)) +
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) +
  facet_grid(cols=vars(variable)) +  facet_wrap(~variable,nrow=3)+
  scale_fill_manual(values = c(ActinobacteriotaCol,ProteobacteriaCol,UnclassifiedCol,VerrucomicrobiotaCol,PlanctomycetotaCol,
                               BacteroidotaCol,AcidobacteriotaCol,CyanobacteriaCol,DesulfobacterotaCol,
                               NitrospirotaCol, CrenarchaeotaCol,NitrospinotaCol,
                               DadabacteriaCol,MyxococcotaCol,DeferrisomatotaCol,CalditrichotaCol,SpirochaetotaCol,
                               ChloroflexiCol,DeinococcotaCol,NB1jCol,EntotheonellaeotaCol,SAR324cladeCol,LatescibacterotaCol,GemmatimonadotaCol)) + theme_bw() +guides(fill=guide_legend(ncol=1))+
  theme(axis.text.x=element_text(angle=90,hjust=1,size=3),axis.title.x = element_text(vjust=-2),legend.position="right",legend.direction="vertical")+
  labs(x="OTU",y="Relative abundance",title="OTUs < 0.05 Relative abundance")
dev.off()
??theme
#continues into AS_NLbubbleplots.R

