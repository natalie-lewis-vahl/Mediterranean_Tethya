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

plotsPath<-"./Figures/18splots"

######RElative abundance
#prep data set for calculations
oat<-bind_cols(otus[order(otus$X.OTU.ID),], taxa[order(taxa$sequence_identifier),])

#Filtering for 95 percent
#countSum<-apply(oat[2:45],1,sum)#sum the rows

#pct<-countSum/sum(countSum)

#oatsum<-cbind(oat, countSum, pct)
#oatsum[desc(oatsum$pct),]
#oatx<-oatsum%>%
#  arrange(desc(pct))
#oat<-oatx%>%
#  mutate(cumpct=cumsum(pct))
#dim(oatx)
#dim(oat)

#oat<-oat[oat$cumpct < 0.95,]

##OR get rid of reads less than 3
oat<-oat[oat$size>4,]

#FIlter out bacteria
oat<-oat%>%filter(Domain!="Bacteria")
unique(oat$Phylum)
head(oat)
#remove otu 1 and 2 because of co-amplification
#FInd which rows to delete
which(oat$X.OTU.ID=="OTU_1")
which(oat$X.OTU.ID=="OTU_2")
dim(xotu_and_taxa)
#Delete row 1 and 95
oat<-oat[-c(1,85),]
oat[oat==""]<-"Unclassified"
#Find last sample column
which(colnames(oat)=="GW1984")

#Melt table to make phylum ID
oatmelt<-melt(oat[-c(46:50,52:56)])
#Make each variable (sample names) a level
levels(oatmelt$variable)<-c("GW1941","GW1942" ,"GW1943","GW1944", "GW1945", "GW1946", "GW1947", "GW1948","GW1949","GW1950","GW1951",
                            "GW1952", "GW1953", "GW1954", "GW1955","GW1956", "GW1957", "GW1958", "GW1959","GW1960","GW1961","GW1962",
                            "GW1963","GW1964","GW1965","GW1966", "GW1967", "GW1968", "GW1969", "GW1970", "GW1971","GW1972","GW1973","GW1974", "GW1975","GW1976", "GW1977", "GW1978","GW1979", "GW1980", "GW1981", "GW1982", "GW1983", "GW1984")
#Assign species names to each sample via the levels
levels(oatmelt$variable)<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau",
                            "Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme",
                            "Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci")

#Group by species and calc relative abundance of each OTU for each species
zoatmelt<-oatmelt%>%
  group_by(variable)%>%
  summarise(relative_abundance=value/sum(value))
dim(zoatmelt)
#Bind the "relative_abundance" col (for each sp seperately) to the main data set
dataset<-cbind(oatmelt,zoatmelt[2])
unique(oat$Phylum)
##Set colours
Unclassified<-"yellowgreen"
Ciliophora<-"green4"
Platyhelminthes<-"wheat3"
Rhodophyta<-"navyblue"
Porifera<-"yellow2"
Annelida<-"violetred4"
Nematoda<-"violetred2"
Cnidaria<-"violetred"
Euglenoza<-"violet"
Bacillariophyta<-"turquoise2"
Magnoliophyta<-"turquoise"
Mollusca<-"tomato3"                 
Bryoza<-"tomato1"   
Echinodermata<-"purple"   
Kinorhyncha<-"thistle4"                    
Chlorophyta<-"tan4"               
Arthropoda<-"tan3"
Dinoflagellata<-"steelblue4"
Brachiopoda<-"steelblue2"
Chordata<-"slategray"

####

ggplot(dataset, aes(variable, y=relative_abundance)) + geom_col(aes(fill=Phylum),position="stack") +
  scale_fill_manual(values = c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,Chordata,
                               Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,Kinorhyncha,
                               Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                               Rhodophyta,Unclassified))+
  labs(x="Species",y="Relative abundance")+scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citroni)))+theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
ggsave("relative_abundances_18s.jpg",path=plotsPath)

#REpeat but without unclassified and porifera OTUS 
boatmelt<-oatmelt[oatmelt$Phylum!="Unclassified"&oatmelt$Phylum!="Porifera",]
dim(boatmelt)
#Repeat grouping and abundance calc
#Group by species and calc relative abundance of each OTU for each species
bboatmelt<-boatmelt%>%
  group_by(variable)%>%
  summarise(relative_abundance=value/sum(value))
#Bind the "relative_abundance" col (for each sp seperately) to the main data set
bdataset<-cbind(boatmelt,bboatmelt[2])
#
ggplot(bdataset, aes(variable, y=relative_abundance)) + geom_col(aes(fill=Phylum),position="stack") +
  scale_fill_manual(values = c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,Chordata,
                               Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,Kinorhyncha,
                               Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                               Rhodophyta,Unclassified))+
  labs(x="Species",y="Relative abundance")+scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citroni)))+theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
ggsave("relative_abundances_18s_wo_unclassandporifera.jpg",path=plotsPath)


unique(oat$Domain)
head(oat)
unique(oat$Class)

oatmelt<-melt(oat[-c(46:50,53:56)])
otucount<-oatmelt%>%
  group_by(Phylum)%>%
  summarise(Numb_otu=n())

ggplot(otucount, aes(x=Numb_otu,y=Phylum))+geom_col(fill=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,Chordata,
                                                                                  Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,Kinorhyncha,
                                                                                  Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                                                                                  Rhodophyta,Unclassified))+labs(x="Number of OTUs",y="Phylum")+theme_bw()
ggsave("Freq_of_otu.jpg",path=plotsPath)
#WIhtout unclassified

botu_count<-otucount[(otucount$Phylum!="Unclassified"),]

ggplot(botu_count, aes(x=Numb_otu,y=Phylum))+geom_col(fill=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,Chordata,
                                                           Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,Kinorhyncha,
                                                           Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                                                           Rhodophyta))+labs(x="Number of OTUs",y="Phylum")+theme_bw()
ggsave("Freq_of_otu_wo_unclass.jpg",path=plotsPath)

#BUbble plot for comparing, richness and abundance of phylums across species
unique(bdataset$Phylum)
otucount2<-oatmelt%>%
  group_by(Phylum,variable)%>%
  summarise(Numb_otu=n())
dim(otucount2)
relabunds<-dataset%>%
  group_by(Phylum,variable)%>%
  summarise(abundance=sum(relative_abundance))
dim(relabunds)
bubbledata<-cbind(otucount2,relabunds[3])
head(bubbledata)
#SHowing between phylums
bubbledata[bubbledata==0]<-NA
ggplot(bubbledata,aes(x=Phylum,y=Numb_otu,colour=variable,size=abundance))+geom_point()
#Showing between species
ggplot(bubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,
                                         Chordata,Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,
                                         Kinorhyncha,Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                                         Rhodophyta,Unclassified))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.7,20))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethia citrina","Tethya meloni"))+ theme(axis.text.x = element_text(face = "italic"))
ggsave("phylum_richness_and_abundance_bubble.jpg",path=plotsPath)
#WIthout unclassified, porifera and ciliophora
bubble2<-bubbledata[(bubbledata$Phylum!="Unclassified"&bubbledata$Phylum!="Porifera"),]
bubble2<-bubble2[(bubble2$Phylum!="Ciliophora"),]
ggplot(bubble2,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,
                              Chordata, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,
                              Kinorhyncha,Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,
                              Rhodophyta))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.7,20))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethia citrina","Tethya meloni"))+ theme(axis.text.x = element_text(face = "italic"))
ggsave("phylum_rich_and_abund_bubble_withoutunclass-cill-porifera.jpg",path=plotsPath)
#TO LOOK AT RICHNESS AS THE SIZE Instead and ABUNDANCE on y axis
#ggplot(bubble2,aes(x=variable,y=abundance,color=Phylum,size=Numb_otu))+
#  scale_color_manual(values=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,
 #                             Chordata,Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,
  #                            Kinorhyncha,Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
   #                           Rhodophyta,Unclassified))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  #scale_size(range = c(.7,20))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)
