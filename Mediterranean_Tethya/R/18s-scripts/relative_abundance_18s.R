library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(vegan)
library(reshape2)
library(tidyr)
library(grDevices)
#ENter data and order the data frames
taxa<-read.csv("./Data/18s/taxa_fixed18s.csv",sep=";")
unique(taxa$Phylum)
otus<-read.csv("./Data/18s/all.otutab.csv", sep=";")
otus<-arrange(otus, X.OTU.ID)
taxa<-arrange(taxa,sequence_identifier)

plotsPath<-"./Figures/18splots"

######RElative abundance
#prep data set for calculations
oat<-bind_cols(otus[order(otus$X.OTU.ID),], taxa[order(taxa$sequence_identifier),])


countSum<-apply(oat[2:45],1,sum)#sum the rows
oatsum<-cbind(oat,countSum)
#Filter for OTUs with more than 5 reads
oat<-oatsum[oatsum$countSum>=5,]

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
unique(oat$Phylum)
head(oat)
#remove otu 1 and 2 because of co-amplification
#Delete row 1 and 2
oat<-oat[-c(1,2),]
oat[oat==""]<-"Unclassified"
#Find last sample column
which(colnames(oat)=="GW1984")
head(oat)
dim(oat)
#Melt table to make phylum ID
#sumarize sum per species across all samples
oatmelt1<-oat%>%
  rowwise() %>%
  mutate(
    sumtau = sum(c_across(c(2:12,45))),
    sumtme = sum(c_across(c(13:23,43,44))),
    sumtci = sum(c_across(24:42))
  )
dim(oatmelt1)
oatmelt1<-oatmelt1[,c(46,51,57:59)]
oatmelt<-melt(oatmelt1)
head(oatmelt)
#Melt data set and remove unnecassary rows
#Make each variable (sample names) a level
levels(oatmelt$variable)<-c("sumtau","sumtme","sumtci")
#Assign species names to each sample via the levels
levels(oatmelt$variable)<-c("Tau","Tme","Tci")
##########
zoatmelt<-oatmelt%>%
  group_by(variable)%>%
  summarise(relative_abundance=value/sum(value))
head(zoatmelt)
#Check relative abundances add up to 1 for each sp
sum(zoatmelt$relative_abundance[zoatmelt$variable=="Tci"])
#Bind the "relative_abundance" col (for each sp seperately) to the main data set
dataset<-cbind(oatmelt,zoatmelt[2])
unique(oat$Phylum)

####Core dataset
#DAta with just the core OTUs and also with just non core OTUs
##core_community: load from core 18s script
head(core_community)
coredata<-oat[oat$X.OTU.ID %in% core_community$X.OTU.ID,]
head(coredata)
######
coredatamelt1<-coredata%>%
  rowwise() %>%
  mutate(
    tau = sum(c_across(c(2:12,45))),
    tme = sum(c_across(c(13:23,43,44))),
    tci = sum(c_across(24:42))
  )
coredatamelt2<-melt(coredatamelt1[,c(1,46,51,57:59)])
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
names(cclimitmelt)[names(cclimitmelt) == 'X.OTU.ID'] <- 'otu'

coredatacomb<-cbind(coredatamelt2,cclimitmelt)
coredatamelt<-coredatacomb%>%
  filter(limit>=0.9)

coredatamelt<-coredatamelt[,-c(6:8)]
head(coredatamelt)
dim(coredatamelt)

#Make each variable (sample names) a level
levels(coredatamelt$variable)<-c("tau","tme","tci")

######

######
##for core dataset
#Group by species and calc relative abundance of each OTU for each species
xcoredatamelt<-coredatamelt%>%
  group_by(variable)%>%
  summarise(relative_abundance=value/sum(value))
#Bind the "relative_abundance" col (for each sp seperately) to the main data set
coredataset<-cbind(coredatamelt,xcoredatamelt[2])
###########
##Set colours
Unclassified<-"slategray"
Ciliophora<-"green4"
Platyhelminthes<-"wheat3"
Rhodophyta<-"navyblue"
Porifera<-"yellow2"
Annelida<-"violetred4"
Nematoda<-"violetred2"
Cnidaria<-"violetred"
Euglenoza<-"violet"
Bacillariophyta<-"turquoise2"
Magnoliophyta<-"darkgreen"
Mollusca<-"tomato3"                 
Bryoza<-"tomato1"   
Echinodermata<-"purple"   
Kinorhyncha<-"springgreen2"                    
Chlorophyta<-"tan4"               
Arthropoda<-"tan3"
Dinoflagellata<-"steelblue4"
Brachiopoda<-"steelblue2"
Chordata<-"yellowgreen"
####
head(dataset)
ggplot(dataset, aes(variable, y=value,fill=Phylum)) + geom_bar(position="fill",stat="identity")+
  scale_fill_manual(values = c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,Chordata,
                               Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,Kinorhyncha,
                               Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                               Rhodophyta,Unclassified))+
  labs(x="Species",y="Relative abundance")+ theme(legend.title=element_blank())+
scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
ggsave("relative_abundances_18s.png",path=plotsPath,dpi=300,units="cm",width=30,height=20)

#REpeat but without unclassified and porifera OTUS 
boatmelt<-oatmelt[oatmelt$Phylum!="Unclassified"&oatmelt$Phylum!="Porifera",]
dim(boatmelt)
#Repeat grouping and abundance calc
#Group by species and calc relative abundance of each OTU for each species
bboatmelt<-boatmelt%>%
  group_by(variable)%>%
  summarise(relative_abundance=value/sum(value))
#Bind the "relative_abundance" col (for each sp seperately) to the main data set
xbdataset<-cbind(boatmelt,bboatmelt[2])
#
ggplot(xbdataset, aes(variable, y=relative_abundance)) + geom_col(aes(fill=Phylum),position="fill") +
  scale_fill_manual(values = c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,Chordata,
                               Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,Kinorhyncha,
                               Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,
                               Rhodophyta))+
  labs(x="Species",y="Relative abundance")+scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
ggsave("2relative_abundances_18s_wo_unclassandporifera.png",path=plotsPath,dpi=300,units="cm",width=30,height=20)
###
#REpeat but without unclassified OTUS only
boatmelt<-oatmelt[oatmelt$Phylum!="Unclassified",]
dim(boatmelt)

#Repeat grouping and abundance calc
#Group by species and calc relative abundance of each OTU for each species
bboatmelt<-boatmelt%>%
  group_by(variable)%>%
  summarise(relative_abundance=value/sum(value))
#Bind the "relative_abundance" col (for each sp seperately) to the main data set
bdataset<-cbind(boatmelt,bboatmelt[2])
#
ggplot(bdataset, aes(variable, y=relative_abundance)) + geom_col(aes(fill=Phylum),position="fill") +
  scale_fill_manual(values = c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,Chordata,
                               Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,Kinorhyncha,
                               Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                               Rhodophyta))+
  labs(x="Species",y="Relative abundance")+scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
ggsave("2relative_abundances_18s_wo_unclass.png",path=plotsPath,dpi=300,units="cm",width=30,height=20)
###################
#coreOTUs
png("./Figures/18splots/relativeabundance-corephyla.png",height=10,width=20,units="cm",res=300)
ggplot(coredataset, aes(variable, y=relative_abundance)) + geom_bar(aes(fill=Phylum), stat="identity", position="fill") + 
  scale_fill_manual(values = c(Chlorophyta,Chordata,Cnidaria,Porifera,Unclassified)) +
  labs(x="Species",y="Relative abundance")+
  scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
dev.off()

##Without unclassified
bcoredatamelt<-coredatamelt[coredatamelt$Phylum!="Unclassified",]
#Repeat grouping and abundance calc
#Group by species and calc relative abundance of each OTU for each species
coredatamelt2<-bcoredatamelt%>%
  group_by(variable)%>%
  summarise(relative_abundance=value/sum(value))
#Bind the "relative_abundance" col (for each sp seperately) to the main data set
bcoredataset<-cbind(bcoredatamelt,coredatamelt2[2])
bcoredataset
png("./Figures/18splots/relativeabundance-corephyla_wounclass.png",height=10,width=20,units="cm",res=300)
ggplot(bcoredataset, aes(variable, y=relative_abundance)) + geom_bar(aes(fill=Phylum), stat="identity", position="fill") + 
  scale_fill_manual(values = c(Chlorophyta,Chordata,Cnidaria,Porifera)) +
  labs(x="Species",y="Relative abundance")+
  scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9))
dev.off()
####################
###gridlayout for stacked bar w and w/o unclass and porifera
legend<-get_legend(ggplot(dataset, aes(variable, y=relative_abundance)) + geom_col(aes(fill=Phylum),position="fill") +
  scale_fill_manual(values = c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,Chordata,
                               Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,Kinorhyncha,
                               Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                               Rhodophyta,Unclassified))+
  labs(x="Species",y="Relative abundance")+scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9)))

plot1<-ggplot(dataset, aes(variable, y=value,fill=Phylum)) + geom_bar(position="fill",stat="identity")+
  scale_fill_manual(values = c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,Chordata,
                               Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,Kinorhyncha,
                               Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                               Rhodophyta,Unclassified))+ggtitle("Including all OTUs")+
  labs(x="Species",y="Relative abundance")+ theme(legend.title=element_blank())+
  scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9),legend.position="none")


plot2<-ggplot(xbdataset, aes(variable, y=relative_abundance)) + geom_col(aes(fill=Phylum),position="fill") +
  scale_fill_manual(values = c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,Chordata,
                               Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,Kinorhyncha,
                               Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,
                               Rhodophyta))+ggtitle("Excluding porifera and unclassified OTUs")+
  labs(x="Species",y="Relative abundance")+scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5,size=9),legend.position="none")
###
png("Figures/18splots/stackedbar18s_wandwo_sponge.png",width=14,height=7,units="in",res=300)
plot_grid(plot1,plot2,legend,nrow=1,rel_widths=c(2/5,2/5,1/5))
dev.off()

#####
otucount<-oatmelt%>%
  group_by(Phylum)%>%
  summarise(Numb_otu=n())

ggplot(otucount, aes(x=Numb_otu,y=Phylum))+geom_col(fill=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,Chordata,
                                                                                  Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,Kinorhyncha,
                                                                                  Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                                                                                  Rhodophyta,Unclassified))+labs(x="Number of OTUs",y="Phylum")+theme_bw()
ggsave("Freq_of_otu.png",path=plotsPath,dpi=300,units="cm",width=30,height=20)
#Withtout unclassified

botu_count<-otucount[(otucount$Phylum!="Unclassified"),]
svg(filename="Figures/18splots/Freq_of_otu_wo_unclass.svg",width=30,height=10)

ggplot(botu_count, aes(y=Numb_otu,x=Phylum))+geom_col(fill=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,Chordata,
                                                           Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,Kinorhyncha,
                                                           Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                                                           Rhodophyta))+labs(x="Phylum",y="Number of OTUs")+theme_bw()

dev.off()
head(oatmelt)
###
####Box plots to show OTU differences
#high abund for OTUs
head(oat)
oat2<-oat[,-c(46:49)]
dataprop<-oat2
dataprop[,2:45]<-apply(oat2[,2:45],2,function(x) x/sum(x))

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
View(nm_otu)

data_high<-m_otu%>%
  arrange(desc(value))%>%
  filter(value>0.05)
unique(data_high$Phylum)
ndata_high<-data_high
ndata_high$Phylum<-factor(ndata_high$Phylum,levels=c(unique(ndata_high$Phylum)))
head(ndata_high)
ndata_high<-ndata_high[!ndata_high$value == 0,]
png("./Figures/18splots/highabundance_otusgraph.png",height=20,width=30,units="cm",res=300)

ggplot(ndata_high, aes(X.OTU.ID, y=value)) +
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) +
  facet_grid(cols=vars(variable)) +  facet_wrap(~variable,nrow=3)+
  scale_fill_manual(values = c(Porifera,Unclassified,Rhodophyta,Annelida,Cnidaria,Bryoza,Mollusca,Ciliophora)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,),axis.title.x = element_text(vjust=-2),legend.position="right",legend.direction="vertical")+
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
png("./Figures/18splots/lowabundance_otusgraphlongwzero.png",height=30,width=100,units="cm",res=300)

ggplot(ndata_low, aes(X.OTU.ID, y=value)) +
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) +
  facet_grid(cols=vars(variable)) +  facet_wrap(~variable,nrow=3)+
  scale_fill_manual(values = c(Unclassified,Mollusca,Cnidaria,Bryoza,Porifera,Ciliophora,Annelida,Dinoflagellata,
                               Kinorhyncha,Rhodophyta,Chlorophyta,Bacillariophyta,
                               Chordata,Arthropoda,Brachiopoda,Euglenoza,Magnoliophyta,
                               Nematoda,Platyhelminthes,Echinodermata)) + theme_bw() +guides(fill=guide_legend(ncol=1))+
  theme(axis.text.x=element_text(angle=90,hjust=1,size=5),axis.title.x = element_text(vjust=-2),legend.position="right",legend.direction="vertical")+
  labs(x="OTU",y="Relative abundance",title="OTUs < 0.05 Relative abundance")
dev.off()

png("./Figures/18splots/lowabundance_otusgraphwzero.png",height=30,width=60,units="cm",res=300)

ggplot(ndata_low, aes(X.OTU.ID, y=value)) +
  geom_boxplot(aes(fill=Phylum)) + geom_jitter(colour="gray80", alpha=0.65) +
  facet_grid(cols=vars(variable)) +  facet_wrap(~variable,nrow=3)+
  scale_fill_manual(values = c(Unclassified,Mollusca,Cnidaria,Bryoza,Porifera,Ciliophora,Annelida,Dinoflagellata,
                               Kinorhyncha,Rhodophyta,Chlorophyta,Bacillariophyta,
                               Chordata,Arthropoda,Brachiopoda,Euglenoza,Magnoliophyta,
                               Nematoda,Platyhelminthes,Echinodermata)) + theme_bw() +guides(fill=guide_legend(ncol=1))+
  theme(axis.text.x=element_text(angle=90,hjust=1,size=5),axis.title.x = element_text(vjust=-2),legend.position="right",legend.direction="vertical")+
  labs(x="OTU",y="Relative abundance",title="OTUs < 0.05 Relative abundance")
dev.off()
??theme
####
###bubble plots

####Group by phylum and also species
botucount<-dataset%>%
  filter(value>0)%>%
  group_by(Phylum,variable)%>%
  dplyr::summarize(Numb_otu=n())
dim(otucount)
relabunds<-dataset%>%
  filter(value>0)%>%
  group_by(Phylum,variable)%>%
  dplyr::summarise(abundance=sum(relative_abundance))

dim(relabunds)
bubbledata<-cbind(botucount,relabunds[3])
head(bubbledata)

#SHowing between phylums
bubbledata[bubbledata==0]<-NA
#####

botu_count<-botucount[(botucount$Phylum!="Unclassified"),]
#BUbble plot for comparing, richness and abundance of phylums across species
unique(bdataset$Phylum)
head(oat)
otucount2<-oatmelt%>%
  filter(value>0)%>%
  group_by(Phylum,variable)%>%
  dplyr::summarise(Numb_otu=n())
dim(otucount2)
relabunds<-dataset%>%
  filter(value>0)%>%
  group_by(Phylum,variable)%>%
  summarise(abundance=sum(relative_abundance))
dim(relabunds)
bubbledata<-cbind(otucount2,relabunds[3])
#SHowing between phylums
bubbledata[bubbledata==0]<-NA
ggplot(bubbledata,aes(x=Phylum,y=Numb_otu,colour=variable,size=abundance))+geom_point()
#Showing between species
ggplot(bubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,
                                         Chordata,Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,
                                         Kinorhyncha,Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                                         Rhodophyta,Unclassified))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethya meloni","Tethia citrina"))+
  theme(axis.text.x = element_text(face = "italic"),legend.margin = margin(0,0,0,-2))+
  guides(color = guide_legend(override.aes = list(size = 5),ncol=2 ),size=guide_legend(nrow=1) )

ggsave("phylum_richness_and_abundance_bubble.png",path=plotsPath,dpi=300,units="cm",width=30,height=20)
##################
#Without unclassified
dim(bubbledata)
xbubbledata<-bubbledata[(!bubbledata$Phylum=="Unclassified"),]
dim(xbubbledata)
ggplot(xbubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,
                              Chordata,Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,
                              Kinorhyncha,Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                              Rhodophyta))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethya meloni","Tethia citrina"))+
  theme(axis.text.x = element_text(face = "italic"),legend.margin = margin(0,0,0,-2))+
  guides(color = guide_legend(override.aes = list(size = 5),ncol=2 ),size=guide_legend(nrow=1) )

ggsave("phylum_richness_and_abundance_bubble_wo_unclass.png",path=plotsPath,dpi=300,units="cm",width=30,height=20)

############################
###For core Otus
#Only core OTUs
head(coredataset)
head(dataset)
coreotucount<-coredataset%>%
  filter(value>0)%>%
  group_by(Phylum,variable)%>%
  dplyr::summarize(Numb_otu=n())

dim(otucount)
corerelabunds<-coredataset%>%
  filter(value>0)%>%
  group_by(Phylum,variable)%>%
  dplyr::summarise(abundance=sum(relative_abundance))

dim(relabunds)
corebubbledata<-cbind(coreotucount,corerelabunds[3])
#SHowing between phylums
corebubbledata[corebubbledata==0]<-NA
##plot
head(corebubbledata)
unique(corebubbledata$Phylum)
png("./Figures/18splots/18bubble_phylum_richness_and_abundanceCORE.png",height=20,width=30,units="cm",res=300)
ggplot(corebubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(Porifera,Unclassified))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethya meloni","Tethia citrina"))+
  theme(axis.text.x = element_text(face = "italic"),legend.margin = margin(0,0,0,-2))+
  guides(color = guide_legend(override.aes = list(size = 5),ncol=1 ),size=guide_legend(nrow=1) )

dev.off()
############3
#Without unclassified
xcorebubbledata<-corebubbledata[(!corebubbledata$Phylum=="Unclassified"),]
png("./Figures/18splots/18bubble_phylum_richness_and_abundanceCORE_wounclass.png",height=20,width=30,units="cm",res=300)
ggplot(xcorebubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(Porifera))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethya meloni","Tethia citrina"))+
  theme(axis.text.x = element_text(face = "italic"),legend.margin = margin(0,0,0,-2))+
  guides(color = guide_legend(override.aes = list(size = 5),ncol=1 ),size=guide_legend(nrow=1) )

dev.off()
#Dateset without unclassified and porifera
ycorebubbledata<-bubbledata[bubbledata$Phylum!="Unclassified"&bubbledata$Phylum!="Porifera",]

###############
library(cowplot)
library(gridExtra)
legend18<-get_legend(ggplot(bubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,
                              Chordata,Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,
                              Kinorhyncha,Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                              Rhodophyta,Unclassified))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethya meloni","Tethia citrina"))+
  theme(axis.text.x = element_text(face = "italic"),legend.margin = margin(0,0,0,-2))+
  guides(color = guide_legend(override.aes = list(size = 5),ncol=2),size=guide_legend(nrow=1) ))

fullbubble18<-ggplot(bubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,
                              Chordata,Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,
                              Kinorhyncha,Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                              Rhodophyta,Unclassified))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethya meloni","Tethia citrina"))+ theme(axis.text.x = element_text(face = "italic"),legend.position = ("none"))+
  ggtitle("Including all OTUs")

exclbubble18<-ggplot(ycorebubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,
                              Chordata,Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,
                              Kinorhyncha,Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Rhodophyta))+
  labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethya meloni","Tethia citrina"))+ theme(axis.text.x = element_text(face = "italic"),legend.position = "none")+
  ggtitle("Excluding porifera and unclassified OTUs")+
  guides(color = guide_legend(override.aes = list(size = 5) ) )

png("Figures/18splots/bubbleplots18s_2.png",width=14,height=7,units="in",res=300)
plot_grid(fullbubble18,exclbubble18,legend18,nrow=1,rel_widths=c(7/20,7/20,6/20))
dev.off()

#WIthout unclassified & porifera 
bubble2<-bubbledata[(bubbledata$Phylum!="Unclassified"&bubbledata$Phylum!="Porifera"),]
#bubble2<-bubble2[(bubble2$Phylum!="Ciliophora"),]
library(tidyverse)
ggplot(bubble2,aes(x=variable,y=Numb_otu,colour=Phylum,size=abundance))+
  scale_color_manual(values=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,
                              Chordata, Ciliophora,Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,
                              Kinorhyncha,Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,
                              Rhodophyta))+
  labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+
  scale_x_discrete(labels=c("Tethya aurantium","Tethya meloni","Tethia citrina"))+ 
  theme(axis.text.x = element_text(face = "italic"),legend.margin=margin(0,0,0,-2))+
  guides(color = guide_legend(override.aes = list(size = 5),ncol=2),size=guide_legend(nrow=1) )

ggsave("Figures/18splots/phylum_rich_and_abund_bubble_withoutunclass-porifera.png",dpi=300,units="cm",width=30,height=20)

#WIthout unclassified graphs together
###Grid arrangement w/o unclass
legend18<-get_legend(ggplot(xbubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
                       scale_color_manual(values=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,
                                                   Chordata,Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,
                                                   Kinorhyncha,Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                                                   Rhodophyta))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
                       scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethya meloni","Tethia citrina"))+
                       theme(axis.text.x = element_text(face = "italic"),legend.margin = margin(0,0,0,-5))+
                       guides(color = guide_legend(override.aes = list(size = 5),ncol=1),size=guide_legend(nrow=1) ))

fullbubble18<-ggplot(xbubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,
                              Chordata,Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,
                              Kinorhyncha,Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
                              Rhodophyta))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethya meloni","Tethia citrina"))+ theme(axis.text.x = element_text(face = "italic"),legend.position = ("none"))+
  ggtitle("Including all OTUs")

corebubble18<-ggplot(xcorebubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(Porifera))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethya meloni","Tethia citrina"))+ theme(axis.text.x = element_text(face = "italic"),legend.position = "none")+
  ggtitle("Including only Core OTUs")+
  guides(color = guide_legend(override.aes = list(size = 5) ) )

png("Figures/18splots/bubbleplots18s_wo_unclass.png",width=14,height=7,units="in",res=300)
plot_grid(fullbubble18,corebubble18,legend18,nrow=1,rel_widths=c(2/5,2/5,1/5))
dev.off()

#geom_text(aes(label=ifelse(as.numeric(variable)>20,as.character(Phylum),'')),vjust=0.3,hjust=0.3,position = position_stack(vjust = 0.3))
  
#TO LOOK AT RICHNESS AS THE SIZE Instead and ABUNDANCE on y axis
#ggplot(bubble2,aes(x=variable,y=abundance,color=Phylum,size=Numb_otu))+
#  scale_color_manual(values=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,
 #                             Chordata,Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,
  #                            Kinorhyncha,Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
   #                           Rhodophyta,Unclassified))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  #scale_size(range = c(.7,20))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)

