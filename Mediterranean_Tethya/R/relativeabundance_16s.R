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
CalditrichotaCol<-"maroon"
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

#BUbble plots
#BUbble plot for comparing, richness and abundance of phylums across species
head(dataset)
unique(dataset$Phylum)
otucount<-dataset%>%
  group_by(Phylum,variable)%>%
  summarise(Numb_otu=n())
dim(otucount)
relabunds<-dataset%>%
  group_by(Phylum,variable)%>%
  summarise(abundance=sum(relative_abundance))
dim(relabunds)
bubbledata<-cbind(otucount,relabunds[3])
head(bubbledata)
#SHowing between phylums
bubbledata[bubbledata==0]<-NA
#Showing between species
ggplot(bubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,BdellovibrionotaCol,
                              CalditrichotaCol,CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
                              DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DependentiaeCol,DesulfobacterotaCol,
                              EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,GemmatimonadotaCol,HydrogenedentesCol,
                              LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,PAUC34fCol,
                              PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,SpirochaetotaCol,SumerlaeotaCol,
                              UnclassifiedCol,VerrucomicrobiotaCol,WPS2Col,WS1Col))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.7,20))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethia citrina","Tethya meloni"))+ theme(axis.text.x = element_text(face = "italic"))
ggsave("phylum_richness_and_abundance_bubble.jpg",path=plotsPath)
#WIthout phylums of high abundance (more than 0.1 relative abundance)and up to 3000 richness


bubble2<-bubbledata%>%
  filter(abundance<=0.1)
bubble2<-bubble2%>%
  filter(Numb_otu<=3500)

unique(bubble2$Phylum)
ggplot(bubble2,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,BdellovibrionotaCol,
                              CalditrichotaCol,CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
                              DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DependentiaeCol,DesulfobacterotaCol,
                              EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,GemmatimonadotaCol,HydrogenedentesCol,
                              LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,PAUC34fCol,
                              PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,SpirochaetotaCol,SumerlaeotaCol,
                              UnclassifiedCol,VerrucomicrobiotaCol,WPS2Col,WS1Col))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.7,20))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethia citrina","Tethya meloni"))+ theme(axis.text.x = element_text(face = "italic"))
ggsave("phylum_rich_and_abund_bubble_withoutunclass-cill-porifera.jpg",path=plotsPath)
#TO LOOK AT RICHNESS AS THE SIZE Instead and ABUNDANCE on y axis
#ggplot(bubble2,aes(x=variable,y=abundance,color=Phylum,size=Numb_otu))+
#  scale_color_manual(values=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,
#                             Chordata,Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,
#                            Kinorhyncha,Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
#                           Rhodophyta,Unclassified))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
#scale_size(range = c(.7,20))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)



