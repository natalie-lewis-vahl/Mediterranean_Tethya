#continuation from "AS_NLcoreCommunity.R" script
#BUbble plots
#BUbble plot for comparing, richness and abundance of phylums across species
head(dataset)
str(dataset)
is.factor(dataset$Phylum)
unique(dataset$Phylum)

otucount<-dataset%>%
  filter(value>0)%>%
  group_by(Phylum,variable)%>%
  dplyr::summarize(Numb_otu=n())

dim(otucount)
relabunds<-dataset%>%
  filter(value>0)%>%
  group_by(Phylum,variable)%>%
  dplyr::summarise(abundance=sum(relative_abundance))

dim(relabunds)
bubbledata<-cbind(otucount,relabunds[3])
head(bubbledata)
#SHowing between phylums
bubbledata[bubbledata==0]<-NA
#Showing between species
png("./Figures/AS_16splots/bubble_phylum_richness_and_abundance.png",height=20,width=30,units="cm",res=300)
ggplot(bubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,
                              CalditrichotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
                              DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DesulfobacterotaCol,
                              EntotheonellaeotaCol,GemmatimonadotaCol,
                              LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,
                              PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,SpirochaetotaCol,
                              UnclassifiedCol,VerrucomicrobiotaCol))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethia citrina","Tethya meloni"))+ theme(axis.text.x = element_text(face = "italic"))
dev.off()

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
png("./Figures/AS_16splots/bubble_phylum_richness_and_abundanceCORE.png",height=20,width=30,units="cm",res=300)
ggplot(corebubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,
                              CrenarchaeotaCol,CyanobacteriaCol,
                              DadabacteriaCol,DesulfobacterotaCol,
                              NitrospirotaCol,
                              PlanctomycetotaCol,ProteobacteriaCol,
                              UnclassifiedCol,VerrucomicrobiotaCol))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethia citrina","Tethya meloni"))+ theme(axis.text.x = element_text(face = "italic"))
dev.off()

library(cowplot)
library(gridExtra)
fullbubble<-ggplot(bubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,
                              CalditrichotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
                              DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DesulfobacterotaCol,
                              EntotheonellaeotaCol,GemmatimonadotaCol,
                              LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,
                              PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,SpirochaetotaCol,
                              UnclassifiedCol,VerrucomicrobiotaCol))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethia citrina","Tethya meloni"))+ theme(axis.text.x = element_text(face = "italic"))
legend<-get_legend(fullbubble)
fullbubble2<-ggplot(bubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,
                              CalditrichotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
                              DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DesulfobacterotaCol,
                              EntotheonellaeotaCol,GemmatimonadotaCol,
                              LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,
                              PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,SpirochaetotaCol,
                              UnclassifiedCol,VerrucomicrobiotaCol))+
  labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+
  scale_x_discrete(labels=c("Tethya aurantium","Tethia citrina","Tethya meloni"))+ 
  theme(axis.text.x = element_text(face = "italic"),legend.position = ("none"))+
  ggtitle("Including all OTUs")

corebubble2<-ggplot(corebubbledata,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
  scale_color_manual(values=c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,
                              CrenarchaeotaCol,CyanobacteriaCol,
                              DadabacteriaCol,DesulfobacterotaCol,
                              NitrospirotaCol,
                              PlanctomycetotaCol,ProteobacteriaCol,
                              UnclassifiedCol,VerrucomicrobiotaCol))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
  scale_size(range = c(.9,15))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethia citrina","Tethya meloni"))+ theme(axis.text.x = element_text(face = "italic"),legend.position = "none")+
  ggtitle("Including only Core OTUs")
png("Figures/AS_16splots/bubbleplots16s.png",width=11,height=5,units="in",res=300)
grid.arrange(fullbubble2,corebubble2,legend,nrow=1)
dev.off()
??png()
#WIthout phylums of high abundance (more than 0.1 relative abundance)and up to 3000 richness

#bubble2<-bubbledata%>%
#  filter(abundance<=0.1)
#bubble2<-bubble2%>%
#  filter(Numb_otu<=3500)

#unique(bubble2$Phylum)
#ggplot(bubble2,aes(x=variable,y=Numb_otu,color=Phylum,size=abundance))+
#  scale_color_manual(values=c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,BdellovibrionotaCol,
#                              CalditrichotaCol,CampylobacterotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
#                              DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DependentiaeCol,DesulfobacterotaCol,
#                              EntotheonellaeotaCol,FirmicutesCol,FusobacteriotaCol,GemmatimonadotaCol,HydrogenedentesCol,
#                              LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,PAUC34fCol,
#                              PlanctomycetotaCol,SAR324cladeCol,SpirochaetotaCol,SumerlaeotaCol,
#                              VerrucomicrobiotaCol,WPS2Col,WS1Col))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
#  scale_size(range = c(.7,20))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethia citrina","Tethya meloni"))+ theme(axis.text.x = element_text(face = "italic"))
#ggsave("phylum_rich_and_abund_bubble_withoutunclass-cill-porifera.jpg",path=plotsPath)

#TO LOOK AT RICHNESS AS THE SIZE Instead and ABUNDANCE on y axis
#ggplot(bubble2,aes(x=variable,y=abundance,color=Phylum,size=Numb_otu))+
#  scale_color_manual(values=c(Annelida,Arthropoda,Bacillariophyta,Brachiopoda,Bryoza,Chlorophyta,
#                             Chordata,Ciliophora, Cnidaria,Dinoflagellata,Echinodermata,Euglenoza,
#                            Kinorhyncha,Magnoliophyta,Mollusca,Nematoda,Platyhelminthes,Porifera,
#                           Rhodophyta,Unclassified))+labs(x="Species",y="Phylum Richness",size="Relative abundance")+
#scale_size(range = c(.7,20))+theme_bw()+geom_jitter(width=0.3,alpha=0.9)

#Visualize Otus seperately- compare abundance across phylas per sp
head(dataset)
dataset[dataset==0]<-NA
dataset<-dataset%>%
  filter(relative_abundance>0)
###
png("./Figures/AS_16splots/bubble_otu_abundance.png",height=30,width=30,units="cm",res=300)

ggplot(dataset,aes(x=variable,y=relative_abundance,color=Phylum))+
  scale_color_manual(values=c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,
                              CalditrichotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
                              DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DesulfobacterotaCol,
                              EntotheonellaeotaCol,GemmatimonadotaCol,
                              LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,
                              PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,SpirochaetotaCol,
                              UnclassifiedCol,VerrucomicrobiotaCol))+labs(x="Species",y="Relative abundance")+
  theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethia citrina","Tethya meloni"))+ theme(axis.text.x = element_text(face = "italic"))
dev.off()

#Abundance divided by richness
png("./Figures/AS_16splots/bubble_phyla_abundance_dividedby_richness.png",height=20,width=30,units="cm",res=300)

ggplot(bubbledata,aes(x=variable,y=abundance/Numb_otu,color=Phylum))+
  scale_color_manual(values=c(AcidobacteriotaCol,ActinobacteriotaCol,BacteroidotaCol,
                              CalditrichotaCol,ChloroflexiCol,CrenarchaeotaCol,CyanobacteriaCol,
                              DadabacteriaCol,DeferrisomatotaCol,DeinococcotaCol,DesulfobacterotaCol,
                              EntotheonellaeotaCol,GemmatimonadotaCol,
                              LatescibacterotaCol,MyxococcotaCol,NB1jCol,NitrospinotaCol,NitrospirotaCol,
                              PlanctomycetotaCol,ProteobacteriaCol,SAR324cladeCol,SpirochaetotaCol,
                              UnclassifiedCol,VerrucomicrobiotaCol))+labs(x="Species",y="Relative abundance/Phylum richness")+
  theme_bw()+geom_jitter(width=0.3,alpha=0.9)+scale_x_discrete(labels=c("Tethya aurantium","Tethia citrina","Tethya meloni"))+ theme(axis.text.x = element_text(face = "italic"))
dev.off()
