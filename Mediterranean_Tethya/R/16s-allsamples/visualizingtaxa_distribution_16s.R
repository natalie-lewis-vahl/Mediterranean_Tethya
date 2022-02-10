library(dplyr)
library(ggplot2)
taxa<-read.table("./Data/16s_allsamples_taxa.csv",header=TRUE,sep=";")
otus<-read.table("./Data/16s_allsamples_otu.csv",sep=";")
names(taxa)
names(otus)
reorder_id<-match(taxa$sequence_identifier,otus$V1)
otus<-otus[reorder_id,]

merged<-merge(taxa,otus,by.x="sequence_identifier",by.y="V1")
str(merged)
dim(merged)

xmerged<-(merged[-c(1:8)])
xmerged<-lapply(xmerged,as.numeric)
head(xmerged)
count<-apply(as.data.frame(xmerged),1,sum)
??apply
merged<-cbind(merged,count)
head(merged)
#Remove OTUS with less than 50 reads
filtered<-merged[merged$count>=50,]
filtered[filtered==""]<-"Unclassified"

unique(filtered$Phylum)
??count
filtered$Phylum<-str_replace_all(filtered$Phylum,"SAR324 clade(Marine group B)","SAR324 clade")
filtered[filtered=="SAR324 clade(Marine group B)"]<-"SAR324 clade"
unique(filtered$Phylum)

plot1data<-filtered%>%
  dplyr::count(Phylum)

data<-plot1data$n
names(data)<-plot1data$Phylum
barplot(data)
View(data)
#The large amount of otus belonging to unclassified makes other phyla harder
#to differentiate - could represent same phyla or lots of different ones
#remove 9215 unknown OTUs optional
#plot1data<-plot1data[!(plot1data$Phylum=="Unclassified"),]
data2<-plot1data$n
names(data2)<-plot1data$Phylum
barplot(data2)
#delete row with unclassified phyla
#filtered<-filtered[!(filtered$Phylum=="Unclassified"),]
counts<-table(filtered$Class,filtered$Phylum)
head(counts)

jpeg("./Figures/AS_16splots/overall16s_taxa_distribution.jpeg",units="in", width=35, height=10, res=300)
barplot(counts,xlab="Phyla",ylab="Number of OTUs",cex.names= 0.5,ylim=c(0,160))
dev.off()
??barplot

#using ggplot
#Arrange table for this and order to have nice plot
datax<-plot1data
datax<-datax %>%
  arrange(desc(n)) %>%    
  mutate(Phylum=factor(Phylum, levels=Phylum))
View(datax)
#colours
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

tiff("./Figures/AS_16splots/overall16s_taxa_distribution2.tiff",units="in", width=28, height=11, res=300)
ggplot(datax,aes(x=Phylum,y=n,fill=Phylum))+
  geom_bar(stat="identity")+scale_y_continuous(limits=c(0,175),expand=c(0,0),minor_breaks = seq(0 , 175, 5), breaks = seq(0, 175, 50))+
  scale_x_discrete(expand=c(0,0))+
  scale_fill_manual(values = c(ProteobacteriaCol,UnclassifiedCol,PlanctomycetotaCol,BacteroidotaCol,
                               DesulfobacterotaCol,ActinobacteriotaCol,AcidobacteriotaCol,VerrucomicrobiotaCol,
                               MyxococcotaCol,NitrospirotaCol,CyanobacteriaCol,DeferrisomatotaCol,NB1jCol,        
                               ChloroflexiCol,CrenarchaeotaCol,DadabacteriaCol,EntotheonellaeotaCol,
                               GemmatimonadotaCol,CalditrichotaCol,DeinococcotaCol,LatescibacterotaCol,
                               NitrospinotaCol,SAR324cladeCol,SpirochaetotaCol)) +
  labs(x="Phylum",y="Number of OTUs")+ggtitle("All 3 tethya species combined")+
  theme_bw() + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=10),legend.position="none",plot.title=element_text(hjust=0.5))+
  theme(text=element_text(size=14),panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey30" ) )

dev.off()
#species seperately
#Tau V2 till v12+v45 (col 9-19,52)
dim(filtered)
filteredtau<-filtered[,-c(20:51)]
dim(filteredtau)
taux<-filteredtau[!apply(filteredtau[,9:20]==0,1,all),]
dim(taux)
#Tme v13 till v23 +v43tillv44 (col 20-30,50,51)
filteredtme<-filtered[,-c(9:19,31:49,52)]
head(filteredtme)
tmex<-filteredtme[!apply(filteredtme[,9:21]==0,1,all),]

#Tci v35 till v42 (col 42-49)
filteredtci<-filtered[-c(9:30,50:52)]
tcix<-filteredtci[!apply(filteredtci[,9:27]==0,1,all),]

#scale_x_discrete(labels = expression(italic(T.aurantium),italic(T.meloni),italic(T.citrina)))+

plot2data<-taux%>%
  dplyr::count(Phylum)
datatau<-plot2data %>%
  arrange(desc(n)) %>%    
  mutate(Phylum=factor(Phylum, levels=Phylum))
##########
tiff("./Figures/AS_16splots/tau16s_taxa_distribution2.tiff",units="in", width=10, height=8, res=300)
ggplot(datatau,aes(x=Phylum,y=n,fill=Phylum))+
  geom_bar(stat="identity")+scale_y_continuous(limits=c(0,150),expand=c(0,0),minor_breaks = seq(0 , 175, 5), breaks = seq(0, 175, 50))+
  scale_x_discrete(expand=c(0,0))+
  scale_fill_manual(values = c(ProteobacteriaCol,UnclassifiedCol,PlanctomycetotaCol,BacteroidotaCol,
                               ActinobacteriotaCol,DesulfobacterotaCol,VerrucomicrobiotaCol,AcidobacteriotaCol,
                               NitrospirotaCol,CyanobacteriaCol,MyxococcotaCol,NB1jCol,        
                               CrenarchaeotaCol,ChloroflexiCol,DadabacteriaCol,EntotheonellaeotaCol,
                               GemmatimonadotaCol,DeferrisomatotaCol,DeinococcotaCol,
                               NitrospinotaCol,SAR324cladeCol,SpirochaetotaCol)) +
  labs(x="Phylum",y="Number of OTUs")+
  ggtitle(expression(italic("Tethya aurantium")))+
  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=0.8,size=9),legend.position="none")+
  theme(text=element_text(size=10),axis.title.x=element_text(vjust=3),panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey30" ))
dev.off()

####Tme
plot3data<-tmex%>%
  dplyr::count(Phylum)
datatme<-plot3data %>%
  arrange(desc(n)) %>%    
  mutate(Phylum=factor(Phylum, levels=Phylum))
##########
tiff("./Figures/AS_16splots/tme16s_taxa_distribution2.tiff",units="in", width=10, height=8, res=300)
ggplot(datatme,aes(x=Phylum,y=n,fill=Phylum))+
  geom_bar(stat="identity")+scale_y_continuous(limits=c(0,150),expand=c(0,0),minor_breaks = seq(0 , 175, 5), breaks = seq(0, 175, 50))+
  scale_fill_manual(values = c(ProteobacteriaCol,UnclassifiedCol,PlanctomycetotaCol,BacteroidotaCol,
                               ActinobacteriotaCol,VerrucomicrobiotaCol,AcidobacteriotaCol,DesulfobacterotaCol,
                               NitrospirotaCol,CyanobacteriaCol,MyxococcotaCol,CrenarchaeotaCol,NB1jCol,        
                               DadabacteriaCol,EntotheonellaeotaCol,GemmatimonadotaCol,LatescibacterotaCol,
                               NitrospinotaCol,SAR324cladeCol,SpirochaetotaCol)) +
  labs(x="Phylum",y="Number of OTUs")+
  ggtitle(expression(italic("Tethya meloni")))+
  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=0.8,size=9),legend.position="none")+
  theme(text=element_text(size=10),axis.title.x=element_text(vjust=3),panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey30" ))

dev.off()

####Tci
plot4data<-tcix%>%
  dplyr::count(Phylum)
datatci<-plot4data %>%
  arrange(desc(n)) %>%    
  mutate(Phylum=factor(Phylum, levels=Phylum))

##########
tiff("./Figures/AS_16splots/tci16s_taxa_distribution2.tiff",units="in", width=10, height=8, res=300)
ggplot(datatci,aes(x=Phylum,y=n,fill=Phylum))+
  geom_bar(stat="identity")+scale_y_continuous(limits=c(0,150),expand=c(0,0),minor_breaks = seq(0 , 175, 5), breaks = seq(0, 175, 50))+
  scale_fill_manual(values = c(ProteobacteriaCol,UnclassifiedCol,PlanctomycetotaCol,BacteroidotaCol,
                               DesulfobacterotaCol,ActinobacteriotaCol,AcidobacteriotaCol,VerrucomicrobiotaCol,
                               MyxococcotaCol,NitrospirotaCol,CyanobacteriaCol,DeferrisomatotaCol,NB1jCol,        
                               ChloroflexiCol,CrenarchaeotaCol,DadabacteriaCol,EntotheonellaeotaCol,
                               GemmatimonadotaCol,CalditrichotaCol,DeinococcotaCol,LatescibacterotaCol,
                               NitrospinotaCol,SAR324cladeCol)) +
  labs(x="Phylum",y="Number of OTUs")+
  ggtitle(expression(italic("Tethya citrina")))+
  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=0.8,size=9),legend.position="none")+
  theme(text=element_text(size=10),axis.title.x=element_text(vjust=3),panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey30" ))


dev.off()

