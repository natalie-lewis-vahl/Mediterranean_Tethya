library(propr)
library(dplyr)
library(tibble)
library(ggdendro)
library(corrgram)
library(ggplot2)
library(bootnet)
#just 16s
otuTaxonomyPath<-"./Data/16s_allsamples_taxa.csv"
#Need to remove sample GW1956 because of low sampling -> the 16th sample
#belonging to Tme
otutable<-read.table("./Data/16s_allsamples_otu.txt", header=TRUE)
taxatable<-read.csv(otuTaxonomyPath, sep=";")
ls()
#reorder otu to match taxa
#get OTU counts and add them to DF
countSum<-apply(otutable[-1],1,sum)
countsDF<-cbind(otutable, countSum)

#bind the taxonomy columns to the otu counts DF
countsWithTaxonomy<-bind_cols(countsDF[order(countsDF$X.OTU.ID),], taxatable[order(taxatable$sequence_identifier),])
head(countsWithTaxonomy)
#Filter from 50 counts up
countsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$countSum >=50,]

#Make a separate otu table without taxa
otu_table=select(countsWithTaxonomy,-c(46:54))
otut<-otu_table[,-1]
#OTU Id as row name
rownames(otut)<-otu_table[,1]
#Samples as row
otu_table<-as.data.frame(t(otut))

#calculate the total number of reads per sample
#remove bad sample
which(rownames(otu_table)=="GW1956")
dim(otu_table)
otu_table<-otu_table[-16,]
otu_table$total <- rowSums(otu_table)
otu_table$sponge <- c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")
#GET NAMES OF ALL OUTS
otus <- names(otu_table)[-c(412,413)]
#PUT SAMPLES BACK AS A ROW
otu_table<-tibble::rownames_to_column(otu_table, "Sample")
xmatrix <- otu_table[,grep('OTU', names(otu_table))] 
fullmatrix<-xmatrix
#For all species
#Need to filter out OTUs not part of the network
#how many samples is each otu in
xmatrix[xmatrix>=1]<-1
otu_samples <- colSums(xmatrix)
active_otus <- names(otu_samples[otu_samples >38])
counts<- otu_table[ ,which((names(otu_table) %in% active_otus)==TRUE)]
species<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")
sample<-c("GW1941","GW1942" ,"GW1943","GW1944", "GW1945", "GW1946", "GW1947", "GW1948","GW1949","GW1950","GW1951",
          "GW1952","GW1953", "GW1954", "GW1955", "GW1957", "GW1958", "GW1959","GW1960","GW1961","GW1962",
          "GW1963","GW1964","GW1965", "GW1966", "GW1967", "GW1968", "GW1969","GW1970","GW1971","GW1972", "GW1973","GW1974", "GW1975", "GW1976", "GW1977", "GW1978", "GW1979","GW1980","GW1981", "GW1982","GW1983", "GW1984")
#Note that the log-ratio transformation, by its nature, fails if the input data contain any zero values. 
#By default, this function replaces all zero values with 1. Alternatively, the user may set the parameter `alpha` greater than zero to approximate log-ratios in the presence of zeros (via the Box-Cox transformation). However, the topic of zero replacement is controversial. Proceed carefully when analyzing data that contain zero values.
rho <- propd(counts, species, metric="rho",ivar="clr",p=100,alpha=NA)
rho <- propr(counts, metric = "rho","clr", alpha=NA,p=100)
best <- rho[">", .50]

getAdj(rho)
rhomatrixall<-getMatrix(rho)
head(rhomatrixall)
tiff("Figures/Networks/allspecies16s_rho.tiff",res=300,units="cm",width=20,height=15)
corrgram(rhomatrixall,order=TRUE,main="Rho proportionality between OTUs present across all 16s samples ",
         lower.panel=panel.shade, upper.panel=panel.cor, text.panel=panel.txt)
dev.off()
 
??corrgram

########################
###
##Tethya citrina
head(xmatrix)
dim(fullmatrix)
tcimatrix_abund<-fullmatrix[22:40,]
tcimatrix<-xmatrix[22:40,]
otu_samplestci <- colSums(tcimatrix)
#present in more than half-(=N.samples/2 +1)
active_otustci <- names(otu_samplestci[otu_samplestci >17])
counts_tci<- tcimatrix_abund[ ,which((names(tcimatrix) %in% active_otustci)==TRUE)]
###
#Note that the log-ratio transformation, by its nature, fails if the input data contain any zero values. 
#By default, this function replaces all zero values with 1. Alternatively, the user may set the parameter `alpha` greater than zero to approximate log-ratios in the presence of zeros (via the Box-Cox transformation). However, the topic of zero replacement is controversial. Proceed carefully when analyzing data that contain zero values.
rho <- propr(counts_tci, metric = "rho","clr", alpha=NA,p=100)
best <- rho[">", .70]
plot(best)
mat<-getMatrix(best)
pca(best)
rhomatrixbest<-getMatrix(best)
best<-propr::simplify(best)

#
getAdj(rho)
rhomatrix<-getMatrix(rho)
dim(rhomatrix)
max(rhomatrix)
min(rhomatrix)

tiff("Figures/Networks/Tcitrina16s_rho.tiff",res=300,units="cm",width=45,height=30)
corrgram(rhomatrix,order=TRUE,main="Rho proportionality between OTUs present in T citrina 16s samples ",
         lower.panel=panel.shade, upper.panel=panel.cor, text.panel=panel.txt)
dev.off()
########################
##Tethy aurantium
dim(xmatrix)
taumatrix<-xmatrix[c(1:11,43),]
taumatrix_abund<-fullmatrix[c(1:11,43),]

otu_samplestau <- colSums(taumatrix)
#present in more than half-(=N.samples/2 +1)
active_otustau <- names(otu_samplestau[otu_samplestau > 10])
filteredmat_tau<- taumatrix_abund[ ,which((names(taumatrix) %in% active_otustau)==TRUE)]
dim(filteredmat_tau) 
##
rho <- propr(filteredmat_tau, metric = "rho","clr", alpha=NA,p=100)
best <- rho[">", .60]
plot(best)
#
rhomatrix<-getMatrix(rho)

tiff("Figures/Networks/Taurantium16s_rho.tiff",res=300,units="cm",width=20,height=15)
corrgram(rhomatrix,order=TRUE,main="Rho proportionality between OTUs present in T aurantium 16s samples ",
         lower.panel=panel.shade, upper.panel=panel.cor, text.panel=panel.txt)
dev.off()
#####
##Tethya meloni
head(xmatrix)
tmematrix<-xmatrix[c(12:21,42,43),]
tmematrix_abund<-fullmatrix[c(12:21,42,43),]

otu_samplestme <- colSums(tmematrix)
#present in more than half-(=N.samples/2 +1)
active_otustme <- names(otu_samplestme[otu_samplestme > 10])
filteredmat_tme<- tmematrix_abund[ ,which((names(tmematrix) %in% active_otustme)==TRUE)]
dim(filteredmat_tme) 

##
rho <- propr(filteredmat_tme, metric = "rho","clr", alpha=NA,p=100)
rhomatrix<-getMatrix(rho)

tiff("Figures/Networks/Tmeloni16s_rho.tiff",res=300,units="cm",width=20,height=15)
corrgram(rhomatrix,order=TRUE,main="Rho proportionality between OTUs present in T meloni 16s samples ",
         lower.panel=panel.shade, upper.panel=panel.cor, text.panel=panel.txt)
dev.off()
#no results
##########################

##############
###############
#Network for 18s data

CountsPath18<-"./Data/18s/all.otutab.csv"
otuTaxonomyPath18<-"./Data/18s/taxa_fixed18s.csv"
#Need to remove sample GW1956 because of low sampling -> the 16th sample
#belonging to Tme
otutable<-read.csv(CountsPath18, sep=";")
taxatable<-read.csv(otuTaxonomyPath18, sep=";")
#reorder otu to match taxa
#get OTU counts and add them to DF
countSum<-apply(otutable[-1],1,sum)
countsDF<-cbind(otutable, countSum)

#bind the taxonomy columns to the otu counts DF
countsWithTaxonomy<-bind_cols(countsDF[order(countsDF$X.OTU.ID),], taxatable[order(taxatable$sequence_identifier),])
head(countsWithTaxonomy)
countsWithTaxonomy<-countsWithTaxonomy%>%
  filter(Domain!="Bacteria")%>%
  filter(Class!="Demospongiae")

#Filter from 50 counts up
countsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$countSum >=5,]

#Make a separate otu table without taxa
otu_table=select(countsWithTaxonomy,-c(46:54))
otut<-otu_table[,-1]
#OTU Id as row name
rownames(otut)<-otu_table[,1]
#Samples as row
otu_table<-as.data.frame(t(otut))

#calculate the total number of reads per sample
#remove bad sample which is different for this dataset
#GW1968 belonging to T citrina
head(otu_table)
which(colnames(otu_table)=="OTU_4076")
which(colnames(otu_table)=="OTU_4077")
which(rownames(otu_table)=="GW1968")
dim(otu_table)
otu_table<-otu_table[-28,-c(1,2)]
dim(otu_table)
otu_table$total <- rowSums(otu_table)
otu_table$sponge <- c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")
#GET NAMES OF ALL OUTS
dim(otu_table)
otus <- names(otu_table)[-c(296,297)]
#PUT SAMPLES BACK AS A ROW
otu_table<-tibble::rownames_to_column(otu_table, "Sample")
xmatrix <- otu_table[,grep('OTU', names(otu_table))] 
fullmatrix<-xmatrix
#For all species
#Need to filter out OTUs not part of the network
#how many samples is each otu in
xmatrix[xmatrix>=1]<-1
otu_samples <- colSums(xmatrix)
active_otus <- names(otu_samples[otu_samples >38])
#none

########################
###Already know T citrina has no core OTUs in 18s so skip

##Tethya aurantium
dim(xmatrix)
taumatrix<-xmatrix[c(1:11,43),]
taumatrix_abund<-fullmatrix[c(1:11,43),]

otu_samplestau <- colSums(taumatrix)
#present in more than half-(=N.samples/2 +1)
active_otustau <- names(otu_samplestau[otu_samplestau > 10])
filteredmat_tau<- taumatrix_abund[ ,which((names(taumatrix) %in% active_otustau)==TRUE)]
dim(filteredmat_tau) 
##
rho <- propr(filteredmat_tau, metric = "rho","clr", alpha=NA,p=100)
best <- rho[">", .60]
plot(best)
#
rhomatrix<-getMatrix(rho)
head(rhomatrix)
tiff("Figures/Networks/Taurantium18s_rho.tiff",res=300,units="cm",width=20,height=15)
corrgram(rhomatrix,order=TRUE,main="Rho proportionality between OTUs present in T aurantium 18s samples ",
         lower.panel=panel.shade, upper.panel=panel.cor, text.panel=panel.txtp.mat = rhomatrix$p, insig = 'p-value')
dev.off()
#####
##Tethya meloni
head(xmatrix)
tmematrix<-xmatrix[c(12:22,42,43),]
tmematrix_abund<-fullmatrix[c(12:22,42,43),]

otu_samplestme <- colSums(tmematrix)
#present in more than half-(=N.samples/2 +1)
active_otustme <- names(otu_samplestme[otu_samplestme > 11])
filteredmat_tme<- tmematrix_abund[ ,which((names(tmematrix) %in% active_otustme)==TRUE)]
dim(filteredmat_tme) 

##
rho <- propr(filteredmat_tme, metric = "rho","clr", alpha=NA,p=100)
rhomatrix<-getMatrix(rho)

tiff("Figures/Networks/Tmeloni18s_rho.tiff",res=300,units="cm",width=20,height=15)
corrgram(rhomatrix,order=TRUE,main="Rho proportionality between OTUs present in T meloni 18s samples ",
         lower.panel=panel.shade, upper.panel=panel.cor, text.panel=panel.txt)
dev.off()

#Repeat
#16 and 18
taotus<-read.csv("./Data/16and18s_otu.csv", sep=";")
otus<-arrange(taotus, X.OTU.ID)
#BInd taxonomy to remove bacterial OTUs from 18s dataset 
taxa<-read.csv("./Data/16ans18s_taxa.csv",sep=";")
head(otus)
taxa<-arrange(taxa,sequence_identifier)
otu_and_taxa<-bind_cols(otus,taxa)
otu_and_taxa<-otu_and_taxa%>%filter(!(Domain=="Bacteria" & dataset=="18s"))
otu_and_taxa<-otu_and_taxa%>%filter(!(Class=="Demospongiae" & dataset=="18s"))

#and remove two otus fromcoamplified sponges for 18s
#Take out OTU for sponge amplification
head(otu_and_taxa)
otu_and_taxa<-otu_and_taxa[-c(47)]
countSum<-apply(otu_and_taxa[3:46],1,sum)
otus_and_taxa<-cbind(otu_and_taxa,countSum)
coamplif_sponge<-otus_and_taxa[otus_and_taxa$dataset=="18s",]%>%
  arrange(desc(countSum))
#
#already removed comaplified sponge
#coamplif_sponge$sequence_identifier[1]#OTU_4076
#otus_and_taxa<-otus_and_taxa[!(otus_and_taxa$X.OTU.ID=="OTU_4076"),]
#coamplif_sponge$sequence_identifier[2]#OTU_4077
#otus_and_taxa<-otus_and_taxa[!(otus_and_taxa$X.OTU.ID=="OTU_4077"),]

#filtering: should be done seperate for both 16 and 18s data set
#Keep rows that have more than 50 counds and are from 16s data set OR have more than 5 counts and belong to 18s ddataset
otus_and_taxab<-otus_and_taxa[otus_and_taxa$countSum >= 50 & otus_and_taxa$dataset=="16s"|otus_and_taxa$countSum >= 5 & otus_and_taxa$dataset=="18s",]
head(otus_and_taxab)
data_all<-otus_and_taxab
which(colnames(data_all)=="GW1984")
data_all<-data_all[,2:46]

otut<-data_all[,-1]
#OTU Id as row name
rownames(otut)<-data_all[,1]
#Samples as row
otu_table<-as.data.frame(t(otut))
head(otu_table)
#calculate the total number of reads per sample
#remove both bad samples from 16s and 18s
#GW1968 and GW1956
head(otu_table)
which(rownames(otu_table)=="GW1956")
which(rownames(otu_table)=="GW1968")
dim(otu_table)
otu_table<-otu_table[-c(16,28),]
dim(otu_table)
otu_table$total <- rowSums(otu_table)
otu_table$sponge <- c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")
#GET NAMES OF ALL OUTS
otu_table
otus <- names(otu_table)
head(otu_table)
#PUT SAMPLES BACK AS A ROW
otu_table<-tibble::rownames_to_column(otu_table, "Sample")
xmatrix <- otu_table[,grep('OTU', names(otu_table))] 
fullmatrix<-xmatrix
#Need to filter out OTUs not part of the network
#how many samples is each otu in
xmatrix[xmatrix>=1]<-1
otu_samples <- colSums(xmatrix)
#For T aurantium
##Tethya aurantium
dim(xmatrix)
taumatrix<-xmatrix[c(1:11,42),]
taumatrix_abund<-fullmatrix[c(1:11,42),]

otu_samplestau <- colSums(taumatrix)
#present in more than 90%
active_otustau <- names(otu_samplestau[otu_samplestau > 10])
filteredmat_tau<- taumatrix_abund[ ,which((names(taumatrix) %in% active_otustau)==TRUE)]
dim(filteredmat_tau) 
##
rho <- propr(filteredmat_tau, metric = "rho","clr", alpha=NA,p=100)
best <- rho[">", .60]
plot(best)
#
rhomatrixtau<-getMatrix(rho)
png("Figures/Networks/Taurantium16sand18s_rho.png",res=300,units="cm",width=30,height=20)
corrgram(rhomatrixtau,order=TRUE,main="Rho proportionality between OTUs present in T aurantium 16s and 18s samples ",
         lower.panel=panel.shade, upper.panel=panel.cor, text.panel=panel.txt)
dev.off()
tiff("Figures/Networks/Taurantium16sand18s_rho.tiff",res=300,units="cm",width=30,height=20)
corrgram(rhomatrixtau,order=TRUE,main="Rho proportionality between OTUs present in T aurantium 16s and 18s samples ",
         lower.panel=panel.shade, upper.panel=panel.cor, text.panel=panel.txt)
dev.off()
#####
#For T meloni
dim(xmatrix)
tmematrix<-xmatrix[c(12:21,40,41),]
tmematrix_abund<-fullmatrix[c(12:21,40,41),]

otu_samplestme <- colSums(tmematrix)
#present in more than 90%
active_otustme <- names(otu_samplestme[otu_samplestme > 10])
filteredmat_tme<- tmematrix_abund[ ,which((names(tmematrix) %in% active_otustme)==TRUE)]
dim(filteredmat_tme) 
##
rho <- propr(filteredmat_tme, metric = "rho","clr", alpha=NA,p=100)
best <- rho[">", .60]
plot(best)
#
rhomatrixtme<-getMatrix(rho)
png("Figures/Networks/Tme16sand18s_rho.png",res=300,units="cm",width=20,height=20)
corrgram(rhomatrixtme,order=TRUE,main="Rho proportionality between OTUs present in T meloni 16s and 18s samples ",
         lower.panel=panel.shade, upper.panel=panel.cor, text.panel=panel.txt)
dev.off()
tiff("Figures/Networks/Tme16sand18s_rho.svg",res=300,units="cm",width=20,height=20)
corrgram(rhomatrixtme,order=TRUE,main="Rho proportionality between OTUs present in T meloni 16s and 18s samples ",
         lower.panel=panel.shade, upper.panel=panel.cor, text.panel=panel.txt)
dev.off()

##For T citrina
dim(xmatrix)
tcimatrix<-xmatrix[c(22:39),]
tcimatrix_abund<-fullmatrix[c(22:39),]

otu_samplestci <- colSums(tcimatrix)
#present in more than 90%
active_otustci <- names(otu_samplestci[otu_samplestci > 16])
filteredmat_tci<- tcimatrix_abund[ ,which((names(tcimatrix) %in% active_otustci)==TRUE)]
head(filteredmat_tci) 
dim(filteredmat_tci)
##
rhotci <- propr(filteredmat_tci, metric = "rho","clr", alpha=NA,p=100)
best <- rho[">", .60]
plot(best)
#
rhomatrixtci<-getMatrix(rhotci)

tiff("Figures/Networks/Tci16sand18s_rho.tiff",res=300,units="cm",width=40,height=30)
corrgram(rhomatrixtci,order=TRUE,main="Rho proportionality between OTUs present in T citrina 16s and 18s samples ",
         lower.panel=panel.shade, upper.panel=panel.cor, text.panel=panel.txt)
dev.off()
png("Figures/Networks/Tci16sand18s_rho.png",res=300,units="cm",width=50,height=40)
corrgram(rhomatrixtci,order=TRUE,main="Rho proportionality between OTUs present in T citrina 16s and 18s samples ",
         lower.panel=panel.shade, upper.panel=panel.cor, text.panel=panel.txt)
dev.off()

## Visual Network
library(igraph)
library(tidygraph)
links_tci<-as_data_frame(graph_from_adjacency_matrix(rhomatrixtci,weighted=TRUE))
otus_tci<-names(filteredmat_tci)
str(otus_tci)
nodes_tci<-names(otus_tci)
net<-graph.data.frame(links_tci,nodes_tci,directed=F)
plot(net)
net <-simplify(net, remove.multiple = T, remove.loops = T,)
delete.vertices(simplify(net), degree(simplify(net))==0)
head(links_tci)
links_tci_filter1<-links_tci[links_tci$weight>=0.7,]
links_tci_filter2<-links_tci[links_tci$weight<=(-0.7),]
links_tcix<-rbind(links_tci_filter1,links_tci_filter2)
dim(links_tcix)
nodes_tcix<-unique(links_tcix$from)
str(nodes_tcix)#13 names
links_tcix2$weight<-abs(links_tcix$weight)
#links_tcix2<-cbind(links_tcix2,interactiontype)
net<-graph_from_data_frame(links_tcix2,nodes_tcix,directed=F)
net
net <-simplify(net, remove.multiple = T, remove.loops = T,)
#deg <-degree(net, mode="all")
head(net)
#head(nodes_tci)
colrs <-c(rep("blue",6), rep("red",7)) 
E(net)$width <-E(net)$weight
l<-layout_nicely(net)
plot(net)
E(net)
View(links_tcix3)
links_tcix3<-links_tcix[c(1,11,12,2,4,5,8,22,23,24,13,16),]
links_tcilabel<-links_tcix3%>%
  mutate_if(is.numeric, round, 2)
colrs <-c("blue","red","red","blue","blue","blue","blue","red","red","red","red","red") 
png("./Figures/Networks/networkdiagramtci2.png",units="in",width=6,height=6,res=300)
plot(net,vertex.size=25,edge.color=colrs,vertex.label.color="black",vertex.label.cex=0.3,vertex.color="lightgray",
     layout=l,edge.label=links_tcilabel$weight,edge.label.cex=0.4,edge.label.color=colrs,vertex.label.family="Arial",
     edge.lty=3,edge.width=0.5,edge.label.font=2,edge.label.family="Arial",vertex.label.font=2)#pt.bg=colrs ,edge.width=E(net)$weight
dev.off()
??vertex.label.family
###for tme
####
links_tme<-as_data_frame(graph_from_adjacency_matrix(rhomatrixtme,weighted=TRUE))
otus_tme<-names(filteredmat_tme)
str(otus_tme)
#dataset_tme<-c("dataset16s","dataset16s","dataset16s","dataset16s","dataset18s","dataset18s")
#nodes_tci<-cbind(otus_tci,dataset_tci)
nodes_tme<-names(otus_tme)
net<-graph.data.frame(links_tme,nodes_tme,directed=F)
plot(net)
net <-simplify(net, remove.multiple = T, remove.loops = T,)
delete.vertices(simplify(net), degree(simplify(net))==0)
head(links_tme)
links_tme_filter1<-links_tme[links_tme$weight>=0.7,]
links_tme_filter2<-links_tme[links_tme$weight<=(-0.7),]
links_tmex<-rbind(links_tme_filter1,links_tme_filter2)
dim(links_tmex)
nodes_tmex<-unique(links_tmex$from)
str(nodes_tmex)#13 names
links_tmex2<-links_tmex
links_tmex2$weight<-abs(links_tmex$weight)
#links_tcix2<-cbind(links_tcix2,interactiontype)
net<-graph_from_data_frame(links_tmex2,nodes_tmex,directed=F)
net <-simplify(net, remove.multiple = T, remove.loops = T,)
#deg <-degree(net, mode="all")
#head(nodes_tci)
colrs <-c(rep("blue",2)) 
E(net)$width <-E(net)$weight
l<-layout_nicely(net)
plot(net)

#set decimal points to 2
E(net)
view(links_tmex)
links_tmex3<-links_tmex[c(1,2),]
links_tmelabel<-links_tmex3%>%
  mutate_if(is.numeric, round, 2)
png("./Figures/Networks/networkdiagramtme2.png",units="in",width=6,height=6,res=300)
plot(net,layout=l,vertex.size=25,edge.color=colrs,vertex.label.color="black",vertex.label.cex=0.4,vertex.color="lightgray",
     edge.label=links_tmelabel$weight,edge.label.cex=0.5,edge.label.color=colrs,vertex.dist=0.3,vertex.label.family="Arial",
     edge.lty=3,edge.width=0.7,edge.label.font=2,edge.label.family="Arial",vertex.label.font=2)#pt.bg=colrs ,edge.width=E(net)$weight
dev.off()


#####
##Tau
links_tau<-as_data_frame(graph_from_adjacency_matrix(rhomatrixtau,weighted=TRUE))
otus_tau<-names(filteredmat_tau)
str(otus_tau)
#dataset_tme<-c("dataset16s","dataset16s","dataset16s","dataset16s","dataset18s","dataset18s")
#nodes_tau<-cbind(otus_tau,dataset_tau)
nodes_tau<-names(otus_tau)
net<-graph.data.frame(links_tau,nodes_tau,directed=F)
plot(net)
net <-simplify(net, remove.multiple = T, remove.loops = T,)
delete.vertices(simplify(net), degree(simplify(net))==0)
head(links_tau)
links_tau_filter1<-links_tau[links_tau$weight>=0.7,]
links_tau_filter2<-links_tau[links_tau$weight<=(-0.7),]
links_taux<-rbind(links_tau_filter1,links_tau_filter2)
dim(links_taux)
nodes_taux<-unique(links_taux$from)
str(nodes_taux)#13 names
links_taux2<-links_taux
links_taux2$weight<-abs(links_taux2$weight)
#links_taux2<-cbind(links_taux2,interactiontype)
net<-graph_from_data_frame(links_taux2,nodes_taux,directed=F)
net <-simplify(net, remove.multiple = T, remove.loops = T,)
#deg <-degree(net, mode="all")
#head(nodes_tau)
colrs <-c(rep("blue",9), rep("red",4)) 
E(net)$width <-E(net)$weight
l<-layout_nicely(net)
plot(net)

View(links_taux2)
#remove edges taken out bc of loop
links_taux3<-links_taux[c(1,22,21,2,3,4,24,23,7,8,10,14,16),]
links_tauxlabel<-links_taux3%>%
  mutate_if(is.numeric, round, 2)
colrs <-c("blue","red","red","blue","blue","blue","red","red",
          "blue","blue","blue","blue","blue") 

png("./Figures/Networks/networkdiagramtau3.png",units="in",width=6,height=6,res=300)
plot(net,vertex.size=25,edge.color=colrs,vertex.label.color="black",vertex.label.cex=0.4,vertex.color="lightgray",
     layout=l,edge.label.cex=0.5,edge.label.color=colrs,edge.label=links_tauxlabel$weight,vertex.label.family="Arial",
     edge.lty=3,edge.width=0.7,edge.label.font=2,edge.label.family="Arial",vertex.label.font=2)#pt.bg=colrs ,edge.width=E(net)$weight
dev.off()

