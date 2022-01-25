library(propr)
library(dplyr)
library(tibble)
library(igraph)
library(network)
library(bipartite)
install.packages("uft8")
library
CountsPath<-"./Data/16s_allsamples_otu.csv"
otuTaxonomyPath<-"./Data/16s_allsamples_taxa.csv"

otutable<-read.csv(CountsPath, sep="\t")
taxatable<-read.csv(otuTaxonomyPath, sep=";")
#reorder otu to match taxa
#get OTU counts and add them to DF
countSum<-apply(otutable[-1],1,sum)
countsDF<-cbind(otutable, countSum)

#bind the taxonomy columns to the otu counts DF
countsWithTaxonomy<-bind_cols(countsDF[order(countsDF$X.OTU.ID),], taxatable[order(taxatable$sequence_identifier),])
head(countsWithTaxonomy)
#Filter from 50 counts up
countsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$countSum >=50,]
dim(countsWithTaxonomy)
#Make a separate otu table for later analysis
otu_table=select(countsWithTaxonomy,-c(46:54))
otut<-otu_table[,-1]
rownames(otut)<-otu_table[,1]
otu_table<-as.data.frame(t(otut))

#calculate the total number of reads per sample
otu_table$total <- rowSums(otu_table)
otu_table$sponge <- c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")
#GET NAMES OF ALL OUTS
otus <- names(otu_table)[-c(412,413)]
#PUT SAMPLES BACK AS A ROW
otu_table<-tibble::rownames_to_column(otu_table, "Sample")
xmatrix <- otu_table[,grep('OTU', names(otu_table))] 
fullmatrix<-xmatrix
# plot network
png("Figures/Networks/bipartite16s.png",units="cm",width=50,height=30,res=300)
plotweb(fullmatrix)
dev.off()
#Need to filter out OTUs not part of the network
#how many samples is each otu in
xmatrix[xmatrix>=1]<-1
otu_samples <- colSums(xmatrix)
active_otus <- names(otu_samples[otu_samples > 23])
counts<- otu_table[ ,which((names(otu_table) %in% active_otus)==TRUE)]
group<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")

rho <- propd(counts, group,p=100)
??propr
getAdj(filteredmat,cutoff=1)
getMatrix(relabund3)
#Matrix now filtered for those that are in at least 10 samples
#of all species however uses absolute abundance from counts
plotweb(filteredmat)
#use presence or absence data
yesorno<-filteredmat
yesorno[yesorno>=1]<-1
is.matrix(yesorno)
visweb(as.matrix(yesorno))

visweb(yesorno,type="none",labsize=2,ylabel="Sponge sample")
#use relative abundance
relabund<-apply(filteredmat,1,function(i) i/sum(i))

relabund1<-as.matrix(relabund)
relabund2<-t(relabund1)
relabund3<-as.data.frame(relabund)

visweb(relabund2,type="none",labsize=2)
dim(relabund2)#102 OTUs
#bipartite networks per species
taumatrix<-xmatrix[1:7,]
otu_samplestau <- colSums(taumatrix)
#present in more than half-(=N.samples/2 +1)
active_otustau <- names(otu_samplestau[otu_samplestau >= 6])
filteredmat_tau<- otu_table[ ,which((names(otu_table) %in% active_otustau)==TRUE)]
dim(filteredmat_tau) #46 OTUs
visweb(filteredmat_tau,type="diagonal")

#use presence or absence data
filteredmat_tau[filteredmat_tau>=1]<-1
visweb(filteredmat_tau,type="diagonal")
#use relative abundance
relabund_tau<-apply(filteredmat_tau,1,function(i) i/sum(i))
relabund2_tau<-t(as.matrix(relabund_tau))
visweb(relabund2_tau,type="none")


#Do for the other 2 species
#repeat but replicate as relative abundance in all samples and compare for the 3 species