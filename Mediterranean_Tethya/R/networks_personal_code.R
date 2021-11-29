library(dplyr)
library(tibble)
library(igraph)
library(network)
library(bipartite)
CountsPath<-"./Data/all.otutab_raw.csv"
otuTaxonomyPath<-"./Data/all_taxa.csv"

otutable<-read.csv(CountsPath, sep="\t")
taxatable<-read.csv(otuTaxonomyPath, sep=";")

#get OTU counts and add them to DF
countSum<-apply(otutable[-1],1,sum)
countsDF<-cbind(otutable, countSum)

#Make cumulative percentage to filter out OTUs which make up less than the 5% cummulative sample count
pct<-countSum/sum(countSum)
countsDF<-cbind(countsDF, pct)
countsDF[desc(countsDF$pct),]
countsDF<-countsDF%>%
  mutate(cumpct=cumsum(pct))

#bind the taxonomy columns to the otu counts DF
countsWithTaxonomy<-bind_cols(countsDF[order(countsDF$X.OTU.ID),], taxatable[order(taxatable$sequence_identifier),])

#Filter up to 95 % additive abundance
countsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$cumpct < 0.95,]
dim(countsWithTaxonomy)
#delete both columns again
countsWithTaxonomy=select(countsWithTaxonomy,-c(cumpct,pct))
#Make a separate otu table for later analysis
otu_table=select(countsWithTaxonomy,-c(23:30))
otu_table<-column_to_rownames(otu_table, "X.OTU.ID")
otu_table<-as.data.frame(t(otu_table))

#calculate the total number of reads per sample
otu_table$total <- rowSums(otu_table)
otu_table$sponge <- c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci")
#GET NAMES OF ALL OUTS
otus <- names(otu_table)[-c(2391,2392)]
#PUT SAMPLES BACK AS A ROW
otu_table<-tibble::rownames_to_column(otu_table, "Sample")
xmatrix <- otu_table[,grep('OTU', names(otu_table))] 
fullmatrix<-xmatrix
#try to plot network
plotweb(xmatrix)
#Need to filter out OTUs not part of the network
#how many samples is each otu in
xmatrix[xmatrix>=1]<-1
otu_samples <- colSums(xmatrix)
active_otus <- names(otu_samples[otu_samples > 10])
filteredmat<- otu_table[ ,which((names(otu_table) %in% active_otus)==TRUE)]
#Matrix now filtered for those that are in at least 10 samples
#of all species however uses absolute abundance from counts
plotweb(filteredmat)
visweb(filteredmat,type="none")
#use presence or absence data
yesorno<-filteredmat
yesorno[yesorno>=1]<-1
tiff("./Figures/networks/bipartite_allsp_presorabs.tiff", units="in", width=10, height=5, res=300)
visweb(yesorno,type="none",labsize=2,ylabel="Sponge sample")
dev.off()
#use relative abundance
relabund<-apply(filteredmat,1,function(i) i/sum(i))
??apply
as.matrix(relabund)
relabund2<-t(relabund)
tiff("./Figures/networks/bipartite_allsp_relabund.tiff", units="in", width=10, height=5, res=300)
visweb(relabund2,type="none",labsize=2)
dev.off()
#bipartite networks per species
taumatrix<-xmatrix[1:7,]
otu_samplestau <- colSums(taumatrix)
#present in at least 3
active_otustau <- names(otu_samplestau[otu_samplestau > 3])
filteredmat_tau<- otu_table[ ,which((names(otu_table) %in% active_otustau)==TRUE)]
plotweb(filteredmat_tau)
visweb(filteredmat_tau,type="diagonal")

#use presence or absence data
filteredmat_tau[filteredmat_tau>=1]<-1
visweb(filteredmat_tau,type="diagonal")
#use relative abundance
relabund_tau<-apply(filteredmat_tau,1,function(i) i/sum(i))
??apply
as.matrix(relabund_tau)
relabund2_tau<-t(relabund_tau)
visweb(relabund2_tau,type="none")

#Do for the other 2 species
#repeat but replicate as relative abundance in all samples and compare for the 3 species