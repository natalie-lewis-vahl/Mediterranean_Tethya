library(propr)
library(dplyr)
library(tibble)
library(ggdendro)
CountsPath<-"./Data/16s_allsamples_otu.csv"
otuTaxonomyPath<-"./Data/16s_allsamples_taxa.csv"
#Need to remove sample GW1956 because of low sampling -> the 16th sample
#belonging to Tme
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
tcicountsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$countSum >=50,]

#Make a separate otu table for later analysis
otu_table=select(countsWithTaxonomy,-c(46:54))
otut<-otu_table[,-1]
rownames(otut)<-otu_table[,1]
otu_table<-as.data.frame(t(otut))

#calculate the total number of reads per sample
#remove bad sample
which(otu_table=="GW1956")
head(otu_table)
otu_table<-otu_table[-16,]
otu_table$total <- rowSums(otu_table)
otu_table$sponge <- c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")
#GET NAMES OF ALL OUTS
otus <- names(otu_table)[-c(412,413)]
#PUT SAMPLES BACK AS A ROW
otu_table<-tibble::rownames_to_column(otu_table, "Sample")
xmatrix <- otu_table[,grep('OTU', names(otu_table))] 
fullmatrix<-xmatrix
str(xmatrix)
#Need to filter out OTUs not part of the network
#how many samples is each otu in
xmatrix[xmatrix>=1]<-1
otu_samples <- colSums(xmatrix)
active_otus <- names(otu_samples[otu_samples >26])
counts<- otu_table[ ,which((names(otu_table) %in% active_otus)==TRUE)]
species<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")
sample<-c("GW1941","GW1942" ,"GW1943","GW1944", "GW1945", "GW1946", "GW1947", "GW1948","GW1949","GW1950","GW1951",
           "GW1952","GW1953", "GW1954", "GW1955", "GW1957", "GW1958", "GW1959","GW1960","GW1961","GW1962",
           "GW1963","GW1964","GW1965", "GW1966", "GW1967", "GW1968", "GW1969","GW1970","GW1971","GW1972", "GW1973","GW1974", "GW1975", "GW1976", "GW1977", "GW1978", "GW1979","GW1980","GW1981", "GW1982","GW1983", "GW1984")
#Note that the log-ratio transformation, by its nature, fails if the input data contain any zero values. 
#By default, this function replaces all zero values with 1. Alternatively, the user may set the parameter `alpha` greater than zero to approximate log-ratios in the presence of zeros (via the Box-Cox transformation). However, the topic of zero replacement is controversial. Proceed carefully when analyzing data that contain zero values.
rho <- propr(counts, sample, metric="rho",ivar=0,p=100)
rho <- propr(counts, metric = "rho", alpha=NA)
best <- rho[">", .90]
plot(best)
dendrogram(best)
??propr

best<-propr::simplify(best)
pca(best, group = sample)
#Second, we look at `snapshot`, a function for visualizing the intensity of the log-ratio transformed data across samples. Heatmap intensity is not scaled.

snapshot(best)

#Finally, we look at the `prism` and `bokeh` functions for visualizing the co-clustering of proportional features. We mention these plots together because they share some key similarities. First, they are all "index-naive". Second, they identify the feature pairs where both constituents co-cluster (with the total number of clusters toggled by `k`). Third, they return a vector of cluster memberships for all features in the `propr` object.

#The `prism` function plots the variance of the ratio of the log-ratio transformed feature pair (VLR) versus the sum of the individual variances of each log-ratio transformed feature (VLS). The ratio of the VLR to the VLS equals $1 - \rho$. As such, we use here seven rainbow colored lines to indicate where $\rho = [.01, .05, .50, 0, 1.50, 1.95, 1.99]$, going from red to violet. A low VLR with a high VLS suggests that the feature pair remains in an equilibrium despite high variability among the individual features.

#```{r, dpi = 66, fig.width = 8, fig.height = 8, results = "hide", message = FALSE}
clusts <- prism(best, k = 4)

#The `bokeh` function plots pairs across the individual variances of the constituent log-ratio transformed features. For clarity of visualization, this figure projects the data on a log-fold scale. Therefore, the highly variable co-clusters appear in the top-right of the figure while the lowly variable co-clusters appear in the bottom-left. Meanwhile, highly proportional pairs tend to aggregate around the $y=x$ diagonal.

#```{r, dpi = 66, fig.width = 8, fig.height = 8, results = "hide", message = FALSE}
clusts <- bokeh(best, k = 4)
```

#These plots can help us conceptualize high-dimensional data and select a highly proportional module for further analysis. In this example, we choose co-cluster 2 for further analysis because it shows high proportionality in the setting of high individual feature variance.

## Down-stream

#We can extract co-cluster 2 from the `propr` object using the `subset` method.

#```{r, results = "hide"}
sub <- subset(best, select = (clusts == 1))
#We can use the `pca` function to see how well this cluster differentiates the two experimental groups. We see here that the first two principal components of this highly proportional module results in good separation between the experimental groups. This matches the separation achieved in the source publication which used `edgeR` for feature selection (Rollins 2015), suggesting that our unsupervised method identified a feature subset that is relevant to the experimental condition.
#{r, dpi = 66, fig.width = 8, fig.height = 8, results = "hide", fig.keep = "last"}
pca(sub, group = sample)
#Having identified a cluster that separates the experimental groups, the next step in this pipeline might involve a gene set enrichment analysis (GSEA) of the transcripts participating in this highly proportional module. The application of GSEA is beyond the scope of this vignette, although we can easily extract the names of the transcripts that belong to this cluster.

transcripts <- colnames(sub@logratio)

#
getAdj(rho)
getMatrix(rho)



########################
###
##Tethya citrina
head(xmatrix)
head(tcimatrix)
tcimatrix_abund<-fullmatrix[23:41,]
tcimatrix<-xmatrix[23:41,]
otu_samplestci <- colSums(tcimatrix)
#present in more than half-(=N.samples/2 +1)
active_otustci <- names(otu_samplestci[otu_samplestci >17])
dim()
counts_tci<- tcimatrix_abund[ ,which((names(tcimatrix) %in% active_otustci)==TRUE)]
View(counts_tci) 
###
#Note that the log-ratio transformation, by its nature, fails if the input data contain any zero values. 
#By default, this function replaces all zero values with 1. Alternatively, the user may set the parameter `alpha` greater than zero to approximate log-ratios in the presence of zeros (via the Box-Cox transformation). However, the topic of zero replacement is controversial. Proceed carefully when analyzing data that contain zero values.
rho <- propr(counts_tci, metric = "rho","clr", alpha=NA,p=100)
best <- rho[">", .5]
plot(best)
dendrogram(best)

best<-propr::simplify(best)
pca(best)
#Second, we look at `snapshot`, a function for visualizing the intensity of the log-ratio transformed data across samples. Heatmap intensity is not scaled.

snapshot(best)

#Finally, we look at the `prism` and `bokeh` functions for visualizing the co-clustering of proportional features. We mention these plots together because they share some key similarities. First, they are all "index-naive". Second, they identify the feature pairs where both constituents co-cluster (with the total number of clusters toggled by `k`). Third, they return a vector of cluster memberships for all features in the `propr` object.

#The `prism` function plots the variance of the ratio of the log-ratio transformed feature pair (VLR) versus the sum of the individual variances of each log-ratio transformed feature (VLS). The ratio of the VLR to the VLS equals $1 - \rho$. As such, we use here seven rainbow colored lines to indicate where $\rho = [.01, .05, .50, 0, 1.50, 1.95, 1.99]$, going from red to violet. A low VLR with a high VLS suggests that the feature pair remains in an equilibrium despite high variability among the individual features.

#```{r, dpi = 66, fig.width = 8, fig.height = 8, results = "hide", message = FALSE}
clusts <- prism(best, k = 3)

#The `bokeh` function plots pairs across the individual variances of the constituent log-ratio transformed features. For clarity of visualization, this figure projects the data on a log-fold scale. Therefore, the highly variable co-clusters appear in the top-right of the figure while the lowly variable co-clusters appear in the bottom-left. Meanwhile, highly proportional pairs tend to aggregate around the $y=x$ diagonal.

#```{r, dpi = 66, fig.width = 8, fig.height = 8, results = "hide", message = FALSE}
clusts <- bokeh(best, k = 3)
```

#These plots can help us conceptualize high-dimensional data and select a highly proportional module for further analysis. In this example, we choose co-cluster 2 for further analysis because it shows high proportionality in the setting of high individual feature variance.

## Down-stream

#We can extract co-cluster 2 from the `propr` object using the `subset` method.

#```{r, results = "hide"}
sub <- subset(best, select = (clusts == 1))
#We can use the `pca` function to see how well this cluster differentiates the two experimental groups. We see here that the first two principal components of this highly proportional module results in good separation between the experimental groups. This matches the separation achieved in the source publication which used `edgeR` for feature selection (Rollins 2015), suggesting that our unsupervised method identified a feature subset that is relevant to the experimental condition.
#{r, dpi = 66, fig.width = 8, fig.height = 8, results = "hide", fig.keep = "last"}
pca(sub)
#Having identified a cluster that separates the experimental groups, the next step in this pipeline might involve a gene set enrichment analysis (GSEA) of the transcripts participating in this highly proportional module. The application of GSEA is beyond the scope of this vignette, although we can easily extract the names of the transcripts that belong to this cluster.

transcripts <- colnames(sub@logratio)

#
getAdj(rho)
rhomatrix<-getMatrix(rho)
dim(rhomatrix)
View(rhomatrix)
max(rhomatrix)
min(rhomatrix)
??heatmap
library(corrgram)
??ggcorr
tiff("Figures/Networks/Tcitrina16s_rho.tiff",res=300,units="cm",width=30,height=25)
corrgram(rhomatrix,order=TRUE,main="Rho proportionality between bacteria present in T citrina samples ",
         lower.panel=panel.shade, upper.panel=panel.cor, text.panel=panel.txt)
dev.off()
########################


##Tethy aurantium
head(xmatrix)
taumatrix<-xmatrix[c(1:11,44),]
otu_samplestau <- colSums(taumatrix)
#present in more than half-(=N.samples/2 +1)
active_otustci <- names(otu_samplestci[otu_samplestau = 12])
filteredmat_tau<- otu_table[ ,which((names(otu_table) %in% active_otustau)==TRUE)]
dim(filteredmat_tau) 
