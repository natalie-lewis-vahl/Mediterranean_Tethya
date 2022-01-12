library(vegan)
library(matrixStats)
library(dplyr)
matrix_file<-"./Data/18s/all.otutab.csv"

otus<-read.csv("./Data/18s/all.otutab.csv", sep="\t")
#BInd taxonomy to remove bacteria
taxa<-read.csv("./Data/18s/taxa_fixed18s.csv",sep=",")
otus<-arrange(otus, X.OTU.ID)
taxa<-arrange(taxa,sequence_identifier)
otu_and_taxa<-bind_cols(otus,taxa)
otu_and_taxa<-otu_and_taxa%>%filter(Domain!="Bacteria")
#Take out OTU 1 for coral amplification
otus_and_taxa<-otu_and_taxa[!(otu_and_taxa$X.OTU.ID=="OTU_1"),]
#Uncomment to also take out OTU 2
otus_and_taxa<-otu_and_taxa[!(otu_and_taxa$X.OTU.ID=="OTU_2"),]
head(otus_and_taxa)
#Remove OTUS with less than 5 reads
countSum<-apply(otus_and_taxa[2:45],1,sum)
otus_and_taxa<-cbind(otus_and_taxa,countSum)
#same thing as using size column from taxa dataset
otus_and_taxa<-otus_and_taxa[otus_and_taxa$countSum>=5,]

d18s_matrix<-otus_and_taxa[,-c(46:57)]
#rEMOVE SEQUENCE id AS A COLUMN AND MAKE THEM THE ROW NAMES
row.names(d18s_matrix)<-d18s_matrix$X.OTU.ID
#remove column with id names and taxa columns again
d18s_matrix<-d18s_matrix[,-1]

#if filtering is wanted at a more thorough level of 95% 
#countSum<-apply(d18s_matrix,1,sum)#sum the rows
#pct<-countSum/sum(countSum)

#d18s_matrix<-cbind(d18s_matrix, countSum, pct)
#d18s_matrix[desc(d18s_matrix$pct),]
#d18s_matrix<-d18s_matrix%>%
#  mutate(cumpct=cumsum(pct))

#d18s_matrix<-d18s_matrix%>%
#  arrange(desc(pct))
#d18s_matrix<-d18s_matrix[d18s_matrix$cumpct < 0.95,]

#358 rows with only OTU1 removed and 499 rows with OTU2 also removed

#groups 
#T aurantium: GW1941 till GW1951 inclusive + GW1984 (12 SAMPLES)
#T meloni: GW1952 till GW1962 +GW1982 +GW1983 (13 SAMPLES)
# T citroni: GW1963 till GW1981 inclusive (19 SAMPLES)

#group vectors for non bactloadCorrected datasets
num_groups<-c(1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,1)
label_groups<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")

head(d18s_matrix)
dim(d18s_matrix)
nmds<-metaMDS(t(d18s_matrix[,-c(45,46,47)]), try=50)

png("Figures/18splots/ordi_nmds_18s2.png",units="cm",height=30,width=30,res=300)
plot(nmds, type = "t", display="sites")

points(nmds, col=num_groups, cex=1.5, pch=16)

ordispider(nmds, label_groups, label=T)

dev.off()
#Check for significant structure
x<-anosim(t(d18s_matrix[,-c(45,46,47)]), num_groups)
hist(x)
x
x$statistic
hist(x$perm)
abline(v=x$statistic)

#vector 
species_vector<-data.frame(Species=c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau"))

#wFreqs_more50Counts % more1PctCounts

#adonis(t(bact_matrix)~Condition, data=temp_treatments, strata=temp_treatments$Sponge_type)
bact_cca<-cca(t(d18s_matrix[,-c(45,46,47)])~Species, data=species_vector)
#bact_cca<-cca(t(bact_matrix[-c(17,18)])~Condition+Sponge_type+Condition:Sponge_type, data=temp_treatments)

png("Figures/18splots/ordi_cca_18s.png",units="cm",height=15,width=15,res=300)

plot(bact_cca, type="t", display="sites")
points(bact_cca, col=num_groups, cex=1.5, pch=16)
ordispider(bact_cca, label_groups, label=T)
dev.off()
summary(bact_cca)
anova(bact_cca)
anova(bact_cca, by="term")

anova(bact_cca, by="margin")

#highly variable otus

is.matrix(d18s_matrix)
xd18s_matrix<-as.matrix(d18s_matrix)
head(xd18s_matrix)
matrix_rMeans<-rowMeans(xd18s_matrix[,-c(45,46,47)])
matrix_rSds<-rowSds(xd18s_matrix[,-c(45,46,47)])
matrix_zscores<-((xd18s_matrix[,-c(45,46,47)]-matrix_rMeans)/matrix_rSds)


OTUVars<-sort(rowVars(xd18s_matrix[,-c(45,46,47)],useNames = T),decreasing = T)

png("Figures/18splots/heatmap_all_18s.png",units="cm",height=30,width=20,res=300)
heatmap(as.matrix(matrix_zscores),ColSideColors = c("green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","red","red","red","red","red","red","red","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","red","red","green"))
dev.off()
variability<-rowVars(xd18s_matrix[,-c(45,46,47)])
OTUVars<-cbind(d18s_matrix[,-c(45,46,47)],variability)
OTUVars<-arrange(OTUVars,by="variability")
variableOTUs<-rownames(head(OTUVars, n=25L))

png("Figures/18splots/heatmap_25mostvariable_18s.png",units="cm",height=20,width=20,res=300)
heatmap(as.matrix(matrix_zscores[rownames(matrix_zscores) %in% variableOTUs,]),ColSideColors = c("green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","red","red","red","red","red","red","red","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue"))
dev.off()

####Vegemite
data(dune)
sponges<-d18s_matrix[,-c(45,46,47)]

dis<-vegdist(t(sponges))
clus<-hclust(dis,"ward.D2")
plot(clus)
#single sites to be joined into large clusters
#or use complete for compact clusters or average -neutral grouping
den<-as.dendrogram(clus)
#ord<-cca(sponges)
wa<-scores(nmds,display="sites",choices=1)
oden<-reorder(den,wa,mean)
vegemite(t(sponges),use=oden,zero="-",scale="log")
tabasco(t(sponges),use=oden,zero="-",scale="Hill") #HIll is the scale publication standard

tabasco(t(sponges),use=oden, zero="-",col = heat.colors(5), scale="Hill")
