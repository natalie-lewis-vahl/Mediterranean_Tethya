library(vegan)
library(matrixStats)
library(dplyr)
matrix_file<-"./Data/16s_allsamples_otu.csv"

bact_matrix<-read.csv(matrix_file, head=T, sep=";", row.names=1)
countSum<-apply(bact_matrix,1,sum)#sum the rows
pct<-countSum/sum(countSum)

bact_matrix<-cbind(bact_matrix, countSum, pct)
bact_matrix<-bact_matrix%>%
  arrange(desc(pct))
bact_matrix<-bact_matrix%>%
  mutate(cumpct=cumsum(pct))
dim(bact_matrix)

#if filtering is wanted up to 95% cumulative sum of abundance can
#use cumpct column
#bact_matrix<-bact_matrix[bact_matrix$cumpct < 95,] 
head(bact_matrix)

bact_matrix<-bact_matrix[bact_matrix$countSum >=50,]


dim(bact_matrix)


#groups 
#Tau:"GW1941" "GW1942" "GW1943" "GW1944" "GW1945" "GW1946" "GW1947" "GW1948" "GW1949" "GW1950" "GW1951" "GW1984"
#Tme:"GW1952" "GW1953" "GW1954" "GW1955" "GW1956" "GW1957" "GW1958" "GW1959" "GW1960" "GW1961" "GW1962" "GW1982" "GW1983" 
#Tci:"GW1963" "GW1964" "GW1965" "GW1966" "GW1967" "GW1968" "GW1969" "GW1970" "GW1971" "GW1972" "GW1973""GW1974" "GW1975" "GW1976" "GW1977" "GW1978" "GW1979""GW1980" "GW1981"

#group vectors for non bactloadCorrected datasets
num_groups<-c(1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,1)
label_groups<-c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau")

#head(bact_matrix)
#creeate plot and delete columns with cumpct, pct and cumsum
nmds<-metaMDS(t(bact_matrix[,-c(45:47)]), try=50)

#With previous incomplete data set got filtered data get Warning message: stress is (nearly) zero: you may have insufficient data
png("./Figures/AS_16splots/nmdsplot_16s.png",height=8,width=8,units="in",res=300)

plot(nmds, type = "t", display="sites")

points(nmds, col=c("green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","red","red","red","red","red","red","red","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","red","red","green"), cex=1.5, pch=16)

ordispider(nmds, label_groups, label=T)
dev.off()
#anosim(t(bact_matrix), num_groups)

#vector 
species_vector<-data.frame(Species=c("Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tau","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tme","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tci","Tme","Tme","Tau"))

#wFreqs_more50Counts % more1PctCounts

#adonis(t(bact_matrix)~Condition, data=temp_treatments, strata=temp_treatments$Sponge_type)
??cca
bact_cca<-cca(t(bact_matrix[,-c(45:47)])~Species, data=species_vector)
#bact_cca<-cca(t(bact_matrix[-c(17,18)])~Condition+Sponge_type+Condition:Sponge_type, data=temp_treatments)
png("./Figures/AS_16splots/ccaplot_16s.png",height=8,width=8,units="in",res=300)

plot(bact_cca, type="t", display="sites")
points(bact_cca, col=c("green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","red","red","red","red","red","red","red","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","red","red","green"), cex=1.5, pch=16)
ordispider(bact_cca, label_groups, label=T)
dev.off()
summary(bact_cca)
anova(bact_cca)
anova(bact_cca, by="term")

anova(bact_cca, by="margin")

#cca with all core species inclusive (defined in AS_NLcoreCommunity script as present
#in 90% of samples of any one of the species).
coreOTUs<-c("OTU_11",  "OTU_3" , "OTU_15",  "OTU_21" ,"OTU_9",  "OTU_25", "OTU_27", "OTU_19" , "OTU_20", "OTU_52", "OTU_14", "OTU_36",
            "OTU_28", "OTU_43", "OTU_114", "OTU_66",  "OTU_93", "OTU_46", "OTU_73", "OTU_65", "OTU_47", "OTU_55", "OTU_103", "OTU_54", 
            "OTU_102", "OTU_119", "OTU_279", "OTU_90","OTU_96", "OTU_63", "OTU_132", "OTU_219", "OTU_360", "OTU_124", "OTU_171", "OTU_215",
            "OTU_210", "OTU_146", "OTU_531", "OTU_182", "OTU_91",  "OTU_153", "OTU_341", "OTU_129","OTU_242", "OTU_888", "OTU_312", 
            "OTU_466","OTU_331", "OTU_455", "OTU_858", "OTU_899")
coreOnlyCounts<-bact_matrix[rownames(bact_matrix[,-c(45:47)]) %in% as.character(coreOTUs),-c(45:47)]
#creeate plot and delete columns with cumpct, pct and cumsum
corenmds<-metaMDS(t(bact_matrix[rownames(bact_matrix[,-c(45:47)]) %in% as.character(coreOTUs),-c(45:47)]))

#With previous incomplete data set got filtered data get Warning message: stress is (nearly) zero: you may have insufficient data
png("./Figures/AS_16splots/nmdsplot_core16s.png",height=8,width=8,units="in",res=300)

plot(corenmds, type = "t", display="sites")

points(corenmds, col=c("green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","red","red","red","red","red","red","red","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","red","red","green"), cex=1.5, pch=16)

ordispider(corenmds, label_groups, label=T)
dev.off()
#cca plot
core_cca<-cca(t(coreOnlyCounts)~Species, data=species_vector)#correspondence analysis,

png("./Figures/AS_16splots/ccaplot_core16s.png",height=8,width=8,units="in",res=300)
plot(core_cca, type="t", display="sites")
points(core_cca, col=c("green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","red","red","red","red","red","red","red","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","red","red","green"), cex=1.5, pch=16)
ordispider(core_cca, label_groups, label=T)
dev.off()
anova(core_cca)

#non-core community OTUs (those not included in the all core species in any species)
#nmds
noncorenmds<-metaMDS(t(bact_matrix[!(rownames(bact_matrix[,-c(45:47)]) %in% as.character(coreOTUs)),-c(45:47)]))

#With previous incomplete data set got filtered data get Warning message: stress is (nearly) zero: you may have insufficient data
png("./Figures/AS_16splots/nmdsplot_noncore16s.png",height=8,width=8,units="in",res=300)

plot(noncorenmds, type = "t", display="sites")

points(noncorenmds, col=c("green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","red","red","red","red","red","red","red","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","red","red","green"), cex=1.5, pch=16)

ordispider(noncorenmds, label_groups, label=T)
dev.off()
#cca
nonCoreCounts<-bact_matrix[!(rownames(bact_matrix[,-c(45:47)]) %in% as.character(coreOTUs)),-c(45:47)]
nonCore_cca<-cca(t(nonCoreCounts)~Species, data=species_vector)
png("./Figures/AS_16splots/2ccaplot_noncore16s.png",height=8,width=8,units="in",res=300)
plot(nonCore_cca, type="t", display="sites")
points(nonCore_cca, col=c("green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","red","red","red","red","red","red","red","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","red","red","green"), cex=1.5, pch=16)
ordispider(nonCore_cca, label_groups, label=T)
dev.off()
anova(nonCore_cca)



#vegemite(t(coreOnlyCounts), core_cca, scale = "log")


#highly variable otus
head(bact_matrix)
is.matrix(bact_matrix)
xbact_matrix<-as.matrix(bact_matrix[,-c(45:47)])
bact_matrix_rMeans<-rowMeans(xbact_matrix)
bact_matrix_rSds<-rowSds(xbact_matrix)
bact_matrix_zscores<-((xbact_matrix-bact_matrix_rMeans)/bact_matrix_rSds)


png("./Figures/AS_16splots/heatmap_variability_16s_allotus.png",height=8,width=8,units="in",res=300)
heatmap(as.matrix(bact_matrix_zscores),ColSideColors =c("green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","red","red","red","red","red","red","red","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","red","red","green"))
dev.off()
#Display large plot with all OTUs visible-to be able to zoom
png("./Figures/AS_16splots/heatmap_variability_16s_allotusfull.png",height=30,width=30,units="in",res=300)
heatmap(as.matrix(bact_matrix_zscores),cexRow=0.5,ColSideColors =c("green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","red","red","red","red","red","red","red","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","red","red","green"))
dev.off()
??heatmap
#Heat map with 25 most variable OTUs
variability<-rowVars(xbact_matrix)
OTUVars<-cbind(bact_matrix,variability)
OTUVars<-arrange(OTUVars,by="variability")
variableOTUs<-rownames(head(OTUVars, n=25L))

png("./Figures/AS_16splots/heatmap_25mostvariable_otus.png",height=8,width=8,units="in",res=300)
heatmap(as.matrix(bact_matrix_zscores[rownames(bact_matrix_zscores) %in% variableOTUs,]),ColSideColors =c("green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","red","red","red","red","red","red","red","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","darkblue","red","red","green"))
dev.off()
