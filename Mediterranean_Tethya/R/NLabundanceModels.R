#!/usr/bin/env Rscript

###########################################
#
# abundanceModels.R
#
# Copyright (C) 2019 Sergio Vargas
#
# Contact: sergio.vargas@lmu.de
#
# This program is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation,
# either version 3 of the License, or any later version.
#
# This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program.
#
# If not, see <https://www.gnu.org/licenses/>
#
################################################################


###########################################
#
# Usage:
#
# This script is a companion for Vargas and Leiva.
# 
# Set the wd to the R folder provided with the package
# once this has been done sourcing should work
#
#
################################################################

################
#Load required libraries
###########################

library(ggplot2)
library(dplyr)
library(vegan)
library(reshape2)

######################
#Private functions
#
#
# with minor modifications from Analytical Methods Using R and Excel
# by: Mark Gardener
# http://www.gardenersown.co.uk/Education/Lectures/Community%20Ecology%20Support%20Files.htm
#
#
####################################################
#significance test of differences in deviance of RAD models
rad_test<-function (x, conf.level = 0.95, log = TRUE) 
{
  if (inherits(x, what = "radfit.frame") == FALSE) 
    stop("\nInput must be the result of a radfit command\n")
  d <- sapply(x, function(y) unlist(lapply(y$models, deviance)))
  a <- sapply(x, function(y) unlist(lapply(y$models, AIC)))
  dev <- t(d)
  dev <- as.data.frame(dev)
  dev <- stack(dev)
  names(dev) <- c("deviance", "model")
  #levels(dev$model) <- c("Log", "Mb", "BS", "Pre", "Zip")
  if (log == TRUE) 
    dev.aov <- aov(log(deviance) ~ model, data = dev)
  else dev.aov <- aov(deviance ~ model, data = dev)
  dev.hsd <- TukeyHSD(dev.aov, ordered = TRUE, conf.level = conf.level)
  result <- list(original = dev, aov = dev.aov, post.hoc = dev.hsd, 
                 aic = a, deviance = d, log = log, conf.level = conf.level)
  class(result) <- "rad.htest"
  invisible(result)
}

print.rad.test<-function (x, ...) 
{
  devmean <- apply(x$deviance, MARGIN = 1, FUN = mean)
  devstdv <- apply(x$deviance, MARGIN = 1, FUN = sd)
  devSE <- devstdv/sqrt(ncol(x$deviance))
  aicmean <- apply(x$aic, MARGIN = 1, FUN = mean)
  result <- rbind(devmean, devSE, aicmean)
  rownames(result) <- c("Mean deviance", "SE of deviance", 
                        "Mean AIC")
  print(result, ...)
}

plot.rad.test<-function (x, which = "original", las = 1, ...) 
{
  WHICH <- c("original", "post.hoc")
  which <- tolower(which)
  which <- match.arg(which, WHICH)
  bpm <- x$original$model
  #bpm <- reorder(x$original$model, -x$original$deviance, FUN = mean)
  if (which == "original") {
    if (x$log == TRUE) 
      boxplot(log(deviance) ~ bpm, data = x$original, 
              las = las, ...)
    else boxplot(deviance ~ bpm, data = x$original, las = las, 
                 ...)
  }
  if (which == "post.hoc") {
    plot(x$post.hoc, las = las, ...)
  }
}

summary.rad.htest<-function (x, ...) 
{
  saov <- summary(x$aov)
  phoc <- x$post.hoc$model
  result <- list(ANOVA = saov, Tukey = phoc, Logarithmic = x$log, 
                 Confidence.Level = x$conf.level)
  print(result, ...)
}

################
#Input/Output paths
###########################

#bactLoadCorrectedCountsPath<-"../Data/cbas_tempVSctrl.otutab.bactloadCorrected_valuesOnly.csv"
uncorrectedCountsPath<-"./Data/all.otutab_raw.csv"
otuTaxonomyPath<-"./Data/all_taxa.csv"

plotsPath<-"./Figures/"

################
#User defined variables
###########################

plotSVG<-FALSE
plotPNG<-TRUE

################
#Main Source
###########################

countsDF<-read.csv(uncorrectedCountsPath, sep="\t")
otuTaxonomy<-read.csv(otuTaxonomyPath, sep=";")

#get OTU counts and add them to DF
countSum<-apply(countsDF[-1],1,sum)
countsDF<-cbind(countsDF, countSum)

#Make cumulative percentage to filter out OTUs which make up less than the 5% cummulative sample count
pct<-countSum/sum(countSum)

countsDF<-cbind(countsDF, pct)
countsDF<-countsDF%>%
  mutate(cumpct=cumsum(pct))
dim(countsDF)
#bind the taxonomy columns to the otu counts DF
countsWithTaxonomy<-bind_cols(countsDF[order(countsDF$X.OTU.ID),], otuTaxonomy[order(otuTaxonomy$sequence_identifier),])


#Filter up to 95 % additive abundance
countsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$cumpct < 0.95,]
dim(countsWithTaxonomy)
#delete both columns again
countsWithTaxonomy=select(countsWithTaxonomy,-c(cumpct,pct))

#remove OTUs with sum zero if any. Uncomment dim lines if check wanted.
#dim(countsWithTaxonomy)
countsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$countSum != 0,]
#dim(countsWithTaxonomy)
head(countsWithTaxonomy)
# Phylum abundance per sample
# AND NOT Getting rid of "unclassified" phylums
#AND Getting rid of "unclassified" phylums if they were left as blanks
abundanceByPhylumBySample<-countsWithTaxonomy %>% group_by(Phylum) %>% filter(Phylum != "") %>% mutate_at(colnames(countsWithTaxonomy)[2:22],sum) %>% distinct_at(vars(colnames(countsWithTaxonomy)[c(2:22,26)])) 

#########################
#Fit rank-abundance dominance models with phylum data
###########################################
library(lattice)

phylumRAD<-radfit(t(as.data.frame(abundanceByPhylumBySample)[,1:21]))
summary(phylumRAD)
phylumRADPlot<-plot(phylumRAD)

if(plotPNG){
  
  phylumRADPlotPath<-paste(plotsPath,"phylumRADPlot.png",sep="")
  png(phylumRADPlotPath)
  print(phylumRADPlot)
  dev.off()
  
}
plotSVG<-TRUE
if(plotSVG){
  
  phylumRADPlotPath<-paste(plotsPath,"phylumRADPlot.svg",sep="")
  svg(phylumRADPlotPath)
  print(phylumRADPlot)
  dev.off()
  
}
par(mfrow=c(3,7))
phylumRADPlot

#For specie T. auratis
s1tauphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,1]))
s2tauphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,2]))
s3tauphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,3]))
s4tauphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,4]))
s5tauphylumRAD<-rad.lognormal(t(as.data.frame(abundanceByPhylumBySample)[,5]))
s6tauphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,6]))
s7tauphylumRAD<-rad.zipf(t(as.data.frame(abundanceByPhylumBySample)[,7]))
par(mfrow=c(1,2))

max(s1tauphylumRAD$y)
plot(s1tauphylumRAD,pch=1,lty=1,col=2,main="Zipf",ylim=c(1,22800))

points(s2tauphylumRAD)
lines(s2tauphylumRAD)
points(s3tauphylumRAD)
lines(s3tauphylumRAD)
points(s4tauphylumRAD)
lines(s4tauphylumRAD)
points(s6tauphylumRAD)
lines(s6tauphylumRAD)
points(s7tauphylumRAD)
lines(s7tauphylumRAD)
plot(s5tauphylumRAD,pch=1,lty=1,col=11,main="Lognormal")

par(mfrow=c(3,7))
phylumRADPlot
#For specie T.meloni
s1tmephylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,8]))
s2tmephylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,9]))
s3tmephylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,10]))
s4tmephylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,11]))
s5tmephylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,12]))
s6tmephylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,13]))
s7tmephylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,14]))
par(mfrow=c(1,2))
plot(s1tmephylumRAD,pch=1,lty=1,col=7,main="Zipf-mandelbrot")
points(s2tmephylumRAD)
lines(s2tmephylumRAD)
points(s3tmephylumRAD)
lines(s3tmephylumRAD)
points(s4tmephylumRAD)
lines(s4tmephylumRAD)
points(s5tmephylumRAD)
lines(s5tmephylumRAD)
points(s6tmephylumRAD)
lines(s6tmephylumRAD)
points(s7tmephylumRAD)
lines(s7tmephylumRAD)
###########
#For species T.citroni
s1tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,15]))
s2tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,16]))
s3tciphylumRAD<-rad.lognormal(t(as.data.frame(abundanceByPhylumBySample)[,17]))
s4tciphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,18]))
s5tciphylumRAD<-rad.zipfbrot(t(as.data.frame(abundanceByPhylumBySample)[,19]))
s6tciphylumRAD<-rad.preempt(t(as.data.frame(abundanceByPhylumBySample)[,20]))
s7tciphylumRAD<-rad.lognormal(t(as.data.frame(abundanceByPhylumBySample)[,21]))
par(mfrow=c(1,3))
plot(s1tciphylumRAD,main="Zipf-mandelbrot")
points(s2tciphylumRAD)
lines(s2tciphylumRAD)
points(s5tciphylumRAD)
lines(s5tciphylumRAD)

plot(s3tciphylumRAD,main="Lognormal")
points(s7tciphylumRAD)
lines(s7tciphylumRAD)
plot(s4tciphylumRAD,main="Preemption")
points(s6tciphylumRAD)
lines(s6tciphylumRAD)
###############
#############
#Plot everything together
png("./Figures/16splots_2/phylumrad_grouped.png")
par(mfrow=c(3,3))
#For sp 1
plot(s1tauphylumRAD,pch=1,lty=1,col=2,ylim=c(1,22800),xlim=c(0,30))
mtext(side=3,"Zipf",cex=1)
points(s2tauphylumRAD,col=2)
lines(s2tauphylumRAD,col=2)
points(s3tauphylumRAD,col=2)
lines(s3tauphylumRAD,col=2)
points(s4tauphylumRAD,col=2)
lines(s4tauphylumRAD,col=2)
points(s6tauphylumRAD,col=2)
lines(s6tauphylumRAD,col=2)
points(s7tauphylumRAD,col=2)
lines(s7tauphylumRAD,col=2)
plot(s5tauphylumRAD,pch=1,lty=1,col=11,ylim=c(1,22800),xlim=c(0,30))
mtext(side=3,"Lognormal",cex=1)
plot.new()
#For sp 2
plot(s1tmephylumRAD,pch=1,lty=1,col=7,ylim=c(1,22800),xlim=c(0,30))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(s2tmephylumRAD,col=7)
lines(s2tmephylumRAD,col=7)
points(s3tmephylumRAD,col=7)
lines(s3tmephylumRAD,col=7)
points(s4tmephylumRAD,col=7)
lines(s4tmephylumRAD,col=7)
points(s5tmephylumRAD,col=7)
lines(s5tmephylumRAD,col=7)
points(s6tmephylumRAD,col=7)
lines(s6tmephylumRAD,col=7)
points(s7tmephylumRAD,col=7)
lines(s7tmephylumRAD,col=7)
plot.new()
plot.new()
#########################
#For species 3
plot(s4tciphylumRAD,col=6)
mtext(side=3,"Preemption",cex=1)
points(s6tciphylumRAD,col=6)
lines(s6tciphylumRAD,col=6)
plot(s3tciphylumRAD,col=11,ylim=c(1,22800),xlim=c(0,30))
mtext(side=3,"Lognormal",cex=1)
points(s7tciphylumRAD,col=11)
lines(s7tciphylumRAD,col=11)
plot(s1tciphylumRAD,col=7,ylim=c(1,22800),xlim=c(0,30))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(s2tciphylumRAD,col=7)
lines(s2tciphylumRAD,col=7)
points(s5tciphylumRAD,col=7)
lines(s5tciphylumRAD,col=7)
dev.off()


phylumRAD_test<-rad_test(phylumRAD, conf.level = 0.95, log = T)
summary.rad.htest(phylumRAD_test)
print.rad.test(phylumRAD_test)
plot.rad.test(phylumRAD_test)
plot.rad.test(phylumRAD_test, "post.hoc")
##########################
#rank-abundance dominance with OTU data
##########################################################
#Without unclassified phyla
otuRAD<-radfit(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(2:22)]))

summary(otuRAD)
otuRADPlot<-plot(otuRAD)

if(plotPNG){
  
  otuRADPlotPath<-paste(plotsPath,"otuRADPlot.png",sep="")
  png(otuRADPlotPath)
  print(otuRADPlot)
  dev.off()
  
}

if(plotSVG){
  
  otuRADPlotPath<-paste(plotsPath,"otuRADPlot.svg",sep="")
  svg(otuRADPlotPath)
  print(otuRADPlot)
  dev.off()
  
}

#To group the same model lines for each sp
par(mfrow=c(3,7))
otuRADPlot
s1otuRAD<-rad.zipf(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(2)]))
s2otuRAD<-rad.zipf(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(3)]))
s3otuRAD<-rad.zipf(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(4)]))
s4otuRAD<-rad.zipf(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(5)]))
s5otuRAD<-rad.zipf(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(6)]))
s6otuRAD<-rad.zipf(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(7)]))
s7otuRAD<-rad.zipf(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(8)]))
tme_s1otuRAD<-rad.zipfbrot(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(9)]))
tme_s2otuRAD<-rad.zipfbrot(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(10)]))
tme_s3otuRAD<-rad.preempt(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(11)]))
tme_s4otuRAD<-rad.zipfbrot(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(12)]))
tme_s5otuRAD<-rad.zipfbrot(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(13)]))
tme_s6otuRAD<-rad.zipfbrot(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(14)]))
tme_s7otuRAD<-rad.zipfbrot(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(15)]))
tci_s1otuRAD<-rad.zipf(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(16)]))
tci_s2otuRAD<-rad.zipfbrot(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(17)]))
tci_s3otuRAD<-rad.zipf(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(18)]))
tci_s4otuRAD<-rad.zipf(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(19)]))
tci_s5otuRAD<-rad.zipf(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(20)]))
tci_s6otuRAD<-rad.zipf(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(21)]))
tci_s7otuRAD<-rad.zipf(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(22)]))

max(otuRAD$y) #HIGHEST VAL IS 21367

#Plot everything together
png("./Figures/16splots_2/oturad_grouped.png")
par(mfrow=c(3,3))
#For sp 1
plot(s1otuRAD,pch=1,lty=1,col=2,ylim=c(1,21300),xlim=c(0,1050))
mtext(side=3,"Zipf",cex=1)
points(s2otuRAD,col=2)
lines(s2otuRAD,col=2)
points(s3otuRAD,col=2)
lines(s3otuRAD,col=2)
points(s4otuRAD,col=2)
lines(s4otuRAD,col=2)
points(s5otuRAD,col=2)
lines(s5otuRAD,col=2)
points(s6otuRAD,col=2)
lines(s6otuRAD,col=2)
points(s7otuRAD,col=2)
lines(s7otuRAD,col=2)
plot.new()
plot.new()
#colours: Log normal is col 2, zipf is 4, mandel is 3, preempt is 5

#For sp 2
plot(tme_s3otuRAD,pch=1,lty=1,col=6,ylim=c(1,21300),xlim=c(0,1050))
mtext(side=3,"Preemption",cex=1)
plot(tme_s1otuRAD,pch=1,lty=1,col=7,ylim=c(1,21300),xlim=c(0,1050))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(tme_s2otuRAD,col=7)
lines(tme_s2otuRAD,col=7)
points(tme_s4otuRAD,col=7)
lines(tme_s4otuRAD,col=7)
points(tme_s5otuRAD,col=7)
lines(tme_s5otuRAD,col=7)
points(tme_s6otuRAD,col=7)
lines(tme_s6otuRAD,col=7)
points(tme_s7otuRAD,col=7)
lines(tme_s7otuRAD,col=7)


plot.new()
#########################
#For species 3
plot(tci_s1otuRAD,col=2,ylim=c(1,22800),xlim=c(0,30))
mtext(side=3,"Zipf",cex=1)
points(tci_s3otuRAD,col=2)
lines(tci_s3otuRAD,col=2)
points(tci_s4otuRAD,col=2)
lines(tci_s4otuRAD,col=2)
points(tci_s5otuRAD,col=2)
lines(tci_s5otuRAD,col=2)
points(tci_s6otuRAD,col=2)
lines(tci_s6otuRAD,col=2)
points(tci_s7otuRAD,col=2)
lines(tci_s7otuRAD,col=2)
plot(tci_s2otuRAD,col=7,ylim=c(1,22800),xlim=c(0,30))
mtext(side=3,"Zipf-mandelbrot",cex=1)
plot.new()

dev.off()


otuRAD_test<-rad_test(otuRAD, conf.level = 0.95, log = T)
summary.rad.htest(otuRAD_test)
print.rad.test(otuRAD_test)
par(mfrow=c(1,1))
plot.rad.test(otuRAD_test)
plot.rad.test(otuRAD_test, "post.hoc")
#######

# AND
#########################
#Fit rank-abundance dominance models  for NOT Getting rid of OTUs belonging to "unclassified" phyla
###########################################
countsWithTaxonomy[countsWithTaxonomy==""]<-"Unclassified"
head(countsWithTaxonomy)
##########################
#rank-abundance dominance with OTU data
#With unclassified phyla
xotuRAD<-radfit(t(countsWithTaxonomy[c(2:22)]))
summary(xotuRAD)
xotuRADPlot<-plot(xotuRAD)
xotuRADPlot

if(plotPNG){
  
  xotuRADPlotPath<-paste(plotsPath,"xotuRADPlot.png",sep="")
  png(xotuRADPlotPath)
  print(xotuRADPlot)
  dev.off()
  
}

if(plotSVG){
  
  xotuRADPlotPath<-paste(plotsPath,"xotuRADPlot.svg",sep="")
  svg(xotuRADPlotPath)
  print(xotuRADPlot)
  dev.off()
  
}

par(mfrow=c(3,7))
xotuRADPlot
xs1otuRAD<-rad.zipf(t(countsWithTaxonomy[c(2)]))
xs2otuRAD<-rad.zipf(t(countsWithTaxonomy[c(3)]))
xs3otuRAD<-rad.zipf(t(countsWithTaxonomy[c(4)]))
xs4otuRAD<-rad.zipf(t(countsWithTaxonomy[c(5)]))
xs5otuRAD<-rad.zipf(t(countsWithTaxonomy[c(6)]))
xs6otuRAD<-rad.zipf(t(countsWithTaxonomy[c(7)]))
xs7otuRAD<-rad.zipf(t(countsWithTaxonomy[c(8)]))
xtme_s1otuRAD<-rad.lognormal(t(countsWithTaxonomy[c(9)]))
xtme_s2otuRAD<-rad.zipfbrot(t(countsWithTaxonomy[c(10)]))
xtme_s3otuRAD<-rad.zipfbrot(t(countsWithTaxonomy[c(11)]))
xtme_s4otuRAD<-rad.lognormal(t(countsWithTaxonomy[c(12)]))
xtme_s5otuRAD<-rad.preempt(t(countsWithTaxonomy[c(13)]))
xtme_s6otuRAD<-rad.lognormal(t(countsWithTaxonomy[c(14)]))
xtme_s7otuRAD<-rad.zipfbrot(t(countsWithTaxonomy[c(15)]))
xtci_s1otuRAD<-rad.zipfbrot(t(countsWithTaxonomy[c(16)]))
xtci_s2otuRAD<-rad.zipfbrot(t(countsWithTaxonomy[c(17)]))
xtci_s3otuRAD<-rad.zipf(t(countsWithTaxonomy[c(18)]))
xtci_s4otuRAD<-rad.zipf(t(countsWithTaxonomy[c(19)]))
xtci_s5otuRAD<-rad.zipf(t(countsWithTaxonomy[c(20)]))
xtci_s6otuRAD<-rad.zipf(t(countsWithTaxonomy[c(21)]))
xtci_s7otuRAD<-rad.zipf(t(countsWithTaxonomy[c(22)]))

#Plot everything together

png("./Figures/16splots_2/xoturad_grouped.png")
par(mfrow=c(3,3))
#For sp 1
plot(xs1otuRAD,pch=1,lty=1,col=2,ylim=c(1,21300),xlim=c(0,1050))
mtext(side=3,"Zipf",cex=1)
points(xs2otuRAD,col=2)
lines(xs2otuRAD,col=2)
points(xs3otuRAD,col=2)
lines(xs3otuRAD,col=2)
points(xs4otuRAD,col=2)
lines(xs4otuRAD,col=2)
points(xs5otuRAD,col=2)
lines(xs5otuRAD,col=2)
points(xs6otuRAD,col=2)
lines(xs6otuRAD,col=2)
points(xs7otuRAD,col=2)
lines(xs7otuRAD,col=2)
plot.new()
plot.new()


#For sp 2
plot(xtme_s5otuRAD,pch=1,lty=1,col=6,ylim=c(1,21300),xlim=c(0,1050))
mtext(side=3,"Preemption",cex=1)
plot(xtme_s1otuRAD,pch=1,lty=1,col=11,ylim=c(1,21300),xlim=c(0,1050))
mtext(side=3,"Lognormal",cex=1)
points(xtme_s4otuRAD,col=11)
lines(xtme_s4otuRAD,col=11)
points(xtme_s6otuRAD,col=11)
lines(xtme_s6otuRAD,col=11)

plot(xtme_s2otuRAD,pch=1,lty=1,col=7,ylim=c(1,21300),xlim=c(0,1050))
mtext(side=3,"Zipf-Mandelbrot",cex=1)
points(xtme_s3otuRAD,col=7)
lines(xtme_s3otuRAD,col=7)
points(xtme_s7otuRAD,col=7)
lines(xtme_s7otuRAD,col=7)

#########################
#For species 3

plot(xtci_s3otuRAD,col=2,ylim=c(1,22800),xlim=c(0,30))
mtext(side=3,"Zipf",cex=1)
points(xtci_s4otuRAD,col=2)
lines(xtci_s4otuRAD,col=2)
points(xtci_s5otuRAD,col=2)
lines(xtci_s5otuRAD,col=2)
points(xtci_s6otuRAD,col=2)
lines(xtci_s6otuRAD,col=2)
points(xtci_s7otuRAD,col=2)
lines(xtci_s7otuRAD,col=2)
plot(xtci_s1otuRAD,col=7,ylim=c(1,22800),xlim=c(0,30))
mtext(side=3,"Zipf-mandelbrot",cex=1)
points(xtci_s2otuRAD,col=7)
lines(xtci_s2otuRAD,col=7)
plot.new()

dev.off()

xotuRAD_test<-rad_test(xotuRAD, conf.level = 0.95, log = T)
summary.rad.htest(xotuRAD_test)
print.rad.test(xotuRAD_test)
par(mfrow=c(1,2),mai=c(1,1,1,1))
plot.rad.test(xotuRAD_test)
plot.rad.test(xotuRAD_test, "post.hoc")


#EOF


