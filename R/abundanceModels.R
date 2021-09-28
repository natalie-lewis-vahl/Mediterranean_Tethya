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
uncorrectedCountsPath<-"../Data/cbas_tempVSctrl.otutab.csv"
otuTaxonomyPath<-"../Data/cbas_otu_taxonomy.csv"

plotsPath<-"../Plots/"

################
#User defined variables
###########################

plotSVG<-TRUE
plotPDF<-FALSE

################
#Main Source
###########################

countsDF<-read.csv(uncorrectedCountsPath, sep="\t")
otuTaxonomy<-read.csv(otuTaxonomyPath, sep="\t")

#get OTU counts and add them to DF
countSum<-apply(countsDF[-1],1,sum)
countsDF<-cbind(countsDF, countSum)

#bind the taxonomy columns to the otu counts DF
countsWithTaxonomy<-bind_cols(countsDF[order(countsDF$X.OTU.ID),], otuTaxonomy[order(otuTaxonomy$sequence_identifier),])

#remove OTUs with sum zero if any. Uncomment dim lines if check wanted.
#dim(countsWithTaxonomy)
countsWithTaxonomy<-countsWithTaxonomy[countsWithTaxonomy$countSum != 0,]
#dim(countsWithTaxonomy)


# Phylum abundance per sample
abundanceByPhylumBySample<-countsWithTaxonomy %>% group_by(Phylum) %>% filter(Phylum != "") %>% mutate_at(colnames(countsWithTaxonomy)[2:21],sum) %>% distinct_at(vars(colnames(countsWithTaxonomy)[c(2:21,25)])) 

#########################
#Fit rank-abundance dominance models with phylum data
###########################################

phylumRAD<-radfit(t(as.data.frame(abundanceByPhylumBySample)[,1:20]))
summary(phylumRAD)
phylumRADPlot<-plot(phylumRAD)

if(plotPDF){
  
  phylumRADPlotPath<-paste(plotsPath,"phylumRADPlot.pdf",sep="")
  pdf(phylumRADPlotPath)
  print(phylumRADPlot)
  dev.off()
  
}

if(plotSVG){
  
  phylumRADPlotPath<-paste(plotsPath,"phylumRADPlot.svg",sep="")
  svg(phylumRADPlotPath)
  print(phylumRADPlot)
  dev.off()
  
}

phylumRAD_test<-rad_test(phylumRAD, conf.level = 0.95, log = T)
summary.rad.htest(phylumRAD_test)
print.rad.test(phylumRAD_test)
plot.rad.test(phylumRAD_test)
plot.rad.test(phylumRAD_test, "post.hoc")


#####################################
#rank-abundance dominance with OTU data
##########################################################
otuRAD<-radfit(t(countsWithTaxonomy[(countsWithTaxonomy$Phylum != ""),c(2:21)]))

summary(otuRAD)
otuRADPlot<-plot(otuRAD)

if(plotPDF){
  
  otuRADPlotPath<-paste(plotsPath,"otuRADPlot.pdf",sep="")
  pdf(otuRADPlotPath)
  print(otuRADPlot)
  dev.off()
  
}

if(plotSVG){
  
  otuRADPlotPath<-paste(plotsPath,"otuRADPlot.svg",sep="")
  svg(otuRADPlotPath)
  print(otuRADPlot)
  dev.off()
  
}


otuRAD_test<-rad_test(otuRAD, conf.level = 0.95, log = T)
summary.rad.htest(otuRAD_test)
print.rad.test(otuRAD_test)
plot.rad.test(otuRAD_test)
plot.rad.test(otuRAD_test, "post.hoc")

#EOF


