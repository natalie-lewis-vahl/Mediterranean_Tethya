library(dplyr)
library(tibble)
library(igraph)
library(network)
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
countsDF<-countsDF%>%
  mutate(cumpct=cumsum(pct))
dim(countsDF)
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
int_matrix <- otu_table[,grep('OTU', names(otu_table))] 
#Make into zeros and 1s
int_matrix[int_matrix >= 1] <- 1

#how many samples is each otu in
otu_samples <- colSums(int_matrix)
active_otus <- names(otu_samples[otu_samples > 10])


require(igraph)
vs <- unique(otu_table$Sample)
g_whole <- graph.empty() + vertices(vs)


for(v in vs){
  print(v)
  cur_samps <- which(otu_table$Sample == v)
  replicates <- otu_table[cur_samps, ]
  replicates <- replicates[,grep('OTU', names(replicates))]
  #   replicates[replicates >= 1] <- 1
  
  possible_otus <- colSums(replicates)
  possible_otus <- names(possible_otus[possible_otus >= 1])
  
  ##############
  
  ##### comment/uncomment if you are basing OTU selection on the number of samples
  ##### the are present in (active otus)
  possible_otus <- intersect(possible_otus, active_otus)
  
  ##############
  
  if(length(possible_otus) == 0){
    print(paste('this sponge does not have OTUs:', v))
  }
  
  replicates <- replicates[,which(names(replicates) %in% possible_otus)]
  
  ##################################################################
  ####### this bit is for relative abundances
  
  rel_abs <- replicates/rowSums(replicates)
  rel_abs <- apply(rel_abs, 2, mean)
  
  ##################################################################
  
  ###### whereas this is for mean abundances....
  abs <- colMeans(replicates)
  ##################################################################
  
  #   replicates <- replicates[,which(colnames(replicates) %in% possible_otus)] 
  #   subset <- replicates[subset,]
  #   n <- names(subset)
  
  #   n <- names(replicates)
  
  n <- possible_otus
  
  g_whole <- g_whole + vertices(n[!n %in% get.vertex.attribute(g_whole, 'name')])  
  g_whole[v,n] <- 1
  
  g_whole <- set.edge.attribute(g_whole, 'abs', which(E(g_whole) %in% E(g_whole)[v %--% n]) , as.matrix(abs))
  
  g_whole <- set.edge.attribute(g_whole, 'rel_abs', which(E(g_whole) %in% E(g_whole)[v %--% n]) , as.matrix(rel_abs))
  
}

otus_in_network <- V(g_whole)$name[! V(g_whole)$name %in% vs]
x <- degree(g_whole, v=otus_in_network)
mean_y <- c()
max_y <- c()

for(v in names(x)){
  #print(v)
  edges <- neighbors(g_whole, v, mode='in')
  abundances <- E(g_whole)[v %--% edges]$rel_abs
  mean_y <- append(mean_y, mean(abundances))
  max_y <- append(max_y, max(abundances))
}


new_y <- as.vector(log10((mean_y)*100))

plot(x,new_y, xlab='OTU degree', ylab='Log10(Mean relative abundance)')
lo <- smooth.spline(x, new_y) #, control = loess.control(surface = "direct"))
lines(predict(lo), col='green', lwd=2)

linear <- lm(new_y ~ x)
abline(linear, col='blue', lwd=2)

e_edges <- E(g)[from(x)]
e_top_weights <- order(e_edges$weight, decreasing=TRUE)[1:3]
E(g)[ as.vector(e_edges)[e_top_weights] ]



# for(i in 1:dim(int_matrix)[1]){
#   
#   if(max(int_matrix[i,]) > 36000){
#     print(i)  
#   }
#   
# }
 

#####Next part of script



######### Script for the results presented in the Nature Comms paper ##########

otus_in_network <- V(g_whole)$name[! V(g_whole)$name %in% vs]
x <- degree(g_whole, v=otus_in_network)

specialists <- names(x[x<5])
opportunists <- names(x[which(x>=2 & x<=2)])
generalists <- names(x[x>5])

total_abs <- colSums(replicates)

total_abs <- log10((total_abs+1))
relative_abs <- total_abs/sum(total_abs)

per_groups <- data.frame(relative_abs, group=NA)

per_groups[which(row.names(per_groups) %in% specialists),]$group <- 'specialists'
per_groups[which(row.names(per_groups) %in% opportunists),]$group <- 'opportunists'
per_groups[which(row.names(per_groups) %in% generalists),]$group <- 'generalists'

per_groups$group <- factor(per_groups$group, levels=c('specialists', 'opportunists', 'generalists'))

pdf("relative-abundances-vs-group.pdf", width=21/2.54, height=21/2.54)
par(mar=c(4.5,5.5,4,4), bty='n')

boxplot(relative_abs~group, data=per_groups, outline=F, ylab='Relative log(abundance)',cex.lab=1.5, cex.axis=1.5, pch=16, tck=0)

means <- tapply(per_groups$relative_abs,per_groups$group,mean)
points(means,col="red",pch=18, cex=1.5)

dev.off()

########################## BEGIN OF FIGURE OF THE CUMULATIVE DEGREE DISTRIBUTIONS ##########################


###### This piece of code produces the figure included in the manuscript showing the cumulative
###### degree distributions of sponges and OTUs

otus_deg <- degree(g_whole, mode='in')[which(degree(g_whole, mode='in') >= 1)]
sponges_deg <- degree(g_whole, mode='out')[which(degree(g_whole, mode='out') >= 1)]

######## we want to show the degree distribution of a random bipartite graph alongside the original network
draw_random <- TRUE

if(draw_random){
  n_otus <- length(which(degree(g_whole, mode='in') != 0))
  n_sponges <- length(which(degree(g_whole, mode='out') != 0))
  
  random_net <- bipartite.random.game(n_sponges, n_otus, typ='gnm', m=ecount(g_whole), directed=TRUE)
  
  otus_deg_rand <- degree(random_net, mode='in')[which(degree(random_net, mode='in') >= 1)]
  sponges_deg_rand <- degree(random_net, mode='out')[which(degree(random_net, mode='out') >= 1)]
}

pdf("dds-both-axis-fits.pdf", width=21/2.54, height=21/2.54)
par(mar=c(4.5,5.5,4,4), bty='n')
occur = as.vector(table(otus_deg))
occur = occur/sum(occur)
p = occur/sum(occur)
y = rev(cumsum(rev(p)))
x = as.numeric(names(table(otus_deg)))

plot(x, y, log="xy", xlab ='log k', ylab='', type='p', cex.lab=1.5, cex.axis=1.5, col='black', pch=16, tck=0, yaxt='n', xaxt='n', xlim=c(1,100))
title(ylab='log Pc(k)', cex.lab = 1.5, line = 4)

box(lwd=1.5)
temp <- data.frame(x,y)

mod1 <- nls(y ~ (x^-a*exp(-x/b)), data = temp, start = list(a = 1, b = 1))

#m <- lm(y ~ x, data = temp)

# add fitted curve
lines(temp$x, predict(mod1, list(x = temp$x)))

if(draw_random){
  occur_temp = as.vector(table(otus_deg_rand))
  occur_temp = occur_temp/sum(occur_temp)
  p_temp = occur_temp/sum(occur_temp)
  y_temp = rev(cumsum(rev(p_temp)))
  x_temp = as.numeric(names(table(otus_deg_rand)))
  
  points(x_temp, y_temp, log="xy", xlab ='', ylab='', type='p', cex.lab=1.5, cex.axis=1.5, col='blue', pch=16, tck=0, yaxt='n', xaxt='n')
}


y1 <- floor(log10(range(y)))
pow <- seq(y1[1], y1[2]+1)
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, 10^pow, tck=0.02, lwd=1.5, lwd.ticks=1.5, cex.axis=1.5, labels=c(expression(10^{-5}), expression(10^{-4}), expression(10^{-3}), expression(10^{-2}), expression(10^{-1}), expression(10^{0}), expression(10^{1})), las=1)
axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1, tck=0.01)


x1 <- floor(log10(range(x)))
pow <- seq(x1[1]-1, x1[2]+2)

ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, 10^pow, tck=0.02, lwd=1.5, lwd.ticks=1.5, cex.axis=1.5, labels=c(expression(10^{-1}), expression(10^{0}), expression(10^{1}), expression(10^{2}), expression(10^{3})))
axis(1, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1, tck=0.01)


occur = as.vector(table(sponges_deg))
occur = occur/sum(occur)
p = occur/sum(occur)
y = rev(cumsum(rev(p)))
x = as.numeric(names(table(sponges_deg)))

par(new=TRUE)
plot(x,y, type='p', log='xy', col='red', pch=16, cex.lab=1.5, cex.axis=1.5, ylim=c(10^-2,10^0), xlim=c(10^1, 10^4), tck=0, yaxt='n', xaxt='n', xlab='', ylab='')

temp <- data.frame(x,y)

mod2 <- nls(y ~ (exp(-x/b)), data = temp, start = list(b = 8))

# add fitted curve
lines(temp$x, predict(mod2, list(x = temp$x)), col='red')

if(draw_random){
  occur_temp = as.vector(table(sponges_deg_rand))
  occur_temp = occur_temp/sum(occur_temp)
  p_temp = occur_temp/sum(occur_temp)
  y_temp = rev(cumsum(rev(p_temp)))
  x_temp = as.numeric(names(table(sponges_deg_rand)))
  
  points(x_temp, y_temp, log='xy', col='orange', pch=16, xlab='', ylab='')
}  


y1 <- floor(log10(range(y)))
pow <- seq(y1[1]-1, y1[2]+1)
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))

axis(4, 10^pow, col='red', tck=0.02, lwd=1.5, lwd.ticks=1.5, cex.axis=1.5, labels=c(expression(10^{-3}), expression(10^{-2}), expression(10^{-1}), expression(10^{0}), expression(10^{1})), las=1)
axis(4, ticksat, labels=NA, tcl=-0.25, col='red', lwd=0, lwd.ticks=1, tck=0.01)

x1 <- floor(log10(range(x)))
pow <- seq(x1[1]-1, x1[2]+2)
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))

axis(3, 10^pow, col='red', tck=0.02, lwd=1.5, lwd.ticks=1.5, cex.axis=1.5, labels=c(expression(10^{0}), expression(10^{1}), expression(10^{2}), expression(10^{3}), expression(10^{4}), expression(10^{5})))
axis(3, ticksat, col='red', labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1, tck=0.01)



dev.off()

########################## END OF FIGURE OF THE CUMULATIVE DEGREE DISTRIBUTIONS ##########################


#### BEGIN OF FIGURE FOR THE FRACTION OF SAMPLES FOUND ACROSS HOST SPECIES VS OTU DEGREE (NUMBER OF HOSTS) ####

###### plot of degree vs. fraction of appearance
###### for the mean fraction of appearance across all the sponge species the OTU appears

otu_indexes <- grep('Otu', get.vertex.attribute(g_whole, 'name'))
otus <- get.vertex.attribute(g_whole, 'name')

fractions_found <- c()

for(otu in otu_indexes){
  
  sponges_for_otu <- get.vertex.attribute(g_whole, 'name')[neighbors(g_whole, otu, mode='all')]
  
  mean_fraction <- c()
  for(s in sponges_for_otu){
    sponge_samples <- otu_table[which(otu_table$sponge == s),which(colnames(otu_table) == otus[otu])]
    sponge_samples[sponge_samples > 1] <- 1
    fraction_found_in_s <- mean(sponge_samples)
    
    if(is.nan(fraction_found_in_s)){
      next
    }
    
    mean_fraction <- append(mean_fraction, fraction_found_in_s)
  }
  
  #print(mean_fraction)
  if(length(mean_fraction) > 1){
    mean_fraction <- mean(mean_fraction)
  }else{
    mean_fraction <- mean_fraction[1]
  }
  fractions_found <- append(fractions_found, mean_fraction)
}


lo <-  smooth.spline(degree(g_whole, otu_indexes), fractions_found, spar=0.35) # loess(fractions_found~degree(g_whole, otu_indexes)) #

pdf("degree-vs-fraction-average.pdf", width=21/2.54, height=21/2.54)
par(mar=c(4.5,5.5,4,4), bty='n')

plot(degree(g_whole, otu_indexes), fractions_found, ylab='Mean fraction of samples found across species', xlab='Degree (number of hosts)', main='Number of hosts vs. fraction of samples across species')
lines(predict(lo), col='red', lwd=2)

dev.off()

#### END OF FIGURE FOR THE FRACTION OF SAMPLES FOUND ACROSS HOST SPECIES VS OTU DEGREE (NUMBER OF HOSTS) ####

