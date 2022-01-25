
xcounts <- matrix(rpois(20*50, 100), 20, 50)
group <- sample(c("A", "B"), size = 20, replace = TRUE)
library(propr)

pr <- propr(counts, # rows as samples, like it should be
            metric = "rho", # or "phi", "phs", "cor", "vlr"
            ivar = "clr", # or can use "iqlr" instead
            alpha = NA, # use to handle zeros
            p = 100) # used by updateCutoffs
??propr
data(iris)
keep<-iris$Species%in% c("setosa","versicolor")
counts<-iris[keep,1:4]*10
group<-ifelse(iris[keep,"Species"]=="setosa","A","B")
dim(counts)
pd<-propd(counts,group,alpha=NA,p=100)

