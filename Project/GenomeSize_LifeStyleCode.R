## PCM Project ##
## This code is to look at how genome size effects plant life style (annual/perinnal)

## That data that will be anaylsis is genome sizes from the Solanaceae family from Kew's database and other websites

# load needed libraries 
library(phytools)
library(ape)
library(nlme)
library(geiger)

# set working directory 
setwd("/Users/Chelsea/Documents/PhD/PCM2018/PCM2018/Project")

## First import tree ## 
soltree <- read.tree(file="smoothprunedtree.tre") ## saved within the same folder
plot(soltree) # plot tree to make sure it's working
taxa <- soltree$tip.label #pull taxon list

# Upload dataset with life form
data <- read.csv("Data_LF1.csv")

#### Prune tree to only data I have ####
colnames(data)
data <- transform(data, newcol=paste(Genus, Species, sep="_"))# need to make genus and species into one row  
length(colnames(data)) # find new column 
states <- structure(data[[6]], names = data[[9]]) #Pull out rows/species that have life form data 
length(states) # this tell you how many states there are
as.character(names(states)) ## need to change vector from factor to character.. 

# Prune the trees to only keep species with known state 
to.drop <- setdiff(taxa, names(states)) #finds which names in tree NOT also in names of states
length(to.drop) - length(states) # how many tips should be left 
prunedtree <- drop.tip(soltree, to.drop) #removes those taxa, new tree has 

smoothprunedtree<-chronopl(prunedtree, lambda=1, age.min=50, CV=TRUE) #make ultrametric, necessary for PCMs
# that step is really slow so only repeat if you need to change the taxon composition
plot(smoothprunedtree, cex = 0.3) #looks good
write.tree(smoothprunedtree, "LifeForm_smoothprunedtree.tre")
 

########
# still need to remove taxa from dataset that are not in the tree
taxa2 <- smoothprunedtree$tip.label #pull taxon list for pruned tree 
states2 <- structure(data[[5]], names = data[[9]]) # take only the genome size from the pruned tree

to.drop2 <- setdiff(taxa2, names(states2)) #finds which names in tree NOT also in names of states
length(states2) - length(to.drop2) # how many tips should be left
#### This line isn't working
tr <- drop.tip(smoothprunedtree, to.drop2) #removes those taxa, new tree has 186 tips
#write.tree(tr) # this should be the same as smoothprunetree

# have chromosome number is a state
states3 <- structure(data[[3]], names = data[[9]])


## make data and tree data match 
df <- data.frame(states, states2, states3, row.names = data[,9])
rname <- names(states2)
row.names(df) <- rname
match <- match(tr$tip.label, rownames(df)) # this matches the data and with tree tips
sorteddata <- df[,][match,] # this get cleans the data


########
# PIC
# source("pic_pcmFunc.R") ## call in function from Hipp's class
#pic.pcm35300(tr, sorteddata)
# x.ic <- pic.pcm35300(tr, data[,2], 'contrasts')
# y.ic <- pic.pcm35300(tr, data[, 1], 'contrasts')
# plot(x.ic, y.ic)

# summary(lm(y.ic ~ x.ic + 0))  # p-value 0.19
# summary(lm(data[,1] ~ data[,2], as.data.frame(data))) # p-value 0.03

#### When consistering PIC there is no corralation ######

## Different PIC -- can't do with discrete traits 
pic.x <- pic(sorteddata$states2, tr, scaled = TRUE, var.contrasts = FALSE,
    rescaled.tree = FALSE)
pic.y <- pic(sorteddata$states, tr, scaled = TRUE, var.contrasts = FALSE,
             rescaled.tree = FALSE)
cor.test(pic.x, pic.y)
fit1 <- lm(pic.y ~ pic.x - 1)
fit2 <- lm(pic.x ~ pic.y - 1)
summary(fit1)
summary(fit2)
anova(fit1)


# Phylo-Anova 

#phylANOVA(tr, x = sorteddata$states, y = sorteddata$states2, nsim=5000) # looking at lf to genome size

## Anova test from Stacey -- looking at genome size ~ life form
fit<-gls(states2 ~ states, sorteddata) #gls of the two varibles 
summary(fit) #p= 2e-04, lnLik= -194.88, r^2 = 1.330237
anova(fit)

# chromosome number ~ life form
fit1<-gls(states3 ~ states, sorteddata) #gls of the two varibles 
summary(fit1) #p= 2e-04, lnLik= -194.88, r^2 = 1.330237
anova(fit1)

# Genome size ~ chromosome number
fit2<-gls(states2 ~ states3, sorteddata) #gls of the two varibles 
summary(fit2) #p= 2e-04, lnLik= -194.88, r^2 = 1.330237
anova(fit2)

# Genome size ~ chromosome number + life form
fit3<-gls(states2 ~ states + states3, sorteddata) #gls of the two varibles 
summary(fit3) #p= 2e-04, lnLik= -194.88, r^2 = 1.330237
anova(fit3)

## have vectors of just the states
statesV <- structure(sorteddata[[1]], names = row.names(sorteddata))
statesV2 <- structure(sorteddata[[2]], names = row.names(sorteddata))
statesV3 <- structure(sorteddata[[3]], names = row.names(sorteddata))
 
 ############
obj <- contMap(tr,  statesV3) #from http://www.phytools.org/eqg2015/asr.html
plot(obj, type = "fan", legend = 0.7 * max(nodeHeights(tr)), fsize = c(0.6, 0.9)) #fan tree
plot(obj, lwd = 1.25, legend = 0.7 * max(nodeHeights(tr)), fsize = c(0.3, 0.7)) #rectangular

#phenogram(tr, sorteddata$states2, fsize = 0.3, method = "pic") # ugh, that's busy!

#plotTree(tr, sorteddata$states2, tiplabels=FALSE)

# Maybe do a reversible and irreverible tree? 

###  Pagel's lambda using the "fitDiscrete()" function in the R package geiger
fitDiscrete(tr, sorteddata, model=c("ER"), treeTransform=c("none", "lambda", "kappa", "delta", "linearChange", "exponentialChange", "twoRate"), data.names=NULL, plotlnl=F, qLimits=c(0.0001, 1000), pLimits=c(0.00001, 10))

### Phylogenetic Eigenvectors Regression.
obj <- contMap(tree, lfmat, plot = FALSE)
