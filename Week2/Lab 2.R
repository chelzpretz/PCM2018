## Lab 2 


library(phytools)
library(geiger)

## partice 

tr <- pbtree(n=20) # a pure birth tree, because it's easy
# a correlated pair of characters on a phylogeny:
xy <- sim.char(tr, matrix(c(20,15,15,20), 2, byrow = T))[,,1]
dimnames(xy)[[2]] <- c('x', 'y')
layout(matrix(1:2, 1, 2))
fancyTree(tr, "contmap", x = xy[, 1])
phylomorphospace(tr, xy, label = 'horizontal', cex = 0.5)

plot(xy)
abline(lm(y ~ x, as.data.frame(xy)))
print(summary(lm(y ~ x, as.data.frame(xy))))


tr <- reorder(tr, 'postorder')
tr.orig <- tr # hanging onto this for now b/c we'll be mucking up the tree later
plot(tr)
nodelabels()
edgelabels()
# the following binds tip labels and node labels as a third column to the edge matrix
tr.edge <- cbind(tr$edge, endlabel = tr$tip.label[tr$edge[, 2]])
tr.edge[is.na(tr.edge[, 'endlabel']), 'endlabel'] <- 
  as.character(tr.edge[is.na(tr.edge[, 'endlabel']), 2])
print(head(tr.edge, 20)) # first 20 rows of the tr$edge, with tip labels



nodesToDo <- unique(tr$edge[, 1])
ind.contrasts <- rep(NA, length(nodesToDo))
names(ind.contrasts) <- sort(as.character(nodesToDo))


contrasts.sd <- ind.contrasts


x <- xy[, 1]
means <- c(x, ind.contrasts)
names(means)[which(names(means) %in% tr$tip.label)] <- 
  match(names(means)[which(names(means) %in% tr$tip.label)], tr$tip.label)
print(means) # make sure it looks right, and it does


for(working.node in nodesToDo) {
  # do a little book-keeping
  edges <- which(tr$edge[, 1] == working.node)
  desc.values <- means[as.character(tr$edge[edges, 2])]
  desc.lengths <- tr$edge.length[edges]
  branch.to.rescale <- which(tr$edge[, 2] == working.node)
  
  
  # 1. make and store the weighted average
  # means[as.character(working.node)] <- weighted.mean(desc.values, 1/desc.lengths)
  means[as.character(working.node)] <- sum(desc.values * (1/desc.lengths))/sum(1/desc.lengths)
  
  # 2. do the contrast and its sd
  ind.contrasts[as.character(working.node)] <- diff(desc.values)
  contrasts.sd[as.character(working.node)] <- sqrt(sum(desc.lengths))
  
  # 3. rescale the remaining branch
  tr$edge.length[branch.to.rescale] <-
    tr$edge.length[branch.to.rescale] +
    (desc.lengths[1] * desc.lengths[2]) / sum(desc.lengths)
}


ancStatesCompare <- cbind(ape = ace(x, tr.orig, method = 'pic')$ace[as.character(nodesToDo)],
                          ourFct = means[as.character(nodesToDo)])
picCompare <- cbind(ape = pic(x, tr.orig)[as.character(nodesToDo)], 
                    ourFct = ind.contrasts[as.character(nodesToDo)] / contrasts.sd[as.character(nodesToDo)])
layout(matrix(1:2, 1))
plot(ancStatesCompare, main = "comparing node states")
abline(0, 1)
plot(picCompare, main = "comparing normalized contrasts")
abline(0, -1)


pic.pcm35300 <- function(tr, dat, out = c('ancStates', 'contrasts'), normalized = TRUE) {
  ## set up our variables
  tr <- reorder(tr, 'postorder')
  nodesToDo <- unique(tr$edge[, 1])
  ic <- ic.sd <- rep(NA, length(nodesToDo))
  names(ic) <- names(ic.sd) <- sort(as.character(nodesToDo))
  means <- c(dat, ic)
  names(means)[which(names(means) %in% tr$tip.label)] <- 
    match(names(means)[which(names(means) %in% tr$tip.label)], tr$tip.label)
  
  for(working.node in nodesToDo) {
    # do a little book-keeping
    edges <- which(tr$edge[, 1] == working.node)
    desc.values <- means[as.character(tr$edge[edges, 2])]
    desc.lengths <- tr$edge.length[edges]
    branch.to.rescale <- which(tr$edge[, 2] == working.node)
    
    
    # 1. make and store the weighted average
    means[as.character(working.node)] <- weighted.mean(desc.values, 1/desc.lengths)
    
    # 2. do the contrast and its sd
    ic[as.character(working.node)] <- diff(desc.values)
    ic.sd[as.character(working.node)] <- sqrt(sum(desc.lengths))
    
    # 3. rescale the remaining branch
    tr$edge.length[branch.to.rescale] <-
      tr$edge.length[branch.to.rescale] +
      (desc.lengths[1] * desc.lengths[2]) / sum(desc.lengths)
  } # close working.node
  
  if(normalized) ic <- ic / ic.sd
  out <- cbind(ancStates = means[as.character(nodesToDo)],
               contrasts = ic[as.character(nodesToDo)])[, out]
  return(out)
}

pic.pcm35300(tr,x)
x.ic <- pic.pcm35300(tr, xy[, 1], 'contrasts')
y.ic <- pic.pcm35300(tr, xy[, 2], 'contrasts')
plot(x.ic, y.ic)

summary(lm(y.ic ~ x.ic + 0))
summary(lm(y ~ x, as.data.frame(xy)))


#To practice
#On your own or with a partner, try out some of these exercises, ordered approximately by complexity. If the first exercises seem too simple, skip them and jump on to the next.

#Simulate one tree of 50 tips and two datasets on that tree (x and y), using the simulation methods we used above or other R packages. Refer to the "Trait simulations" and "Tree simulations" paragraphs on the CRAN phylogenetics task view (https://cran.r-project.org/web/views/Phylogenetics.html). Generate independent contrasts for x and y, then calculate the correlation coefficient (cor.test in R) or least-squares regression coefficient of determination (lm in R) and p-values in R. Compare your analysis of the independent contrasts against your analysis of the untransformed data.




#Simulate 1000 trees of 50 tips each and two uncorrelated traits on each tree under a Brownian motion process. For each dataset, calculate the mean, the Pearson product-moment (cor.test), and the p-value for the correlation using both standard methods and independent contrasts. How often should you expect to find a P < 0.05? How often do you for both the PIC and standard methods? and which yields a better estimate of the mean? Is either biased? is there a lower variance to one or the other?

#Modify our pic.35300 function to return Garland et al.'s confidence intervals on the root state of the tree.