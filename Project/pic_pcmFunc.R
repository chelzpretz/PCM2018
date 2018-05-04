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