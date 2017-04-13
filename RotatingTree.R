require(ape)
require(phytools)
require(adephylo)
library(geiger)

###################### Block calculation ####################
blockLevel <- function(listX) {
  nblocks <- 0
  previous <- "randomstring"
  for (i in listX) {
    if (i != previous) {
      nblocks <- nblocks + 1
      previous <- i
    }
  }
  return(nblocks)
}

blockLevelCluster <- function(listX) {
  allTypes<- unique(listX)
  differences <- c()
  for (i in allTypes) {
    differences <- c(differences, (diff(grep(paste0("^",i,"$"), listX))-1))
  }
  diffMean <- sum(differences)
  return(diffMean)
}

############################# Get Trait list function ################################

getTraitList <- function(metadata = metadata, type = type, treeTips = treeTips) {
  traitListType <- sapply(treeTips, function(i) {
    traitsTemp <- metadata[,type][grep(paste0("^",i,"$"), rownames(metadata))]
    return(traitsTemp)
  })
  traitListColor <- sapply(traitListType, function(y) getColor(y))
  names(traitListColor) <- treeTips
  return(traitListColor)
}

###################### Get position from traits ####################

getPositionFromTraits <- function(x,traits) {
  pos <- grep(paste0("^",x,"$"), names(traits))
  return(pos)
}

###################### switch node test T/F ####################

switchNodeTest <- function(tree, position, traits) {
  # partialTree <- subtrees(tree)
  descendants <- getDescendants(tree, node = position)
  twoD <- descendants[0:2]
  tipsOne <- tips(tree, twoD[1])
  tipsTwo <- tips(tree, twoD[2])
  # partialTree <- extract.clade(tree, node = position)
  # get colours for tree
  coloursSubOne <- traits[sapply(tipsOne, function(x) getPositionFromTraits(x, traits))]
  coloursSubTwo <- traits[sapply(tipsTwo, function(x) getPositionFromTraits(x, traits))]
  # get most recent ancestor
  # mrcaTemp <- findMRCA(partialTree, tips = partialTree$tip.label)
  # rotTree <- rotateNodes(partialTree, node = mrcaTemp)
  # coloursSubRot <- traits[sapply(rotTree$tip.label, function(x) getPositionFromTraits(x, traits))]
  if (blockLevel(c(coloursSubOne,coloursSubTwo)) > blockLevel(c(coloursSubTwo,coloursSubOne))) {
    # the original is more fragmented than block
    return(TRUE)
  } else {
    # the same or original is better
    return(FALSE)
  }
}

###################### switch node test global ####################

switchNodeTestGlobal <- function(tree, position, traits) {
  allTips <- tree$tip.label
  # Get descendants
  descendants <- getDescendants(tree, node = position)
  tipsOne <- tips(tree, descendants[1])
  tipsOneStart <- grep(paste0("^",tipsOne[1],"$"), allTips)
  tipsTwo <- tips(tree, descendants[2])
  tipsTwoEnd <- grep(paste0("^",tipsTwo[length(tipsTwo)],"$"), allTips)
  #grep and replace
  allTipsREV <- allTips
  allTipsREV[tipsOneStart:tipsTwoEnd] <- c(tipsTwo, tipsOne)
  coloursOriginal <- traits[sapply(allTips, function(x) getPositionFromTraits(x, traits))]
  coloursReverse <- traits[sapply(allTipsREV, function(x) getPositionFromTraits(x, traits))]
  if (blockLevel(coloursOriginal) > blockLevel(coloursReverse)) {
    # the original is more fragmented than block
    return(TRUE)
  } else {
    # the same or original is better
    return(FALSE)
  }
}
###################### switch node test cluster ####################

switchNodeTestCluster <- function(tree, position, traits) {
  # partialTree <- subtrees(tree)
  # partialTree <- extract.clade(tree, node = position)
  # Alltips
  allTips <- tree$tip.label
  # Get descendants
  descendants <- getDescendants(tree, node = position)
  tipsOne <- tips(tree, descendants[1])
  tipsOneStart <- grep(paste0("^",tipsOne[1],"$"), allTips)
  tipsTwo <- tips(tree, descendants[2])
  tipsTwoEnd <- grep(paste0("^",tipsTwo[length(tipsTwo)],"$"), allTips)
  #grep and replace
  allTipsREV <- allTips
  allTipsREV[tipsOneStart:tipsTwoEnd] <- c(tipsTwo, tipsOne)
  # tipsOrderORIGINAL <- paste0(c(tipsOne,tipsTwo), collapse = ",")
  # tipsOrderREPLACE <- paste0(c(tipsTwo,tipsOne),collapse = ",")
  # allTipsCollREV <- gsub(tipsOrderORIGINAL,tipsOrderREPLACE, allTipsColl)
  # allTipsREV <- unlist(strsplit(allTipsCollREV, split = ",", fixed = TRUE))
  # get colours for tree
  coloursOriginal <- traits[sapply(allTips, function(x) getPositionFromTraits(x, traits))]
  coloursReverse <- traits[sapply(allTipsREV, function(x) getPositionFromTraits(x, traits))]
  
  # coloursSub <- traits[sapply(tree$tip.label, function(x) getPositionFromTraits(x, traits))]
  # get most recent ancestor
  # mrcaTemp <- findMRCA(tree, tips = partialTree$tip.label)
  # rotTree <- rotateNodes(tree, node = mrcaTemp)
  # coloursSubRot <- traits[sapply(rotTree$tip.label, function(x) getPositionFromTraits(x, traits))]
  if (blockLevelCluster(coloursOriginal) < blockLevelCluster(coloursReverse)) {
    # the original is more fragmented than block
    return(TRUE)
  } else {
    # the same or original is better
    return(FALSE)
  }
}

###################### reorder tree ####################

reorderTree <- function(treeOriginal = treeOriginal, traits = traits) {
  tree <- compute.brlen(treeOriginal)
  internalNodes <- (Ntip(tree)+1):(Ntip(tree)+Nnode(tree))
  for (i in rev(1:Nnode(tree))) {
    print(paste0("Testing node #", i, "/", Nnode(tree), " of nodes in tree (Local)"))
    output <- switchNodeTest(tree, internalNodes[i], traits)
    print("Switching nodes...")
    if (TRUE) {
      tree <- rotateNodes(tree, node = internalNodes[i])
    } else {
      next
    }
  }
  for (i in rev(1:Nnode(tree))) {
    print(paste0("Testing node #", i, "/", Nnode(tree), " of nodes in tree (Global)"))
    output <- switchNodeTestGlobal(tree, internalNodes[i], traits)
    print("Switching nodes...")
    if (output[[1]]) {
      tree <- rotateNodes(tree, node = internalNodes[i])
    } else {
      next
    }
  }
  for (i in 1:Nnode(tree)) {
    print(paste0("Testing node #", i, "/", Nnode(tree), " of nodes in tree (Cluster)"))
    output <- switchNodeTestCluster(tree, internalNodes[i], traits)
    print("Switching nodes...")
    if (output[[1]]) {
      tree <- rotateNodes(tree, node = internalNodes[i])
    } else {
      next
    }
  }
  return(tree)
}


# # 
# reorderTreeORIGINAL <- function(treeOriginal = treeOriginal, traits = traits) {
#   tree <- compute.brlen(treeOriginal)
#   internalNodes <- (Ntip(tree)+1):(Ntip(tree)+Nnode(tree))
#   for (i in rev(1:Nnode(tree))) {
#     print(paste0("Testing node #", i, "/", Nnode(tree), " of nodes in tree"))
#     output <- switchNodeTest(tree, i, traits)
#     print("Switching nodes...")
#     if (TRUE) {
#       tree <- rotateNodes(tree, node = internalNodes[i])
#     } else {
#       next
#     }
#   }
#   return(tree)
# }
