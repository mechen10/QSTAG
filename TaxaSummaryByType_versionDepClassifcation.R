#!/bin/bash

### This script is to make taxa summaries of different "categories" (aka salinity types)


########## OPT PARSE LATER##########


############### FOR TESTING #################

setwd('/Users/melissachen/Documents/Masters/Project_Environmental/DUMMYDATA/16sBaltic_frombotaclust/')
otuByTypePWD <- '/Users/melissachen/Documents/Masters/Project_Environmental/DUMMYDATA/16sBaltic_frombotaclust/CLASSIFICATION/TaxaPercentage/Species.txt'
otuByTaxaPWD <- '/Users/melissachen/Documents/Masters/Project_Environmental/DUMMYDATA/16sBaltic_frombotaclust/CLASSIFICATION/TypeCounts/Species.txt'
taxaIDsPWD <- '/Users/melissachen/Documents/Masters/Project_Environmental/DUMMYDATA/16sBaltic_frombotaclust/taxaIDLegend.txt'
level <- 'Species'
# collOTUTablePWD <- '/Users/melissachen/Documents/Masters/Project_Environmental/DUMMYDATA/taxa_sum_fraser16s/OTU_Table_nochpmito_nocon_COLL_L5.txt'
# otuPresAbsPWD <- '/Users/melissachen/Documents/Masters/Project_Environmental/DUMMYDATA/16sSALBIN/CLASSIFICATION/OTUTablebyType.txt'
# collOTUTableListPWD <- unlist(strsplit(collOTUTablePWD, split = ","))

############### LOAD DATA ################

otuByType <- read.delim(paste0(otuByTypePWD)
                         , header = TRUE
                         , row.names = 1
                         , stringsAsFactors = FALSE
                         , na.strings = c("","NA","na")
                        , check.names = FALSE)
otuByTaxa <- read.delim(paste0(otuByTaxaPWD)
                        , header = TRUE
                        , row.names = 1
                        , stringsAsFactors = FALSE
                        , na.strings = c("","NA","na")
                        ,check.names = FALSE)
taxaIDs <- read.delim(paste0(taxaIDsPWD)
                      , header = FALSE
                      , stringsAsFactors = FALSE)


############### START ###################
# Generate color gradient; red/blue and purple
bluered <- colorRampPalette(c("blue","white","red"))

# Collapse marine and fresh (from type to typeSimple)
# Each row (taxa) sums to 1; so each column represents % of OTUs in that group that represents Fresh/marine/brack
freshCollapsed <-rowSums(as.matrix(otuByType[rownames(otuByType),grep("^.*fresh.*$", colnames(otuByType))]))
marineCollapsed <-rowSums(as.matrix(otuByType[rownames(otuByType),grep("^.*marine.*$", colnames(otuByType))]))
brackishCollapsed <-rowSums(as.matrix(otuByType[rownames(otuByType),grep("^.*brack.*$", colnames(otuByType))]))

# Make into one matrix
collapsedTaxa <- cbind(freshCollapsed, marineCollapsed, brackishCollapsed)

# Make a gradient of red-> blue; make uneven so that 'centre' is white.
# number of steps in color gradient
ngrad <- 10
fremarCol <- bluered(2*ngrad+1)

# Make a colorVector where each taxa is judged as either fresh or marine; or brackish
colorVector <- c()
for (r in 1:nrow(collapsedTaxa)) {
  # If there are more brackish reps in taxa group, then it is called 'brackish'
  if (collapsedTaxa[r,"brackishCollapsed"] > collapsedTaxa[r,"freshCollapsed"] & collapsedTaxa[r,"brackishCollapsed"] > collapsedTaxa[r,"marineCollapsed"]) {
    colTemp <- "purple"
  } else if ((collapsedTaxa[r,"freshCollapsed"] == 0) & (collapsedTaxa[r,"marineCollapsed"] == 0)) {
    colTemp <- "black"
  } else if (sum(collapsedTaxa[r,]) < 0.10) {
    colTemp <- "black"
  } else {
    # For everything else, color it by the proportion of marine and fresh representatives
    # Need to adjust for cases of '0'--> make 0.001
    ratio <- (collapsedTaxa[r,"marineCollapsed"]+0.0001)/(collapsedTaxa[r,"freshCollapsed"]+0.0001)
    # Maximum is 2 because color position in fremarCol must range from 0-2
    ratio <- min(2,ratio)
    # multiple ratio by ngrad and round to nearest whole number
    colTemp <- fremarCol[round(ratio*ngrad)+1]
  }
  colorVector <- c(colorVector, colTemp)
}
names(colorVector) <- rownames(collapsedTaxa)

# Collapse otusByTaxa so that it is 'fresh', 'marine', or 'brack'
# Done by COUNTS because summing percentages will not make sense. 
freshTaxaSummary <- colSums(otuByTaxa[grep("fresh", rownames(otuByTaxa)),])
marineTaxaSummary <- colSums(otuByTaxa[grep("marine", rownames(otuByTaxa)),])
brackTaxaSummary <- colSums(otuByTaxa[grep("brack", rownames(otuByTaxa)),])

# Combine into single matrix
allTaxaSummarybyType <- cbind(freshTaxaSummary, brackTaxaSummary, marineTaxaSummary)
toDelete <- grep("__$|Unassigned|unassigned", rownames(allTaxaSummarybyType))
if (length(toDelete) > 0) {
  allTaxaSummarybyType <- allTaxaSummarybyType[-toDelete,]
}

# Make allTaxaSummarybyType same order as colorVector
colorVectorOrdered <- colorVector[match(rownames(allTaxaSummarybyType), names(colorVector))]



pdf(paste0(level,"_barplot_representatives.pdf"), pointsize = 14)
par(fig = c(0,0.7,0,1))
barplot(prop.table(allTaxaSummarybyType, margin = 2)
        , col = colorVectorOrdered
        , border = NA
        , )
par(fig = c(0.7,1,0,1), mar = c(5.1,0,4.1,2.1))
dev.off()

# Figure out what brackish species are

BrackishSpecies <- names(colorVectorOrdered)[which(colorVectorOrdered == "purple")]

capture.output(BrackishSpecies, file = paste0(level,"_BrackishSpecies.txt"))
