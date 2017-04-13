#! /usr/bin/Rscript
########################## OPTPARSE ####################################

library(optparse)

option_list = list(
  make_option(c("-m", "--metadataShareOTUs"), type="character",
              help="metadataShareOTUsoutput"),
  make_option(c("-N", "--gradientNames"), type="character",
              help="gradient file, must be same between sets"),
  make_option(c("-d", "--dataNames"), type="character",
              help="Comma separated list of file names; must match those in metadataShareOTUs"),
  make_option(c("-t", "--tablesList"),
              help="comma separated OTU table files", type="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

metadataShareOTUsPWD = opt$metadataShareOTUs
gradientNamesPWD = opt$gradientNames
dataNames = opt$dataNames
tablesListPWD = opt$tables

########################## LOAD FILES ####################################

# Load files
print("Loading files for graphing...")

metadataShareOTUs <- read.delim(paste0(metadataShareOTUsPWD)
                                , header = TRUE
                                , row.names = 1
                                , stringsAsFactors = FALSE
                                )

gradientNames <- read.delim(paste0(gradientNamesPWD), header = FALSE, strip.white = TRUE, sep = ",", stringsAsFactors = FALSE)
gradientNames <- unlist(strsplit(paste0(gradientNames), split = ",", fixed = TRUE))

dataNames <- strsplit(dataNames, split = ",", fixed = TRUE)


tablesListPWD <- unlist(strsplit(tablesListPWD, split = ",", fixed = TRUE))



firstSet <- read.delim(paste0(tablesListPWD[1])
                       , skip = 1
                       , header = TRUE
                       , stringsAsFactors = FALSE
                       , row.names = 1)
firstSet <- firstSet[-ncol(firstSet)]
secondSet <- read.delim(paste0(tablesListPWD[2])
                       , skip = 1
                       , header = TRUE
                       , stringsAsFactors = FALSE
                       , row.names = 1)
secondSet <- secondSet[-ncol(secondSet)]

########################## FUNCTIONS ####################################

getColor <- function(x) {
  if (x == as.character(gradientNames[1])) {
    return("blue")
  } else if (x == as.character(gradientNames[2])) {
    return("purple")
  } else if (x == as.character(gradientNames[3])) {
    return("red")
  } else if (x == paste0("hi",gradientNames[2])) {
    return("pink")
  } else if (x == paste0("lo", gradientNames[2])) {
    return("magenta")
  } else if (x == paste0("half", gradientNames[1])) {
    return("cyan")
  } else if (x == paste0("half", gradientNames[2])) {
    return("pink4")
  } else if (x == paste0("half", gradientNames[3])) {
    return("orange")
  } else if (x == paste0("bloom", gradientNames[1])) {
    return("darkblue")
  } else if (x == paste0("bloom", gradientNames[2])) {
    return("purple4")
  } else if (x == paste0("bloom", gradientNames[3])) {
    return("darkred")
  } else {
    return("NA")
  }
}

makeRel <- function(x) {
  headerNames <- colnames(x)
  rowNames <- rownames(x)
  newMatrix <- matrix()
  first <- TRUE
  for (sampleCol in 1:ncol(x)) {
    total <- sum(x[sampleCol])
    tempCol <- x[sampleCol]/total
    if (first) {
      newMatrix <- tempCol
      first <- FALSE
    } else {
      newMatrix <- cbind(newMatrix, tempCol)
    }
  }
  colnames(newMatrix) <- headerNames
  rownames(newMatrix) <- rowNames
  return(newMatrix)
}

########################## START ####################################

colorsListFirst <- c()
for (i in 1:nrow(firstSet)) {
  if (rownames(firstSet)[i] %in% rownames(metadataShareOTUs)) {
    typeTemp <- metadataShareOTUs[,1][i]
    if (is.na(typeTemp)) {
      typeTemp <- 'noclass'
    }
    colorsListFirst <- c(colorsListFirst, getColor(typeTemp))
  } else {
    colorsListFirst <- c(colorsListFirst, "white")
  }
}

colorsListSecond <- c()
for (i in 1:nrow(secondSet)) {
  if (rownames(secondSet)[i] %in% rownames(metadataShareOTUs)) {
    typeTemp <- metadataShareOTUs[,2][i]
    if (is.na(typeTemp)) {
      typeTemp <- 'noclass'
    }
    colorsListSecond <- c(colorsListSecond, getColor(typeTemp))
  } else {
    colorsListSecond <- c(colorsListSecond, "white")
  }
}


firstSetRel <- makeRel(firstSet)
secondSetRel <- makeRel(secondSet)

jpeg("FirstSetShared.jpeg")
barplot(as.matrix(firstSetRel)
        , col = colorsListFirst
        , border = NA)
dev.off()

jpeg("SecondSetShared.jpeg")
barplot(as.matrix(secondSetRel)
        , col = colorsListSecond
        , border = NA)
dev.off()
