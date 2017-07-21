#! /usr/bin/Rscript
# Gradient Binning Output- graphing
# setwd("~")
# setwd("/Volumes/MCDataDrive/Users/melissachen/Documents/Masters_data/Undergrad_Environmental/Roughwork/ALLDATASETS_7march2017_16s/2_analysis_ALL/SALBIN_16s_FRASER_13march2017/")
########################## OPTPARSE ####################################
library(optparse)

option_list = list(
  make_option(c("-b", "--boundaries"), type="character",
              help="Boundaries.txt output from binning script", default = "boundaries.txt"),
  make_option(c("-M", "--model_boundaries"),
              help="Model Boundaries output from binning script", type="character", default = 'modelBoundaries_type.txt'),
  make_option(c("-t", "--taxa_abund"), 
              help="Taxa abundances across gradient file output from binning script", type="character", default = 'taxa_abundances_across_gradient.txt'),
  make_option(c("-T", "--types"),
              help="Types across gradient output file from binning script", type="character", default = "types_across_gradient"),
  make_option(c("-N","--Namesgradient"), 
              help= "Gradient text output from binning script", type = "character", default = "gradient.txt"),
  make_option(c("-a","--allBins"), 
              help= "If set as true, will graph 'all' of bins. Default = TRUE", action = "store_false", default = "TRUE"),
  make_option(c("-c","--condensedBins"), 
              help= "If set as true, will graph condensed version of bins. Default = TRUE", action = "store_false", default = "TRUE"),
  make_option(c("-B", "--biom"), type="character",
              help="biom-- already collapsed, but in biom format.", metavar="character", default = "OTUTableText.txt"),
  make_option(c("-n", "--minthreshold"), type="character",
              help="Minimum relabundance to be plotted", metavar="numeric", default = 0.05),
  make_option(c("-O", "--OTUTable"), type = "character",
              help="Type OTU Table output from Classification script", metavar = "character", default = "OTUTablebyType.txt"),
  make_option(c("-L","--LOG"), 
              help= "LOG output from binning script", type = "character", default = "LOG.txt")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

boundariesPWD = opt$boundaries
modelboundariesPWD = opt$model_boundaries
taxaAbundPWD = opt$taxa_abund
typesPWD = opt$types
gradientNamesPWD = opt$Namesgradient
allBins = opt$allBins
condensedBins = opt$condensedBins
biomPWD = opt$biom
minthreshold = opt$minthreshold
LOGPWD = opt$LOG
typetablePWD = opt$OTUTable


########################## For Testing ####################################
# setwd('/Users/melissachen/Documents/Masters/Project_Environmental/DUMMYDATA/TESTING')
# wd <- getwd()
# boundariesPWD = unlist(paste0(wd, "/boundaries.txt"))
# modelboundariesPWD = paste0(wd, "/modelBoundaries_type.txt")
# taxaAbundPWD = paste0(wd, "/taxa_abundances_across_gradient.txt")
# typesPWD = paste0(wd, "/types_across_gradient")
# gradientNamesPWD = paste0(wd, "/gradient.txt")
# allBins = TRUE
# condensedBins = TRUE
# biomPWD = paste0(wd, "/OTUTableText.txt")
# minthreshold = 0.05
# LOGPWD = paste0(wd, "/LOG.txt")
# typetablePWD = paste0(wd, "/OTUTablebyType.txt")


########################## LOAD FILES ####################################

# Load files
print("Loading files for graphing...")

boundaries <- read.delim(paste0(boundariesPWD), header = FALSE, strip.white = TRUE)
histBoundaries <- as.numeric(boundaries$V1)

logvalues <- read.delim(paste0(LOGPWD), header = FALSE, strip.white = TRUE)
minGrad <- as.character(logvalues[grep("minVal", logvalues$V1),])
minGrad <- floor(as.numeric(gsub("minVal =","", minGrad)))
maxGrad <- as.character(logvalues[grep("maxVal", logvalues$V1),])
maxGrad <- ceiling(as.numeric(gsub("maxVal =","", maxGrad)))
treePWD <- as.character(logvalues[grep("tree", logvalues$V1),])
treePWD <- gsub("tree = ","",treePWD)
metadataPWD <- as.character(logvalues[grep("metadata", logvalues$V1),])
metadataPWD <- gsub("metadata = ","",metadataPWD)
threshold <- as.character(logvalues[grep("threshold", logvalues$V1),])
threshold <- as.numeric(gsub("threshold =", "", threshold))

modelBoundaries <- read.delim(paste0(modelboundariesPWD), header = TRUE, row.names = 1, strip.white = TRUE, stringsAsFactors = FALSE, na.strings = c('','na','NA'))

taxaAbundances <- read.delim(paste0(taxaAbundPWD), header = FALSE, row.names = 1, strip.white = TRUE, stringsAsFactors = FALSE, na.strings = c('','na','NA'))
gradient <- as.numeric(taxaAbundances[1,])
taxaAbundances <- taxaAbundances[-1,]

typesAll <- read.delim(paste0(typesPWD,"_all.txt"), header = TRUE, strip.white = TRUE, stringsAsFactors = FALSE, row.names = 1, na.strings = c('','na','NA'))
colnames(typesAll) <- gsub("X", "", colnames(typesAll))
typesCondensed <- read.delim(paste0(typesPWD,"_condensed.txt"), header = TRUE, strip.white = TRUE, stringsAsFactors = FALSE, row.names = 1, na.strings = c('','na','NA'))
colnames(typesCondensed) <- gsub("X", "", colnames(typesCondensed))

gradientNames <- read.delim(paste0(gradientNamesPWD), header = FALSE, strip.white = TRUE, sep = ",", stringsAsFactors = FALSE)


########################## HISTOGRAM (Single) ####################################

# Plot histogram
print("Plotting single histogram")
meanBoundaries <- mean(histBoundaries)
total <- length(histBoundaries)
# expProportion <- dpois(seq(minGrad,maxGrad,2),meanBoundaries)
# expFrequency <- expProportion*total*2

pdf(file = "BoundaryDistributions.pdf")
hist(histBoundaries
     , main = "Boundary distribution"
     , xlab = "Gradient"
     # , breaks = seq(minGrad,maxGrad,2)
     , col = 'grey'
     , xlim = c(minGrad,maxGrad))
# lines(expFrequency ~ seq(0,36,2), lwd = 2, col = "red")
abline(v = meanBoundaries
       , col= "red"
       , lwd = 2
       )
dev.off()

########################## TAXA ABUND (Individual) ############

# Taxa abundances for each taxa
print("Plotting individual taxa abund")
system("mkdir INDIV_PLOTS")

for (x in (1:nrow(taxaAbundances))) {
  taxa <- rownames(taxaAbundances)[x]
  Avalue <- modelBoundaries[taxa,'A']
  Bvalue <- modelBoundaries[taxa,'B']
  testB <- modelBoundaries[taxa,'sigAB']
  meanA <- modelBoundaries[taxa,'meanA']
  meanB <- modelBoundaries[taxa,'meanB']
  meanC <- modelBoundaries[taxa,'meanC']
  typeMB <- modelBoundaries[taxa,'type']
  bloom <- modelBoundaries[taxa, 'bloom']
  pdf(file = paste0("INDIV_PLOTS/",typeMB,'_Indiv',taxa,'.pdf'))
  plot(as.numeric(taxaAbundances[x,]) ~gradient
       , xlab = 'Gradient'
       , ylab = 'Relative Abundance'
       , main = paste(taxa)
       , sub = paste(typeMB,bloom,signif(modelBoundaries[taxa,'sigAB'],3)," ",signif(modelBoundaries[taxa,'sigBC'],3), " ",signif(modelBoundaries[taxa,'sigAC'],1)))
  if (testB != 'None') {
  lines(c(minGrad,Avalue), c(meanA,meanA)
        , lwd = 3
        , col = 'blue')
  lines(c(Avalue,Bvalue), c(meanB,meanB)
        , lwd = 3
        , col = 'purple')
  lines(c(Bvalue,maxGrad), c(meanC,meanC)
        , lwd = 3
        , col = 'red')
  } else if (testB == 'None') {
    lines(c(minGrad,Avalue), c(meanA,meanA)
          , lwd = 3
          , col = 'blue')
    lines(c(Avalue,Bvalue), c(meanA,meanC)
          , lwd = 3
          , col = 'purple')
    lines(c(Bvalue,maxGrad), c(meanC,meanC)
          , lwd = 3
          , col = 'red')
        }

  dev.off()
}


########################## BARPLOT: Taxa types across gradient ####################################

# Taxa types across gradient
print("Plotting barplots for taxa types")

colorsGradient <- c("blue","darkblue","lightblue","cyan","purple","green","magenta","pink","orange","red","darkred","black", "grey","white")
colorsGradientCondensed <- c("blue","purple","red","black")


if (allBins == TRUE) {
  typesAll <- as.matrix(typesAll)
  pdf(file = "taxatypesacrossgradient_all.pdf")
  par(mar=c(5.1,4.1,4.1,8.1), xpd = TRUE)
  barplot(typesAll, beside = FALSE
          , col = colorsGradient
          , main = 'Taxa types across gradient'
          , xlab = 'Gradient'
          , ylab = 'Relative Abundance')
  legend('right', inset = c(-0.35,0), legend = rev(rownames(typesAll)), pch = 19, col = rev(colorsGradient))
  dev.off()
}

if (condensedBins == TRUE) {
  typesCondensed <- as.matrix(typesCondensed)
  pdf(file = "taxatypesacrossgradient_condensed.pdf")
  par(mar=c(5.1,4.1,4.1,8.1), xpd = TRUE)
  barplot(typesCondensed, beside = FALSE
          , col = colorsGradientCondensed
          , main = 'Taxa types across gradient'
          , xlab = 'Gradient'
          , ylab = 'Relative Abundance')
  legend('right', inset = c(-0.35,0), legend = rev(rownames(typesCondensed)), pch = 19, col = rev(colorsGradientCondensed))
  dev.off()
}
  
  ########################## TAXA ABUND (by type) ####################################
system('mkdir TaxaAbund_bytype_all')
system('mkdir TaxaAbund_bytype_condensed')
print("Making taxa abundance plots by type")

if (allBins == TRUE) {
  for (i in 1:nrow(typesAll)){
    pdf(file = paste0('./TaxaAbund_bytype_all/Indiv_', rownames(typesAll)[i], '_all.pdf'))
    plot(typesAll[i,] ~ colnames(typesAll)
         , type = 'h'
         , lwd = 10
         , col = 'darkgrey'
         , main = paste0(rownames(typesAll)[i],' abundance across gradient')
         , xlab = 'Gradient'
         , ylab = 'Relative abundance')
    dev.off()
  }
}

if (condensedBins == TRUE) {
  for (i in 1:nrow(typesCondensed)){
    pdf(file = paste0('./TaxaAbund_bytype_condensed/Indiv_', rownames(typesCondensed)[i], '.pdf'))
    plot(typesCondensed[i,] ~ colnames(typesCondensed)
         , type = 'h'
         , lwd = 10
         , col = 'darkgrey'
         , main = paste0(rownames(typesCondensed)[i],' abundance across gradient')
         , xlab = 'Gradient'
         , ylab = 'Relative abundance')
    dev.off()
  }
}
# 
# # Make one with fresh and marine combined to see if it dips in middle
# typesCondensedCombo <- typesCondensed[paste0(gradientNames[1]),] + typesCondensed[paste0(gradientNames[3]),]
# pdf(file = paste0('./TaxaAbund_bytype_condensed/Indiv_', gradientNames[1],gradientNames[3], '.pdf'))
# plot(unlist(typesCondensedCombo) ~ colnames(typesCondensed)
#      , type = 'h'
#      , lwd = 10
#      , col = 'darkgrey'
#      , main = paste0(gradientNames[1],gradientNames[3],' abundance across gradient')
#      , xlab = 'Gradient'
#      , ylab = 'Relative abundance')
# dev.off()

########################## HISTOGRAM (by type) ####################################
# Histogram but separated by 'type'
print ("Making histograms, by type")

lowOnly <- modelBoundaries[grep(paste0(gradientNames[1]), modelBoundaries$typeSimple),]
lowOnlyTrue <- modelBoundaries[grep(paste0("^",gradientNames[1],"Restricted$"), modelBoundaries$type),]
lowOnlyHalf <- modelBoundaries[grep(paste0(gradientNames[1],"Peak"), modelBoundaries$type),]
lowOnlyBloom <- modelBoundaries[grep(paste0(gradientNames[1],"Bloom"), modelBoundaries$type),]


highOnly <- modelBoundaries[grep(paste0(gradientNames[3]), modelBoundaries$typeSimple),]
highOnlyTrue <- modelBoundaries[grep(paste0("^",gradientNames[3],"Restricted$"), modelBoundaries$type),]
highOnlyHalf <- modelBoundaries[grep(paste0(gradientNames[3], "Peak"), modelBoundaries$type),]
highOnlyBloom <- modelBoundaries[grep(paste0(gradientNames[3], "Bloom"), modelBoundaries$type),]

interOnly <- modelBoundaries[grep(paste0(gradientNames[2]), modelBoundaries$typeSimple),]
interOnlyTrue <- modelBoundaries[grep(paste0("^",gradientNames[2],"Restricted$"), modelBoundaries$type),]
interOnlyHalf <- modelBoundaries[grep(paste0(gradientNames[2], "Peak"), modelBoundaries$type),]
interOnlyBloom <- modelBoundaries[grep(paste0(gradientNames[2], "Bloom"), modelBoundaries$type),]
interOnlyHi <- modelBoundaries[grep(paste0(gradientNames[2], "PeakHiToler"), modelBoundaries$type),]
interOnlylo <- modelBoundaries[grep(paste0(gradientNames[2], "PeakLoToler"), modelBoundaries$type),]

# Set even bin widths 
is.even <- function(x) { x%%2 == 0}
if (is.even(maxGrad-minGrad)){
  binBreaks <- seq(minGrad,maxGrad, by = 2)
} else {
  binBreaks <- seq(minGrad,maxGrad+1, by = 2)
}

system("mkdir Histograms")

########################## Low ####################################

allLow <- c(lowOnly$boundaries,lowOnly$boundaries)
TrueLow <- c(lowOnlyTrue$boundaries)
HalfLow <- c(lowOnlyHalf$boundaries)
BloomLow <- c(lowOnlyBloom$boundaries)
nonBloomLow <- c(HalfLow, TrueLow)



if (condensedBins == TRUE) {
  lowBoundaries <- c()
  if (length(allLow) > 0){
    for (i in 1:length(allLow)) {
      if (!is.na(allLow[i])) {
        lowBoundaries <- c(lowBoundaries, allLow[i])
      } 
    }
    meanLow <- round(mean(lowBoundaries),2)
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[1],"Boundaries_all.pdf"))
    hist(as.numeric(lowBoundaries)
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , breaks = binBreaks
         , col = "grey"
         , main = paste0("Boundary distribution for ALL ", gradientNames[1],": u=",as.character(meanLow)))
    abline(v=meanLow
           , col = "red"
           , lwd = 2)
    dev.off()
  }

    lowBoundaries <- c()
  if (length(nonBloomLow) > 0){
    for (i in 1:length(nonBloomLow)) {
      if (!is.na(nonBloomLow[i])) {
        lowBoundaries <- c(lowBoundaries, nonBloomLow[i])
      } 
    }
    meanLow <- round(mean(lowBoundaries),2)
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[1],"Boundaries_nonBloom.pdf"))
    hist(lowBoundaries
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , col = "grey"
         , breaks = binBreaks
         , main = paste0("Boundary distribution for nonBLOOM ", gradientNames[1],": u=",as.character(meanLow)))
    abline(v=meanLow
           , col = "red"
           , lwd = 2)
    dev.off()
  }
}

if (allBins == TRUE) {
  # True
  lowBoundaries <- c()
  if (length(TrueLow) > 0){
    for (i in 1:length(TrueLow)) {
      if (!is.na(TrueLow[i])) {
        lowBoundaries <- c(lowBoundaries, TrueLow[i])
      } 
    }
    meanLow <- round(mean(lowBoundaries),2)
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[1],"Boundaries_True.pdf"))
    hist(lowBoundaries
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , breaks = binBreaks
         , col = "grey"
         , main = paste0("Boundary distribution for TRUE ", gradientNames[1],": u=",as.character(meanLow)))
    abline(v=meanLow
           , col = "red"
           , lwd = 2)
    dev.off()
  }

    # Half
  if (length(HalfLow) > 0){
    for (i in 1:length(HalfLow)) {
      if (!is.na(HalfLow[i])) {
        lowBoundaries <- c(lowBoundaries, HalfLow[i])
      } 
    }
    meanLow <- round(mean(lowBoundaries),2)
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[1],"Boundaries_Half.pdf"))
    hist(as.numeric(lowBoundaries)
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , col = "grey"
         , breaks = binBreaks
         , main = paste0("Boundary distribution for HALF ", gradientNames[1],": u=",as.character(meanLow)))
    abline(v=meanLow
           , col = "red"
           , lwd = 2)
    dev.off()
  }
  
  # Bloom
  lowBoundaries <- c()
  if (length(BloomLow) > 0){
    for (i in 1:length(BloomLow)) {
      if (!is.na(BloomLow[i])) {
        lowBoundaries <- c(lowBoundaries, BloomLow[i])
      } 
    }
    meanLow <- round(mean(lowBoundaries),2)
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[1],"Boundaries_Bloom.pdf"))
    hist(lowBoundaries
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , breaks = binBreaks
         , col = "grey"
         , main = paste0("Boundary distribution for BLOOM ", gradientNames[1],": u=",as.character(meanLow)))
    abline(v=meanLow
           , col = "red"
           , lwd = 2)
    dev.off()
  }
}

########################## High ####################################

# Now pull out all boundaries for fresh
allHigh <- c(highOnly$boundaries, highOnly$boundariestwo)
TrueHigh <- c(highOnlyTrue$boundaries)
HalfHigh <- c(highOnlyHalf$boundaries)
BloomHigh <- c(highOnlyBloom$boundaries)
nonBloomHigh <- c(HalfHigh,TrueHigh)

if (condensedBins == TRUE) {
  HighBoundaries <- c()
  if (length(allHigh) > 0){
    for (i in 1:length(allHigh)) {
      if (!is.na(allHigh[i])) {
        HighBoundaries <- c(HighBoundaries, allHigh[i])
      } 
    }
    meanHigh <- round(mean(HighBoundaries),2)
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[3],"Boundaries_all.pdf"))
    hist(HighBoundaries
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , breaks = binBreaks
         , col = "grey"
         , main = paste0("Boundary distribution for ALL ", gradientNames[3],": u=",as.character(meanHigh)))
    abline(v=meanHigh
           , col = "red"
           , lwd = 2)
    dev.off()
  }
  HighBoundaries <- c()
  if (length(nonBloomHigh) > 0){
    for (i in 1:length(nonBloomHigh)) {
      if (!is.na(nonBloomHigh[i])) {
        HighBoundaries <- c(HighBoundaries, nonBloomHigh[i])
      } 
    }
    meanHigh <- round(mean(HighBoundaries),2)
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[3],"Boundaries_nonBloom.pdf"))
    hist(HighBoundaries
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , breaks = binBreaks
         , col = "grey"
         , main = paste0("Boundary distribution for nonBlOOM ", gradientNames[3],": u=",as.character(meanHigh)))
    abline(v=meanHigh
           , col = "red"
           , lwd = 2)
    dev.off()
  }
}

if (allBins == TRUE) {
  # True
  HighBoundaries <- c()
  if (length(TrueHigh) > 0){
    for (i in 1:length(TrueHigh)) {
      if (!is.na(TrueHigh[i])) {
        HighBoundaries <- c(HighBoundaries, TrueHigh[i])
      } 
    }
    meanHigh <- round(mean(HighBoundaries),2)
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[3],"Boundaries_True.pdf"))
    hist(HighBoundaries
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , breaks = binBreaks
         , col = "grey"
         , main = paste0("Boundary distribution for TRUE ", gradientNames[3],": u=",as.character(meanHigh)))
    abline(v=meanHigh
           , col = "red"
           , lwd = 2)
    dev.off()
  }
  
  # Half
  HighBoundaries <- c()
  if (length(HalfHigh) > 0) {
    for (i in 1:length(HalfHigh)) {
      if (!is.na(HalfHigh[i])) {
        HighBoundaries <- c(HighBoundaries, HalfHigh[i])
      } 
    }
    meanHigh <- round(mean(HighBoundaries),2)
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[3],"Boundaries_Half.pdf"))
    hist(HighBoundaries
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , breaks = binBreaks
         , col = "grey"
         , main = paste0("Boundary distribution for HALF ", gradientNames[3],": u=",as.character(meanHigh)))
    abline(v=meanHigh
           , col = "red"
           , lwd = 2)
    dev.off()
  }
  
  # Bloom
  HighBoundaries <- c()
  if (length(BloomHigh) > 0) {
    for (i in 1:length(BloomHigh)) {
      if (!is.na(BloomHigh[i])) {
        HighBoundaries <- c(HighBoundaries, BloomHigh[i])
      } 
    }
    meanHigh <- round(mean(HighBoundaries),2)
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[3],"Boundaries_Bloom.pdf"))
    hist(HighBoundaries
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , breaks = binBreaks
         , col = "grey"
         , main = paste0("Boundary distribution for BLOOM ", gradientNames[3],": u=",as.character(meanHigh)))
    abline(v=meanHigh
           , col = "red"
           , lwd = 2)
    dev.off()
  }
}

# 
# ########################## Low and High ####################################
# 
# # Truelow and Truehigh combined
# combinedLH <- c(TrueLow, TrueHigh)
# LHBoundaries <- c()
# if (length(combinedLH) > 0) {
#   for (i in 1:length(combinedLH)) {
#     if (!is.na(combinedLH[i])) {
#       LHBoundaries <- c(LHBoundaries, combinedLH[i])
#     } 
#   }
#   meancomb <- round(mean(LHBoundaries),2)
#   jpeg(file = paste0("./Histograms/Histogram_",gradientNames[1],gradientNames[3],"Boundaries_True.jpeg"))
#   hist(LHBoundaries
#        , xlim = range(minGrad,maxGrad)
#        , xlab = "Gradient"
#        , breaks = binBreaks
#        , col = "grey"
#        , main = paste0("Boundary distribution for True COMBINED ", gradientNames[1],gradientNames[3],": u=",as.character(meancomb)))
#   abline(v=meancomb
#          , col = "red"
#          , lwd = 2)
#   dev.off()
#   
# }
########################## Inter ####################################

# Now pull out all boundaries for fresh
allInter <- cbind(interOnly$boundaries, interOnly$boundariestwo)
TrueInter <- cbind(interOnlyTrue$boundaries,interOnlyTrue$boundariestwo)
HalfInter <- cbind(interOnlyHalf$boundaries,interOnlyHalf$boundariestwo)
BloomInter <- cbind(interOnlyBloom$boundaries, interOnlyBloom$boundariestwo)
HiInter <- cbind(interOnlyHi$boundaries, interOnlyHi$boundariestwo)
LoInter <- cbind(interOnlylo$boundaries, interOnlylo$boundariestwo)
nonBloomInter <- rbind(HalfInter,TrueInter)

if (condensedBins == TRUE) {
  InterBoundaries <- c()
  if (length(allInter) > 0){
    InterBoundaries <- allInter
    meanInterOne <- round(mean(InterBoundaries[,1]),2)
    meanInterTwo <- round(mean(InterBoundaries[,2]),2)
    maximumY <- max(max(table(InterBoundaries[,1])),max(table(InterBoundaries[,2])))*1.5
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[2],"Boundaries_all.pdf"))
    hist(InterBoundaries[,1]
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , breaks = binBreaks
         , ylim = c(0,maximumY)
         , col = rgb(0,0,1,0.5)
         , main = paste0("Boundary distribution for ALL ", gradientNames[2],": u=",as.character(meanInterOne),",",as.character(meanInterTwo)))
    abline(v=meanInterOne
           , col = "blue"
           , lwd = 2)
    hist(InterBoundaries[,2]
         , xlim = range(minGrad,maxGrad)
         , breaks = binBreaks
         , ylim = c(0,maximumY)
         , col = rgb(1,0,0,0.5)
         , add = TRUE)
    abline(v=meanInterTwo
           , col = "red"
           , lwd = 2)
    dev.off()
  }
  InterBoundaries <- c()
  if (length(nonBloomInter) > 0){
    for (i in 1:length(nonBloomInter)) {
    InterBoundaries <- nonBloomInter
    meanInterOne <- round(mean(InterBoundaries[,1]),2)
    meanInterTwo <- round(mean(InterBoundaries[,2]),2)
    maximumY <- max(max(table(InterBoundaries[,1])),max(table(InterBoundaries[,2])))*1.5
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[2],"Boundaries_nonBloom.pdf"))
    hist(InterBoundaries[,1]
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , ylim = c(0,maximumY)
         , col = rgb(0,0,1,0.5)
         , breaks = binBreaks
         , main = paste0("Boundary distribution for nonBLOOM ", gradientNames[2],": u=",as.character(meanInterOne),",",as.character(meanInterTwo)))
    abline(v=meanInterOne
           , col = "blue"
           , lwd = 2)
    hist(InterBoundaries[,2]
         , xlim = range(minGrad,maxGrad)
         , ylim = c(0,maximumY)
         , breaks = binBreaks
         , col = rgb(1,0,0,0.5)
         , add = TRUE)
    abline(v=meanInterTwo
           , col = "red"
           , lwd = 2)
    dev.off()
  }
  }
}

if (allBins == TRUE) {
  # True
  InterBoundaries <- c()
  if (length(TrueInter) > 0){
    InterBoundaries <- TrueInter
    meanInterOne <- round(mean(InterBoundaries[,1]),2)
    meanInterTwo <- round(mean(InterBoundaries[,2]),2)
    maximumY <- max(max(table(InterBoundaries[,1])),max(table(InterBoundaries[,2])))*1.5
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[2],"Boundaries_True.pdf"))
    hist(InterBoundaries[,1]
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , ylim = c(0,maximumY)
         , breaks = binBreaks
         , col = rgb(0,0,1,0.5)
         , main = paste0("Boundary distribution for TRUE ", gradientNames[2],": u=",as.character(meanInterOne),",",as.character(meanInterTwo)))
    abline(v=meanInterOne
           , col = "blue"
           , lwd = 2)
    hist(InterBoundaries[,2]
         , xlim = range(minGrad,maxGrad)
         , ylim = c(0,maximumY)
         , breaks = binBreaks
         , col = rgb(1,0,0,0.5)
         , add = TRUE)
    abline(v=meanInterTwo
           , col = "red"
           , lwd = 2)
    dev.off()
  }
  
  # Half
  InterBoundaries <- c()
  if (length(HalfInter) > 0){
    InterBoundaries <- HalfInter
    meanInterOne <- round(mean(InterBoundaries[,1]),2)
    meanInterTwo <- round(mean(InterBoundaries[,2]),2)
    maximumY <- max(max(table(InterBoundaries[,1])),max(table(InterBoundaries[,2])))*1.5
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[2],"Boundaries_Half.pdf"))
    hist(InterBoundaries[,1]
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , ylim = c(0,maximumY)
         , breaks = binBreaks
         , col = rgb(0,0,1,0.5)
         , main = paste0("Boundary distribution for Half ", gradientNames[2],": u=",as.character(meanInterOne),",",as.character(meanInterTwo)))
    abline(v=meanInterOne
           , col = "blue"
           , lwd = 2)
    hist(InterBoundaries[,2]
         , xlim = range(minGrad,maxGrad)
         , ylim = c(0,maximumY)
         , breaks = binBreaks
         , col = rgb(1,0,0,0.5)
         , add = TRUE)
    abline(v=meanInterTwo
           , col = "red"
           , lwd = 2)
    dev.off()
  }
  
  # Bloom
  InterBoundaries <- c()
  if (length(BloomInter) > 0){
    InterBoundaries <- BloomInter
    meanInterOne <- round(mean(InterBoundaries[,1]),2)
    meanInterTwo <- round(mean(InterBoundaries[,2]),2)
    maximumY <- max(max(table(InterBoundaries[,1])),max(table(InterBoundaries[,2])))*1.5
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[2],"Boundaries_Bloom.pdf"))
    hist(InterBoundaries[,1]
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , ylim = c(0,maximumY)
         , breaks = binBreaks
         , col = rgb(0,0,1,1)
         , main = paste0("Boundary distribution for BLOOM ", gradientNames[2],": u=",as.character(meanInterOne),",",as.character(meanInterTwo)))
    abline(v=meanInterOne
           , col = "blue"
           , lwd = 2)
    hist(InterBoundaries[,2]
         , xlim = range(minGrad,maxGrad)
         , ylim = c(0,maximumY)
         , breaks = binBreaks
         , col = rgb(1,0,0,0.5)
         , add = TRUE)
    abline(v=meanInterTwo
           , col = "red"
           , lwd = 2)
    dev.off()
  }
  
  # Hi
  InterBoundaries <- c()
  if (length(HiInter) > 0){
    InterBoundaries <- HiInter
    meanInterOne <- round(mean(InterBoundaries[,1]),2)
    meanInterTwo <- round(mean(InterBoundaries[,2]),2)
    maximumY <- max(max(table(InterBoundaries[,1])),max(table(InterBoundaries[,2])))*1.5
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[2],"Boundaries_Hi.pdf"))
    hist(InterBoundaries[,1]
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , ylim = c(0,maximumY)
         , breaks = binBreaks
         , col = rgb(0,0,1,0.5)
         , main = paste0("Boundary distribution for HI ", gradientNames[2],": u=",as.character(meanInterOne),",",as.character(meanInterTwo)))
    abline(v=meanInterOne
           , col = "blue"
           , lwd = 2)
    hist(InterBoundaries[,2]
         , xlim = range(minGrad,maxGrad)
         , ylim = c(0,maximumY)
         , breaks = binBreaks
         , col = rgb(1,0,0,0.5)
         , add = TRUE)
    abline(v=meanInterTwo
           , col = "red"
           , lwd = 2)
    dev.off()
  }
  
  # Lo
  InterBoundaries <- c()
  if (length(LoInter) > 0){
    InterBoundaries <- LoInter
    meanInterOne <- round(mean(InterBoundaries[,1]),2)
    meanInterTwo <- round(mean(InterBoundaries[,2]),2)
    maximumY <- max(max(table(InterBoundaries[,1])),max(table(InterBoundaries[,2])))*1.5
    pdf(file = paste0("./Histograms/Histogram_",gradientNames[2],"Boundaries_Lo.pdf"))
    hist(InterBoundaries[,1]
         , xlim = range(minGrad,maxGrad)
         , xlab = "Gradient"
         , ylim = c(0,maximumY)
         , breaks = binBreaks
         , col = rgb(0,0,1,0.5)
         , main = paste0("Boundary distribution for LO ", gradientNames[2],": u=",as.character(meanInterOne),",",as.character(meanInterTwo)))
    abline(v=meanInterOne
           , col = "blue"
           , lwd = 2)
    hist(InterBoundaries[,2]
         , xlim = range(minGrad,maxGrad)
         , ylim = c(0,maximumY)
         , breaks = binBreaks
         , col = rgb(1,0,0,0.5)
         , add = TRUE)
    abline(v=meanInterTwo
           , col = "red"
           , lwd = 2)
    dev.off()
  }
  
}



########################## TAXASUMMARIES SECTION ####################################
print("In taxa summaries section")

library(outliers)
library(stats)
output <- "TaxaSummaries"

system(paste0("mkdir ",output))
system(paste0("mkdir ./",output,"/IndividualPlots"))
# system(paste0("biom convert -i ",biom," -o ./",output,"/OTUTable_text.txt --to-tsv --header-key taxonomy"))

############## Pre-amble; loading files, formatting, relative abund ###################
# gradientNames <- read.delim(paste0(gradient)
#                             , strip.white = TRUE
#                             , header = FALSE
#                             , sep = ","
#                             , stringsAsFactors = FALSE)

# Load taxa summaries and strip taxonomy
# taxa <- read.delim(paste0("./",output,"/OTUTable_text.txt"), strip.white = TRUE, stringsAsFactors = FALSE, header = TRUE, row.names = 1, skip = 1)
taxa <- read.delim(paste0(biomPWD)
                   , strip.white = TRUE
                   , stringsAsFactors = FALSE
                   , header = TRUE
                   , row.names = 1
                   , skip = 1)
colnames(taxa) <- gsub(".","-", colnames(taxa), fixed = TRUE)

# Make taxonomy ref
taxonomyRef <- data.frame(taxa[,ncol(taxa)])
rownames(taxonomyRef) <- rownames(taxa)

# Make taxa summary separate
taxa <- taxa[,-(ncol(taxa))]

# Make taxa summary relative abundance
relAbund <- function(x) {
  total <- sum(x)
  newx <- x/total
  return(newx)
}

# Apply function to taxa
taxaAbund <- data.frame(apply(taxa, 2, FUN = relAbund))
names(taxaAbund) <- names(taxa)

#--

# Load metadata
metadata <- read.delim(paste0(metadataPWD), stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""), row.names = 1)
rownames(metadata) <- gsub(".","-",rownames(metadata), fixed = TRUE)

# Make sure only samples in metadata are ones in taxa summaries

metadata <- metadata[rownames(metadata) %in% names(taxaAbund),]

# Sort taxa summaries sites and metadata so they are increasing in gradient
# Use re-order with gradientNames[4] to make in order of gradient
metadata <- metadata[order(metadata[[paste0(gradientNames[4])]]),]
# then make taxa the same
taxaAbund <- taxaAbund[,rownames(metadata)]

# 
# # Split ferdous and I's taxas up by metadata$Person
# 
# taxa.F <- taxa[,grep("F", metadata$Person)]
# taxa.M <- taxa[,grep("M", metadata$Person)]

#-- 

# Load modelBoundaries
modelBoundaries <- read.delim(paste0(modelboundariesPWD)
                              , header = TRUE
                              , row.names = 1
                              , strip.white = TRUE
                              , stringsAsFactors = FALSE
                              , na.strings = c('','na','NA'))

#--
# # OTU table by type
# otuType <- read.delim(paste0(typetable)
#                           , header = TRUE
#                           , row.names = 1
#                           , stringsAsFactors = FALSE)
# colnames(otuType) <- gsub('X', '', colnames(otuType))
# # Phylogeny
# phyloTree <- read.tree(paste0(phylotree))
# outgrouptree <- phyloTree$tip.label[length(phyloTree$tip.label)]
# phyloTree.rooted <- root(phyloTree, outgrouptree, resolve.root = TRUE)

#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~
# 
# #### Collapsing replicates
# 
# # Making new stuff
# taxa.F.col <- taxa.F
# taxa.M.col <- taxa.M
# 
# # Collapsing reps by deleting ends and doing means
# names(taxa.F.col) <- gsub(".r[0-9]","",names(taxa.F.col))
# taxa.F.col <- data.frame(do.call(cbind, by (t(taxa.F.col), INDICES = names(taxa.F.col), FUN = colMeans)))
# 
# # Doing same as above to my data set
# names(taxa.M.col) <- gsub("[.][0-9][.]",".", names(taxa.M.col))
# names(taxa.M.col) <- gsub("[0-9][A-C]","", names(taxa.M.col))
# taxa.M.col <- data.frame(do.call(cbind, by(t(taxa.M.col), INDICES = names(taxa.M.col), FUN = colMeans)))
# 
# # Now, collapse metadata file by rownames option
# metadata$ColRep <- rownames(metadata)
# 
# # Collapsing reps by deleting ends and doing median
# metadata$ColRep <- gsub(".r[0-9]","", metadata$ColRep)
# metadata$ColRep <- gsub("[.][0-9][.]Jul",".Jul", metadata$ColRep)
# metadata$ColRep <- gsub("[.][0-9][.]Aug",".Aug", metadata$ColRep)
# metadata$ColRep <- gsub("[0-9][A-C]","", metadata$ColRep)
# 
# # Aggregate data by collapsing via ColRep: 
# # If values are numeric, then find the mean. 
# # But if they are not, take the first value if they are all the same.
# # Finally, they are non-numeric and not identical, just put "NA"
# 
# metadata.col <- aggregate(metadata, by = list(metadata$ColRep), function(x) {
#   if (is.numeric(x) == TRUE) {
#     mean(x)
#   } else{
#     if(length(unique(x)) == 1){
#       x[1]
#     } else{
#       NA
#     }
#   }
# })
# 
# # Make into data frame and change rownames to column 1
# 
# metadata.col <- data.frame(metadata.col, row.names = 1)
# 
# # Delete last 3 rows because they are controls
# metadata.col$Description
# delete <- c(grep("Control", metadata.col$Description))
# metadata.col <- metadata.col[-(delete),]
# 
# #~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~

###### Now for the exciting part (kind of). Just adjusting orders and things before plotting #####
# 
# # First, I want to sort my metadata so that it is increasing in gradient.
# metadata.col <- metadata.col[order(metadata.col$SalinityEnviron),]
# 
# # Now, I use this to sort my data. same names to a vector and check for matches. Save these matches as the 'ordered.sites' and then extract things in that order, and re-save it.
# 
# tst.F.col <- c(names(taxa.F.col), row.names(metadata.col))
# ordered.sites.F.col <- tst.F.col[duplicated(tst.F.col)]
# 
# taxa.F.col.ordered <- taxa.F.col[,c(ordered.sites.F.col)]
# 
# 
# # Now for "M"
# tst.M.col <- c(names(taxa.M.col), row.names(metadata.col))
# ordered.sites.M.col <- tst.M.col[duplicated(tst.M.col)]
# 
# taxa.M.col.ordered <- taxa.M.col[,c(ordered.sites.M.col)]


## Now, I can make plots based on taxonomy.

# I want to space it out so it is relevant to how saline it is.

# # Recall the ordered sites is called: ordered.sites
# # match the sites to the rows names in metadata
# rows.tokeep.col <- match(row.names(metadata), ordered.sites.F.col, nomatch = FALSE)
# rows.tokeep.F.col <- rows.tokeep.F.col > 0
# Then, extract gradient values for these
ordered.gradient <- metadata[[paste0(gradientNames[4])]]

# Now, I want to replace the headers of the euktaxa.ordered with these values.
taxaAbund.grad <- taxaAbund
names(taxaAbund.grad) <- ordered.gradient


# For M
# I want to space it out so it is relevant to how saline it is.

# # Recall the ordered sites is called: ordered.sites
# # match the sites to the rows names in metadata
# rows.tokeep.M.col <- match(row.names(metadata.col), ordered.sites.M.col, nomatch = FALSE)
# rows.tokeep.M.col <- rows.tokeep.M.col > 0
# # Then, extract gradient values for these
# ordered.gradient.M.col <- metadata.col$SalinityEnviron[rows.tokeep.M.col]
# 
# # Now, I want to replace the headers of the euktaxa.ordered with these values.
# 
# taxa.M.col.ordered.grad <- as.data.frame(taxa.M.col.ordered)
# names(taxa.M.col.ordered.grad) <- ordered.gradient.M.col



# Now, I want to take these 'trend' lines and plot them all on single plots so I can see general trends.
color.random <- sample(colors()[grep('gr(a|e)y|white', colors(), invert = TRUE)], 200)

# Changing the dataframe into a matrix so that I can plug it directly into barplot. 
taxaAbund.matrix <- as.matrix(taxaAbund.grad, stringsAsFactors = FALSE, header = TRUE, rownames.force = TRUE)

# Want to look at ONLY things that are grthan the 0.05 threshold
# Create an empty vector
grthanTF <- vector()

# Now, for each row in taxa.matrix, I want to test whether ANY value in that line is greater than 0.05.
# I take this vector (T/F,T/F,T/F) and ask whether any of them are true; if one of them are, then that line is saved to the vector grthan.
for (i in 1:nrow(taxaAbund.matrix)) {
  line <- as.factor(taxaAbund.matrix[i,] >= as.numeric(minthreshold))
  if (any(line == "TRUE")) {
    grthanTF[i] <- TRUE
  } else{
    grthanTF[i] <- FALSE
  }
}


#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~

### INDIVIDUAL PLOTS
print("Making individual plots of taxa distributions")
#### USE THIS ONE IF YOU WANT TO LOOK AT INDIVDIUAL PLOTS TO SEE IF FIT OF LINES ARE GOOD OR NOT
# I've commented out the loess lines so they look cleaner. you can include them if you want by un-commenting



# Plot only major taxa
taxaAbund.grad.filt <- taxaAbund.grad[grthanTF,]



# Find maximum relative abundance
maxAbund <- max(taxaAbund.grad.filt)

positionsTypes <- unlist(lapply(rownames(taxaAbund.grad.filt), FUN = function(x) grep(paste0("^",x,"$"), rownames(modelBoundaries))))
typesforIndiv <- modelBoundaries$type[positionsTypes]

# fits <- data.frame()
for (i in (1:nrow(taxaAbund.grad.filt))) {
  pdf(paste0('./',output,'/IndividualPlots/Indiv-',row.names(taxaAbund.grad.filt)[i],typesforIndiv[i],'.pdf'))
  # quartz(paste(row.names(taxaAbund.grad.filt)[i]),5,5)
  plot(as.numeric(taxaAbund.grad.filt[i,]) ~ ordered.gradient
       , las = 2
       , ylim = c(0,maxAbund)
       , type = "p"
       , ylab = "Relative Abundance"
       , lwd = 2
       , xlab = "Gradient"
       , pch = 21
       , col = "grey"
       , bg = "grey"
       , main = rownames(taxaAbund.grad.filt[i,])
       , cex.main = 0.8
       , xlim = c(minGrad,maxGrad)
  )

  # Fit a polnomial model to the data
  # First, make function that gets rid of 'outliers'
  outliertest <- scores(as.numeric(taxaAbund.grad.filt[i,]), type = "chisq", prob = 0.99)
  if (anyNA(outliertest)) {
    outliertestTF <- rep(TRUE,length(outliertest[,1]))
  } else {
    outliertestTF <- !outliertest
  }
  
  taxaAbundNooutliers <- taxaAbund.grad.filt[i,][outliertestTF]
  gradientNooutliers <- ordered.gradient[outliertestTF]
  
  assign(sub("^.*__.*__","", paste0(row.names(taxaAbund.grad.filt)[i],".lo")), loess(as.numeric(taxaAbundNooutliers) ~ gradientNooutliers))
  xl <- seq(min(gradientNooutliers, na.rm = TRUE),max(gradientNooutliers, na.rm = TRUE), (max(gradientNooutliers, na.rm = TRUE)-min(gradientNooutliers, na.rm = TRUE))/1000)
  lines(xl
        , predict(get(sub("^.*__.*__"
                          ,""
                          , paste0(row.names(taxaAbund.grad.filt)[i],".lo"))),xl)
        , col = "blue"
        , lwd = 2)

  dev.off()
}


#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~

######## Line plot of taxa distributions #############
print("Making Line plots of taxa distributions")

# First, look at taxa and set colours.
# rownames(taxaAbund.grad)
# [1] "k__Archaea;p__Euryarchaeota;c__Thermoplasmata"         
# [2] "k__Bacteria;p__Actinobacteria;c__Acidimicrobiia"       
# [3] "k__Bacteria;p__Actinobacteria;c__Actinobacteria"       
# [4] "k__Bacteria;p__Actinobacteria;c__Thermoleophilia"      
# [5] "k__Bacteria;p__Armatimonadetes;c__Armatimonadia"       
# [6] "k__Bacteria;p__Bacteroidetes;c__Cytophagia"            
# [7] "k__Bacteria;p__Bacteroidetes;c__Flavobacteriia"        
# [8] "k__Bacteria;p__Bacteroidetes;c__Sphingobacteriia"      
# [9] "k__Bacteria;p__Bacteroidetes;c__[Saprospirae]"         
# [10] "k__Bacteria;p__Cyanobacteria;c__Chloroplast"           
# [11] "k__Bacteria;p__Cyanobacteria;c__Synechococcophycideae" 
# [12] "k__Bacteria;p__Firmicutes;c__Bacilli"                  
# [13] "k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria"  
# [14] "k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria"  
# [15] "k__Bacteria;p__Proteobacteria;c__Epsilonproteobacteria"
# [16] "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria"  
# [17] "k__Bacteria;p__Verrucomicrobia;c__Verrucomicrobiae"    
# [18] "k__Bacteria;p__Verrucomicrobia;c__[Spartobacteria]"

# colors.lines <- c("firebrick3" # thermoplasmata
#                   ,"goldenrod4" # Acidimicrobiia
#                   ,"gold" # Aactinobacteria
#                   ,"coral3" # thermoleophilia
#                   ,"deeppink4" # Armatimonadia
#                   ,"blue" # Cytophagia
#                   ,"orange" # Flavobacteriia
#                   ,"darkorange3" # Sphingobacteriia
#                   ,"darkgoldenrod4" # Saprospirae
#                   ,"chartreuse3" # chloroplast
#                   ,"aquamarine3" # Snechococcophycideae
#                   ,"brown" # bacilli
#                   ,"pink" # alphaproteobacteria
#                   , "beige" # Deltaprteobacteria
#                   ,"forestgreen" # Epsilonproteobacteria
#                   ,"bisque4" # Gammaproteobacteria
#                   ,"purple" # Verrucomicrobiae
#                   ,"cornsilk4" # Spartobacteria
# )
# For this, we want to get rid of extremely low abundance taxa I think...

# First, make an empty plot with the limits set as the maximum found in my data.frame
pdf(paste0('./',output,'/Taxadistributions.pdf'))
# quartz("Taxa distributions",10,5)
# par(mar = c(4,4,1,1), oma = c(2,2,2,2))
plot(x = NULL, y = NULL
     , xlim = c(1,max(ordered.gradient, na.rm = TRUE))
     , ylim = c(0,max(taxaAbund.grad.filt)*1.25)
     , ylab = "Abundance (%)"
     , xlab = "Salinity"
     , main = "Loess curve taxa distributions across salinitiy"
)
# par(fig=c(0,0.5,0,1), new = TRUE)


# Now, plot loess prediction lines on a single plot for all taxa in table. This is to see whether there are trends in OTU abundance across gradient.
# We save the summary() residuals for each 'taxa' in the "fits" file so you can look at it later and see how well it fits.

# Get random colours
colors.lines <- sample(colors()[grep('gr(a|e)y|white', colors(), invert = TRUE)], nrow(taxaAbund.grad.filt), replace = TRUE)

fits <- data.frame()
for (i in (1:nrow(taxaAbund.grad.filt))) {# First, make function that gets rid of 'outliers'
  outliertest <- scores(as.numeric(taxaAbund.grad.filt[i,]), type = "chisq", prob = 0.99)
  if (anyNA(outliertest)) {
    outliertestTF <- rep(TRUE,length(outliertest[,1]))
  } else {
    outliertestTF <- !outliertest
  }
  
  taxaAbundNooutliers <- taxaAbund.grad.filt[i,][outliertestTF]
  gradientNooutliers <- ordered.gradient[outliertestTF]
  
  assign(sub("^.*__.*__","", paste0(row.names(taxaAbund.grad.filt)[i],".lo")), loess(as.numeric(taxaAbundNooutliers) ~ gradientNooutliers))
  xl <- seq(min(gradientNooutliers, na.rm = TRUE),max(gradientNooutliers, na.rm = TRUE), (max(gradientNooutliers, na.rm = TRUE)-min(gradientNooutliers, na.rm = TRUE))/1000)
  lines(xl, predict(get(sub("^.*__.*__","", paste0(row.names(taxaAbund.grad.filt)[i],".lo"))),xl)
        , col = colors.lines[i]
        , lwd = 2)
  
}
dev.off()

# names(fits) <- c("Taxa","ResStErr")
# par(fig = c(0,1,0,1), mar = c(0,0,0,0), oma = c(0,0,0,2), new = TRUE)
# plot(0,0, type = "n", bty="n", xaxt="n", yaxt = "n")
# 
# legend("right",
#        legend = row.names(taxaAbund.grad.filt),
#        col = colors.lines[c(1:nrow(taxaAbund.grad.filt))],
#        lty = c(1,1,1),
#        lwd = 4,
#        cex = 0.8,
#        bty = "n",
#        y.intersp = 1.2,
# )

#~~~#~~~~#~~~~#~~~~#~~~~#~~~~#~~~~#~~~~#~~~~#~~~~#~~~~#~~~#

######## Line plot, by type (all) ###########
print("Making line plots, by type")

modelBoundariesEdited <- modelBoundaries
# for (i in 1:length(modelBoundariesEdited$type)) {
#   if (modelBoundariesEdited$type[i] == 'noclass'){
#     if (modelBoundariesEdited$bloom[i] != 'No')
#       modelBoundariesEdited$type[i] <- paste0('bloom',modelBoundariesEdited$bloom[i])
#   }
# }

typeUnique <- unique(modelBoundariesEdited$type)
taxaAbundTypes <- data.frame()
metaTypes <- data.frame()
for (i in 1:length(typeUnique)){
  assign(paste0(typeUnique[i],'OTUs'), rownames(modelBoundariesEdited)[grep(paste0("^",typeUnique[i],"$"), modelBoundariesEdited$type)])
  get(paste0(typeUnique[i],'OTUs'))
  for (j in 1:length(get(paste0(typeUnique[i],'OTUs')))){
    tempAbund <- taxaAbund.grad[grep(paste0("^",get(paste0(typeUnique[i],'OTUs'))[j],"$"), rownames(taxaAbund.grad)),]
    tempAbund
    taxaAbundTypes <- rbind(taxaAbundTypes,tempAbund)
  }
  tempMeta <- cbind(get(paste0(typeUnique[i],'OTUs')), rep_len(typeUnique[i],length(paste0(typeUnique[i],'OTUs'))))
  metaTypes <- rbind(metaTypes,tempMeta)
}

# Make the taxaAbundTypes equal in value for each OTU (not location)

# Make taxa summary relative abundance
# Apply function to taxa

taxaAbundTypesRelative <- t(apply(taxaAbundTypes, 1, FUN = relAbund))
# is.matrix(taxaAbundTypesRelative)
# taxaAbundTypesRelative <- data.frame(apply(taxaAbundTypes, 1, FUN = relAbundOTU))
names(taxaAbund) <- names(taxa)

metaTypes <- data.frame(metaTypes,row.names = 1)
colnames(metaTypes) <- "type"

# Make colors for number of types
typesforApprovedColors <- c(paste0(gradientNames[1], "Restricted")
                            , paste0(gradientNames[1], "Bloom")
                            , paste0(gradientNames[1], "Peak")
                            , paste0(gradientNames[2], "Restricted")
                            , paste0(gradientNames[2], "Bloom")
                            , paste0(gradientNames[2], "PeakLoToler")
                            , paste0(gradientNames[2], "Peak")
                            , paste0(gradientNames[2], "Peak,HiToler")
                            , paste0(gradientNames[3], "Peak")
                            , paste0(gradientNames[3], "Restricted")
                            , paste0(gradientNames[3], "Bloom")
                            ,"noclass"
                            , "ubiquitous"
                            , paste0('inv',gradientNames[2])
)
colorsApproved <- c("blue","yellow","cyan","purple","lightgreen","magenta","lightblue","pink","orange","red", "green","black","gray", "white")

# typesforApprovedColors <- c(paste0(gradientNames[1])
#                             , paste0("half",gradientNames[1])
#                             , paste0(gradientNames[2])
#                             , paste0('lo',gradientNames[2])
#                             , paste0('half',gradientNames[2])
#                             , paste0('hi',gradientNames[2])
#                             , paste0('half',gradientNames[3])
#                             , paste0(gradientNames[3])
#                             ,"noclass"
#                             , paste0('inv',gradientNames[2]))
# colorsApproved <- c("blue","cyan","purple","magenta","lightblue","pink","orange","red","green","darkgreen")
positionColors <- unlist(lapply(typeUnique,FUN = function(x) grep(paste0('^',x,'$'), typesforApprovedColors)))
colors.multiline <- unlist(lapply(positionColors, FUN = function(x) colorsApproved[x]))
# colors.multiline <- sample(colorsApproved, length(typeUnique), replace = FALSE)
colors.multilineplot <- colors.multiline[factor(metaTypes$type)]

# Get maxheight of stuff for graph
# maxValue <- max(taxaAbundTypesRelative)


# First, make an empty plot with the limits set as the maximum found in my data.frame
pdf(paste0('./',output,'/TaxadistributionsRelative.pdf'))
# par(mar = c(4,4,1,1), oma = c(2,2,2,2))
par(fig=c(0,0.67,0,1))
plot(x = NULL, y = NULL
     , xlim = c(1,max(ordered.gradient, na.rm = TRUE))
     , ylim = c(0,max(taxaAbundTypesRelative)*1.25)
     , ylab = "Abundance (%)"
     , xlab = "Salinity"
     , main = "Loess curves of taxa distributions (Relative)"
)

# Now, plot loess prediction lines on a single plot for all taxa in table. This is to see whether there are trends in OTU abundance across gradient.
# We save the summary() residuals for each 'taxa' in the "fits" file so you can look at it later and see how well it fits.

fits <- data.frame()
for (i in (1:nrow(taxaAbundTypesRelative))) {
  # First, make function that gets rid of 'outliers'
  outliertest <- scores(as.numeric(taxaAbundTypesRelative[i,]), type = "chisq", prob = 0.99)
  if (anyNA(outliertest)) {
    outliertestTF <- rep(TRUE,length(outliertest))
  } else {
    outliertestTF <- !outliertest
  }
  
  taxaAbundNooutliers <- as.numeric(taxaAbundTypesRelative[i,])[which(outliertestTF)]
  gradientNooutliers <- ordered.gradient[outliertestTF]
  
  assign(paste0(row.names(taxaAbundTypesRelative)[i],".lo"), loess(as.numeric(taxaAbundNooutliers) ~ gradientNooutliers))
  xl <- seq(min(gradientNooutliers, na.rm = TRUE),max(gradientNooutliers, na.rm = TRUE), (max(gradientNooutliers, na.rm = TRUE)-min(gradientNooutliers, na.rm = TRUE))/1000)
  lines(xl, predict(get(paste0(row.names(taxaAbundTypesRelative)[i],".lo")),xl)
        , col = adjustcolor(colors.multilineplot[i],alpha.f = 0.2)
        , lwd = 2)
}
par(fig = c(0.5,1,0,1), new = TRUE)
plot(0,0
     , bty = "n"
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , pch = "")
legend("center"
       , legend = levels(as.factor(metaTypes$type))
       , pch = 19
       , col = colors.multiline
)
dev.off()

# names(fits) <- c("Taxa","ResStErr")
# par(fig = c(0,1,0,1), mar = c(0,0,0,0), oma = c(0,0,0,2), new = TRUE)
# plot(0,0, type = "n", bty="n", xaxt="n", yaxt = "n")
# 
# legend("right",
#        legend = row.names(taxaAbund.grad.filt),
#        col = colors.lines[c(1:nrow(taxaAbund.grad.filt))],
#        lty = c(1,1,1),
#        lwd = 4,
#        cex = 0.8,
#        bty = "n",
#        y.intersp = 1.2,
# )

######## Line plot, by type (condensed) ###########
print("Making line plots, by type (condensed)")

modelBoundariesEdited <- modelBoundaries

typeUnique <- unique(modelBoundariesEdited$typeSimple)
taxaAbundTypes <- data.frame()
metaTypes <- data.frame()
for (i in 1:length(typeUnique)){
  assign(paste0(typeUnique[i],'OTUs'), rownames(modelBoundariesEdited)[grep(paste0("^",typeUnique[i],"$"), modelBoundariesEdited$typeSimple)])
  get(paste0(typeUnique[i],'OTUs'))
  # tempAbund <- data.frame()
  for (j in 1:length(get(paste0(typeUnique[i],'OTUs')))){
    tempAbund <- taxaAbund.grad[grep(paste0("^",get(paste0(typeUnique[i],'OTUs'))[j],"$"), rownames(taxaAbund.grad)),]
      taxaAbundTypes <- rbind(taxaAbundTypes,tempAbund)
  }
 
  tempMeta <- cbind(get(paste0(typeUnique[i],'OTUs')), rep_len(typeUnique[i],length(paste0(typeUnique[i],'OTUs'))))
  metaTypes <- rbind(metaTypes,tempMeta)
}

# Make the taxaAbundTypes equal in value for each OTU (not location)

# Make taxa summary relative abundance
# Apply function to taxa

taxaAbundTypesRelative <- t(apply(taxaAbundTypes, 1, FUN = relAbund))
# is.matrix(taxaAbundTypesRelative)
# taxaAbundTypesRelative <- data.frame(apply(taxaAbundTypes, 1, FUN = relAbundOTU))
names(taxaAbund) <- names(taxa)

metaTypes <- data.frame(metaTypes,row.names = 1)
colnames(metaTypes) <- "type"

# Make colors for number of types
typesforApprovedColors <- c(paste0(gradientNames[1], "Restricted")
                            # , paste0('bloom',gradientNames[1])
                            # , paste0("half",gradientNames[1])
                            , paste0(gradientNames[2], "Restricted")
                            # , paste0('bloom',gradientNames[2])
                            # , paste0('lo',gradientNames[2])
                            # , paste0('half',gradientNames[2])
                            # , paste0('hi',gradientNames[2])
                            # , paste0('half',gradientNames[3])
                            , paste0(gradientNames[3], "Restricted")
                            # , paste0('bloom',gradientNames[3])
                            ,"noclass"
                            # , paste0('inv',gradientNames[2])
)
colorsApproved <- c("blue","purple","red", "black")

# typesforApprovedColors <- c(paste0(gradientNames[1])
#                             , paste0("half",gradientNames[1])
#                             , paste0(gradientNames[2])
#                             , paste0('lo',gradientNames[2])
#                             , paste0('half',gradientNames[2])
#                             , paste0('hi',gradientNames[2])
#                             , paste0('half',gradientNames[3])
#                             , paste0(gradientNames[3])
#                             ,"noclass"
#                             , paste0('inv',gradientNames[2]))
# colorsApproved <- c("blue","cyan","purple","magenta","lightblue","pink","orange","red","green","darkgreen")
positionColors <- unlist(lapply(typeUnique,FUN = function(x) grep(paste0('^',x,'$'), typesforApprovedColors)))
colors.multiline <- unlist(lapply(positionColors, FUN = function(x) colorsApproved[x]))
# colors.multiline <- sample(colorsApproved, length(typeUnique), replace = FALSE)
colors.multilineplot <- colors.multiline[factor(metaTypes$type)]

# Get maxheight of stuff for graph
# maxValue <- max(taxaAbundTypesRelative)


# First, make an empty plot with the limits set as the maximum found in my data.frame
pdf(paste0('./',output,'/TaxadistributionsRelative_condensed.pdf'))
# par(mar = c(4,4,1,1), oma = c(2,2,2,2))
par(fig=c(0,0.67,0,1))
plot(x = NULL, y = NULL
     , xlim = c(1,max(ordered.gradient, na.rm = TRUE))
     , ylim = c(0, max(taxaAbundTypesRelative)*1.25)
     , ylab = "Abundance (%)"
     , xlab = "Salinity"
     , main = "Loess curves of taxa distributions (Relative)"
)

# Now, plot loess prediction lines on a single plot for all taxa in table. This is to see whether there are trends in OTU abundance across gradient.
# We save the summary() residuals for each 'taxa' in the "fits" file so you can look at it later and see how well it fits.

fits <- data.frame()
for (i in (1:nrow(taxaAbundTypesRelative))) {
  # First, make function that gets rid of 'outliers'
  outliertest <- scores(as.numeric(taxaAbundTypesRelative[i,]), type = "chisq", prob = 0.99)
  if (anyNA(outliertest)) {
    outliertestTF <- rep(TRUE,length(outliertest))
  } else {
    outliertestTF <- !outliertest
  }
  
  taxaAbundNooutliers <- as.numeric(taxaAbundTypesRelative[i,])[which(outliertestTF)]
  gradientNooutliers <- ordered.gradient[outliertestTF]
  
  assign(paste0(row.names(taxaAbundTypesRelative)[i],".lo"), loess(as.numeric(taxaAbundNooutliers) ~ gradientNooutliers))
  xl <- seq(min(gradientNooutliers, na.rm = TRUE),max(gradientNooutliers, na.rm = TRUE), (max(gradientNooutliers, na.rm = TRUE)-min(gradientNooutliers, na.rm = TRUE))/1000)
  lines(xl, predict(get(paste0(row.names(taxaAbundTypesRelative)[i],".lo")),xl)
        , col = adjustcolor(colors.multilineplot[i],alpha.f = 0.2)
        , lwd = 2)
  
}
par(fig = c(0.5,1,0,1), new = TRUE)
plot(0,0
     , bty = "n"
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , pch = "")
legend("center"
       , legend = levels(as.factor(metaTypes$type))
       , pch = 19
       , col = colors.multiline
)
dev.off()

#~~~#~~~~#~~~~#~~~~#~~~~#~~~~#~~~~#~~~~#~~~~#~~~~#~~~~#~~~#

######## Single Line plot, by type ###########

print("Single line plot by type")

aggregateRelativeAbund.df <- taxaAbundTypesRelative

aggregateRelativeAbund.df <- data.frame(aggregate(aggregateRelativeAbund.df, list(metaTypes$type), FUN = mean), row.names = 1 )
colnames(aggregateRelativeAbund.df) <- gsub("X","",colnames(aggregateRelativeAbund.df))

# First, make an empty plot with the limits set as the maximum found in my data.frame
pdf(paste0('./',output,'/TaxadistributionsRelativeAggregate.pdf'), width = 7, height = 5)
# par(mar = c(4,4,1,1), oma = c(2,2,2,2))
par(fig=c(0,0.67,0,1))
plot(x = NULL, y = NULL
     , xlim = c(1,max(ordered.gradient, na.rm = TRUE))
     , ylim = c(0,max(aggregateRelativeAbund.df)*1.25)
     , ylab = "Abundance (%)"
     , xlab = "Salinity"
     , main = "Aggregate loess curves of taxa distributions (Relative)"
)

# Now, plot loess prediction lines on a single plot for all taxa in table. This is to see whether there are trends in OTU abundance across gradient.
# We save the summary() residuals for each 'taxa' in the "fits" file so you can look at it later and see how well it fits.

fits <- data.frame()
for (i in (1:nrow(aggregateRelativeAbund.df))) {
  
  assign(paste0(row.names(aggregateRelativeAbund.df)[i],".lo"), loess(as.numeric(aggregateRelativeAbund.df[i,])~ ordered.gradient))
  xl <- seq(min(ordered.gradient, na.rm = TRUE),max(ordered.gradient, na.rm = TRUE), (max(ordered.gradient, na.rm = TRUE)-min(ordered.gradient, na.rm = TRUE))/1000)
  lines(xl, predict(get(paste0(row.names(aggregateRelativeAbund.df)[i],".lo")),xl)
        , col = adjustcolor(unique(colors.multilineplot)[i],alpha.f = 0.99)
        , lwd = 2)
  
}
par(fig = c(0.5,1,0,1), new = TRUE)
plot(0,0
     , bty = "n"
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , pch = "")
legend("center"
       , legend = levels(as.factor(metaTypes$type))
       , pch = 19
       , col = unique(colors.multiline)
)
dev.off()




#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~

######## Making barplot; taxa summaries ##########
print("Making barplots")

# For the data matrix, I want to colour the bars so that all taxanomic groups that represent more than 5% of the total OTUs found in any sample are coloured. Note that this is slightly different than in the paper; they claim to use a threshold of 1%, but when I looked at the data it seemed they were wrong. There were samples with more than 1% (but less than 5%) representation that they did not label. Thus, I will be using the 5% threshold.

# Then, I use the vector 'grthan' to do some things:

# If I don't want to just have random colors, I can do this:
# levels(coloredbars)
# [1] "k__Archaea;p__Euryarchaeota;c__Thermoplasmata"        
# [2] "k__Bacteria;p__Actinobacteria;c__Actinobacteria"      
# [3] "k__Bacteria;p__Armatimonadetes;c__Armatimonadia"      
# [4] "k__Bacteria;p__Bacteroidetes;c__[Saprospirae]"        
# [5] "k__Bacteria;p__Bacteroidetes;c__Cytophagia"           
# [6] "k__Bacteria;p__Bacteroidetes;c__Flavobacteriia"       
# [7] "k__Bacteria;p__Bacteroidetes;c__Sphingobacteriia"     
# [8] "k__Bacteria;p__Cyanobacteria;c__Chloroplast"          
# [9] "k__Bacteria;p__Cyanobacteria;c__Synechococcophycideae"
# [10] "k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria" 
# [11] "k__Bacteria;p__Verrucomicrobia;c__[Spartobacteria]"   
# [12] "k__Bacteria;p__Verrucomicrobia;c__Verrucomicrobiae" 
# color.purpose <- c("firebrick3" # Thermoplasmata
#                    ,"gold" # Actinobacteria
#                    ,"deeppink4" # Armatimonadia
#                    ,"darkgoldenrod4" # Saprospirae
#                    ,"blue" # Cytophagia
#                    ,"orange" # Flavobacteria
#                    ,"darkorange3" # Sphingobacteriia
#                    ,"chartreuse3" # Chloroplast
#                    ,"aquamarine3" # Synechococcophycideae
#                    ,"pink" # Alphaproteobacteria
#                    ,"cornsilk4" # Spartobacteria
#                    ,"purple" # Verrucomicrobiae
#                    )


# Try to do a calculation-based way to space out taxa so they are spread out by gradient.
# First, list all unique salinities
unique.num <- unique(names(taxaAbund.grad.filt))

# Then, count how many of each unique gradient there is. 
# This loop counts how many reps there are of each 'unique' gradient in the taxa file
spacing <- data.frame()
for (i in 1:length(unique.num)) {
  spacing[i,1] <- as.numeric(unique.num[i])
  spacing[i,2] <- sum(unique.num[i] == names(taxaAbund.grad.filt))
}
names(spacing) <- c("Sal","count")

# Let us imagine that each "gradient" is worth 1 unit. each separation of gradient, then, is also worth 1 unit. 
# Therefore, the barwidth should be calculated as 1/#ofsites@thesamegradient, then multiplied by a constant to make the average bar width 1
# We need the average bar width to be 1 because the 'space' in the barplot function is calculated as a multiple of the average bar width. Thus, we want a '1' 'to actually equal '1'
# The distance between clusters should be calculated as 1*difference, where the differences are scaled to the average bar width.


# Calculating proportional band widths.
spacing.width <- vector()
for (i in 1:length(spacing[,2])) {
  
  spacing.width <- c(spacing.width, rep((1/spacing[i,2]), spacing[i,2]))
  
}
# And now multiplying by 1/barwidth so the mean bar width is one

spacing.wid <- spacing.width*(1/mean(spacing.width))

# Calculating distances
spacing.num <- vector()
for (i in 1:length(spacing[,2])) {
  if (i == 1) {
    spacing.num <- c(0, rep(0, spacing[i,2]-1))
  } else {
    spacing.num <- c(spacing.num, (spacing[i,1]-spacing[i-1,1])*1 , rep(0, spacing[i,2]-1))
  }
}

# COLORS
grthan <- gsub(FALSE,NA,grthanTF)
grthan <- as.logical(grthan)

# I make a factor that tells me which row names should have coloured bars. Bars with less than 5% representation across all samples are listed as "NA".
coloredbars <- as.factor(rownames((taxaAbund.matrix[grthan,])))

# I make a factor that tells me how many different coloured bars there are. Since the factor above has NA values as well, I count the ones that are NOT "NA".
ncoloredbars <-sum(!is.na(coloredbars))

# # I make a factor that samples random colours (that are not grey-- because I'm using white for <5% samples and I don't want it to get ambiguous) to be used in graphing. It has a length of ncoloredbars.
# color.random <- sample(colors()[grep('gr(a|e)y', colors(), invert = TRUE)], ncoloredbars)

pdf(paste0("./",output,"/TaxaSummariesPartitioned.pdf"), width = 10, height = 5)
# quartz("Taxa summaries", 10,5)
# par(mar = c(4,4,0,0), oma = c(5,1,2,1))
# par(fig=c(0,0.7,0,1), new= TRUE)
barplot(taxaAbund.matrix 
        , col = c(rev(colors.lines))[coloredbars] # Using the generated colours from above, I colour each bar that has >5% representation a different colour.
        , space = spacing.num # This lines the bar plot up with the boxplot in spacing
        , cex.axis = 0.6 # scaling axis
        , cex.lab = 1 # scaling label
        , ylab = "% Representation" # y axis label
        , xlab = "Salinity"
        , xaxt = 'n'
        , width = spacing.wid
        , xlim = c(0,sum(spacing.wid,spacing.num))
        , border = NA
        , main = 'Taxa Summaries; partitioned'
        , sub = paste0('Minimum max count = ', minthreshold)
        
)

# Add manual x axis

axis(side = 1
     , at = seq(0,sum(spacing.wid,spacing.num), by = sum(spacing.wid,spacing.num)/(30/5))
     , labels = seq(0,30, by = 5))
dev.off()

# Type out legend labels because I want them a certain way
# legendlabels <- c("Class Thermoplasmata", "Class Actinobacteria", "Class Armatimonadia","Cytophagia","Flavobacteriia","Sphingobacteriia","Saprospirae","Chloroplast (Cyanobacteria)","Synechococcophycidae","Alphaproteobacteria","Verrumicrobiae","Spartobacteria")

# This makes an 'empty' plot with no labels or axis.
# par(fig=c(0.5,1,0,1), new = TRUE)
# plot(0,0, type = "n", bty="n", xaxt="n", yaxt = "n", ylab = "", xlab = "", pch = "")
# legend("right"
#        , legend = rev(c(legendlabels, "Other"))
#        , pch = 22 # squares
#        , col= "black" # black outline
#        , pt.bg = rev(c(as.character(color.purpose), "white")) #fill squares with colours
#        , bty = "n" # Get rid of outer box
#        , title = "Taxa"
#        , cex = 0.95
#        , xpd = FALSE
#        
# )

####### OTU TOLERANCE RANGES ##########
# use taxaAbund to get 'extreme' points as well.

# Order so that increasing in salinity (no necessary, but nice)
taxaAbundances.ordered <- taxaAbundances[,order(gradient)]
gradient.ordered <- gradient[order(gradient)]
# Make gradient a dataframe with same column names as taxaAbundances
gradient.ordered <- t(as.data.frame(gradient.ordered))
colnames(gradient.ordered) <- colnames(taxaAbundances.ordered)

getextreme <- function(x,taxaAbundances.ordered, gradient, threshold,editedModelBoundaries) {
  # Check to make sure is not 'noclass'
  # Get ID name
  tempName <- rownames(editedModelBoundaries)[x]
  # Get abundances for that ID name
  totalAbundances.temp <- taxaAbundances.ordered[grep(paste0("^",tempName,"$"), rownames(taxaAbundances.ordered)),]
  # Calculate relative threshold for that ID
  relThreshold <- max(as.vector(totalAbundances.temp[1,]))*threshold
  # Get all abundances where salinity is before and after the START and STOP
  taxaAbundances.ordered.BEFORETEMP <- totalAbundances.temp[as.vector(editedModelBoundaries['START'][x,] >= gradient.ordered)]
  taxaAbundances.ordered.AFTERTEMP <- totalAbundances.temp[as.vector(editedModelBoundaries['STOP'][x,] <= gradient.ordered)]
  
  skipposMin = FALSE
  if (length(taxaAbundances.ordered.BEFORETEMP) == 0) {
    extrStart <- editedModelBoundaries['START'][x,]
    skipposMin = TRUE
  }
  skipposMax = FALSE
  if (length(taxaAbundances.ordered.AFTERTEMP) == 0) {
    extrStop <- editedModelBoundaries['STOP'][x,]
    skipposMax = TRUE
  }
  
  # Get maximum salinity where you see above threshold
  
  if (!skipposMin) {
    whichtoMin <- which(taxaAbundances.ordered.BEFORETEMP >= relThreshold)
    if (length(whichtoMin) > 0) {
      posMin <- colnames(taxaAbundances.ordered.BEFORETEMP)[min(whichtoMin)]
      extrStart <- gradient.ordered[,posMin]
    } else {
      posMin <- colnames(taxaAbundances.ordered.BEFORETEMP)[ncol(taxaAbundances.ordered.BEFORETEMP)]
      extrStart <- gradient.ordered[,posMin]
      }
  }
  if (!skipposMax) {
    whichtoMax <- which(taxaAbundances.ordered.AFTERTEMP >= relThreshold)
    if (length(whichtoMax) > 0) {
      posMax <- colnames(taxaAbundances.ordered.AFTERTEMP)[max(whichtoMax)]
      extrStop <- gradient.ordered[,posMax]
    } else {
      posMax <- colnames(taxaAbundances.ordered.AFTERTEMP)[1]
      extrStop <- gradient.ordered[,posMax]
      }
  }

  extremeValues <- c(extrStart,extrStop)
  return(extremeValues)
}


editedModelBoundaries <- modelBoundaries[,c('type','typeSimple','boundaries','boundariestwo','bloom')]
# for (i in 1:nrow(editedModelBoundaries)) {
#   if (editedModelBoundaries$type[i] == "noclass"){
#     if (editedModelBoundaries$typeSimple[i] != "noclass") {
#       editedModelBoundaries$type[i] <- paste0("bloom",editedModelBoundaries$bloom[i])
#     }
#   }
# }
editedModelBoundaries <- editedModelBoundaries[-grep("noclass",editedModelBoundaries$typeSimple),]
editedModelBoundaries <- editedModelBoundaries[grep("[0-9].*", editedModelBoundaries$boundaries),]

editedModelBoundaries$START <- rep(NA,nrow(editedModelBoundaries))
editedModelBoundaries$STOP <- rep(NA, nrow(editedModelBoundaries))
editedModelBoundaries$EXTRMIN <- rep(NA, nrow(editedModelBoundaries))
editedModelBoundaries$EXTRMAX <- rep(NA, nrow(editedModelBoundaries))
for (i in 1:nrow(editedModelBoundaries)) {
  if (editedModelBoundaries$typeSimple[i] == paste0(gradientNames[1], "Restricted")) {
    editedModelBoundaries$START[i] <- minGrad
    editedModelBoundaries$STOP[i] <- editedModelBoundaries$boundaries[i]
  } else if (editedModelBoundaries$typeSimple[i] == paste0(gradientNames[3], "Restricted")) {
    editedModelBoundaries$START[i] <- editedModelBoundaries$boundaries[i]
    editedModelBoundaries$STOP[i] <- maxGrad
  } else if (editedModelBoundaries$typeSimple[i] == paste0(gradientNames[2], "Restricted")) {
    editedModelBoundaries$START[i] <- editedModelBoundaries$boundaries[i]
    editedModelBoundaries$STOP[i] <- editedModelBoundaries$boundariestwo[i]
  }
  extremeValuesTemp <- getextreme(i,taxaAbundances.ordered, gradient.ordered, threshold, editedModelBoundaries)
  editedModelBoundaries$EXTRMIN[i] <- extremeValuesTemp[1]
  editedModelBoundaries$EXTRMAX[i] <- extremeValuesTemp[2]
}



getColor <- function(x) {
  if (x == paste0(gradientNames[1], "Restricted")) {
    return("blue")
  } else if (x == paste0(gradientNames[2], "Restricted")) {
    return("purple")
  } else if (x == paste0(gradientNames[3], "Restricted")) {
    return("red")
  } else if (x == paste0(gradientNames[2], "PeakHiToler")) {
    return("pink")
  } else if (x == paste0(gradientNames[2], "PeakLoToler")) {
    return("magenta")
  } else if (x == paste0(gradientNames[1], "Peak")) {
    return("cyan")
  } else if (x == paste0(gradientNames[2], "Peak")) {
    return("pink4")
  } else if (x == paste0(gradientNames[3], "Peak")) {
    return("orange")
  } else if (x == paste0(gradientNames[1], "Bloom")) {
    return("darkblue")
  } else if (x == paste0( gradientNames[2], "Bloom")) {
    return("purple4")
  } else if (x == paste0(gradientNames[3], "Bloom")) {
    return("darkred")
  } else{
    return("NA")
  }
}

# PLOT ORDERED

# Get subset only
loLines <- editedModelBoundaries[grep(paste0(gradientNames[1]), editedModelBoundaries$typeSimple),]
interLines <- editedModelBoundaries[grep(paste0(gradientNames[2]), editedModelBoundaries$typeSimple),]
hiLines <- editedModelBoundaries[grep(paste0(gradientNames[3]), editedModelBoundaries$typeSimple),]

# Sorting lo and hi lines are easy
loLines <- loLines[order(loLines$STOP, decreasing = FALSE), ]
hiLines <- hiLines[order(hiLines$START, decreasing = FALSE),]

# But I want inter lines to be sorted by 'middle' I think
sizeForInter <- c()
midForInter <- c()
for (i in 1:nrow(interLines)) {
  size <- (interLines$STOP[i]-interLines$START[i])
  mid <- mean(c(interLines$STOP[i],interLines$START[i]))
  sizeForInter <- c(sizeForInter, size)
  midForInter <- c(midForInter, mid)
}

interLines <- interLines[order(midForInter,sizeForInter),]

# Combine fresh, brackish, marine
comboLinesOrdered <- rbind(loLines,interLines,hiLines)
# Get rid of NAs from fresh and marine dataset
comboLinesOrdered <- comboLinesOrdered[grep("[0-9]",comboLinesOrdered$STOP),]
comboLinesOrdered <- comboLinesOrdered[grep("[0-9]",comboLinesOrdered$START),]
# Get y intersect
sequenceOrderCombo <- seq(0,100,by = (100-0)/nrow(comboLinesOrdered))

# Get 'extreme' points; over threshold.
minPoints <- editedModelBoundaries['EXTRMIN'][unlist(lapply(rownames(comboLinesOrdered), function(x) {grep(paste0("^",x,"$"), rownames(editedModelBoundaries))})),]
maxPoints <- editedModelBoundaries['EXTRMAX'][unlist(lapply(rownames(comboLinesOrdered), function(x) {grep(paste0("^",x,"$"), rownames(editedModelBoundaries))})),]


pdf(file = "OTUrangelimits_Ordered_all.pdf", width = 5, height = 10)
plot(NULL
     , xlim = c(minGrad,maxGrad)
     , ylim = c(0,100)
     , xlab = "Salinity"
     , ylab = ""
     , bty = "n"
     , yaxt = "n"
     , main = "OTU salinity tolerance ranges across salinity")
for (i in 1:nrow(comboLinesOrdered)) {
  colorTempFull <- getColor(comboLinesOrdered$type[i])
  lines(c(comboLinesOrdered$START[i], comboLinesOrdered$STOP[i])
        , c(sequenceOrderCombo[i],sequenceOrderCombo[i])
        , col = colorTempFull)
  colorTempLess <- rbind(col2rgb(colorTempFull), alpha = 25.5)[,1]/255
  lines(c(minPoints[i], maxPoints[i])
        , c(sequenceOrderCombo[i],sequenceOrderCombo[i])
        , cex = 0.1
        , col = rgb(colorTempLess[1], colorTempLess[2], colorTempLess[3], colorTempLess[4]))
}
dev.off()


pdf(file = "OTUrangelimits_OrderedCondensed.pdf", width = 5, height = 10)
plot(NULL
     , xlim = c(minGrad,maxGrad)
     , ylim = c(0,100)
     , xlab = "Salinity"
     , ylab = ""
     , bty = "n"
     , yaxt = "n"
     , main = "OTU salinity tolerance ranges across salinity")
for (i in 1:nrow(comboLinesOrdered)) {
  colorTempFull <- getColor(comboLinesOrdered$typeSimple[i])
  lines(c(comboLinesOrdered$START[i], comboLinesOrdered$STOP[i])
        , c(sequenceOrderCombo[i],sequenceOrderCombo[i])
        , col = colorTempFull)
  colorTempLess <- rbind(col2rgb(colorTempFull), alpha = 25.5)[,1]/255
  lines(c(minPoints[i], maxPoints[i])
        , c(sequenceOrderCombo[i],sequenceOrderCombo[i])
         , cex = 0.1
         , col = rgb(colorTempLess[1], colorTempLess[2], colorTempLess[3], colorTempLess[4]))
  # points(maxPoints[i], sequenceOrderCombo[i]
  #        , cex = 0.1
  #        , col = c(col2rgb(getColor(comboLinesOrdered$typeSimple[i])),0.1))
}

dev.off()

#### GET RID OF BLOOM #### TBA********************


# Get subset only
loLinesSimp <- editedModelBoundaries[grep(paste0("^",gradientNames[1],"Restricted$"), editedModelBoundaries$type),]
interLinesSimp <- editedModelBoundaries[grep(paste0("^",gradientNames[2],"Restricted$"), editedModelBoundaries$type),]
hiLinesSimp <- editedModelBoundaries[grep(paste0("^",gradientNames[3],"Restricted$"), editedModelBoundaries$type),]

# Sorting lo and hi lines are easy
loLinesSimp <- loLinesSimp[order(loLinesSimp$STOP, decreasing = FALSE), ]
hiLinesSimp <- hiLinesSimp[order(hiLinesSimp$START, decreasing = FALSE),]

# But I want inter lines to be sorted by 'middle' I think
sizeForInterSimp <- c()
midForInterSimp <- c()
for (i in 1:nrow(interLinesSimp)) {
  size <- (interLinesSimp$STOP[i]-interLinesSimp$START[i])
  mid <- mean(c(interLinesSimp$STOP[i],interLinesSimp$START[i]))
  sizeForInterSimp <- c(sizeForInterSimp, size)
  midForInterSimp <- c(midForInterSimp, mid)
}

interLinesSimp <- interLinesSimp[order(midForInterSimp,sizeForInterSimp),]

# Combine fresh, brackish, marine
comboLinesOrderedSimp <- rbind(loLinesSimp,interLinesSimp,hiLinesSimp)
# Get rid of NAs from fresh and marine dataset
comboLinesOrderedSimp <- comboLinesOrderedSimp[grep("[0-9]",comboLinesOrderedSimp$STOP),]
comboLinesOrderedSimp <- comboLinesOrderedSimp[grep("[0-9]",comboLinesOrderedSimp$START),]
# Get y intersect
sequenceOrderComboSimp <- seq(0,100,by = (100-0)/nrow(comboLinesOrderedSimp))

# Get 'extreme' points; over threshold.
minPointsSimp <- editedModelBoundaries['EXTRMIN'][unlist(lapply(rownames(comboLinesOrderedSimp), function(x) {grep(paste0("^",x,"$"), rownames(editedModelBoundaries))})),]
maxPointsSimp <- editedModelBoundaries['EXTRMAX'][unlist(lapply(rownames(comboLinesOrderedSimp), function(x) {grep(paste0("^",x,"$"), rownames(editedModelBoundaries))})),]



pdf(file = "OTUrangelimits_Ordered_PureOnly.pdf", width = 5, height = 10)
plot(NULL
     , xlim = c(minGrad,maxGrad)
     , ylim = c(0,100)
     , xlab = "Salinity"
     , ylab = ""
     , bty = "n"
     , yaxt = "n"
     , main = "OTU salinity tolerance ranges across salinity")
for (i in 1:nrow(comboLinesOrderedSimp)) {
  colorTempFullSimp <- getColor(comboLinesOrderedSimp$type[i])
  lines(c(comboLinesOrderedSimp$START[i], comboLinesOrderedSimp$STOP[i])
        , c(sequenceOrderComboSimp[i],sequenceOrderComboSimp[i])
        , col = colorTempFullSimp)
  colorTempLessSimp <- rbind(col2rgb(colorTempFullSimp), alpha = 25.5)[,1]/255
  lines(c(minPointsSimp[i], maxPointsSimp[i])
        , c(sequenceOrderComboSimp[i],sequenceOrderComboSimp[i])
        , cex = 0.1
        , col = rgb(colorTempLessSimp[1], colorTempLessSimp[2], colorTempLessSimp[3], colorTempLessSimp[4]))
}
dev.off()



#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~
####### BETA PLOTS OF TYPES ##########
# library(vegan)
# library(ape)
# library(GUniFrac)
# library(picante)
# library(MASS)
# 
# # OTU table by type
# otuType <- read.delim(paste0(typetablePWD)
#                       , header = TRUE
#                       , row.names = 1
#                       , stringsAsFactors = FALSE)
# 
# colnames(otuType) <- gsub('X', '', colnames(otuType))
# 
# # Phylogeny
# phyloTree <- read.tree(paste0(treePWD))
# if (!is.rooted(phyloTree)){
#   outgrouptree <- phyloTree$tip.label[length(phyloTree$tip.label)]
#   phyloTree.rooted <- root(phyloTree, outgrouptree, resolve.root = TRUE)
# } else {
#   phyloTree.rooted <- phyloTree
# }
# 
# otuType.allpresent <- otuType[,(colnames(otuType) %in% phyloTree.rooted$tip.label)]
# 
# # Unifrac and PCOa plot
# 
# ##_______ TESTING BEGIN
# # Try something fancy-- iterate through each combo and calculate Unifrac with multiple rarefications
# # Rarefy does it randomly each time-- good.
# iterations <- 10
# distanceUnifracRare <- matrix(nrow = nrow(otuType.allpresent), ncol = nrow(otuType.allpresent))
# for (i in 1:nrow(otuType.allpresent)) {
#   # print(rownames(otuType.allpresent)[i]) #TBD
#   for (j in i:nrow(otuType.allpresent)) {
#     if (i == j) {
#       distanceUnifracRare[i,j] <- 0
#     } else {
#       print(paste0("Doing: ",rownames(otuType.allpresent)[i], "vs ", rownames(otuType.allpresent)[j])) #TBD
#       tempDistances <- c()
#       for (k in 1:iterations){
#         print(paste0("Iteration: ",k, " out of ",iterations," iterations"))
#         tempRare <- Rarefy(otuType.allpresent[c(i,j),])
#         tempUnifrac <- GUniFrac(tempRare$otu.tab.rff, phyloTree.rooted)
#         tempDistances[length(tempDistances)+1] <- tempUnifrac$unifracs[2,1,'d_UW']
#         # print(tempDistances) #TBD
#       }
#       # print(mean(tempDistances))
#       distanceUnifracRare[i,j] <- mean(tempDistances)
#       distanceUnifracRare[j,i] <- mean(tempDistances)
#     }
#   }
#   # print('--') #TBD
# }
# rownames(distanceUnifracRare) <- rownames(otuType.allpresent)
# colnames(distanceUnifracRare) <- rownames(otuType.allpresent)
# unweightedUnifracs <- distanceUnifracRare
# 
# #_______ TESTING END
# 
# 
# 
# # unifracOTUType <- GUniFrac(otuType.allpresent, phyloTree.rooted)
# # unweightedUnifracs <- unifracOTUType$unifracs[,,'d_UW']
# 
# otuType.UWUF <- pcoa(unweightedUnifracs)
# PC1value <- round(otuType.UWUF$values$Relative_eig[1]*100)
# PC2value <- round(otuType.UWUF$values$Relative_eig[2]*100)
# 
# 
# # Rarefy all types
# # otuType.allpresent.rare <- Rarefy(otuType.allpresent)
# # otuType.allpresent.final <- otuType.allpresent.rare$otu.tab.rff
# # rowSums(otuType.allpresent)
# 
# 
# 
# # For Legend
# 
# # Colours
# 
# # Old random colors
# # nongreycolors <- colors()[-grep('gray|grey|white',colors())]
# # Specific colors
# typesforApprovedColorsBeta <- c(paste0(gradientNames[1])
#                                 , paste0('bloom',gradientNames[1])
#                                 , paste0("half",gradientNames[1])
#                                 , paste0(gradientNames[2])
#                                 , paste0('bloom',gradientNames[2])
#                                 , paste0('lo',gradientNames[2])
#                                 , paste0('half',gradientNames[2])
#                                 , paste0('hi',gradientNames[2])
#                                 , paste0('half',gradientNames[3])
#                                 , paste0(gradientNames[3])
#                                 , paste0('bloom',gradientNames[3])
#                                 ,"noclass"
#                                 , paste0('inv',gradientNames[2])
# )
# colorsApprovedBeta <- c("blue","yellow","cyan","purple","lightgreen","magenta","lightblue","pink","orange","red", "green","black","gray")
# 
# correctedOrder <- otuType.UWUF$vectors[unlist(lapply(typesforApprovedColorsBeta, function(x) {
#   grep(paste0("^",x,"$"), rownames(otuType.UWUF$vectors))
# })),]
# 
# firstaxis <- correctedOrder[,1]
# secondaxis <- correctedOrder[,2]
# 
# labelsType <- rownames(correctedOrder)
# # colorlabels <- sample(nongreycolors, length(labelsType))
# colorlabels <- unlist(lapply(labelsType, function(x) {
#   colorsApprovedBeta[grep(paste0("^",x,"$"),typesforApprovedColorsBeta)]
# }))
# 
# 
# tempList <- otuType.allpresent[1,]
# # Calculate PD value for each group
# PD.values <- sapply(seq(1,nrow(otuType.allpresent)), function(x) {
#   tempList <- otuType.allpresent[x,]
#   pd(tempList, phyloTree.rooted)
# } )
# 
# # Set names and make scaled so that mean value = 1
# colnames(PD.values) <- rownames(otuType.allpresent)
# scaleFactor <- 2/(mean(unlist(PD.values[1,])))
# PD.values.scaled <- PD.values
# PD.values.scaled[1,] <- (unlist(PD.values[1,])*scaleFactor+1)
# # Make PD.values in correct order
# PD.values.scaled <- PD.values.scaled[,unlist(lapply(typesforApprovedColorsBeta, function(x) {
#   grep(paste0("^",x,"$"), colnames(PD.values.scaled))
# }))]
# 
# 
# # PLOT
# jpeg(paste0("./",output,"/UWUnifrac_byTypePCOA.jpg"), width = 700, height = 500)
# par(fig=c(0,0.67,0,1))
# plot(firstaxis, secondaxis
#      , bg = colorlabels
#      , xlab = paste0('PCO1 (',PC1value,'% of variation)')
#      , ylab = paste0('PCO2 (',PC2value,'% of variation)')
#      , main = "Unweighted Unifrac distances between 'Types'"
#      , pch = 21
#      , col = "black"
#      , cex = unlist(PD.values.scaled[1,]))
# par(fig=c(0.5,1,0,1), new = TRUE)
# plot(0,0
#      , bty = "n"
#      , xaxt = "n"
#      , yaxt = "n"
#      , xlab = ""
#      , ylab = ""
#      , pch = "")
# legend("center"
#        , legend = labelsType
#        , pch = 19
#        , col = colorlabels
#        , bty = "n")
# dev.off()
# 
# ##------------------
# # NMDS plot
# NMDS <- isoMDS(unweightedUnifracs, y = cmdscale(unweightedUnifracs, 2), k= 2)
# 
# # Correct order
# correctedOrder.NMDS <- NMDS$points[unlist(lapply(typesforApprovedColorsBeta, function(x) {
#   grep(paste0("^",x,"$"), rownames(NMDS$points))
# })),]
# 
# labelsTypeNMDS <- rownames(correctedOrder.NMDS)
# 
# jpeg(paste0("./",output,"/UWUnifrac_byTypeNMDS.jpg"), width = 700, height = 500)
# par(fig=c(0,0.67,0,1))
# plot(correctedOrder.NMDS
#      , pch = 19
#      , cex = unlist(PD.values.scaled[1,])
#      , col = colorlabels
#      , xlab = "NMDS1"
#      , ylab = "NMDS2"
#      , sub = paste0("Stress: ",round(NMDS$stress/100,2)))
# par(fig=c(0.5,1,0,1), new = TRUE)
# plot(0,0
#      , bty = "n"
#      , xaxt = "n"
#      , yaxt = "n"
#      , xlab = ""
#      , ylab = ""
#      , pch = "")
# legend("center"
#        , legend = labelsTypeNMDS
#        , pch = 19
#        , col = colorlabels)
# dev.off()

# detach("package:picante", unload=TRUE)
# detach("package:GUniFrac", unload = TRUE)
# detach("package:ape", unload = TRUE)
# detach("package:vegan", unload=TRUE)
# 




