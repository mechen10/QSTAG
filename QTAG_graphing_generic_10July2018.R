#! /usr/bin/Rscript
# Gradient Binning Output- graphing
########################## OPTPARSE ####################################
require(optparse)

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
biomPWD = opt$biom
minthreshold = opt$minthreshold
LOGPWD = opt$LOG
typetablePWD = opt$OTUTable

########################## LOAD FILES ####################################

# Load files
print("Loading files for graphing...")

boundaries <- read.delim(paste0(boundariesPWD), header = FALSE, strip.white = TRUE)
histBoundaries <- as.numeric(boundaries$V1)

logvalues <- read.delim(paste0(LOGPWD), header = FALSE, strip.white = TRUE)
minGrad <- as.character(logvalues[grep("minX", logvalues$V1),])
minGrad <- floor(as.numeric(gsub("minX =","", minGrad)))
maxGrad <- as.character(logvalues[grep("maxY", logvalues$V1),])
maxGrad <- ceiling(as.numeric(gsub("maxY =","", maxGrad)))
metadataPWD <- as.character(logvalues[grep("metadata", logvalues$V1),])
metadataPWD <- gsub("metadata = ","",metadataPWD)

threshold <- as.character(logvalues[grep("^threshold:", logvalues$V1),])
threshold <- as.numeric(gsub("threshold:", "", threshold))

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

colorsGradient <- c("blue" #lorestr
                    ,"darkblue" # lobloom
                    ,"cyan" #lopeak
                    ,"slateblue1" #bracklotoler
                    ,"purple" #brackrestr
                    ,"magenta" # brack bloom
                    ,"pink4" #brack peak
                    ,"violetred" #brackpeakhitoler
                    ,"orange" #hipeak
                    ,"red" #hirestr
                    ,"darkred"#hibloom
                    ,"black" #noclass
                    , "grey" #ubiq
                    ,"white" #invbrack
)
colorsGradientCondensed <- c("blue","purple","red","black")


typesAll <- as.matrix(typesAll)
pdf(file = "taxatypesacrossgradient_all.pdf")
par(mar=c(5.1,4.1,4.1,9), xpd = TRUE)
barplot(typesAll, beside = FALSE
        , col = colorsGradient
        , main = 'Taxa types across gradient'
        , xlab = 'Gradient'
        , ylab = 'Relative Abundance')
legend('right', inset = c(-0.4,0), legend = rev(rownames(typesAll)), pch = 21, col = c("black",rev(colorsGradient)[-1]), pt.bg = rev(colorsGradient), bty="n")
dev.off()


typesCondensed <- as.matrix(typesCondensed)
pdf(file = "taxatypesacrossgradient_condensed.pdf")
par(mar=c(5.1,4.1,4.1,9), xpd = TRUE)
barplot(typesCondensed, beside = FALSE
        , col = colorsGradientCondensed
        , main = 'Taxa types across gradient'
        , xlab = 'Gradient'
        , ylab = 'Relative Abundance')
legend('right', inset = c(-0.4,0), legend = rev(rownames(typesCondensed)), pch = 21, pt.bg = rev(colorsGradientCondensed), col = rev(colorsGradientCondensed), bty="n")
dev.off()


########################## TAXA ABUND (by type) ####################################
system('mkdir TaxaAbund_bytype_all')
system('mkdir TaxaAbund_bytype_condensed')
print("Making taxa abundance plots by type")


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


########################## High ####################################

# Now pull out all boundaries for fresh
allHigh <- c(highOnly$boundaries, highOnly$boundariestwo)
TrueHigh <- c(highOnlyTrue$boundaries)
HalfHigh <- c(highOnlyHalf$boundaries)
BloomHigh <- c(highOnlyBloom$boundaries)
nonBloomHigh <- c(HalfHigh,TrueHigh)

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
        return("violetred")
    } else if (x == paste0(gradientNames[2], "PeakLoToler")) {
        return("slateblue1")
    } else if (x == paste0(gradientNames[1], "Peak")) {
        return("cyan")
    } else if (x == paste0(gradientNames[2], "Peak")) {
        return("pink4")
    } else if (x == paste0(gradientNames[3], "Peak")) {
        return("orange")
    } else if (x == paste0(gradientNames[1], "Bloom")) {
        return("darkblue")
    } else if (x == paste0( gradientNames[2], "Bloom")) {
        return("magenta")
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



