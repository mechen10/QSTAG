#!usr/Bin/Rscript

library(optparse)

######################## OPTPARSE ###########################

option_list = list(
  make_option(c("-d","--distanceMatrixFP"), type = "character"
              , help = "distance matrix from Stilian's script", default = "distance_matrix.tsv"),
  make_option(c("-m","--modelBoundariesFP"), type = "character"
                , help = "model Boundaries output from binning script", default = "modelBoundaries.txt"),
  make_option(c("-N", "--gradientNamesFP"), type = "character"
              , help = "gradient file from binnin script", default = "gradient.txt")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

distanceMatrixFP <- opt$distanceMatrixFP
modelBoundariesFP <- opt$modelBoundariesFP
gradientFP <- opt$gradientNamesFP




##################################################################
# Temp 

setwd("/Users/melissachen/Documents/BETATEST_13march2017/16sFraser")

distanceMatrixFP  <- "distance_matrix.tsv"
modelBoundariesFP <- "modelBoundaries_type.txt"
gradientFP <- "gradient.txt"
uniqueOTUs <- "OTUs.txt"



######################## Loading Files ###########################

dir.create("BETAPLOTS")

# Load distance matrix
dm <- read.delim(paste0(distanceMatrixFP)
           , header = TRUE
           , stringsAsFactors = FALSE
           , row.names = 1
           )
# dm[,1]
# which(duplicated(dm[,1]))
# unique(colnames(dm))
# match(dm[,1], colnames(dm))
# 
# dm[2016,1]
# grep(paste0("^",dm[2017,1],"$"), colnames(dm))
# colnames(dm)[2018]
# 
# grep(paste0("^",dm[8,1],"$"), colnames(dm))
# dm[,grep(paste0("^",dm[1,1],"$"), colnames(dm))]


# rownames(dm) <- gsub("-",".", rownames(dm), fixed = TRUE)
# colnames(dm) <- gsub("-",".", colnames(dm), fixed = TRUE)

# Load modelBoundaries

modelBoundaries <- read.delim(paste0(modelBoundariesFP)
                              , header = TRUE
                              , row.names = 1
                              , stringsAsFactors = FALSE
                              )
# Change the 'bloom' to 'bloomfresh'
for (i in 1:nrow(modelBoundaries)) {
  if (modelBoundaries$type[i] == "noclass") {
    if (modelBoundaries$bloom[i] != "No") {
      modelBoundaries$type[i] = paste0("bloom",modelBoundaries$bloom[i])
    }
  }
}

# Load gradientNames
gradientNames <- read.csv(paste0(gradientFP)
                            , header = FALSE
                          , stringsAsFactors = FALSE
                            )

uniqueOTUList <- read.delim(paste0(uniqueOTUs))


######################## PCOA ###########################
require(MASS)
require(ape)

# Make PCOA
dm.PCOA <- pcoa(dm)
PCAvar1 <- round(dm.PCOA$values$Relative_eig[1]*100,2)
PCAvar2 <- round(dm.PCOA$values$Relative_eig[2]*100,2)
PCAvar3 <- round(dm.PCOA$values$Relative_eig[3]*100,2)
dm.PCOA.twoAxis <- dm.PCOA$vectors[,1:3]

dm.PCOA.twoAxis.correctOrder <- dm.PCOA.twoAxis[unlist(lapply(rownames(modelBoundaries), function(x) {
  grep(paste0("^",x,"$"), rownames(dm.PCOA.twoAxis))
})),]

rownames(modelBoundaries) %in% rownames(dm.PCOA.twoAxis)

# # Make NMDS
# dm.NMDS <- isoMDS(as.matrix(dm), y = cmdscale(dm, 2), k= 2)
# dm.NMDS.correctOrder <- dm.NMDS$points[unlist(lapply(rownames(modelBoundaries), function(x) {
#   grep(paste0("^",x,"$"), rownames(dm.NMDS$points))
# }))]

# Make order of colours
modelBoundaries$type <- factor(modelBoundaries$type, levels = c(paste0("bloom",gradientNames[1])
                                   , gradientNames[1]
                                   , paste0("half",gradientNames[1])
                                   , paste0("lo", gradientNames[2])
                                   , gradientNames[2]
                                   , paste0("half", gradientNames[2])
                                   , paste0("hi", gradientNames[2])
                                   , paste0("bloom", gradientNames[2])
                                   , paste0("half", gradientNames[3])
                                   , gradientNames[3]
                                   , paste0("bloom", gradientNames[3])
                                   , "noclass"))
listAllTypes <- as.vector(unlist(c(paste0("bloom",gradientNames[1])
                  , gradientNames[1]
                  , paste0("half",gradientNames[1])
                  , paste0("lo", gradientNames[2])
                  , gradientNames[2]
                  , paste0("half", gradientNames[2])
                  , paste0("hi", gradientNames[2])
                  , paste0("bloom", gradientNames[2])
                  , paste0("half", gradientNames[3])
                  , gradientNames[3]
                  , paste0("bloom", gradientNames[3])
                  , "noclass")))

orderedColorsAll <- cbind(v1 = c(col2rgb("darkblue"), 0.25)
                   , v2 = c(col2rgb("blue"), 0.25)
                   , v3 = c(col2rgb("lightblue"), 0.25)
                   , v4 = c(col2rgb("cyan"), 0.25)
                   , v5 = c(col2rgb("purple"), 0.25)
                   , v6 = c(col2rgb("pink"), 0.25)
                   , v7 = c(col2rgb("magenta"), 0.25)
                   , v8 = c(col2rgb("green"), 0.25)
                   , v9 = c(col2rgb("orange"), 0.25)
                   , v10 = c(col2rgb("red"), 0.25)
                   , v11 = c(col2rgb("darkred"), 0.25)
                   , v12 = c(col2rgb("black"), 0.00001))

orderedColorsCondensed <- cbind(v1 = c(col2rgb("blue"), 0.25)
                          , v2 = c(col2rgb("blue"), 0.25)
                          , v3 = c(col2rgb("blue"), 0.25)
                          , v4 = c(col2rgb("purple"), 0.25)
                          , v5 = c(col2rgb("purple"), 0.25)
                          , v6 = c(col2rgb("purple"), 0.25)
                          , v7 = c(col2rgb("purple"), 0.25)
                          , v8 = c(col2rgb("purple"), 0.25)
                          , v9 = c(col2rgb("red"), 0.25)
                          , v10 = c(col2rgb("red"), 0.25)
                          , v11 = c(col2rgb("red"), 0.25)
                          , v12 = c(col2rgb("white"), 0))
# orderedColorsAll <- cbind("0 0 139 0.1"
#                           , "0 0 255 0.1"
#                           , "173 216 230 0.1"
#                           , "0 255 255 0.1"
#                           , "160 32 240 0.1"
#                           , "255 192 203 0.1"
#                           , "255 0 255 0.1"
#                           , "0 255 0 0.1"
#                           , "255 165 0 0.1"
#                           , "255 0 0 0.1"
#                           , "139 0 0 0.1"
#                           , "0 0 0 0.1")
colnames(orderedColorsAll) <- as.character(listAllTypes)
colnames(orderedColorsCondensed) <- as.character(listAllTypes)

# make nice dataframe
colorsKeep <- data.frame(row.names = c("r","g","b","alpha"))
# colorsKeep <- data.frame(row.names = c("rgb"))
for (i in 1:ncol(orderedColorsAll)) {
  if (colnames(orderedColorsAll)[i] %in% levels(modelBoundaries$type)) {
    colorsKeep <- cbind(colorsKeep,orderedColorsAll[,i])
  }
}

colorsKeepCondensed <- data.frame(row.names = c("r","g","b","alpha"))
# colorsKeep <- data.frame(row.names = c("rgb"))
for (i in 1:ncol(orderedColorsCondensed)) {
  if (colnames(orderedColorsCondensed)[i] %in% levels(modelBoundaries$type)) {
    colorsKeepCondensed <- cbind(colorsKeepCondensed,orderedColorsCondensed[,i])
  }
}


colorsKeepFinal <- unlist(lapply(1:ncol(colorsKeep), function(x) {
  rgb(colorsKeep[1,x]/255, colorsKeep[2,x]/255, colorsKeep[3,x]/255, colorsKeep[4,x])
}))

colorsKeepFinalCondensed <- unlist(lapply(1:ncol(colorsKeepCondensed), function(x) {
  rgb(colorsKeepCondensed[1,x]/255, colorsKeepCondensed[2,x]/255, colorsKeepCondensed[3,x]/255, colorsKeepCondensed[4,x])
}))
    
jpeg("./BETAPLOTS/Beta_ALL_axis1and2.jpeg")
plot(dm.PCOA.twoAxis.correctOrder[,1] ~ dm.PCOA.twoAxis.correctOrder[,2]
     , cex = 0.5
     , pch = 19
     , xlab = paste0(PCAvar2, "% of variation (Axis2)")
     , ylab = paste0(PCAvar1, "% of variation (Axis1)")
     , main = "PCOA plot of phylogenetic distance between OTUs (1&2)"
     , col = colorsKeepFinalCondensed[modelBoundaries$type])
dev.off()

jpeg("./BETAPLOTS/Beta_ALL_axis2and3.jpeg")
plot(dm.PCOA.twoAxis.correctOrder[,2] ~ dm.PCOA.twoAxis.correctOrder[,3]
     , cex = 0.5
     , pch = 19
     , xlab = paste0(PCAvar3, "% of variation (Axis3)")
     , ylab = paste0(PCAvar2, "% of variation (Axis2)")
     , main = "PCOA plot of phylogenetic distance between OTUs (2&3)"
     , col = colorsKeepFinalCondensed[modelBoundaries$type])
dev.off()

# quartz("Axis 1,3",5,5)
# plot(dm.PCOA.twoAxis.correctOrder[,1] ~ dm.PCOA.twoAxis.correctOrder[,3]
#      , cex = 0.5
#      , pch = 19
#      , xlab = paste0(PCAvar3, "% of variation")
#      , ylab = paste0(PCAvar1, "% of variation")
#      , col = colorsKeepFinal[modelBoundaries$type])
# 

# quartz(,5,5)
# plot(dm.PCOA.twoAxis.correctOrder[,1] ~ dm.PCOA.twoAxis.correctOrder[,2]
#      , cex = 0.5
#      , pch = 19
#      , col = colorsKeepFinalCondensed[modelBoundaries$type])

position <- c()
for (i in 1:length(dm.PCOA.twoAxis.correctOrder[,1])) {
  if (dm.PCOA.twoAxis.correctOrder[i,1] <= -3) {
    position <- c(position, i)
  }
}

tempdm.noarch <- dm.PCOA.twoAxis.correctOrder[-position,]


jpeg("./BETAPLOTS/Beta_NoArch_axis2and3.jpeg")
plot(tempdm.noarch[,2] ~ tempdm.noarch[,3]
     , cex = 0.5
     , pch = 19
     , xlab = paste0(PCAvar3, "% of variation (Axis3)")
     , ylab = paste0(PCAvar2, "% of variation (Axis2)")
     , main = "PCOA plot of phylogenetic distance between OTUs, noarch (2&3)"
     , col = colorsKeepFinalCondensed[modelBoundaries$type[-position]])
dev.off()

jpeg("./BETAPLOTS/Beta_NoArch_axis1and2.jpeg")
plot(tempdm.noarch[,1] ~ tempdm.noarch[,2]
     , cex = 0.5
     , pch = 19
     , xlab = paste0(PCAvar2, "% of variation (Axis2)")
     , ylab = paste0(PCAvar1, "% of variation (Axis3)")
     , main = "PCOA plot of phylogenetic distance between OTUs, noarch (1&2)"
     , col = colorsKeepFinalCondensed[modelBoundaries$type[-position]])
dev.off()

################ NMDS ###################
# # 
# # Make NMDS
# dm.MDS <- cmdscale(dm, 2)
# dm.MDS.correctOrder <- dm.MDS[unlist(lapply(rownames(modelBoundaries), function(x) {
#   grep(paste0("^",x,"$"), rownames(dm.MDS))
# })),]
# 
# quartz(,5,5)
# plot(dm.MDS.correctOrder
#      , cex = 0.5
#      , pch = 19
#      , col = colorsKeepFinalCondensed[modelBoundaries$type])
# 
# 
# 
# # dm.NMDS <- isoMDS(as.matrix(dm), y = cmdscale(dm, 2), k= 2)
# dm.NMDS.correctOrder <- dm.NMDS$points[unlist(lapply(rownames(modelBoundaries), function(x) {
#   grep(paste0("^",x,"$"), rownames(dm.NMDS$points))
# }))]
# 
