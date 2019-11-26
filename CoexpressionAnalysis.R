getwd()
workingDir <- "./Coexpression Network/Medicago-truncatula-salt-rhizobia-coexpression-network/"
setwd(workingDir)

library(WGCNA)
options(stringsAsFactors = FALSE)

unfilteredData <- read.csv("All_Unfiltered_50,443 genes.csv", stringsAsFactors = FALSE) #LiverFemale3600.csv")#  /All_Unfiltered_50,443 genes.csv"
#unfilteredData <- read.table(file = "PRJNA524006.ose2-lmin50-mm2.count.tsv", sep = '\t', header = TRUE)
dim(unfilteredData)
names(unfilteredData)

datExpr0 <- as.data.frame(t(unfilteredData))#[, -c(1)]))
colnames(datExpr0) <- datExpr0[1, ]
#datExpr0 = datExpr0[-1, ] #for rna seq from jeanne, only have 1 of these
datExpr0 = datExpr0[-1, ] 

datExpr0[] <- lapply(datExpr0, function(x) as.numeric(as.character(x)))

normalize <- function(x) {
        return ((x - min(x)) / (max(x) - min(x)))
}
datExpr0_normalized <- normalize(datExpr0)

sampleTree <- hclust(dist(datExpr0_normalized), method = "average")

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr0_normalized, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


net <- blockwiseModules(datExpr0_normalized, power = 6, maxBlockSize = 10000, 
                        TOMType = "unsigned", minModuleSize = 200,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "baseData",
                       verbose = 3) #minModuleSize used to be 30

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs;
geneTree <- net$dendrograms[[1]];
save(net, datExpr0, datExpr0_normalized, MEs, moduleLabels, moduleColors, geneTree,
     file = "toLoad.RData")

load("toLoad.RData") 
load("baseData-block.1.RData")
load("baseData-block.2.RData")
load("baseData-block.3.RData")
load("baseData-block.4.RData")


nGenes <- ncol(datExpr0_normalized)
nSamples <- nrow(datExpr0_normalized)

# Open a graphical window
sizeGrWindow(6,6)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(net$dendrograms[[2]], moduleColors[net$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 3
plotDendroAndColors(net$dendrograms[[3]], moduleColors[net$blockGenes[[3]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 3", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 4
plotDendroAndColors(net$dendrograms[[4]], moduleColors[net$blockGenes[[4]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 4", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


geneNames <- colnames(datExpr0_normalized)

for (color in moduleColors)
{
        tempList <- names(datExpr0_normalized)[moduleColors==color]
        write.csv(tempList, paste(color,".csv"))
        
}



#plotNetworkHeatmap(datExpr=datExpr0_normalized, 
#                   plotGenes=colnames(datExpr0_normalized),
#                   power=6,
#                   networkType="unsigned",
#                   main="Network Heatmap")


#restGenes <- (moduleColors != "grey")
#diss1 <- 1-TOMsimilarityFromExpr(datExpr0_normalized[,restGenes], power = 6)
