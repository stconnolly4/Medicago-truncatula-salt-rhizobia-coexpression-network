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

# get adjacency matrix - note, will need to break into smaller blocks for visualization

multiExpr = fixDataStructure(datExpr0_normalized, verbose = 0, indent = 0)


blockwiseIndividualTOMs(
        multiExpr,
        multiWeights = NULL,
        
        # Data checking options
        
        checkMissingData = TRUE,
        
        # Blocking options
        
        blocks = NULL,
        maxBlockSize = 10000,
        blockSizePenaltyPower = 5,
        nPreclusteringCenters = NULL,
        randomSeed = 54321,
        
        # Network construction arguments: correlation options
        
        corType = "pearson",
        maxPOutliers = 1,
        quickCor = 0,
        pearsonFallback = "individual",
        cosineCorrelation = FALSE,
        
        # Adjacency function options
        
        power = 6,
        networkType = "unsigned",
        checkPower = TRUE,
        replaceMissingAdjacencies = FALSE,
        
        # Topological overlap options
        
        TOMType = "unsigned",
        TOMDenom = "min",
        suppressTOMForZeroAdjacencies = FALSE,
        suppressNegativeTOM = FALSE,
        
        # Save individual TOMs? If not, they will be returned in the session.
        
        saveTOMs = TRUE,
        individualTOMFileNames = "individualTOM-Set%s-Block%b.RData",
        
        # General options
        
        nThreads = 0,
        useInternalMatrixAlgebra = FALSE,
        verbose = 2, indent = 0)

load("individualTOM-Set1-Block1.RData")
block1_matrix <- as.matrix(tomDS)
dim(block1_matrix)
write.csv(block1_matrix,file="block1.csv")

rm(list = ls())

load("individualTOM-Set1-Block2.RData")
block2_matrix <- as.matrix(tomDS)
dim(block2_matrix)
write.csv(block1_matrix,file="block2.csv")

rm(list = ls())

load("individualTOM-Set1-Block3.RData")
block3_matrix <- as.matrix(tomDS)
dim(block3_matrix)
write.csv(block3_matrix,file="block3.csv")

rm(list = ls())

load("individualTOM-Set1-Block4.RData")
block4_matrix <- as.matrix(tomDS)
dim(block4_matrix)
write.csv(block4_matrix,file="block4.csv")
