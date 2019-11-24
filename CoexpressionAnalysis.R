getwd()
workingDir <- "./Coexpression Network/Medicago-truncatula-salt-rhizobia-coexpression-network/"
setwd(workingDir)

library(WGCNA)
options(stringsAsFactors = FALSE)

# unfilteredData <- read.csv("PRJNA524006.ose2-lmin50-mm2.count.tsv", stringsAsFactors = FALSE) #LiverFemale3600.csv")#  /All_Unfiltered_50,443 genes.csv"
unfilteredData <- read.table(file = "PRJNA524006.ose2-lmin50-mm2.count.tsv", sep = '\t', header = TRUE)
dim(unfilteredData)
names(unfilteredData)

datExpr0 <- as.data.frame(t(unfilteredData))#[, -c(1)]))
colnames(datExpr0) <- datExpr0[1, ]
datExpr0 = datExpr0[-1, ] #for rna seq from jeanne, only have 1 of these
datExpr0 = datExpr0[-1, ] 

datExpr0[] <- lapply(datExpr0, function(x) as.numeric(as.character(x)))

sampleTree <- hclust(dist(datExpr0), method = "average")

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Choose a set of soft-thre   sholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
 # Plot the results:
sizeGr .   Window(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


net <- blockwiseModules(datExpr0, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs;
geneTree <- net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "try1.RData")