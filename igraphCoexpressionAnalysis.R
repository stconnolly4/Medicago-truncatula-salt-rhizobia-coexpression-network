getwd()
workingDir <- "./Coexpression Network/Medicago-truncatula-salt-rhizobia-coexpression-network/"
setwd(workingDir)

library(igraph)

unfilteredData <- read.csv("All_Unfiltered_50,443 genes.csv", stringsAsFactors = FALSE) #LiverFemale3600.csv")#  /All_Unfiltered_50,443 genes.csv"
dim(unfilteredData)
names(unfilteredData)

datExpr0 <- as.data.frame(t(unfilteredData))#[, -c(1)]))
colnames(datExpr0) <- datExpr0[1, ]
datExpr0 = datExpr0[-1, ] 

datExpr0[] <- lapply(datExpr0, function(x) as.numeric(as.character(x)))

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
datExpr0_normalized <- normalize(datExpr0)


require(doParallel)
cores <- makeCluster(detectCores(), type='PSOCK')

cl <- makeCluster(getOption('cl.cores', cores))
registerDoParallel(cl)
registerDoSEQ()

parLapply(cl, 

g <- graph.adjacency(
  as.matrix(as.dist(cor(datExpr0_normalized, method="pearson"))),
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)

)
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)

# Colour negative correlation edges as blue
E(g)[which(E(g)$weight<0)]$color <- "darkblue"

# Colour positive correlation edges as red
E(g)[which(E(g)$weight>0)]$color <- "darkred"

# Convert edge weights to absolute values
E(g)$weight <- abs(E(g)$weight)


g <- delete_edges(g, E(g)[which(E(g)$weight<0.8)])

# Assign names to the graph vertices (optional)
V(g)$name <- V(g)$name

# Change shape of graph vertices
V(g)$shape <- "sphere"

# Change colour of graph vertices
V(g)$color <- "skyblue"

# Change colour of vertex frames
V(g)$vertex.frame.color <- "white"

# Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
# Multiply scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(datExpr0_normalized, 1, mean)) + 1.0) * 10

# Amplify or decrease the width of the edges
edgeweights <- E(g)$weight * 2.0

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(g, algorithm="prim")

# Plot the tree object
plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="My first graph"
)
