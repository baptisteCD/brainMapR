
#' Identify significant cluster around a (significant) vertex of interest
#'
#' This function reads a formatted brain association map (with X Y Z coordinates).
#' From the target vertex, the algorithm identifies neighboring vertices that are also significant.
#'
#' @param bwas Formatted brain association map
#' @param signifThreshold Significance threshold
#' @param targetVertex Name of the vertex of interest
#' @param nbNearestNeighbours Number of nearest neighbors to add to the cluster at each iterative step
#' @return An updated brain association map with cluster membership
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
significantClusterAroundVertex=function(bwas, signifThreshold, targetVertex, nbNearestNeighbours ){

  if (bwas$p[which(bwas$Probe == targetVertex)]> signifThreshold){
  print(paste("Vertex", targetVertex, "located in", bwas$LABELNAME[which(bwas$Probe==targetVertex)], "not significant" ))
  } else {

bwas$voxelStatus="nonSignificant"
bwas$voxelStatus[which(bwas$p < signifThreshold)]="bwSignificant"
bwas$voxelStatus[which(bwas$Probe == targetVertex)]="truePositive"

# First iter - closest vertices
clost= Rvcg::vcgKDtree( target = as.matrix(bwas[,c("X", "Y", "Z")]) ,  query = as.matrix(bwas[which(bwas$Probe==targetVertex),c("X", "Y", "Z")]), k=nbNearestNeighbours)$index # Get vertices number

# Mark the vertices and exclude NS ones and TP one
bwas$voxelStatus[clost]="belongInCluster"
bwas$voxelStatus[which(bwas$p > signifThreshold)]="nonSignificant"
bwas$voxelStatus[which(bwas$Probe==targetVertex)]="truePositive"
nCluster=length(which(bwas$voxelStatus=="belongInCluster"))

# Loop until all cluster has been identified
nIncrementCluster=nCluster
sizeCluster=0

while(nIncrementCluster>0){
clost=NULL
clost <- Rvcg::vcgKDtree( target = as.matrix(bwas[,c("X", "Y", "Z")]) ,  query = as.matrix(bwas[which(bwas$voxelStatus %in% c("belongInCluster", "truePositive")),c( "X","Y", "Z")]), k=nbNearestNeighbours)$index
bwas$voxelStatus[unique(c(clost))]="belongInCluster"
bwas$voxelStatus[which(bwas$p > signifThreshold)]="nonSignificant"
nIncrementCluster=length(which(bwas$voxelStatus=="belongInCluster"))-sizeCluster
sizeCluster=length(which(bwas$voxelStatus=="belongInCluster"))
}

bwas[,paste0("cluster_", targetVertex)]=ifelse(bwas$voxelStatus %in% c("belongInCluster", "truePositive"), 1 , 0 )
print(paste("Found ", sum(bwas[,paste0("cluster_", targetVertex)]),"significant vertices in cluster including ", targetVertex, "located in", bwas$LABELNAME[which(bwas$Probe==targetVertex)] ))
bwas$voxelStatus=NULL

}

return(bwas)
}


#' Identify and describe all significant clusters from a brain association maps
#'
#' This function reads a formatted brain association map (with X Y Z coordinates).
#'
#'
#' @param outputFolder path to the output/writing directory
#' @param bwas Formatted brain association map
#' @param moda Measurement type ("area" or "thickness" for cortical measurements; "LogJacs" or "thick" for subcortical)
#' @param signifThreshold Significance threshold
#' @param KNearestNeighbours Number of nearest neighbors to add to the cluster at each iterative step
#' @param phenotypeLabel Label of the phenotype for enhanced plotting
#' @return Summary file of significant clusters and their properties
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
clusterIdentification=function(outputFolder, bwas, moda, signifThreshold, KNearestNeighbours, phenotypeLabel){

bwas11=NULL
twoTPinCLuster=fwerProblemIncreaseK=0
bwas$signifVoxel=ifelse(bwas$p < signifThreshold, 1 , 0)

# Depending on modality - loop on hemisphere or subcortical structures
if (moda %in% c("area", "thickness")){ ROIlist=c("lh", "rh")}
if (moda %in% c("LogJacs", "thick")){ ROIlist=c(10, 11, 12, 13, 17, 18, 26, 49, 50, 51, 52, 53, 54, 58)}

# Loop on the different structures and concatenate results
  for (ROI in ROIlist){
  print(paste0("Looking for clusters in ", ROI))
  if (moda %in% c("LogJacs", "thick")){ bwas1=bwas[which(bwas$ROINb==ROI),] }
  if (moda %in% c("area", "thickness")){ bwas1=bwas[which(bwas$hemi==ROI),] }

  # Find significant clusters
  bwas1$inCluster=0
  while (length(which(bwas1$signifVoxel==1 & bwas1$inCluster==0))>0){
  signifVertex=bwas1$Probe[which(bwas1$signifVoxel==1  & bwas1$inCluster==0)][1]
  bwas1=significantClusterAroundVertex(bwas=bwas1, signifThreshold=signifThreshold, targetVertex = signifVertex, nbNearestNeighbours = KNearestNeighbours)
  bwas1$inCluster=rowSums(bwas1[,grep(x = colnames(bwas1), pattern = "cluster_"), drop=FALSE])

    if (length(table(bwas1$inCluster))>2){
 print("WARNING: some vertices have been attributed to several clusters. You may want to visually check the results. In the meantime, the 2 clusters are merged")

    # Merge clusters and resets variable
    TPvertex=bwas1$Probe[which(bwas1$inCluster==2 & bwas1$signifVoxel==1)][1]
    bwas1[which(bwas1[,paste0("cluster_", signifVertex)]==1),paste0("cluster_", TPvertex)]=1
    bwas1[,paste0("cluster_", signifVertex)]=NULL
    fwerProblemIncreaseK=1
    bwas1$inCluster=rowSums(bwas1[,grep(x = colnames(bwas1), pattern = "cluster_"), drop=FALSE])
    }
  }
    # Combine all ROI - hemisphere of subcortical volumes
    bwas11=plyr::rbind.fill(bwas11, bwas1)
  }
    # Extract statistics and information on the full data
    nbClusters=length(grep(x = colnames(bwas11), pattern = "cluster_"))
    maxSizeClusters=max(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_"), drop=FALSE], na.rm = T), na.rm = T)
    medianSizeClusters=median(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_"), drop=FALSE], na.rm = T), na.rm = T)
    minSizeClusters=min(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_"), drop=FALSE], na.rm = T), na.rm = T)

    FWERcluster_ROI=0
    for (iii in grep(x = colnames(bwas11), pattern = "cluster_")){
      if( length(table(bwas11$LABELNAME[which(bwas11[,iii]==1)]))>1 ){ FWERcluster_ROI=1   }   }

     # Format results
    res=as.data.frame(t(c( nbClusters, maxSizeClusters,medianSizeClusters, minSizeClusters, FWERcluster_ROI)))
    colnames(res)=c("NumberClusters", "maxSizeCluster", "medianSizeCluster", "minSizeCluster", "FWERcluster_ROI")
    # Write outputs and bwas table with all the informations (for plots etc..)

    write.table(bwas11, file = paste0(outputFolder, phenotypeLabel , "_", moda, "_clustersAndCoordinates"), quote = F )
    return(res)
}



#' Run all the cluster identification from a raw brain map
#'
#' This function combines, reading and formatting fuctions as well as the clustering ones to perform all the analyses for you.
#'
#'
#' @param inputPath path (folder) to the raw brain association maps (outputs of OSCA)
#' @param bwasFile name of the brain association map (without extension)
#' @param outputFolder path to output folder
#' @param signifThreshold significance threshold
#' @return Summary file of significant clusters and their properties
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices vroom
#' @export
identifyClustersBWAS=function(inputPath, bwasFile, outputFolder, signifThreshold){

# Open BWAS result
bwasResult=vroom(paste0(inputPath, bwasFile), show_col_types = F)
resTot=NULL

# Loop on all modalities
for (moda in c("area", "thickness", "thick", "LogJacs")){
# Open and format bwas summary statistics
if (moda %in% c("thick", "LogJacs")){
bwasFormattedlh=formatBWASsubcortical(BWASsumstat = bwasResult, hemi="lh", mod=moda)
bwasFormattedrh=formatBWASsubcortical(BWASsumstat = bwasResult, hemi="rh", mod=moda)
bwas=rbind(bwasFormattedlh, bwasFormattedrh)
}

if (moda %in% c("thickness", "area")){
# Format left and right hemisphere data
bwasFormattedlh=formatBWAScortical(BWASsumstat = bwasResult, hemi = "lh", mod = moda)
bwasFormattedrh=formatBWAScortical(BWASsumstat = bwasResult, hemi = "rh", mod = moda)
# Rbind left and right
bwas=rbind(bwasFormattedlh, bwasFormattedrh)
}

# Run cluster statistic
res=clusterIdentification(outputFolder =outputFolder, bwas = bwas,  moda = moda, signifThreshold = signifThreshold, KNearestNeighbours = 10, phenotypeLabel = bwasFile)
colnames(res)=paste0(colnames(res), "_", moda)
resTot=c(resTot, res)
}

# Write final output - with all iterations
write.table(resTot, paste0(outputFolder, "/Results_clusterFWER_", bwasFile ,  ".txt" ))


}
