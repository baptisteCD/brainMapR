#' Manhattan plot of a brain association map
#'
#' This function takes in raw brain association maps, formats it and writes plot
#'
#'
#' @param path path (folder) to the raw brain association maps (outputs of OSCA)
#' @param bwasFile name of the brain association map
#' @param yMax maximum of the y-axis
#' @param phenotypeLabel Label of the phenotype - for plotting
#' @param signifThreshold significance threshold
#' @return Manhattan plot as well as formatted brain association map (full and subsetted for significant vertices)
#' @import plyr png qqman readr Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
BrainMapManhattanPlot<-function(path , bwasFile , yMax, phenotypeLabel, signifThreshold){

BWASsignif=NULL
chi2All=NULL
pAll=NULL
lambda=NULL
signifLog10= -log( signifThreshold, base=10 )

# Initialise plot
png(paste0(path, "/Manhathan_", bwasFile, "_simple.png"), width = 50, height = 16, units = "cm", res = 400, type="cairo")
plot.new()
layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,6,7,8), nrow =1 , ncol = 20, byrow = F))

BWASsumstatfilePath = paste0(path ,  bwasFile)
BWASsumstatfile<-  read.table(BWASsumstatfilePath, header=T, stringsAsFactors = F)
BWASsumstatfile$log10p<- -log10(BWASsumstatfile$p)

for (mod in c( "thickness","area", "thick", "LogJacs")){
  for (hemi in c("lh", "rh")){
BWASsumstat=NULL

  if (mod %in% c("area", "thickness")){
BWASsumstatHemi<-formatBWAScortical(BWASsumstat = BWASsumstatfile,hemi = hemi, mod=mod)
BWASsumstat=rbind(BWASsumstat,BWASsumstatHemi)
  }

   if (mod %in% c("LogJacs", "thick")){
BWASsumstat<-formatBWASsubcortical(BWASsumstat = BWASsumstatfile)
BWASsumstat=BWASsumstat[grep(BWASsumstat$Probe, pattern = mod),]
if (hemi == "lh"){
  BWASsumstat=BWASsumstat[which(BWASsumstat$subcv<30),]
} else{
  BWASsumstat=BWASsumstat[which(BWASsumstat$subcv>30),]
  BWASsumstat$xax=1:length(BWASsumstat$ProbeID)
}
}

# Get all chi2 and pvalues
BWASsumstat$CHI2=(BWASsumstat$b/BWASsumstat$se)**2

# Set color
laab<-unique(BWASsumstat$LABEL)
colvect<-rep(c("#56B4E9","#0072B2"),length(laab)+1)
BWASsumstat$col<-0
jjj<-1
for (iii in laab){
  BWASsumstat$col[which(BWASsumstat$LABEL==iii)]<-colvect[jjj]
  jjj<-jjj+1
}

# Get rid of unused vertices
if(length(which(is.na(BWASsumstat$log10p)))>0){
BWASsumstat<-BWASsumstat[-which(is.na(BWASsumstat$log10p)),]
}

# Add hemi-moda to text file
BWASsumstat$hemi=hemi
BWASsumstat$mod=mod

# Store all vertices
BWASsignif=plyr::rbind.fill(BWASsignif, BWASsumstat)

# Manhathan plot of association
if (hemi=="lh" & mod=="thickness"){
  par(mar=c(8,4,2,0.3))
plot(BWASsumstat$xax, BWASsumstat$log10p, col=BWASsumstat$col, pch=20, cex=1,  xlab="", xaxt="n", cex.lab=1, ylim=c(0, yMax+1), ylab="-log10(p-value)", bty="n", main="", cex.main=0.7, xaxs="i", yaxs="i", xlim=c(-10, max(BWASsumstat$xax)+10))
title("Cortical thickness \n left", line = -1, cex.main=0.7)
title(paste0(phenotypeLabel), line = 1, cex.main=0.9, adj=0)
} else{

if(mod=="area"){ moda<- "Cortical area"} else if (mod=="thickness"){ moda<- "Cortical thickness" } else if (mod=="LogJacs"){ moda<- "Subcortical area" } else if (mod=="thick"){ moda<- "Subcortical thickness" }

titl<-paste( moda , "\n",ifelse(hemi=="lh", "Left", "Right"))
par(mar=c(8,0.5,2,0.5))
plot(BWASsumstat$xax, BWASsumstat$log10p, col=BWASsumstat$col, pch=20, cex=1,  xlab="", xaxt="n",  yaxt="n", cex.lab=0.9, ylim=c(0, yMax+1), ylab="", bty="n", main="", cex.main=0.7, xaxs="i", yaxs="i", xlim=c(-10, max(BWASsumstat$xax)+10) )
title(titl, line = -1, cex.main=0.7)
}
# add ticks and labels
# Change axis
laab<-unique(BWASsumstat$LABEL)
idchange<-NULL
for (iii in laab){
 idchange<-c(idchange, BWASsumstat$xax[which(BWASsumstat$LABEL==iii)][1])
}
idchange<-c(idchange, max(BWASsumstat$xax))
tickNb<-idchange[1:(length(idchange)-1)]+(idchange[2:length(idchange)]-idchange[1:(length(idchange)-1)])/2
tickNb[-1]=tickNb[-1]-0.002*max(BWASsumstat$xax)

odd <- function(x) x%%2 != 0
axis(1, at = tickNb[odd(which(tickNb>0))], las=2, labels = unique(BWASsumstat$LABELNAME)[odd(which(tickNb>0))], lwd.ticks = 1.5, lwd = 1.5 ,cex.axis=0.8, col.axis = c("#56B4E9"), xaxs="i")
axis(1, at = tickNb[!odd(which(tickNb>0))], las=2, labels = unique(BWASsumstat$LABELNAME)[!odd(which(tickNb>0))], lwd.ticks = 1.5, lwd = 1.5 ,cex.axis=0.8, col.axis = c("#0072B2"), xaxs="i")

abline(h= signifLog10 , col="darkred", lwd=2)
    }
}
dev.off()

write.csv(BWASsignif, paste0(path, "/BWAS_fullSummary_", phenotypeLabel, ".csv" ))
write.csv(BWASsignif[which(BWASsignif$log10p > signifLog10),], paste0(path, "/BWAS_signif_", phenotypeLabel , ".csv" ))

}


#' QQplot of an association map (or any vector of pvalues)
#'
#' This function takes in a vector of pvalues and creates a qqplot (in console)
#'
#'
#' @param pvector Vector of pvalues
#' @param col colour of plot
#' @param add Add to existing plot (T/F)
#' @param ylim limit of y axis
#' @return QQplot
#' @import plyr png qqman readr Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
qqplot <- function(pvector, col = "black", add = F, ylim = NULL){
expectedP <- -log10(ppoints(length(pvector)))
  observedP <- -log10(sort(pvector, decreasing = F))
  if (add == F) {
    plot(x = expectedP, y = observedP, col = col,
         xlab = expression(Expected ~ ~-log[10](italic(p))),
         ylab = expression(Observed ~ ~-log[10](italic(p))),
         pch = 16, cex = 0.7, ylim = NULL)
    abline(0, 1, col = "red")
  }else{
    points(x = expectedP, y = observedP, col = col,
           xlab = expression(Expected ~ ~-log[10](italic(p))),
           ylab = expression(Observed ~ ~-log[10](italic(p))),
           pch = 16, cex = 0.7, ylim = NULL)
  }
}


#' QQplot of several brain association maps
#'
#' This function reads in a list of brain association maps, which have been formatted and created as part of the Manhattan plot creation (this allows inferring more easily the file names)
#'
#'
#' @param path path to working directory
#' @param scenarioList List of scenarios/folders each containing a brain association map
#' @param legendList List of labels for the different scenarios/folders
#' @param variableLabel label of the variable
#' @param colourList List of colours for the different scenarios/folders
#' @return QQplot with several distributions
#' @export
superimposedQQplot=function(path , scenarioList,  legendList, colourList, variableLabel){

jjj=2
png(paste0(path, "QQplotCombined_assocVertices", variableLabel, ".png"), width = 15, height = 15, units = "cm", res=400)
par(mar=c(4,4,2,1))
bwas=read_csv(paste0(path, scenarioList[1], "/BWAS_fullSummary_", variableLabel, ".csv"))
qqplot(bwas$p, col=colourList[1])
 print(round(median(bwas$CHI2,na.rm=T)/qchisq(0.5,df=1),3))
legend(x = 3.5, y = 0.75*max(bwas$log10p), legend = legendList ,pch=20, pt.cex=1.5,  col = colourList)

 for (scenario in scenarioList[-1] ){
bwas=read_csv(paste0(path, scenarioList[jjj], "/BWAS_fullSummary_", variableLabel, ".csv"))
qqplot(bwas$p, col=colourList[jjj], add=T)
 print(print(round(median(bwas$CHI2,na.rm=T)/qchisq(0.5,df=1),3)))
jjj=jjj+1
 }

dev.off()
}


#' Subcortical brain plot
#'
#' This function reads in an annotated brain association map (cluster and coordinates).
#' It produces snapshots of the brain surfaces.
#' The function opens the phenotype file in order to calculate the variance of the phenotype, used to convert association betas into correlation coefficients.
#'
#'
#' @param path path to working directory (input and output folder)
#' @param bwasFile name of the brain association map (without extension)
#' @param pathPhenotypeFile Path to the variable files (e.g. variable.phen) used get variance of the phenotype
#' @return Outside and Inside snapshots of the surfaces.
#' @import plyr png qqman readr Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
plotSubcortical=function(path, bwasFile, pathPhenotypeFile){
  for (moda in c("thick", "LogJacs")){
 for (hemi in c("lh", "rh")){
bwasPlot=read.table( paste0(path, bwasFile , "_", moda, "_clustersAndCoordinates"), header=T)
bwasPlot=bwasPlot[which(bwasPlot$hemi==hemi),]
bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Y=bwasPlot$Y*(-1)

# Open phenotype file
pheno=read.table(paste0(pathPhenotypeFile))

# Transform betas in correlations
bwasPlot$cor=bwasPlot$b/sqrt(var(pheno$V3, na.rm = T))

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(-0.1, 0.1, len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6],RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=ifelse(  bwasPlot$signifVoxel==1, 1.5 ,0.8 )

# Draw plots and save screenshots
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Z",  "X", "Y")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(path, "/BWAS_", bwasFile, "_", hemi, "_", moda , "_clustersAndCoordinates_inside.png"))
rgl.close()

bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Z=bwasPlot$Z*(-1)
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Z",  "X", "Y")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(path, "/BWAS_", bwasFile,"_", hemi, "_", moda , "_clustersAndCoordinates_outside.png"))
rgl.close()

}}
}

#' Better rendition of the subcortical brain plots
#'
#' This function reads in an annotated brain association map (cluster and coordinates).
#' It produces snapshots of the brain surfaces.
#' The function opens the phenotype file in order to calculate the variance of the phenotype, used to convert association betas into correlation coefficients.
#'
#'
#' @param path path to working directory (input and output folder)
#' @param bwasFile name of the brain association map (without extension)
#' @param pathPhenotypeFile Path to the variable files (e.g. variable.phen) used get variance of the phenotype
#' @return Outside and Inside snapshots of the surfaces - the subcortical structures are spread out to faciliate viewing.
#' @import plyr png qqman readr Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
plotSubcortical_flat=function(path, bwasFile, pathPhenotypeFile){

  for (moda in c("thick", "LogJacs")){
 for (hemi in c("lh", "rh")){
bwasPlot=read.table(paste0(path, bwasFile , "_", moda, "_clustersAndCoordinates"), header=T)
bwasPlot=bwasPlot[which(bwasPlot$hemi==hemi),]
bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Y=bwasPlot$Y*(-1)

# Open phenotype file
pheno=read.table(paste0(pathPhenotypeFile))

# Transform betas in correlations
bwasPlot$cor=bwasPlot$b/sqrt(var(pheno$V3, na.rm = T))

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(-0.4, 0.4, len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6],RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=ifelse(  bwasPlot$signifVoxel==1, 2 ,0.8 )

# Change coordinates for flat plotting
if (hemi=="rh"){

bwasPlot$Xf=bwasPlot$X
bwasPlot$Yf=bwasPlot$Y
bwasPlot$Zf=bwasPlot$Z

  # Hippocampus
bwasPlot$Yf[bwasPlot$subcv %in% c(17,53)]=bwasPlot$Yf[bwasPlot$subcv %in% c(17,53)]-17

# Amygdala
bwasPlot$Yf[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Yf[bwasPlot$subcv %in% c(18,54)]-17
bwasPlot$Zf[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Zf[bwasPlot$subcv %in% c(18,54)]+7

# Thalamus
bwasPlot$Zf[bwasPlot$subcv %in% c(10,49)]=bwasPlot$Zf[bwasPlot$subcv %in% c(10,49)]-20

# Caudate
bwasPlot$Yf[bwasPlot$subcv %in% c(11,50)]=bwasPlot$Yf[bwasPlot$subcv %in% c(11,50)]+20

# Accumbens
bwasPlot$Zf[bwasPlot$subcv %in% c(26,58)]=bwasPlot$Zf[bwasPlot$subcv %in% c(26,58)]+12

# Pallidum
bwasPlot$Yf[bwasPlot$subcv %in% c(13,52)]=bwasPlot$Yf[bwasPlot$subcv %in% c(13,52)]-15

# Second set of coordinates
bwasPlot$Xf2=bwasPlot$Xf*(-1)
bwasPlot$Zf2=bwasPlot$Zf*(-1)
bwasPlot$Yf2=bwasPlot$Yf

# Accumbens
bwasPlot$Zf2[bwasPlot$subcv %in% c(26,58)]=bwasPlot$Zf2[bwasPlot$subcv %in% c(26,58)]-7

# Pallidum
bwasPlot$Yf2[bwasPlot$subcv %in% c(13,52)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(13,52)]-8

# Amygdala
bwasPlot$Yf2[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(18,54)]-5

# Hippocampus
bwasPlot$Yf2[bwasPlot$subcv %in% c(17,53)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(17,53)]-5

# Caudate
bwasPlot$Yf2[bwasPlot$subcv %in% c(11,50)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(11,50)]-8
##########################
} else if (hemi == "lh"){
# Second set of coordinates

bwasPlot$Xf=bwasPlot$X
bwasPlot$Yf=bwasPlot$Y
bwasPlot$Zf=bwasPlot$Z

  # Hippocampus
bwasPlot$Yf[bwasPlot$subcv %in% c(17,53)]=bwasPlot$Yf[bwasPlot$subcv %in% c(17,53)]-21

# Amygdala
bwasPlot$Yf[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Yf[bwasPlot$subcv %in% c(18,54)]-21
bwasPlot$Zf[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Zf[bwasPlot$subcv %in% c(18,54)]+7

# Thalamus
bwasPlot$Zf[bwasPlot$subcv %in% c(10,49)]=bwasPlot$Zf[bwasPlot$subcv %in% c(10,49)]-20

# Caudate
bwasPlot$Yf[bwasPlot$subcv %in% c(11,50)]=bwasPlot$Yf[bwasPlot$subcv %in% c(11,50)]+12

# Accumbens
bwasPlot$Zf[bwasPlot$subcv %in% c(26,58)]=bwasPlot$Zf[bwasPlot$subcv %in% c(26,58)]+19

# Pallidum
bwasPlot$Yf[bwasPlot$subcv %in% c(13,52)]=bwasPlot$Yf[bwasPlot$subcv %in% c(13,52)]-22

# Second set of coordinates
bwasPlot$Xf2=bwasPlot$Xf*(-1)
bwasPlot$Zf2=bwasPlot$Zf*(-1)
bwasPlot$Yf2=bwasPlot$Yf

# Caudate
bwasPlot$Yf2[bwasPlot$subcv %in% c(11,50)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(11,50)]+8

# Pallidum
bwasPlot$Yf2[bwasPlot$subcv %in% c(13,52)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(13,52)]+6

# Amygdala
bwasPlot$Yf2[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(18,54)]+3

# Hippocampus
bwasPlot$Yf2[bwasPlot$subcv %in% c(17,53)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(17,53)]+3

# Accumbens
bwasPlot$Zf2[bwasPlot$subcv %in% c(26,58)]=bwasPlot$Zf2[bwasPlot$subcv %in% c(26,58)]+6

# Thalamus
bwasPlot$Zf2[bwasPlot$subcv %in% c(10,49)]=bwasPlot$Zf2[bwasPlot$subcv %in% c(10,49)]+2

}

# Draw plots and save screenshots
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Zf",  "Xf", "Yf")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(path, "/BWAS_", bwasFile, "_", hemi, "_", moda , "_clustersAndCoordinates_inside.png"))
rgl.close()

par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Zf2",  "Xf2", "Yf2")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(path, "/BWAS_", bwasFile,"_", hemi, "_", moda , "_clustersAndCoordinates_outside.png"))
rgl.close()

}}
}

#' Cortical brain plot
#'
#' This function reads in an annotated brain association map (cluster and coordinates).
#' It produces snapshots of the brain surfaces.
#' The function opens the phenotype file in order to calculate the variance of the phenotype, used to convert association betas into correlation coefficients.
#'
#'
#' @param path path to working directory (input and output folder)
#' @param bwasFile name of the brain association map (without extension)
#' @param pathPhenotypeFile Path to the variable files (e.g. variable.phen) used get variance of the phenotype
#' @return Outside and Inside snapshots of the surfaces.
#' @import plyr png qqman readr Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
plotCortical=function(path, bwasFile, pathPhenotypeFile){

 for (moda in c("area", "thickness")){

 for (hemi in c("rh", "lh")){
bwasPlot=read.table(paste0(path, bwasFile , "_", moda, "_clustersAndCoordinates"), header=T)
bwasPlot=bwasPlot[which(bwasPlot$hemi==hemi),]

# Open phenotype file
pheno=read.table(paste0(pathPhenotypeFile))

# Transform betas in correlations
bwasPlot$cor=bwasPlot$b/sqrt(var(pheno$V3, na.rm = T))

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(-0.4, 0.4, len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6], RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=ifelse(bwasPlot$signifVoxel==1, 2, 0.8)

#bwasPlot2=bwasPlot[sample(1:length(bwasPlot$ProbeID), size = 10000),]
#bwasPlot2=rbind(bwasPlot2, bwasPlot[which(bwasPlot$signifVoxel==1 | bwasPlot$inCluster>0 ),])
bwasPlot2=bwasPlot


# Draw plot
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot2[,c( "Y","X", "Z")]), col=bwasPlot2$color, radius = bwasPlot2$radius)
rgl.snapshot(paste0(path, "/BWAS_", bwasFile,"_", hemi, "_", moda , "_clustersAndCoordinates_inside.png"))
rgl.close()

bwasPlot2$X=bwasPlot2$X*(-1)
bwasPlot2$Y=bwasPlot2$Y*(-1)
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot2[,c( "Y","X", "Z")]), col=bwasPlot2$color, radius = bwasPlot2$radius)
rgl.snapshot(paste0(path, "/BWAS_", bwasFile,"_", hemi, "_", moda , "_clustersAndCoordinates_outside.png"))
rgl.close()

}}
}


#' Combine cortical and subcortical figures to create a multi-panel plot
#'
#' This function uses the cortical ans subcortical figures created using the specific functions
#'
#'
#' @param path path to working directory (input and output folder)
#' @param bwasFile Variable name (used to name the clustersAndCoordinates file)
#' @param pathToLegendBar Path to the manually created legend bar
#' @return Outside and Inside snapshots of the surfaces.
#' @import plyr png qqman readr Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
combineCorticalSubcorticalPlots=function(path, bwasFile, pathToLegendBar){

# List of files
ll=NULL
for (moda in c("thickness", "area","thick", "LogJacs") ){
ll=c(ll, paste0(path, "/BWAS_", bwasFile, "_", "lh", "_", moda, "_clustersAndCoordinates_outside.png"))
ll=c(ll, paste0(path, "/BWAS_", bwasFile, "_","lh", "_", moda, "_clustersAndCoordinates_inside.png"))
ll=c(ll, paste0(path, "/BWAS_", bwasFile, "_", "rh", "_", moda, "_clustersAndCoordinates_inside.png"))
ll=c(ll, paste0(path, "/BWAS_", bwasFile, "_", "rh", "_", moda, "_clustersAndCoordinates_outside.png"))
}
ll=c(ll, paste0(pathToLegendBar))

# Crop and convert format
plots2 <- lapply(ll<-ll,function(x){
  if(x!=paste0(pathToLegendBar)){
   img <- as.raster(readPNG(x)[,100:1100,])} else {
    # img <- as.raster(readPNG(x))} else {
     img <- as.raster(readPNG(x)[,,])
    }
    rasterGrob(img, interpolate = T)

})

# Lay and write png
lay <- rbind(c(1,1,3,3,5,5,7,7,9,9,11,11,13,13,15,15,17),
            c(2,2,4,4,6,6,8,8,10,10,12,12,14,14,16,16, 17))

gs=grid.arrange(grobs = plots2, layout_matrix = lay)
ggsave(paste0(path, "/Plots_Combined", bwasFile, ".png"),width=18, height=4, gs)


}


