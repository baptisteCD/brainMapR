#' QQplot of an association map (or any vector of pvalues)
#'
#' This function takes in a vector of pvalues and creates a qqplot (in console)
#'
#'
#' @param pvector Vector of pvalues
#' @param col colour of points
#' @param add Add to existing plot (T/F)
#' @return QQplot
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
qqplot <- function(pvector, col = "black", add = F){
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

#' Manhattan plot of a brain association map
#'
#' This function takes in raw brain association maps, formats it and writes plot
#'
#'
#' @param inputPath path (folder) to the raw brain association maps (outputs of OSCA)
#' @param bwasFile name of the brain association map
#' @param yMax maximum of the y-axis in the Manhattan plot
#' @param phenotypeLabel Label of the phenotype - used for plotting
#' @param signifThreshold pvalue significance threshold used to account for multiple testing
#' @param outputPath path where the outputs will be written
#' @param qqPlot True/False: should a qqplot also be created?
#' @return Manhattan plot as well as annotated brain association maps (full and subset of significant vertices)
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices vroom
#' @export
BrainMapAnnotAndManhattanPlot<-function(inputPath , bwasFile , yMax, phenotypeLabel, signifThreshold, outputPath, qqPlot=T){

# Initialise objects
BWASsignif=NULL
signifLog10= -log( signifThreshold, base=10 )

# Initialise plot
png(paste0(outputPath, "Manhathan_", bwasFile, "_simple.png"), width = 50, height = 16, units = "cm", res = 400, type="cairo")
plot.new()
layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,6,7,8), nrow =1 , ncol = 20, byrow = F))

# Open summary statistics
BWASsumstatfilePath = paste0(inputPath ,  bwasFile)
BWASsumstatfile<-  vroom(BWASsumstatfilePath, show_col_types = F)
BWASsumstatfile$log10p<- -log10(BWASsumstatfile$p)

# Loop on modality and hemispheres to annotate summary statistics and plot sections of Manhattan plot
for (mod in c( "thickness","area", "thick", "LogJacs")){
  for (hemi in c("lh", "rh")){
BWASsumstat=NULL

  if (mod %in% c("area", "thickness")){
BWASsumstat<-formatBWAScortical(BWASsumstat = BWASsumstatfile, hemi = hemi, mod=mod)
  }

   if (mod %in% c("LogJacs", "thick")){
BWASsumstat<-formatBWASsubcortical(BWASsumstat = BWASsumstatfile , hemi = hemi, mod=mod)
}

# Get all chi2 value
BWASsumstat$CHI2=(BWASsumstat$b/BWASsumstat$se)**2

# Set colors for Manhattan plot
laab<-unique(BWASsumstat$ROINb)
colvect<-rep(c("#56B4E9","#0072B2"),length(laab)+1)
BWASsumstat$col<-0
jjj<-1
for (iii in laab){
  BWASsumstat$col[which(BWASsumstat$ROINb==iii)]<-colvect[jjj]
  jjj<-jjj+1
}

# Get rid of vertices not included in the analysis
if(length(which(is.na(BWASsumstat$log10p)))>0){
BWASsumstat<-BWASsumstat[-which(is.na(BWASsumstat$log10p)),]
}

# rbind all vertices with annotations
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

# Add ticks and labels for the axes
laab<-unique(BWASsumstat$ROINb)
idchange<-NULL
for (iii in laab){
 idchange<-c(idchange, BWASsumstat$xax[which(BWASsumstat$ROINb==iii)][1])
}
idchange<-c(idchange, max(BWASsumstat$xax))
tickNb<-idchange[1:(length(idchange)-1)]+(idchange[2:length(idchange)]-idchange[1:(length(idchange)-1)])/2
tickNb[-1]=tickNb[-1]-0.002*max(BWASsumstat$xax)

odd <- function(x) x%%2 != 0
axis(1, at = tickNb[odd(which(tickNb>0))], las=2, labels = unique(BWASsumstat$ROIlabel)[odd(which(tickNb>0))], lwd.ticks = 1.5, lwd = 1.5 ,cex.axis=0.8, col.axis = c("#56B4E9"), xaxs="i")
axis(1, at = tickNb[!odd(which(tickNb>0))], las=2, labels = unique(BWASsumstat$ROIlabel)[!odd(which(tickNb>0))], lwd.ticks = 1.5, lwd = 1.5 ,cex.axis=0.8, col.axis = c("#0072B2"), xaxs="i")

# Add significance threshold line
abline(h= signifLog10 , col="darkred", lwd=2)
    }
}
dev.off()

# Write annotated summary statistics
write.csv(BWASsignif, paste0(outputPath, "/BWAS_fullSummary_", bwasFile, ".csv" ))
write.csv(BWASsignif[which(BWASsignif$log10p > signifLog10),], paste0(outputPath, "/BWAS_signif_", bwasFile , ".csv" ))

# Produce QQplot
if (qqPlot==T){
    png(paste0(outputPath, "QQplot_", bwasFile, ".png"), width = 10, height = 10, units = "cm", res = 400, type="cairo")
    qqplot(pvector = BWASsumstat$p )
    dev.off()
}

}

#' QQplot of several brain association maps
#'
#' This function reads in a list of brain association maps, which have been formatted and created as part of the Manhattan plot creation (this allows inferring more easily the file names)
#'
#'
#' @param inputPaths paths to working directories containing the files BWAS_fullSummary_....csv, created by BrainMapManhattanPlot()
#' @param bwasFiles list of brain association maps within each inputPath
#' @param colourList List of colours for the different scenarios/folders
#' @param legendList List of labels for the different scenarios/folders
#' @param outputPath path where the output will be written
#' @param phenotypeLabel Label of the phenotype - used for plotting
#' @return QQplot with several distributions
#' @export
superimposedQQplot=function(inputPaths , bwasFiles, legendList, colourList, phenotypeLabel, outputPath){

# Initialise plot with first scenario
png(paste0(outputPath, "QQplotCombined_assocVertices", phenotypeLabel, ".png"), width = 15, height = 15, units = "cm", res=400)
par(mar=c(4,4,2,1))
if (file.exists(paste0(inputPaths[1], "/BWAS_fullSummary_", bwasFiles[1], ".csv"))){
bwas=vroom(paste0(inputPaths[1], "/BWAS_fullSummary_", bwasFiles[1], ".csv"), show_col_types = F)} else {
    print(paste0(inputPaths[1], "/BWAS_fullSummary_", bwasFiles[1], ".csv", "  not found, please check that the path is correct or that you have run BrainMapAnnotAndManhattanPlot() first"))
}
qqplot(bwas$p, col=colourList[1])
 print(round(median(bwas$CHI2,na.rm=T)/qchisq(0.5,df=1),3))
legend(x = 3.5, y = 0.75*max(bwas$log10p), legend = legendList ,pch=20, pt.cex=1.5,  col = colourList)

jjj=2
 for (scenario in legendList[-1] ){
     if (file.exists(paste0(inputPaths[jjj], "/BWAS_fullSummary_", bwasFiles[jjj], ".csv"))) {
bwas=vroom(paste0(inputPaths[jjj], "/BWAS_fullSummary_", bwasFiles[jjj], ".csv"), show_col_types = F) } else {
 print(paste0(inputPaths[jjj], "/BWAS_fullSummary_", bwasFiles[jjj], ".csv", "  not found, please check that the path is correct or that you have run BrainMapAnnotAndManhattanPlot() first"))
}
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
#' @param inputPath path (folder) to the raw brain association maps (outputs of identifyClustersBWAS())
#' @param bwasFile name of the brain association map
#' @param variancePheno Variance of the phenotype (used to standardise the effect sizes into correlations)
#' @param signifThreshold pvalue significance threshold used to account for multiple testing
#' @return Outside and Inside snapshots of the surfaces.
#' @param outputPath path where the outputs will be written
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
plotSubcortical=function(inputPath, bwasFile, variancePheno, outputPath, signifThreshold){

  for (moda in c("thick", "LogJacs")){
 for (hemi in c("lh", "rh")){
bwasPlot=vroom( paste0(inputPath, bwasFile ), show_col_types = F)
bwasPlot=formatBWASsubcortical(BWASsumstat=bwasPlot, hemi=hemi, mod=moda)
bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Y=bwasPlot$Y*(-1)

# Transform betas in correlations (assumes all vertices have been standardised)
bwasPlot$cor=bwasPlot$b/sqrt(variancePheno)

# Identify significant vertices to plot
bwasPlot$signifVoxel=ifelse(bwasPlot$p < signifThreshold, 1 ,0)

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
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile, "_", hemi, "_", moda , "_inside.png"))
rgl.close()

bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Z=bwasPlot$Z*(-1)
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Z",  "X", "Y")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile,"_", hemi, "_", moda , "_outside.png"))
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
#' @param inputPath path (folder) to the raw brain association maps (outputs of identifyClustersBWAS())
#' @param bwasFile name of the brain association map
#' @param variancePheno Variance of the phenotype (used to standardise the effect sizes into correlations)
#' @param signifThreshold pvalue significance threshold used to account for multiple testing
#' @return Outside and Inside snapshots of the surfaces.
#' @param outputPath path where the outputs will be written
#' @return Outside and Inside snapshots of the surfaces.
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices vroom
#' @export
plotSubcortical_flat=function(inputPath, bwasFile, variancePheno, signifThreshold, outputPath){

for (moda in c("thick", "LogJacs")){
 for (hemi in c("lh", "rh")){
bwasPlot=vroom( paste0(inputPath, bwasFile ), show_col_types = F)
bwasPlot=formatBWASsubcortical(BWASsumstat=bwasPlot, hemi=hemi, mod=moda)
bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Y=bwasPlot$Y*(-1)

# Transform betas in correlations (assumes all vertices have been standardised)
bwasPlot$cor=bwasPlot$b/sqrt(variancePheno)

# Identify significant vertices to plot
bwasPlot$signifVoxel=ifelse(bwasPlot$p < signifThreshold, 1 ,0)

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(-0.1, 0.1, len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6],RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=ifelse(  bwasPlot$signifVoxel==1, 1.5 ,0.8 )

# Change coordinates for flat plotting
if (hemi=="rh"){

bwasPlot$Xf=bwasPlot$X
bwasPlot$Yf=bwasPlot$Y
bwasPlot$Zf=bwasPlot$Z

# Hippocampus
bwasPlot$Yf[bwasPlot$ROINb %in% c(17,53)]=bwasPlot$Yf[bwasPlot$ROINb %in% c(17,53)]-17
# Amygdala
bwasPlot$Yf[bwasPlot$ROINb %in% c(18,54)]=bwasPlot$Yf[bwasPlot$ROINb %in% c(18,54)]-17
bwasPlot$Zf[bwasPlot$ROINb %in% c(18,54)]=bwasPlot$Zf[bwasPlot$ROINb %in% c(18,54)]+7
# Thalamus
bwasPlot$Zf[bwasPlot$ROINb %in% c(10,49)]=bwasPlot$Zf[bwasPlot$ROINb %in% c(10,49)]-20
# Caudate
bwasPlot$Yf[bwasPlot$ROINb %in% c(11,50)]=bwasPlot$Yf[bwasPlot$ROINb %in% c(11,50)]+20
# Accumbens
bwasPlot$Zf[bwasPlot$ROINb %in% c(26,58)]=bwasPlot$Zf[bwasPlot$ROINb %in% c(26,58)]+12
# Pallidum
bwasPlot$Yf[bwasPlot$ROINb %in% c(13,52)]=bwasPlot$Yf[bwasPlot$ROINb %in% c(13,52)]-15

# Second set of coordinates
bwasPlot$Xf2=bwasPlot$Xf*(-1)
bwasPlot$Zf2=bwasPlot$Zf*(-1)
bwasPlot$Yf2=bwasPlot$Yf
# Accumbens
bwasPlot$Zf2[bwasPlot$ROINb %in% c(26,58)]=bwasPlot$Zf2[bwasPlot$ROINb %in% c(26,58)]-7
# Pallidum
bwasPlot$Yf2[bwasPlot$ROINb %in% c(13,52)]=bwasPlot$Yf2[bwasPlot$ROINb %in% c(13,52)]-8
# Amygdala
bwasPlot$Yf2[bwasPlot$ROINb %in% c(18,54)]=bwasPlot$Yf2[bwasPlot$ROINb %in% c(18,54)]-5
# Hippocampus
bwasPlot$Yf2[bwasPlot$ROINb %in% c(17,53)]=bwasPlot$Yf2[bwasPlot$ROINb %in% c(17,53)]-5
# Caudate
bwasPlot$Yf2[bwasPlot$ROINb %in% c(11,50)]=bwasPlot$Yf2[bwasPlot$ROINb %in% c(11,50)]-8
##########################
} else if (hemi == "lh"){
# Second set of coordinates

bwasPlot$Xf=bwasPlot$X
bwasPlot$Yf=bwasPlot$Y
bwasPlot$Zf=bwasPlot$Z
# Hippocampus
bwasPlot$Yf[bwasPlot$ROINb %in% c(17,53)]=bwasPlot$Yf[bwasPlot$ROINb %in% c(17,53)]-21
# Amygdala
bwasPlot$Yf[bwasPlot$ROINb %in% c(18,54)]=bwasPlot$Yf[bwasPlot$ROINb %in% c(18,54)]-21
bwasPlot$Zf[bwasPlot$ROINb %in% c(18,54)]=bwasPlot$Zf[bwasPlot$ROINb %in% c(18,54)]+7
# Thalamus
bwasPlot$Zf[bwasPlot$ROINb %in% c(10,49)]=bwasPlot$Zf[bwasPlot$ROINb %in% c(10,49)]-20
# Caudate
bwasPlot$Yf[bwasPlot$ROINb %in% c(11,50)]=bwasPlot$Yf[bwasPlot$ROINb %in% c(11,50)]+12
# Accumbens
bwasPlot$Zf[bwasPlot$ROINb %in% c(26,58)]=bwasPlot$Zf[bwasPlot$ROINb %in% c(26,58)]+19
# Pallidum
bwasPlot$Yf[bwasPlot$ROINb %in% c(13,52)]=bwasPlot$Yf[bwasPlot$ROINb %in% c(13,52)]-22

# Second set of coordinates
bwasPlot$Xf2=bwasPlot$Xf*(-1)
bwasPlot$Zf2=bwasPlot$Zf*(-1)
bwasPlot$Yf2=bwasPlot$Yf
# Caudate
bwasPlot$Yf2[bwasPlot$ROINb %in% c(11,50)]=bwasPlot$Yf2[bwasPlot$ROINb %in% c(11,50)]+8
# Pallidum
bwasPlot$Yf2[bwasPlot$ROINb %in% c(13,52)]=bwasPlot$Yf2[bwasPlot$ROINb %in% c(13,52)]+6
# Amygdala
bwasPlot$Yf2[bwasPlot$ROINb %in% c(18,54)]=bwasPlot$Yf2[bwasPlot$ROINb %in% c(18,54)]+3
# Hippocampus
bwasPlot$Yf2[bwasPlot$ROINb %in% c(17,53)]=bwasPlot$Yf2[bwasPlot$ROINb %in% c(17,53)]+3
# Accumbens
bwasPlot$Zf2[bwasPlot$ROINb %in% c(26,58)]=bwasPlot$Zf2[bwasPlot$ROINb %in% c(26,58)]+6
# Thalamus
bwasPlot$Zf2[bwasPlot$ROINb %in% c(10,49)]=bwasPlot$Zf2[bwasPlot$ROINb %in% c(10,49)]+2

}

# Draw plots and save screenshots
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Zf",  "Xf", "Yf")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile, "_", hemi, "_", moda , "_clustersAndCoordinates_inside.png"))
rgl.close()

par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Zf2",  "Xf2", "Yf2")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile,"_", hemi, "_", moda , "_clustersAndCoordinates_outside.png"))
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
#' @param inputPath path (folder) to the raw brain association maps (outputs of identifyClustersBWAS())
#' @param bwasFile name of the brain association map
#' @param variancePheno Variance of the phenotype (used to standardise the effect sizes into correlations)
#' @param signifThreshold pvalue significance threshold used to account for multiple testing
#' @return Outside and Inside snapshots of the surfaces.
#' @param outputPath path where the outputs will be written
#' @return Outside and Inside snapshots of the surfaces.
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices vroom
#' @export
plotCortical=function(inputPath, bwasFile, variancePheno, signifThreshold, outputPath){

 for (moda in c("area", "thickness")){

 for (hemi in c("rh", "lh")){

bwasPlot=vroom(paste0(inputPath, bwasFile ) , show_col_types = F)
bwasPlot=formatBWAScortical(BWASsumstat = bwasPlot, hemi = hemi, mod = moda)

# Transform betas in correlations
bwasPlot$cor=bwasPlot$b/sqrt(variancePheno)

# Identify significant vertices to plot
bwasPlot$signifVoxel=ifelse(bwasPlot$p < signifThreshold, 1 ,0)

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(-0.4, 0.4, len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6], RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=ifelse(bwasPlot$signifVoxel==1, 2, 0.8)

# Draw plot
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Y","X", "Z")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile,"_", hemi, "_", moda , "_inside.png"))
rgl.close()

bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Y=bwasPlot$Y*(-1)
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Y","X", "Z")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile,"_", hemi, "_", moda , "_outside.png"))
rgl.close()

}}
}


#' Combine cortical and subcortical figures to create a multi-panel plot
#'
#' This function uses the cortical ans subcortical figures created using the specific functions
#'
#'
#' @param inputPath path (folder) to the cortical and subcortical plots
#' @param bwasFile Variable name (used to name the clustersAndCoordinates file)
#' @param pathToLegendBar Path to the manually created legend bar
#' @param outputPath path where the outputs will be written
#' @return Outside and Inside snapshots of the surfaces.
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
combineCorticalSubcorticalPlots=function(inputPath, bwasFile, pathToLegendBar, outputPath){

# List of files
ll=NULL
for (moda in c("thickness", "area","thick", "LogJacs") ){
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_", "lh", "_", moda, "_clustersAndCoordinates_outside.png"))
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_","lh", "_", moda, "_clustersAndCoordinates_inside.png"))
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_", "rh", "_", moda, "_clustersAndCoordinates_inside.png"))
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_", "rh", "_", moda, "_clustersAndCoordinates_outside.png"))
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
ggsave(paste0(outputPath, "/Plots_Combined", bwasFile, ".png"),width=18, height=4, gs)

}


