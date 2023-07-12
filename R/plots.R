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
BWASsumstat$CHI2=(BWASsumstatfile$b/BWASsumstatfile$se)**2

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
    par(mar=c(5,5,2,0.5))
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
png(paste0(outputPath, "QQplotCombined_assocVertices", phenotypeLabel, ".png"), width = 20, height = 15, units = "cm", res=400)
par(mar=c(4,4,2,10))
if (file.exists(paste0(inputPaths[1], "/BWAS_fullSummary_", bwasFiles[1], ".csv"))){
bwas=vroom(paste0(inputPaths[1], "/BWAS_fullSummary_", bwasFiles[1], ".csv"), show_col_types = F)} else {
    print(paste0(inputPaths[1], "/BWAS_fullSummary_", bwasFiles[1], ".csv", "  not found, please check that the path is correct or that you have run BrainMapAnnotAndManhattanPlot() first"))
}
qqplot(bwas$p, col=colourList[1])
coord <- par("usr") # options for legend position
legend(x = coord[2] * 1.05, y = coord[4], legend = legendList ,pch=20, pt.cex=1.5,  col = colourList, xpd = TRUE)
print(paste0("lambda=",round(median(bwas$CHI2,na.rm=T)/qchisq(0.5,df=1),3)))

jjj=2
 for (scenario in legendList[-1] ){
     if (file.exists(paste0(inputPaths[jjj], "/BWAS_fullSummary_", bwasFiles[jjj], ".csv"))) {
bwas=vroom(paste0(inputPaths[jjj], "/BWAS_fullSummary_", bwasFiles[jjj], ".csv"), show_col_types = F) } else {
 print(paste0(inputPaths[jjj], "/BWAS_fullSummary_", bwasFiles[jjj], ".csv", "  not found, please check that the path is correct or that you have run BrainMapAnnotAndManhattanPlot() first"))
}
qqplot(bwas$p, col=colourList[jjj], add=T)
print(paste0("lambda=",round(median(bwas$CHI2,na.rm=T)/qchisq(0.5,df=1),3)))
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
#' @param correlationRange range of the correlation coefficients (for improved colors) - default is (-1; 1)
#' @param outputPath path where the outputs will be written
#' @return Outside and Inside snapshots of the surfaces.
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
plotSubcortical=function(inputPath, bwasFile, variancePheno, outputPath, signifThreshold, correlationRange=c(-1,1)){

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
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(correlationRange[1], correlationRange[2], len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6],RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=ifelse(  bwasPlot$signifVoxel==1, 1.5 ,0.8 )

# Draw plots and save screenshots
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Z",  "X", "Y")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile, "_", hemi, "_", moda , "_inside_bundled.png"))
rgl.close()

bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Z=bwasPlot$Z*(-1)
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Z",  "X", "Y")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile,"_", hemi, "_", moda , "_outside_bundled.png"))
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
#' @param correlationRange range of the correlation coefficients (for improved colors) - default is (-1; 1)
#' @param outputPath path where the outputs will be written
#' @return Outside and Inside snapshots of the surfaces.
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices vroom
#' @export
plotSubcortical_flat=function(inputPath, bwasFile, variancePheno, signifThreshold, outputPath, correlationRange=c(-1,1)){

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
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(correlationRange[1], correlationRange[2], len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6],RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=ifelse(  bwasPlot$signifVoxel==1, 1.5 ,0.8 )

# Add flat coordinates
bwasPlot=addFlatCoordinatesSubcortical(annotBwas = bwasPlot, hemi = hemi)

# Draw plots and save screenshots
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Zf",  "Xf", "Yf")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile, "_", hemi, "_", moda , "_inside.png"))
rgl.close()

par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Zf2",  "Xf2", "Yf2")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile,"_", hemi, "_", moda , "_outside.png"))
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
#' @param correlationRange range of the correlation coefficients (for improved colors) - default is (-1; 1)
#' @param outputPath path where the outputs will be written
#' @param style style for plotting (based on the different cortical fsaverage surface). Possible options: "orig", "pial", "sphere", "inflated",  "inflated_pre",  "pial_semi_inflated",  "smoothwm" (default)
#' @return Outside and Inside snapshots of the surfaces.
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices vroom
#' @export
plotCortical=function(inputPath, bwasFile, variancePheno, signifThreshold, outputPath, correlationRange=c(-1,1), style="smoothwm"
){

 for (moda in c("area", "thickness")){

 for (hemi in c("rh", "lh")){

bwasPlot=vroom(paste0(inputPath, bwasFile ) , show_col_types = F)
bwasPlot=formatBWAScortical(BWASsumstat = bwasPlot, hemi = hemi, mod = moda)

# Transform betas in correlations
bwasPlot$cor=bwasPlot$b/sqrt(variancePheno)

# Identify significant vertices to plot
bwasPlot$signifVoxel=ifelse(bwasPlot$p < signifThreshold, 1 ,0)

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(correlationRange[1], correlationRange[2], len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6], RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=ifelse(bwasPlot$signifVoxel==1, 2, 0.8)

# Draw plot
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[, paste(c("Y","X", "Z"), style, sep = "_")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile,"_", hemi, "_", moda ,"_",style, "_inside.png"))
rgl.close()

bwasPlot[, paste(c("X"), style, sep = "_")]=bwasPlot[, paste(c("X"), style, sep = "_")]*(-1)
bwasPlot[, paste(c("Y"), style, sep = "_")]=bwasPlot[, paste(c("Y"), style, sep = "_")]*(-1)
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[, paste(c("Y","X", "Z"), style, sep = "_")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile,"_", hemi, "_", moda ,"_",style,  "_outside.png"))
rgl.close()

}}
}


#' Scan association file to return the range of the association statistics (for the significant vertices)
#'
#' This function takes in raw brain association maps, and returns the range of association, which can be used to inform plotting
#'
#'
#' @param inputPath path (folder) to the raw brain association maps (outputs of OSCA)
#' @param bwasFile name of the brain association map
#' @param variancePheno variance of the phenotypes (to estimate correlation)
#' @param signifThreshold pvalue significance threshold used to account for multiple testing
#' @return Range of correlation coefficients, useful to set appropriate legend bar and cor range in plotting function
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices vroom
#' @export
getCorrelationRange<-function(inputPath , bwasFile , variancePheno, signifThreshold){

# Open summary statistics
BWASsumstatfile<-  vroom( paste0(inputPath ,  bwasFile), show_col_types = F)
BWASsumstatfile=BWASsumstatfile[which(BWASsumstatfile$p < signifThreshold),]
return(range(BWASsumstatfile$b/sqrt(variancePheno), na.rm = T))

}

#' Combine cortical and subcortical figures to create a multi-panel plot
#'
#' This function uses the cortical and subcortical figures created using the specific functions
#'
#'
#' @param outputPath path where the outputs will be written
#' @param correlationRange range of the correlation coefficients (for improved colors) - default is (-1; 1)
#' @param statisticsName name of the statistics, which is to appear on top of the legend (default is "Correlation")
#' @return The legend bar with manually selected range of effect sizes.
#' @import RColorBrewer
#' @export
createLegendBar=function(outputPath, correlationRange=c(-1,1), statisticsName="Correlation"){
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6],RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours

png(paste0(outputPath, "/legendbar", correlationRange[1],"_",correlationRange[2], ".png"), width = 5, height = 15, units = "cm", res = 400)
par(mar=c(2,4,2,2))
my.colors = colorRampPalette(cols)
z=matrix(1:100,nrow=1)
x=1
y=seq(correlationRange[1],correlationRange[2],len=100)
image(x,y,z,col=my.colors(20),axes=FALSE,xlab="",ylab="", main=statisticsName)
axis(2)
dev.off()
}



#' Combine cortical and subcortical figures to create a multi-panel plot
#'
#' This function uses the cortical and subcortical figures created using the specific functions
#'
#'
#' @param inputPath path (folder) to the cortical and subcortical plots
#' @param bwasFile name of bwas file (used in naming of cortical and subcortical plots)
#' @param pathToLegendBar Path to the legend bar, created with createLegendBar()
#' @param outputPath folder where the outputs will be written
#' @param style style for plotting (based on the different cortical fsaverage surface). Possible options: "orig", "pial", "sphere", "inflated",  "inflated_pre",  "pial_semi_inflated",  "smoothwm" (default)
#' @return A combined plot with cortical and subcortical surface plots and legendbar
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
combineCorticalSubcorticalPlots=function(inputPath, bwasFile, pathToLegendBar, outputPath, style="smoothwm"){

# List of files
ll=NULL
for (moda in c("thickness", "area") ){
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_", "lh", "_", moda,"_",style,  "_outside.png"))
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_","lh", "_", moda,"_",style,  "_inside.png"))
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_", "rh", "_", moda,"_",style,  "_inside.png"))
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_", "rh", "_", moda, "_",style, "_outside.png"))
}
for (moda in c("thick", "LogJacs") ){
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_", "lh", "_", moda,  "_outside.png"))
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_","lh", "_", moda,  "_inside.png"))
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_", "rh", "_", moda,  "_inside.png"))
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_", "rh", "_", moda,  "_outside.png"))
}
ll=c(ll, paste0(pathToLegendBar))

# Crop and convert format
plots2 <- lapply(ll<-ll,function(x){
  if(x!=paste0(pathToLegendBar)){
   img <- as.raster(readPNG(x)[,50:1150,])} else {
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

#' Subcortical GIF
#'
#' This function reads in a brain association map.
#' It produces snapshots of the brain surfaces, with a slight rotation angle, as well as GIF
#' The function needs the variance of the phenotype, in order to convert association betas into correlation coefficients.
#'
#'
#' @param inputPath path (folder) to the raw brain association maps
#' @param bwasFile name of the brain association map
#' @param variancePheno Variance of the phenotype (used to standardise the effect sizes into correlations)
#' @param signifThreshold pvalue significance threshold
#' @param moda modality to plot ("thick" or "LogJacs")
#' @param hemi hemisphere to plot ("lh" or "rh")
#' @param nbImagesForGif Number of png images to generate for the gif. The larger the smoother the gif, but the longer it takes to generate. Please use multiples of 6 if you want the GIF to be created automatically.
#' @param outputPath path where the outputs will be written. Absolute path may be required to produce GIF.
#' @param correlationRange range of the correlation coefficients (for improved colors) - default is (-1; 1)
#' @return Several snapshots of the surfaces.
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices mvMonitoring
#' @export
plotSubcorticalGIF=function(inputPath, bwasFile, variancePheno, outputPath, hemi, moda, signifThreshold, nbImagesForGif, correlationRange=c(-1,1)){

# Open and format
bwasPlot=vroom( paste0(inputPath, bwasFile ), show_col_types = F)
bwasPlot=formatBWASsubcortical(BWASsumstat=bwasPlot, hemi=hemi, mod=moda)
bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Y=bwasPlot$Y*(-1)

# Transform betas in correlations (assumes all vertices have been standardised)
bwasPlot$cor=bwasPlot$b/sqrt(variancePheno)

# Identify significant vertices to plot
bwasPlot$signifVoxel=ifelse(bwasPlot$p < signifThreshold, 1 ,0)

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(correlationRange[1], correlationRange[2], len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6],RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=ifelse( bwasPlot$signifVoxel==1, 1.5 ,0.8 )

# Draw plots and save screenshots
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Z",  "X", "Y")]), col=bwasPlot$color, radius = bwasPlot$radius)
movie3d( spin3d(rpm=10), duration = 6, fps = nbImagesForGif/6 , frames = paste0(outputPath, "/BWAS_", bwasFile,"_", hemi, "_", moda , "_GIF_"), dir = outputPath,  convert=NULL, clean=F, movie = paste0(outputPath, "/BWAS_", bwasFile,"_", hemi, "_", moda , "_GIF") )
rgl.close()

}

#' Add new coordinates to subcortical annotated maps - which allow for flat plotting
#'
#' This function reads in an annotated subcortical brain association map
#' It add additional sets of coordinates for flat plotting
#'
#'
#' @param annotBwas Annotated subcortical brain map
#' @param hemi Hemisphere ("lh" or "rh")
#' @return Annotated subcortical brain map with additional coordiantes
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices mvMonitoring
#' @export
addFlatCoordinatesSubcortical=function(annotBwas, hemi){

# Change coordinates for flat plotting
if (hemi=="rh"){

annotBwas$Xf=annotBwas$X
annotBwas$Yf=annotBwas$Y
annotBwas$Zf=annotBwas$Z

# Hippocampus
annotBwas$Yf[annotBwas$ROINb %in% c(17,53)]=annotBwas$Yf[annotBwas$ROINb %in% c(17,53)]-22
# Amygdala
annotBwas$Yf[annotBwas$ROINb %in% c(18,54)]=annotBwas$Yf[annotBwas$ROINb %in% c(18,54)]-22
annotBwas$Zf[annotBwas$ROINb %in% c(18,54)]=annotBwas$Zf[annotBwas$ROINb %in% c(18,54)]+9
# Thalamus
annotBwas$Zf[annotBwas$ROINb %in% c(10,49)]=annotBwas$Zf[annotBwas$ROINb %in% c(10,49)]-22
# Caudate
annotBwas$Yf[annotBwas$ROINb %in% c(11,50)]=annotBwas$Yf[annotBwas$ROINb %in% c(11,50)]+24
# Accumbens
annotBwas$Zf[annotBwas$ROINb %in% c(26,58)]=annotBwas$Zf[annotBwas$ROINb %in% c(26,58)]+14
# Pallidum
annotBwas$Yf[annotBwas$ROINb %in% c(13,52)]=annotBwas$Yf[annotBwas$ROINb %in% c(13,52)]-17

# Second set of coordinates (outside view)
annotBwas$Xf2=annotBwas$Xf*(-1)
annotBwas$Zf2=annotBwas$Zf*(-1)
annotBwas$Yf2=annotBwas$Yf
# Accumbens
annotBwas$Zf2[annotBwas$ROINb %in% c(26,58)]=annotBwas$Zf2[annotBwas$ROINb %in% c(26,58)]-7
# Pallidum
annotBwas$Yf2[annotBwas$ROINb %in% c(13,52)]=annotBwas$Yf2[annotBwas$ROINb %in% c(13,52)]-8
# Amygdala
annotBwas$Yf2[annotBwas$ROINb %in% c(18,54)]=annotBwas$Yf2[annotBwas$ROINb %in% c(18,54)]-5
# Hippocampus
annotBwas$Yf2[annotBwas$ROINb %in% c(17,53)]=annotBwas$Yf2[annotBwas$ROINb %in% c(17,53)]-5
# Caudate
annotBwas$Yf2[annotBwas$ROINb %in% c(11,50)]=annotBwas$Yf2[annotBwas$ROINb %in% c(11,50)]-8
##########################
} else if (hemi == "lh"){
# Second set of coordinates

annotBwas$Xf=annotBwas$X
annotBwas$Yf=annotBwas$Y
annotBwas$Zf=annotBwas$Z
# Hippocampus
annotBwas$Yf[annotBwas$ROINb %in% c(17,53)]=annotBwas$Yf[annotBwas$ROINb %in% c(17,53)]-25
# Amygdala
annotBwas$Yf[annotBwas$ROINb %in% c(18,54)]=annotBwas$Yf[annotBwas$ROINb %in% c(18,54)]-25
annotBwas$Zf[annotBwas$ROINb %in% c(18,54)]=annotBwas$Zf[annotBwas$ROINb %in% c(18,54)]+9
# Thalamus
annotBwas$Zf[annotBwas$ROINb %in% c(10,49)]=annotBwas$Zf[annotBwas$ROINb %in% c(10,49)]-26
# Caudate
annotBwas$Yf[annotBwas$ROINb %in% c(11,50)]=annotBwas$Yf[annotBwas$ROINb %in% c(11,50)]+14
# Accumbens
annotBwas$Zf[annotBwas$ROINb %in% c(26,58)]=annotBwas$Zf[annotBwas$ROINb %in% c(26,58)]+21
# Pallidum
annotBwas$Yf[annotBwas$ROINb %in% c(13,52)]=annotBwas$Yf[annotBwas$ROINb %in% c(13,52)]-24

# Second set of coordinates (outside view)
annotBwas$Xf2=annotBwas$Xf*(-1)
annotBwas$Zf2=annotBwas$Zf*(-1)
annotBwas$Yf2=annotBwas$Yf
# Caudate
annotBwas$Yf2[annotBwas$ROINb %in% c(11,50)]=annotBwas$Yf2[annotBwas$ROINb %in% c(11,50)]+8
# Pallidum
annotBwas$Yf2[annotBwas$ROINb %in% c(13,52)]=annotBwas$Yf2[annotBwas$ROINb %in% c(13,52)]+6
# Amygdala
annotBwas$Yf2[annotBwas$ROINb %in% c(18,54)]=annotBwas$Yf2[annotBwas$ROINb %in% c(18,54)]+3
# Hippocampus
annotBwas$Yf2[annotBwas$ROINb %in% c(17,53)]=annotBwas$Yf2[annotBwas$ROINb %in% c(17,53)]+3
# Accumbens
annotBwas$Zf2[annotBwas$ROINb %in% c(26,58)]=annotBwas$Zf2[annotBwas$ROINb %in% c(26,58)]+6
# Thalamus
annotBwas$Zf2[annotBwas$ROINb %in% c(10,49)]=annotBwas$Zf2[annotBwas$ROINb %in% c(10,49)]+2

}

return(annotBwas)

}




#' Subcortical GIF - translation from grouped subcortical volumes to flat plotting
#'
#' This function reads in a brain association map
#' It produces snapshots of the brain surfaces, with a slight rotation angle, which can be used to make a GIF
#' The function needs the variance of the phenotype, in order to convert association betas into correlation coefficients.
#'
#'
#' @param inputPath path (folder) to the raw brain association maps
#' @param bwasFile name of the brain association map
#' @param variancePheno Variance of the phenotype (used to standardise the effect sizes into correlations)
#' @param signifThreshold pvalue significance threshold
#' @param moda modality to plot ("thick" or "LogJacs")
#' @param hemi hemisphere to plot ("lh" or "rh")
#' @param nbImagesForGif Number of png images to generate for the gif. The larger the smoother the gif, but the longer it takes to generate.
#' @param leftOrRightView The side from which the snapshots will be taken
#' @param outputPath path where the outputs will be written
#' @param correlationRange range of the correlation coefficients (for improved colors) - default is (-1; 1)
#' @return Several snapshots of the surfaces.
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices mvMonitoring magick
#' @export
plotSubcorticalToFlatGIF=function(inputPath, bwasFile, variancePheno, outputPath, hemi, moda, signifThreshold, nbImagesForGif, leftOrRightView, correlationRange=c(-1,1)){

# Open and format
bwasPlot=vroom( paste0(inputPath, bwasFile ), show_col_types = F)
bwasPlot=formatBWASsubcortical(BWASsumstat=bwasPlot, hemi=hemi, mod=moda)
bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Y=bwasPlot$Y*(-1)

# Transform betas in correlations (assumes all vertices have been standardised)
bwasPlot$cor=bwasPlot$b/sqrt(variancePheno)

# Identify significant vertices to plot
bwasPlot$signifVoxel=ifelse(bwasPlot$p < signifThreshold, 1 ,0)

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(correlationRange[1], correlationRange[2], len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6],RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=ifelse( bwasPlot$signifVoxel==1, 1.5 ,0.8 )

# Add flat coordinates
bwasPlot=addFlatCoordinatesSubcortical(annotBwas = bwasPlot, hemi = hemi)

# Produce plots
if (leftOrRightView=="left"){

for (iii in 0:nbImagesForGif){

plotMax=bwasPlot[,c( "Z",  "X", "Y")]+(bwasPlot[,c( "Zf",  "Xf", "Yf")]-bwasPlot[,c( "Z",  "X", "Y")]) / nbImagesForGif*(iii)
plotMax=as.data.frame(plotMax)
bwasPlot$Zplot=plotMax[,1]
bwasPlot$Xplot=plotMax[,2]
bwasPlot$Yplot=plotMax[,3]

# Draw plots and save screenshots
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Zplot",  "Xplot", "Yplot")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile, "_", hemi, "_", moda , "_ToFlatGIF_", leftOrRightView, "_view_",  sprintf(fmt = "%03d", iii), ".png"))
rgl.close()
}

} else if (leftOrRightView=="right"){

bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Z=bwasPlot$Z*(-1)
bwasPlot$Y=bwasPlot$Y

for (iii in 0:nbImagesForGif){

plotMax=bwasPlot[,c( "Z",  "X", "Y")]+(bwasPlot[,c( "Zf2",  "Xf2", "Yf2")]-bwasPlot[,c( "Z",  "X", "Y")]) / nbImagesForGif*(iii)
plotMax=as.data.frame(plotMax)
bwasPlot$Zplot=plotMax[,1]
bwasPlot$Xplot=plotMax[,2]
bwasPlot$Yplot=plotMax[,3]

# Draw plots and save screenshots
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Zplot",  "Xplot", "Yplot")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile, "_", hemi, "_", moda , "_ToFlatGIF_", leftOrRightView, "_view_",  sprintf(fmt = "%03d", iii), ".png"))
rgl.close()

}
}

# Create GIF using magick functions
## Open files
print("Making GIF")
imgs <- list.files(path = outputPath , pattern = paste0("BWAS_", bwasFile, "_", hemi, "_", moda , "_ToFlatGIF_", leftOrRightView, "_view_"),  full.names = T )
imgs<-c(imgs[1:length(imgs)], imgs[length(imgs):1])

img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 20 frames per second
img_animated <- image_animate(img_joined, fps = 20)
## save to disk
image_write(image = img_animated, path = paste0(outputPath, "/BWAS_", bwasFile, "_", hemi, "_", moda , "_ToFlatGIF_", leftOrRightView, "_view", ".gif"))

}



#' Subcortical GIF flat - rotation of each subcortical structure from the flat format
#'
#' This function reads in a brain association map
#' It produces snapshots of the subcortical brain surfaces, with a slight rotation angle, which can be used to make a GIF
#' The function needs the variance of the phenotype, in order to convert association betas into correlation coefficients.
#'
#'
#' @param inputPath path (folder) to the raw brain association maps
#' @param bwasFile name of the brain association map
#' @param variancePheno Variance of the phenotype (used to standardise the effect sizes into correlations)
#' @param signifThreshold pvalue significance threshold
#' @param moda modality to plot ("thick" or "LogJacs")
#' @param hemi hemisphere to plot ("lh" or "rh")
#' @param nbImagesForGif Number of png images to generate for the gif. The larger the smoother the gif, but the longer it takes to generate.
#' @param leftOrRightView The side from which the snapshots will be taken
#' @param outputPath path where the outputs will be written
#' @param correlationRange range of the correlation coefficients (for improved colors) - default is (-1; 1)
#' @return Several snapshots of the surfaces.
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices mvMonitoring magick
#' @export
plotSubcortical_FlatGIF=function(inputPath, bwasFile, variancePheno, outputPath, hemi, moda, signifThreshold, nbImagesForGif, leftOrRightView, correlationRange=c(-1,1)){

# Open and format
bwasPlot=vroom( paste0(inputPath, bwasFile ), show_col_types = F)
bwasPlot=formatBWASsubcortical(BWASsumstat=bwasPlot, hemi=hemi, mod=moda)
bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Y=bwasPlot$Y*(-1)

# Transform betas in correlations (assumes all vertices have been standardised)
bwasPlot$cor=bwasPlot$b/sqrt(variancePheno)

# Identify significant vertices to plot
bwasPlot$signifVoxel=ifelse(bwasPlot$p < signifThreshold, 1 ,0)

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(correlationRange[1], correlationRange[2], len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6],RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=ifelse( bwasPlot$signifVoxel==1, 1.5 ,0.8 )

# Add flat coordinates
bwasPlot=addFlatCoordinatesSubcortical(annotBwas = bwasPlot, hemi = hemi)

# Produce plots
if (leftOrRightView=="left"){

for (iii in 0:nbImagesForGif){

# Draw plots and save screenshots
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Zf",  "Xf", "Yf")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile, "_", hemi, "_", moda , "_FlatGIF_", leftOrRightView, "_view_", sprintf(fmt = "%03d", iii), ".png"))
rgl.close()

# Each structure is rotated

# First translation to align the rotation axis with origin
for (ROI in unique(bwasPlot$ROINb)){

    # Get coordinates of barycentre
coordBar=c(min(bwasPlot$Zf[which(bwasPlot$ROINb==ROI)]) + (max(bwasPlot$Zf[which(bwasPlot$ROINb==ROI)])-min(bwasPlot$Zf[which(bwasPlot$ROINb==ROI)]))/2, min(bwasPlot$Xf[which(bwasPlot$ROINb==ROI)]) + (max(bwasPlot$Xf[which(bwasPlot$ROINb==ROI)])-min(bwasPlot$Xf[which(bwasPlot$ROINb==ROI)]))/2 , min(bwasPlot$Yf[which(bwasPlot$ROINb==ROI)])+ (max(bwasPlot$Yf[which(bwasPlot$ROINb==ROI)])-min(bwasPlot$Yf[which(bwasPlot$ROINb==ROI)]))/2 )

# Translate shape
bwasPlot$Zf[which(bwasPlot$ROINb==ROI)]=bwasPlot$Zf[which(bwasPlot$ROINb==ROI)]-coordBar[1]
bwasPlot$Xf[which(bwasPlot$ROINb==ROI)]=bwasPlot$Xf[which(bwasPlot$ROINb==ROI)]-coordBar[2]
bwasPlot$Yf[which(bwasPlot$ROINb==ROI)]=bwasPlot$Yf[which(bwasPlot$ROINb==ROI)]-coordBar[3]

# Rotation
plotMax=as.matrix(bwasPlot[which(bwasPlot$ROINb==ROI),c( "Zf",  "Xf", "Yf")]) %*% rotateScale3D(rot_angles = c(360/nbImagesForGif,0,0))
plotMax=as.data.frame(plotMax)
bwasPlot$Zf[which(bwasPlot$ROINb==ROI)]=plotMax[,1]
bwasPlot$Xf[which(bwasPlot$ROINb==ROI)]=plotMax[,2]
bwasPlot$Yf[which(bwasPlot$ROINb==ROI)]=plotMax[,3]

# Back translation
bwasPlot$Zf[which(bwasPlot$ROINb==ROI)]=bwasPlot$Zf[which(bwasPlot$ROINb==ROI)]+coordBar[1]
bwasPlot$Xf[which(bwasPlot$ROINb==ROI)]=bwasPlot$Xf[which(bwasPlot$ROINb==ROI)]+coordBar[2]
bwasPlot$Yf[which(bwasPlot$ROINb==ROI)]=bwasPlot$Yf[which(bwasPlot$ROINb==ROI)]+coordBar[3]

} } } else if (leftOrRightView=="right"){

bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Z=bwasPlot$Z*(-1)
bwasPlot$Y=bwasPlot$Y

for (iii in 0:nbImagesForGif){

# Draw plots and save screenshots
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Zf2",  "Xf2", "Yf2")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(outputPath, "/BWAS_", bwasFile, "_", hemi, "_", moda , "_FlatGIF_", leftOrRightView, "_view_", sprintf(fmt = "%03d", iii), ".png"))
rgl.close()

# Each structure is rotated

# First translation to align the rotation axis with origin
for (ROI in unique(bwasPlot$ROINb)){

    # Get coordinates of barycentre
coordBar=c(min(bwasPlot$Zf2[which(bwasPlot$ROINb==ROI)]) + (max(bwasPlot$Zf2[which(bwasPlot$ROINb==ROI)])-min(bwasPlot$Zf2[which(bwasPlot$ROINb==ROI)]))/2, min(bwasPlot$Xf2[which(bwasPlot$ROINb==ROI)]) + (max(bwasPlot$Xf2[which(bwasPlot$ROINb==ROI)])-min(bwasPlot$Xf2[which(bwasPlot$ROINb==ROI)]))/2 , min(bwasPlot$Yf2[which(bwasPlot$ROINb==ROI)])+ (max(bwasPlot$Yf2[which(bwasPlot$ROINb==ROI)])-min(bwasPlot$Yf2[which(bwasPlot$ROINb==ROI)]))/2 )

# Translate shape
bwasPlot$Zf2[which(bwasPlot$ROINb==ROI)]=bwasPlot$Zf2[which(bwasPlot$ROINb==ROI)]-coordBar[1]
bwasPlot$Xf2[which(bwasPlot$ROINb==ROI)]=bwasPlot$Xf2[which(bwasPlot$ROINb==ROI)]-coordBar[2]
bwasPlot$Yf2[which(bwasPlot$ROINb==ROI)]=bwasPlot$Yf2[which(bwasPlot$ROINb==ROI)]-coordBar[3]

# Rotation
plotMax=as.matrix(bwasPlot[which(bwasPlot$ROINb==ROI),c( "Zf2",  "Xf2", "Yf2")]) %*% rotateScale3D(rot_angles = c(360/nbImagesForGif,0,0))
plotMax=as.data.frame(plotMax)

bwasPlot$Zf2[which(bwasPlot$ROINb==ROI)]=plotMax[,1]
bwasPlot$Xf2[which(bwasPlot$ROINb==ROI)]=plotMax[,2]
bwasPlot$Yf2[which(bwasPlot$ROINb==ROI)]=plotMax[,3]

# Back translation
bwasPlot$Zf2[which(bwasPlot$ROINb==ROI)]=bwasPlot$Zf2[which(bwasPlot$ROINb==ROI)]+coordBar[1]
bwasPlot$Xf2[which(bwasPlot$ROINb==ROI)]=bwasPlot$Xf2[which(bwasPlot$ROINb==ROI)]+coordBar[2]
bwasPlot$Yf2[which(bwasPlot$ROINb==ROI)]=bwasPlot$Yf2[which(bwasPlot$ROINb==ROI)]+coordBar[3]
} }

}

# Create GIF using magick functions
## Open files
print("Making GIF")
imgs <- list.files(path = outputPath , pattern = paste0("BWAS_", bwasFile, "_", hemi, "_", moda , "_FlatGIF_", leftOrRightView, "_view_"),  full.names = TRUE )
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

# Reduce size for gif
img_joined<-image_scale(img_joined, "x600")

## animate at 20 frames per second
img_animated <- image_animate(img_joined, fps = 20, optimize = T)
## save to disk
image_write(image = img_animated, path = paste0(outputPath, "/BWAS_", bwasFile, "_", hemi, "_", moda , "_FlatGIF_", leftOrRightView, "_view", ".gif") )

}


#' Cortical brain plots for GIF
#'
#' This function reads in a brain association map.
#' It produces snapshots of the brain surfaces.
#' The function needs the variance of the phenotype, in order to convert association betas into correlation coefficients.
#'
#'
#' @param inputPath path (folder) to the raw brain association maps
#' @param bwasFile name of the brain association map
#' @param variancePheno Variance of the phenotype (used to standardise the effect sizes into correlations)
#' @param signifThreshold pvalue significance threshold
#' @param moda modality to plot ("thickness" or "area")
#' @param hemi hemisphere to plot ("lh" or "rh")
#' @param nbImagesForGif Number of png images to generate for the gif. The larger the smoother the gif, but the longer it takes to generate. Please use a multiple of 6 if you want the GIF to be created automatically.
#' @param outputPath path where the outputs will be written. Absolute path may be required to produce GIF.
#' @param correlationRange range of the correlation coefficients (for improved colors) - default is (-1; 1)
#' @param faster downsample mesh for faster plotting (e.g. testing)
#' @param style style for plotting (based on the different cortical fsaverage surface). Possible options: "orig", "pial", "sphere", "inflated",  "inflated_pre",  "pial_semi_inflated",  "smoothwm" (default)
#' @return Several snapshots of the rotating cortex and a GIF
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices mvMonitoring
#' @export
plotCorticalGIF=function(inputPath, bwasFile, variancePheno, signifThreshold, moda, hemi,  outputPath, nbImagesForGif, correlationRange=c(-1,1), faster, style){

# Open files and annotate
bwasPlot=vroom(paste0(inputPath, bwasFile ) , show_col_types = F)
bwasPlot=formatBWAScortical(BWASsumstat = bwasPlot, hemi = hemi, mod = moda)

# Transform betas in correlations
bwasPlot$cor=bwasPlot$b/sqrt(variancePheno)

# Identify significant vertices to plot
bwasPlot$signifVoxel=ifelse(bwasPlot$p < signifThreshold, 1 ,0)

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(correlationRange[1], correlationRange[2], len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6], RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=ifelse(bwasPlot$signifVoxel==1, 2, 0.8)

if (faster==T){
set.seed(8)
exctr=sample(x = 1:dim(bwasPlot)[1], size = dim(bwasPlot)[1]/4, replace = F)
bwasPlot=bwasPlot[exctr,]
}

# Add max - min points to set camera
# Add max points to fix camera
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[, paste(c("Y","X", "Z"), style, sep = "_")]), col=bwasPlot$color, radius = bwasPlot$radius )
movie3d( spin3d(rpm=10), duration = 6, fps = nbImagesForGif/6 , frames = paste0(outputPath, "/BWAS_", bwasFile,"_", hemi, "_", moda , "_", style,"_GIF_"), dir = outputPath,  convert=NULL, clean=F, movie = paste0(outputPath, "/BWAS_", bwasFile,"_", hemi, "_", moda ,"_", style, "_GIF") )
rgl.close()

}




#' Combine subcortical figures to create a multi-panel plot
#'
#' This function uses the subcortical figures created using the specific functions
#'
#'
#' @param inputPath path (folder) to the cortical and subcortical plots
#' @param bwasFile name of bwas file (used in naming of cortical and subcortical plots)
#' @param pathToLegendBar Path to the legend bar, created with createLegendBar()
#' @param outputPath folder where the outputs will be written
#' @return A combined plot with cortical and subcortical surface plots and legendbar
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
combineSubcorticalPlots=function(inputPath, bwasFile, pathToLegendBar, outputPath){

# List of files
ll=NULL
for (moda in c("thick", "LogJacs") ){
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_", "lh", "_", moda, "_outside.png"))
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_","lh", "_", moda, "_inside.png"))
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_", "rh", "_", moda, "_inside.png"))
ll=c(ll, paste0(inputPath, "/BWAS_", bwasFile, "_", "rh", "_", moda, "_outside.png"))
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
lay <- rbind(c(1,1,3,3,5,5,7,7,9),
            c(2,2,4,4,6,6,8,8,9))

gs=grid.arrange(grobs = plots2, layout_matrix = lay)
ggsave(paste0(outputPath, "/Plots_SubcorticalCombined", bwasFile, ".png"),width=10, height=4, gs)

}

