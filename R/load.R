
#' Open and format association brain map (Cortical)
#'
#' This function opens and format brain association maps to include coordinates and ROI (Desikan atlas).
#' The brain association map requires a ProbeID column, which contains the name of vertices/voxels.
#' It opens the relevant atlas file to merge with the association map.
#'
#' @param BWASsumstat Brain association map (OSCA format)
#' @param hemi Hemisphere ("lh" or "rh")
#' @param mod modality ("area" or "thickness")
#' @return An enriched brain association maps which contains more vertex information
#' @import plyr png qqman readr Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
formatBWAScortical<-function(BWASsumstat, hemi, mod){
BWASsumstat$ProbeID<-BWASsumstat$Probe

atlas<-read.table(system.file("extdata/atlas", "fsaverage_vertex_labels_names_order2.txt", package = "brainMapR",mustWork = TRUE), header=T, stringsAsFactors = F)
# adapt probe prefix to modality and hemisphere
if (hemi=="lh" & mod=="thickness"){ prefix="lht"}
if (hemi=="lh" & mod=="area"){ prefix="lha" }
if (hemi=="rh" & mod=="thickness"){ prefix="rht" }
if (hemi=="rh" & mod=="area"){ prefix="rha"}
atlas$ProbeID<-paste0(prefix, "_",substr(atlas$VertexNum, 11,nchar(atlas$VertexNum)))
BWASsumstat<-merge(BWASsumstat, atlas, by="ProbeID" )

if (hemi=="lh"){
  if ( length(which(is.na(BWASsumstat$LHLabel)) )>0){
  BWASsumstat<-BWASsumstat[-which(is.na(BWASsumstat$LHLabel)),] }
BWASsumstat<-BWASsumstat[order(BWASsumstat[,"newOrderLeft"]),]
BWASsumstat$LABEL<-BWASsumstat[,"LeftHem_Label"]
BWASsumstat$LABELNAME<-BWASsumstat[,"LHLabel"]
BWASsumstat$rgbx<-BWASsumstat[,"Lrgbx"]
BWASsumstat$rgby<-BWASsumstat[,"Lrgby"]
BWASsumstat$rgbz<-BWASsumstat[,"Lrgbz"]
}
if (hemi=="rh"){
  if ( length(which(is.na(BWASsumstat$RHLabel)) )>0){
  BWASsumstat<-BWASsumstat[-which(is.na(BWASsumstat$RHLabel)),] }
BWASsumstat<-BWASsumstat[order(BWASsumstat[,"newOrderRight"]),]
BWASsumstat$LABEL<-BWASsumstat[,"RightHem_Label"]
BWASsumstat$LABELNAME<-BWASsumstat[,"RHLabel"]
BWASsumstat$rgbx<-BWASsumstat[,"Rrgbx"]
BWASsumstat$rgby<-BWASsumstat[,"Rrgby"]
BWASsumstat$rgbz<-BWASsumstat[,"Rrgbz"]
}
BWASsumstat$xax<-1:length(BWASsumstat[,1])

return(BWASsumstat)
}


#' Open and format association brain map (Subcortical)
#'
#' This function opens and format brain association maps to include coordinates and ROI (ENIGMA-shape atlas).
#' The brain association map requires a ProbeID column, which contains the name of vertices/voxels.
#' It opens the relevant atlas file to merge with the association map.
#'
#' @param BWASsumstat Brain association map (OSCA format)
#' @return An enriched brain association maps which contains more vertex information
#' @import plyr png qqman readr Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
formatBWASsubcortical<-function(BWASsumstat){

BWASsumstat$ProbeID<-gsub(".*[a-z][_]", "", BWASsumstat$Probe)

atlas<-read.table(system.file("extdata/atlas", "fsaverage_subcortical_vertices_label.txt", package = "brainMapR",mustWork = TRUE), header=T, stringsAsFactors = F)
BWASsumstat<-merge(BWASsumstat, atlas, by="ProbeID" )

BWASsumstat<-BWASsumstat[order(BWASsumstat$subcv),]
BWASsumstat$LABEL<-BWASsumstat[,"subcv"]
BWASsumstat$LABELNAME<-BWASsumstat[,"subcvName"]

BWASsumstat$xax<-1:length(BWASsumstat[,1])

return(BWASsumstat)
}


#' Add coordinates to formatted association brain map (Cortical)
#'
#' This function opens and format brain association maps to include coordinates (useful for plotting)
#' The brain association map requires a ProbeID column, which contains the name of vertices/voxels.
#' It opens the relevant atlas file to merge with the association map.
#'
#' @param bwasFormatted Formatted brain association map (created using formatBWAScortical)
#' @param hemi Hemisphere ("lh" or "rh")
#' @param moda modality ("area" or "thickness")
#' @return An enriched brain association maps which contains spatial coordinates
#' @import plyr png qqman readr Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
bwasAddCoordinatesCortical=function(bwasFormatted, hemi, moda){
# Get vertices coordinates
coor=read.table(paste0(system.file("extdata/atlas", "", package = "brainMapR",mustWork = TRUE), hemi, ".vertexCoordinates.asc"))
coor$V1=paste0(hemi, substr(moda, 1,1), "_", coor$V1)
colnames(coor)=c("ProbeID", "X", "Y", "Z", "V4")
# Add coordinates
bwasFormatted=merge(bwasFormatted, coor[,1:4], by="ProbeID")
bwasFormatted$X=bwasFormatted$X*(-1) # coordinates are left/right swapped in asc format
bwasFormatted$hemi=hemi
return(bwasFormatted)
}


#' Add coordinates to formatted association brain map (Subortical)
#'
#' This function opens and format brain association maps to include coordinates (useful for plotting)
#' The brain association map requires a ProbeID column, which contains the name of vertices/voxels.
#' It opens the relevant atlas file to merge with the association map.
#'
#' @param bwasFormatted Formatted brain association map (created using formatBWAScortical)
#' @param hemi Hemisphere ("lh" or "rh")
#' @param moda modality ("LogJacs" or "thick")
#' @return An enriched brain association maps which contains spatial coordinates
#' @import plyr png qqman readr Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices
#' @export
bwasAddCoordinatesSubcortical=function(bwasFormatted, hemi, moda){
# Get vertices coordinates
coor=read.table(system.file("extdata/atlas", "Atlas_coordinates_subortical.txt", package = "brainMapR",mustWork = TRUE), stringsAsFactors = F, header=T)
# Add coordinates
bwasFormatted=merge(bwasFormatted, coor, by="ProbeID")
return(bwasFormatted)
}
