
#' Open and format association brain map (Cortical)
#'
#' This function opens and format brain association maps to include coordinates and ROI (Desikan atlas). All this information is required to produce Manhattan plots.
#' The brain association map requires a ProbeID column, which contains the name of vertices/voxels.
#' It opens the relevant atlas file (incl. in the R package) to annotate the association map.
#'
#' @param BWASsumstat Brain association map (OSCA format)
#' @param hemi Hemisphere ("lh" or "rh")
#' @param mod modality ("area" or "thickness")
#' @return An enriched brain association maps which contains more vertex information, including corresponding grey-matter region, 3D coordinates, and ad hoc information for plotting.
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices vroom
#' @export
formatBWAScortical<-function(BWASsumstat, hemi, mod){

# Open atlas with details about FreeSurfer fsaverage vertices
atlas<-vroom(system.file("extdata/atlas", "Atlas_coordinates_corticalVertices.txt", package = "brainMapR",mustWork = TRUE), show_col_types = FALSE)

# Merge atlas and summary statistics
BWASsumstat<-merge(BWASsumstat, atlas, by.x="Probe", by.y="vertexID" )

# Subset table, sort by ROI and add xax which is the order of vertices for Manhattan plotting
BWASsumstat<-BWASsumstat[which(BWASsumstat$hemi==hemi & BWASsumstat$moda==mod),]
BWASsumstat<-BWASsumstat[order(BWASsumstat$order),]
BWASsumstat$xax<-1:length(BWASsumstat[,1])

# Give some details to user
print(paste(length(BWASsumstat$Probe), "vertices annotated") )

return(BWASsumstat)
}

#' Open and format association brain map (Subcortical)
#'
#' This function opens and format brain association maps to include coordinates and ROI (from ENIGMA-shape atlas). All this information is required to produce Manhattan plots.
#' The brain association map requires a ProbeID column, which contains the name of vertices/voxels.
#' It opens the relevant atlas file (incl. in the R package) to annotate the association map.
#'
#' @param BWASsumstat Brain association map (OSCA format)
#' @param hemi Hemisphere ("lh" or "rh")
#' @param mod modality ("LogJacs" or "thick")
#' @return An enriched brain association maps which contains more vertex information
#' @import plyr png qqman Rvcg rgl RColorBrewer grid gridExtra viridis Morpho ggplot2 utils stats graphics grDevices vroom
#' @export
formatBWASsubcortical<-function(BWASsumstat, hemi, mod){

# Open atlas with vertex information
atlas<-vroom(system.file("extdata/atlas", "Atlas_coordinates_ROI_subortical.txt", package = "brainMapR",mustWork = TRUE), show_col_types = FALSE)

# Merge with summary statistics
BWASsumstat<-merge(BWASsumstat, atlas, by.x="Probe", by.y="ProbeID" )

# Order by ROI and create xax which will serve as the x axis in Manhattan plot.
BWASsumstat<-BWASsumstat[which(BWASsumstat$hemi==hemi & BWASsumstat$moda==mod),]
BWASsumstat<-BWASsumstat[order(BWASsumstat$ROINb),]
BWASsumstat$xax<-1:length(BWASsumstat[,1])

# Give some details to user
print(paste(length(BWASsumstat$Probe), "vertices annotated") )

return(BWASsumstat)
}

