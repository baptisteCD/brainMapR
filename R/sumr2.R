
#' SumR2 regression - estimate morphometricity
#'
#' This function reads a brain association map
#' It annotates the association map to add SumR2 pre-calculated on reference samples
#' Then performs sumR2 regression to estimate morphometricity
#'
#' @param inputPath path (folder) to the raw brain association maps (outputs of OSCA)
#' @param bwasFile Brain association map in (OSCA format)
#' @param refPanel Reference panel for the sumR2 -  (default UK Biobank, options include ADNI1, ADNI2GO3, AIBL, AVERAGE). If several reference panels are listed, results using each sumR2 is returned.
#' @param nblock Number of blocks used in jacknife to estimate SE  - (default 200)
#' @param chi2Threshold Chi2 threshold to exclude vertices with outlying large association - (default 80)
#' @param bwasSampleSize Optional : If character, the name of the columns in bwas file that contains the sample size (default "NMISS"). If numeric, the sample size.
#' @param outputPath path where the outputs will be written
#' @return Morphometricity, SE, confidence intervals and pvalue, SumR2 intercept and SE.
#' @import bigsnpr vroom plyr
#' @export
sumR2_regression_univariate=function(inputPath , bwasFile, refPanel, nblock=200, chi2Threshold=80, bwasSampleSize="NMISS" ){

# Open summary statistics
    BWASsignif=NULL
BWASsumstatfilePath = paste0(inputPath ,  bwasFile)
BWASsumstatfile<-  vroom(BWASsumstatfilePath, show_col_types = F)
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
# rbind all vertices with annotations
BWASsignif=plyr::rbind.fill(BWASsignif, BWASsumstat)
  } }


# Get rid of vertices not included in the analysis
if(length(which(is.na(BWASsignif$p)))>0){
BWASsignif<-BWASsignif[-which(is.na(BWASsignif$p)),]
}
# Remove vertices with no SumR2 (missing from calculations due to missingness or lack of variability in data)
if(length(which(is.na(BWASsignif[,paste0("SumR2_", refPanel)])))>0){
BWASsignif<-BWASsignif[-which(rowSums(is.na(BWASsignif[,grep(pattern = "SumR2_", x = colnames(BWASsignif))]))>0),]
}

if (is.character(bwasSampleSize)){
if (is.null(BWASsignif[1,bwasSampleSize])){
    print(paste0("Sample size (column ", bwasSampleSize, ") not present in the summary statistics, please provide a number of the name of the column"))
} else {
    ssize=BWASsignif[1,bwasSampleSize]
} }
if (is.numeric(bwasSampleSize)){
    ssize=bwasSampleSize
}

ldresOutAll=NULL
for (panel in refPanel){

    print(paste0("Estimating morphometricity using SumR2 regression. Using ", dim(BWASsignif)[1], " vertices and SumR2 from reference panel: ", panel ))

ldres<-bigsnpr::snp_ldsc(ld_score = BWASsignif[,paste0("SumR2_", panel)], ld_size = ldsize, chi2 = BWASsignif$CHI2, sample_size = ssize ,chi2_thr1 = 80)
pv= 2*pnorm(-abs(ldres[3]/ldres[4]))
cil=ldres[3]-1.96*ldres[4]
ciu=ldres[3]+1.96*ldres[4]

ldresOut=t(as.data.frame(c(panel,ldres, pv, cil, ciu  )))
colnames(ldresOut)=c("sumR2RefPanel","int" , "int_se", "m2", "m2_se",  "pvalue",   "m2_CI_lb", "m2_CI_ub"   )
rownames(ldresOut)=NULL
ldresOutAll=rbind(ldresOutAll, ldresOut)
}

write.table(ldresOutAll, paste0(outputPath, "/SumR2_regression_univariate_" , bwasFile, ".rsq"), col.names = T, row.names = F)
return(ldresOutAll)

}
