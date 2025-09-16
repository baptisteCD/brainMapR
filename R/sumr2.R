
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
sumR2_regression_univariate=function(inputPath , bwasFile, refPanel, nblock=200, chi2Threshold=80, bwasSampleSize="NMISS", outputPath ){

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

# Get sample size in right format
if (is.character(bwasSampleSize)){
if (is.null(BWASsignif2[1,bwasSampleSize])){
    print(paste0("Sample size (column ", bwasSampleSize, ") not present in the summary statistics, please provide a number of the name of the column"))
} else {
    ssize=BWASsignif2[1,bwasSampleSize]
} }
if (is.numeric(bwasSampleSize)){
    ssize=bwasSampleSize
}

# SumR2 regression
ldresOutAll=NULL
for (panel in refPanel){

    # Remove vertices with no SumR2 (missing from calculations due to missingness or lack of variability in data)
if(length(which(is.na(BWASsignif[,paste0("SumR2_", refPanel)])))>0){
BWASsignif2<-BWASsignif[-which(is.na(BWASsignif[,paste0("SumR2_", refPanel)])),]
}

    print(paste0("Estimating morphometricity using SumR2 regression. Using ", dim(BWASsignif2)[1], " vertices and SumR2 from reference panel: ", panel ))

ldres<-bigsnpr::snp_ldsc(ld_score = BWASsignif2[,paste0("SumR2_", panel)], ld_size = dim(BWASsignif2)[1], chi2 = BWASsignif2$CHI2, sample_size = ssize ,chi2_thr1 = 80, blocks = nblock )
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



#' SumR2 regression - estimate grey-matter correlation (rGM)
#'
#' This function reads two brain association map
#' It annotates the association maps to add SumR2 pre-calculated on reference samples
#' Then performs sumR2 regression to estimate rGM
#'
#' @param inputPath Paths (folders) to the raw brain association maps (outputs of OSCA)
#' @param bwasFile Brain association maps in (OSCA format)
#' @param refPanel Reference panel for the sumR2 -  (default UK Biobank, options include ADNI1, ADNI2GO3, AIBL). If several reference panels are listed, results using each sumR2 is returned.
#' @param nblock Number of blocks used in jacknife to estimate SE  - (default 200)
#' @param chi2Threshold Chi2 threshold to exclude vertices with outlying large association - (default 80)
#' @param bwasSampleSize Optional : If character, the name of the columns in bwas file that contains the sample size (default "NMISS"). If numeric, the sample sizes.
#' @param varConstrained (TRUE/FALSE, default = T) whether variance components estimates should be constrained to be between 0 and 1
#' @return Morphometricity, SE, confidence intervals and pvalue, SumR2 intercept and SE.
#' @import GFA vroom plyr foreach
#' @export
sumR2_regression_bivariate(inputPath = inputPath , bwasFile = bwasFile , refPanel = c("ADNI1", "ADNI2GO3", "AIBL"), nblock=200, chi2Threshold=80, bwasSampleSize=c("NMISS", "NMISS"), varConstrained=T ){

# Open summary statistics
    BWASsignif=NULL
BWASsumstatfilePath = paste0(inputPath[1] ,  bwasFile[1])
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

BWASsumstatfilePath2 = paste0(inputPath[2] ,  bwasFile[2])
BWASsumstatfile2<-  vroom(BWASsumstatfilePath2, show_col_types = F)

BWASsignif=merge(BWASsignif, BWASsumstatfile2, by="Probe", suffixes = c("_1", "_2"))

BWASsignif$z_1=BWASsignif$b_1/BWASsignif$se_1
BWASsignif$z_2=BWASsignif$b_2/BWASsignif$se_2


# Remove vertices with missing association in either bwas
if (length(which(is.na(BWASsignif$p_1)))>0 | length(which(is.na(BWASsignif$p_2)))>0){
BWASsignif<-BWASsignif[-which(is.na(BWASsignif$p_1) | is.na(BWASsignif$p_2) ),]
}

# Get sample sizes in right format
if (is.na(suppressWarnings(as.numeric(bwasSampleSize[1])))){
if (is.null(BWASsignif[1,bwasSampleSize[1]])){
    print(paste0("Sample size (column ", bwasSampleSize[1], ") not present in the summary statistics, please provide a number of the name of the column"))
} else {
    ssize1=BWASsignif[1,bwasSampleSize[1]]
} } else {
    ssize1=as.numeric(bwasSampleSize[1])
}


# Get sample size in right format
if (is.na(suppressWarnings(as.numeric(bwasSampleSize[2])))){
if (is.null(BWASsignif[1,bwasSampleSize[2]])){
    print(paste0("Sample size (column ", bwasSampleSize[1], ") not present in the summary statistics, please provide a number of the name of the column"))
} else {
    ssize2=BWASsignif[1,bwasSampleSize[2]]
} } else {
    ssize2=as.numeric(bwasSampleSize[2])
}

# SumR2 regression - estimate rGM
ldresOutAll=NULL
for (panel in refPanel){

# Remove vertices with no SumR2 (missing from calculations due to missingness or lack of variability in data)
if(length(which(is.na(BWASsignif[,paste0("SumR2_", panel)])))>0){
BWASsignif2<-BWASsignif[-which(is.na(BWASsignif[,paste0("SumR2_", panel)])),]
}

    print(paste0("Estimating rGM using SumR2 regression. Using ", dim(BWASsignif2)[1], " vertices and SumR2 from reference panel: ", panel ))


ldres=GFA::ldsc_rg(ld_score = BWASsignif2[,paste0("SumR2_", panel)], ld_size = dim(BWASsignif2)[1], z1=BWASsignif2$z_1 ,z2=BWASsignif2$z_2 , sample_size_1 = ssize1, sample_size_2=ssize2, blocks=nblock, chi2_thr2=chi2Threshold)


# Constrained or unconstrained variance estimates
if (varConstrained==T){
   if (ldres[7]>1){
       ldres[7]=1
      ldres[11]=ldres[3]/(sqrt(ldres[7])*sqrt(ldres[9]))
      print("Variance of trait 1 constrained to 1") }
    if (ldres[9]>1){
       ldres[9]=1
      ldres[11]=ldres[3]/(sqrt(ldres[7])*sqrt(ldres[9]))
      print("Variance of trait 2 constrained to 1") }
     if (ldres[7]<0){
       ldres[7]=0.00001
      ldres[11]=ldres[3]/(sqrt(ldres[7])*sqrt(ldres[9]))
      print("Variance of trait 1 constrained to 0") }
    if (ldres[9]<0){
       ldres[9]=0.00001
      ldres[11]=ldres[3]/(sqrt(ldres[7])*sqrt(ldres[9]))
      print("Variance of trait 2 constrained to 0") }
}

pv1= 2*pnorm(-abs(ldres[7]/ldres[8]))
cil1=ldres[7]-1.96*ldres[8]
ciu1=ldres[7]+1.96*ldres[8]

pv2= 2*pnorm(-abs(ldres[9]/ldres[10]))
cil2=ldres[9]-1.96*ldres[10]
ciu2=ldres[9]+1.96*ldres[10]

pvrg= 2*pnorm(-abs(ldres[11]/ldres[12]))
cilrg=ldres[11]-1.96*ldres[12]
ciurg=ldres[11]+1.96*ldres[12]


ldresOut=t(as.data.frame(c(panel,c(ldres[1:8]), pv1, cil1, ciu1, c(ldres[9:10]), pv2, cil2, ciu2, c(ldres[11:12]), pvrg, cilrg, ciurg)))

colnames(ldresOut)=c("sumR2RefPanel", "int",       "int_se",       "braincov",    "braincov_se",       "t1_int",       "t2_int",         "m2_1", "m2_1_se", "pvalue_1",   "m2_1_CI_lb", "m2_1_CI_ub",  "m2_2", "m2_2_se", "pvalue_2",   "m2_2_CI_lb", "m2_2_CI_ub", "rGM",   "rGM_se", "pvalue_rGM",   "rGM_CI_lb", "rGM_CI_ub"  )
rownames(ldresOut)=NULL

ldresOutAll=rbind(ldresOutAll, ldresOut)
}

return(ldresOutAll)
}

