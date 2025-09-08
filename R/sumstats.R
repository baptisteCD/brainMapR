
#' SumR2 regression
#'
#' This function reads a brain association map
#' It annotates the association map to add SumR2 pre-calculated on reference samples
#' Then performs sumR2 regression to estimate morphometricity
#'
#' @param inputPath path (folder) to the raw brain association maps (outputs of OSCA)
#' @param bwas Brain association map in (OSCA format)
#' @param refPanel Reference panel for the sumR2 -  (default UK Biobank, options include ADNI1, ADNI2GO3, AIBL). If several reference panels are listed, results using each sumR2 is returned.
#' @param nblock Number of blocks used in jacknife to estimate SE  - (default 200)
#' @param chi2Threshold Chi2 threshold to exclude vertices with outlying large association - (default 80)
#' @return Morphometricity, SE, confidence intervals and pvalue, SumR2 intercept and SE.
#' @import data.table bigsnpr
#' @export
sumR2_regression_univariate=function(inputPath , bwas, refPanel, nblock=200, chi2Threshold=80 ){

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


ldres<-bigsnpr::snp_ldsc(ld_score = dat[,"SumR2_full"], ld_size = ldsize, chi2 = dat$chi2, sample_size = dat$NMISS[1] ,chi2_thr1 = 80)
pv= 2*pnorm(-abs(ldres[3]/ldres[4]))
cil=ldres[3]-1.96*ldres[4]
ciu=ldres[3]+1.96*ldres[4]


ldresOut=t(as.data.frame(c(ldres, pv, cil, ciu  )))
colnames(ldresOut)=c("int" , "int_se", "m2", "m2_se",  "pvalue",   "m2_CI_lb", "m2_CI_ub"   )
rownames(ldresOut)=NULL

return(ldresOut)
}


### Add pvalue and confidence intervals
### Add SumR2 (multiple ones) in atlas data - to have option
### Add option to add plot (default is yes)
### Add flag about large SE? - based on theory - i.e. SE twice as large as theory - check SumR2 plot and consider REML analysis, or varying chi2Threshold
