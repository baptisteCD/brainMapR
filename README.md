# brainMapR - an R package to analyse and visualise brain association maps
 
<br>

Left cortical thickness (regions associated with age)

![Alt Text](https://github.com/baptisteCD/brainMapR/blob/main/lh_thickness.gif)

<br><br>

Left subcortical thickness (regions associated with age)

![Alt Text](https://github.com/baptisteCD/brainMapR/blob/main/lh_thick_combined.gif)

All surface styles from FreeSurfer are available, and handy one-line gif functions are implemented

![Alt Text](https://github.com/baptisteCD/brainMapR/blob/main/BWAS_BWAS_Age_MRI.linear_lh_thickness_AllStylesPial_leftview.gif)


ROI plots are also possible

![Alt Text](https://github.com/baptisteCD/brainMapR/blob/main/ROIbasedAD_vs._HC_ST_Independent_ROI_smoothwm)

--> See code example and usage 

--> [**Tutorial and examples for analysing BWAS results, annotating maps, making plots and GIFs e**](https://baptistecd.github.io/Brain-Mapping-LMM/RR_9_Extractresults_RealPhenotypes.html).

--> [**Tutorial and examples for SumR2 regression**](https://baptistecd.github.io/Brain-Mapping-LMM/RR_11_SumR2_wiki.html).



If you use brainMapR to plot and annotate BWAS, please cite the relevant publication(s)
[A unified framework for association and prediction from vertex-wise grey-matter structure](https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.25109)
[Parsimonious model for mass-univariate vertexwise analysis](https://doi.org/10.1117/1.JMI.9.5.052404)
[Grey-Matter Structure Markers of Alzheimer's Disease, Alzheimer's Conversion, Functioning and Cognition: A Meta-Analysis Across 11 Cohorts]( https://doi.org/10.1002/hbm.70089)

If you use brainMapR to estimate morphometricity and grey-matter correlation, please cite:
[Fast and Efficient Estimation of Morphometricity and Grey-Matter Correlation from Neuroimaging Summary Statistics](https://hal.science/hal-05356974v1/file/Paper_VF.pdf)

--> Version updates

2025: Version brainMapR_1.1.0.9000 - adding SumR2 regression plots

2025: Version brainMapR_1.0.0.9000 - now includes SumR2 regression for the brain surface, and ROI plots

2023: Version brainMapR_0.8.0.9000 - first release


<br><br>
Overview
--------

brainMapR is a package to analyse and plot brain association maps (BWAS results). It is tailored for brain MRI vertex-wise analyses, extracted from T1w images. brainMapR requires brain MRI processed with [**FreeSurfer**](https://surfer.nmr.mgh.harvard.edu/) (for cortical vertices) and/or [**ENIGMA-shape package**](https://enigma.ini.usc.edu/ongoing/enigma-shape-analysis/) (for subcortical vertices). 

For more details about the processing and the analyses - see my two previous repositories [brain-LMM](https://baptistecd.github.io/Brain-LMM/) and [brain-mapping-LMM](https://baptistecd.github.io/Brain-Mapping-LMM/).
 
This package accompanies the publications:   
[A unified framework for association and prediction from vertex-wise grey-matter structure](https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.25109)

[Parsimonious model for mass-univariate vertexwise analysis](https://doi.org/10.1117/1.JMI.9.5.052404)

[Grey-Matter Structure Markers of Alzheimer's Disease, Alzheimer's Conversion, Functioning and Cognition: A Meta-Analysis Across 11 Cohorts]( https://doi.org/10.1002/hbm.70089)
 

<br><br>
Getting started   
-------- 

- Install the R package:
```
library(devtools)
devtools::install_github("baptisteCD/brainMapR")
```
- Load package
```
library(brainMapR)
```