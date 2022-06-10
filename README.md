# brainMapR - an R package to analyse and visualise brain association maps
 
<br>
<br>

![Cortical plots and GIF](https://github.com/baptisteCD/brainMapR/lh_thickness.gif)

![Subcortical plots and GIF](https://github.com/baptisteCD/brainMapR/lh_thick_combined.gif)


--> See code example and usage [**here**](https://baptistecd.github.io/Brain-Mapping-LMM/RR_9_Extractresults_RealPhenotypes.html).

<br>
Overview
--------

brainMapR is a package to analyse and plot brain association maps (BWAS results). It is tailored for brain MRI vertex-wise analyses, extracted from T1w images. brainMapR requires brain MRI processed with [**FreeSurfer**](https://surfer.nmr.mgh.harvard.edu/) (for cortical vertices) and/or [**ENIGMA-shape package**](https://enigma.ini.usc.edu/ongoing/enigma-shape-analysis/) (for subcortical vertices). 

For more details about the processing and the analyses - see my two previous repositories [brain-LMM](https://baptistecd.github.io/Brain-LMM/) and [brain-mapping-LMM](https://baptistecd.github.io/Brain-Mapping-LMM/).
 
This package accompanies two publications:   
[A unified framework for association and prediction from vertex-wise grey-matter structure](https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.25109)

[Parsimonious model for mass-univariate vertexwise analysis](https://doi.org/10.1117/1.JMI.9.5.052404)
 
<br>
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