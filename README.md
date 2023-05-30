# Estimating the causal relationship between smoking and abdominal adiposity using Mendelian randomization

This repository contains the contain to reproduce the MR analysis to assess the associations between smoking and fat distribution.

## Software used in the analysis:

A total of 4 Mendelian Randomization Methods were used:

Traditional mendelian randomization pipeline using TwoSampleMR package (v4.6: https://mrcieu.github.io/TwoSampleMR/)
CAUSE-MR (v1.0.0, install it here: https://github.com/jean997/cause)
LHC-MR (v1.0.0, install it from here: https://github.com/LizaDarrous/lhcMR)
MVMR (vX, install it from here: https://github.com/WSpiller/MVMR)

Additionally, genetic risk scores of IVs with fat depots and hormonal GWAS summary statistics were performed with the function grs.summary from the package (gtx, v.0.8.0). This function calculates the joint effects of a set of variants on an outcome trait using only GWAS summary statistics from the exposure and the outcome.
The package can be downloaded from here: https://github.com/tobyjohnson/gtx

### Packages requiered:

Importantly, all code for 2SMR codes start by loading several libraries: 

```
All of the following can be downloaded using install.packages() function:

library(ggplot2)
library(ggrepel)
library(dplyr)
library(gridExtra)
library(rmarkdown)
library(data.table)
library(jsonlite)
library(httr)
library(tidyverse)
library(phenoscanner)
```

### Additional information:

2SMR, MVMR and GRS code were run locally in a Windows 10 computer. CAUSE and LHC-MR were run in Computerome C2 (https://www.computerome.dk/).

## GWAS summary statistics used:

In this publication we are testing whether there is any causal association between smoking traits and fat distribution. To do so we used the latest and largest GWAS summary statistics for each of the traits. 

### Smoking traits:

For smoking traits we decided to use:

Lifetime smoking from Wottonn et al 2020 (available here: https://data.bris.ac.uk/data/dataset/10i96zb8gm0j81yz0q6ztei23d)
Smoking initiation from Liu et al 2019 (available here: https://conservancy.umn.edu/handle/11299/201564).
Cigarettes per day (past and current smokers) from Liu et al (available here: https://conservancy.umn.edu/handle/11299/201564).
Cigarettes per day (current smokers) from Elseworht et al (available here: )
Packs-year from Elseworth et al (available here: ).

For detailed information on the curation of the data and selection of IVs for each MR method, please check the Methods and Supplementary information of each paper.

### Fat distribution traits:

As outcome we used the following fat distribution traits: waist-hip ratio (WHR), body mass index adjusted WHR (WHRadjBMI), waist circumference (WC), BMI-adjusted WC (WCadjBMI), hip circumference (HC), BMI-adjusted hip circumference (HCadjBMI), WHRadjBMI statified for smokers and WCadjBMI stratified for smokers. 

The GWAS summary statistics utilized for each analysis might differ since traditional MR analysis are sensitive to sample overlap, while CAUSE and LHC-MR can control for it.

For 2SMR the summary statistics used are:

WHR: Shungin et al 2015 (available at: )
WHRadjBMI: Shunging et al 2015 (available at: )
WC: Shunging et al 2015 (available at: )
HC: Shungin et al 2015 (available at: )
WCadjBMI: Shunging et al 2015 (available at: )
HCadjBMI: Shunging et al 2015 (available at: )
WHRadjBMI stratified for smokers: Justice et al 2017 (available at: )
WCadjBMI stratified for smokers: Justice et al 2017 (available at: )

For CAUSE and LHC-MR the summary statistics we particularly used GWAS summary statistics with bigger sample sizes for HC, WC, WHR, and WHRadjBMI:

WC Elseworth et al 2018 (available at: )
HC: Elseworth et al 2018 (available at: )
WHR : Pulit et al 2019: (available: https://zenodo.org/record/1251813/files/whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1)
WHRadjBMI: Pulit et al 2019: (availbale: https://zenodo.org/record/1251813/files/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1)
