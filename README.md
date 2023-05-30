# Estimating the causal relationship between smoking and abdominal adiposity using Mendelian randomization

The repository consists of two sections:

1) The folder structure and code to reproduce our analysis in the /R folder.
2) The variants and summary statistics used to produce our results, plus a code to easily navigate through them in /data_&_results

## GWAS summary statistics used:

In this publication we are testing whether there is any causal association between smoking traits and fat distribution. To do so we used the latest and largest GWAS summary statistics for each of the traits. 

### Smoking traits:

For smoking traits we decided to use:

Lifetime smoking from Wottonn et al 2020 (available here: )
Smoking initiation from Liu et al 2019 (available here: ).
Cigarettes per day (past and current smokers) from Liu et al (available here: ).
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

### Curation of GWAS summary statistics:

For each analysis, we curated the data according to the needs of each software. While we won't be adding the curated data in the github due to their big file sizes, you can either:
1) replicate our analysis by having all GWAS summary statistics in the RAW_DATA folder.
2) ask for a zenodo link where all the curated data will be stored. We will only provide this is the replication through step 1 is not successful **(NOTE: ideally this won't happen, but sources like PhenoScanner, SNPNexus or dbSNP are constantly being updated, so results might differ if a long time passes between the original analysis and the replication).**

**All code in this github relies on having the raw data from these sources in the RAW_DATA folder**

## Running CAUSE:

The analysis for this publication started in 03/2020 and, thus, we used the first version of CAUSE, which the authors used to publish in Nature Communication: (link.) Thus, to reproduce our results you should have the versions of the following packages installed:

```
devtools::install_version("mixsqp", version = "0.1-97", repos = "http://cran.us.r-project.org")
devtools::install_version("ashr", version = "2.2-32", repos = "http://cran.us.r-project.org")
devtools::install_github("jean997/cause@v1.0.0")

## Running 2SMR:

The 2SMR analysis are sensitive to sample overlap, so we decided to use GIANT BMI summary statistics, since moderate and vigorous physical activity and sedentary time GWAS are from UK Biobank, exclusively. You can find the original GIANT BMI summary statistics here: https://portals.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz

To run the analysis following the code in the folder /R/2SMR/running_2SMR, you will need to curate the data using the code in /R/2SMR/curating_data_4_2SMR/.

### Packages requiered

All 2SMR codes start by loading several libraries. 

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

*Though TwoSampleMR needs remotes:*
```
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR") #analysis were performed with version 4.26
```
