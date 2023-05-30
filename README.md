# Estimating the causal relationship between smoking and abdominal adiposity using Mendelian randomization

The repository consists of two sections:

1) The folder structure and code to reproduce our analysis in the /R folder.
2) The variants and summary statistics used to produce our results, plus a code to easily navigate through them in /data_&_results

## GWAS summary statistics used:

In this publication we are testing whether there is any causal association between smoking traits and fat distribution. To do so we used the latest and largest GWAS summary statistics for each of the traits. 

### Smoking traits:

For smoking traits we decided to use:

Lifetime smoking from Wottonn et al ().
Smoking initiation from Liu et al ().
Cigarettes per Day from Liu et al ().
Age of Smoking from Liu et al ().

### Anthropometric traits:

As explained before, the main analysis are for BMI as the main adiposity trait, though we also ran analysis for other traits that inform about central obesity or abdominal fat accumualation: body fat percentage, waist circumference adjusted for BMI and waist-to-hip ratio adjusted for BMI. 

BMI GWAS summary statistics from Pulit et al can be found here: (https://zenodo.org/record/1251813/files/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1)

WHRadjBMI GWAS summary statistics from Pulit et al can be found here: (https://zenodo.org/record/1251813/files/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1)

WCadjBMI GWAS summary statistics from Shungin et al 2015 can be found here: (https://portals.broadinstitute.org/collaboration/giant/images/7/73/GIANT_2015_WCadjBMI_COMBINED_EUR.txt.gz)

BFP GWAS summary statistics from Elsworth et al 2018 can be found here: (https://gwas.mrcieu.ac.uk/datasets/ukb-b-8909/)

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
