##############
#INTRODUCTION#
##############

#This is a code to do some checks on the TwoSampleMR thingie with both: TDI and Cigarettes per Day.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(TwoSampleMR)

##############
#Loading data#
##############

CigDay_source <- fread("C:/Users/zlc436/Desktop/SMK_WHRAdjBMI/CigarettesPerDay.txt.gz")

CigDay_check <- extract_outcome_data(
  snps = CigDay_source$RSID[seq(1,3)],
  outcomes = 'ieu-b-25',
  proxies = TRUE
)

#The data matches perfectly. I don't think that the analysis we did were wrong.
#In any case, let's continue to do some good ol' checks. But this time,
#we are gonna do it with TDI. 
#We are gonna download the data in VCF format. 
#Then we are gonna compare the data to that obtained with extract_outomce_data.


