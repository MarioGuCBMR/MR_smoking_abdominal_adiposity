##############
#INTRODUCTION#
##############

#In this code we are going to prepare the data for the loci we want to co-localize with SHBG, cortisol, current smokers and WHRadjBMI for current smokers.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

#We are going to load the data for all the traits.
#the SNP that we found that the causal effect is significant: rs1317286.

whradjbmi <- fread("C:/Users/zlc436/Desktop/PhD/copy_data/1_curated_data/whradjbmi_curated.txt")

######################
#STEP 2: SAVING DATA:#
######################

fwrite(whradjbmi, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Colocalization/Colocalization/info_4_loci/1_curated_data_4_locus_zoom/whradjbmi_4_locus_zoom.txt")
