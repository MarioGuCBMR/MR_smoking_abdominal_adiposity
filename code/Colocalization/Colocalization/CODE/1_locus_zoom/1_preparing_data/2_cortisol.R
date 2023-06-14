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

cortisol <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Colocalization/Curating_data/Curated_data/cortisol_curated.txt")

##################################
#STEP 1: obtain the right columns#
##################################

cortisol$eaf.outcome <- NA

cortisol_corrected <- cortisol %>%
  select("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "chr_pos")

colnames(cortisol_corrected) <-  c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")

################################################################
#STEP 2: be careful that we need the alleles in capital letters#
################################################################

cortisol_corrected$effect_allele <- toupper(cortisol_corrected$effect_allele)
cortisol_corrected$other_allele <- toupper(cortisol_corrected$other_allele)

######################
#STEP 2: SAVING DATA:#
######################

fwrite(cortisol_corrected, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Colocalization/Colocalization/info_4_loci/1_curated_data_4_locus_zoom/cortisol_4_locus_zoom.txt")
