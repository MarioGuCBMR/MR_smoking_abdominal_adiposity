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

cigday <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/RAW_Data/CigarettesPerDay.txt.gz")

##################################
#STEP 1: obtain the right columns#
##################################

cigday$SNP <- paste("chr", cigday$CHROM, ":", cigday$POS, sep = "")

cigday_corrected <- cigday %>%
  select("RSID", "CHROM", "POS", "ALT", "REF", "AF", "BETA", "SE", "PVALUE", "N", "SNP")

colnames(cigday_corrected) <-  c("variant","chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")

######################
#STEP 2: SAVING DATA:#
######################

fwrite(cigday_corrected, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Colocalization/Colocalization/info_4_loci/1_curated_data_4_locus_zoom/cigday-combined_4_locus_zoom.txt")
