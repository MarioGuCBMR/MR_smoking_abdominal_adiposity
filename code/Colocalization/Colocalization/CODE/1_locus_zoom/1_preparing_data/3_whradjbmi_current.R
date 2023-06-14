##############
#INTRODUCTION#
##############

#In this code we are going to prepare the data for the loci we want to co-localize with SHBG, cortisol, current smokers and WHRadjBMI for current smokers.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

chr_parser <- function(chr_pos){
  
  tmp <- strsplit(chr_pos, ":")[[1]][1]
  chr_ <- strsplit(tmp, "chr")[[1]][2]
  
  return(chr_)
  
}

pos_parser <- function(chr_pos){
  
  pos <- strsplit(chr_pos, ":")[[1]][2]
  
  return(pos)
  
}

##############
#Loading data#
##############

#We are going to load the data for all the traits.
#the SNP that we found that the causal effect is significant: rs1317286.

whradjbmi <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMIsmk/WHRAdjBMI_Smk_Strat_Curated.txt")

##################################
#STEP 1: obtain the right columns#
##################################

whradjbmi$pos <- as.numeric(as.character(unlist(sapply(whradjbmi$chr_pos_37, pos_parser))))

whradjbmi_corrected <- whradjbmi %>%
  select("rs_id", "chromosome", "pos", "Effect_allele", "Other_allele", "EAF_HapMapCEU", "Effect_SMK", "StdErr_SMK", "P_value_SMK", "N_SMK", "chr_pos_37")

colnames(whradjbmi_corrected) <-  c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")

######################
#STEP 2: SAVING DATA:#
######################

fwrite(whradjbmi_corrected, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Colocalization/Colocalization/info_4_loci/1_curated_data_4_locus_zoom/whradjbmi-current_4_locus_zoom.txt")
