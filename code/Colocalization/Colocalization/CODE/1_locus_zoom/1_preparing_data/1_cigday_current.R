##############
#INTRODUCTION#
##############

#In this code we are going to prepare the data for the loci we want to co-localize with SHBG, cortisol, current smokers and WHRadjBMI for current smokers.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
#library(TwoSampleMR)

##################
#Loading function#
##################

loci_aligner <- function(cigday_ss, other_ss){
  
  #We are gonna put the example so that we can know what happens:
  
  #cigday_ss <- cigday_loci_aligned
  #other_ss <- whradjbmi_clean
  
  #STEP 0: let's take into account we do not have EAF for SHBG:
  
  check <- which(colnames(other_ss) == "eaf.outcome")
  
  if(is_empty(check)){
    
    other_ss$eaf.outcome <- NA
    
  }
  
  #STEP 1: match the data. This will probably fail with tri-allelic SNPs. Here I think we are gonna be OK...
  
  harm_data <- harmonise_data(cigday_ss, other_ss, action=1) #action = 1 because we do not wanna remove any kind of SNP.
  
  harm_data_clean <- harm_data[which(harm_data$remove == FALSE),]
  
  #And we remove the duplicates:
  
  harm_data_clean <- harm_data_clean[which(duplicated(harm_data_clean$SNP) == FALSE),]

  #STEP 2: reorganize the dataframe, cuz we need something clean:
  
  harm_data_end <- harm_data_clean %>%
    select("rsid.x", "effect_allele.outcome", "other_allele.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "SNP")
  
  colnames(harm_data_end) <- c("variant","effect_allele", "other_allele", "beta", "standard_error", "p_value", "sample_size", "chr_pos")
  
  return(harm_data_end)
  
}


##############
#Loading data#
##############

#We are going to load the data for all the traits.
#the SNP that we found that the causal effect is significant: rs1317286.

cigday <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/CigDayCurrent/Cig_Day_Current_Curated_FULL.txt")

##################################
#STEP 1: obtain the right columns#
##################################

cigday$SNP <- paste("chr", cigday$chr.exposure, ":", cigday$pos.exposure, sep = "")

cigday$samplesize.exposure <- 33229

cigday_corrected <- cigday %>%
  select("rsid", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "SNP")

colnames(cigday_corrected) <-  c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")

######################
#STEP 2: SAVING DATA:#
######################

fwrite(cigday_corrected, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Colocalization/Colocalization/info_4_loci/1_curated_data_4_locus_zoom/cigday_4_locus_zoom.txt")
