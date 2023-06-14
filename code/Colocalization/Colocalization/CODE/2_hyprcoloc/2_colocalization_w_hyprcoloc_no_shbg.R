##############
#INTRODUCTION#
##############

#In this code we are going to prepare the data for the loci we want to co-localize with SHBG, cortisol, current smokers and WHRadjBMI for current smokers.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(TwoSampleMR)

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

#We are going to load the data for all the traits for the loci that we are interested in:

coloc_df <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Colocalization/Colocalization/info_4_loci/1_ss_of_all_traits/15_78396129_79396129_no_shbg.txt")

########################################################
#STEP 1: let's define the loci that we are going to use#
########################################################

trait.names <- c("cigday", "whradjbmi",  "cortisol")

betas <- coloc_df %>%
  select(BETA.cigday, BETA.whradjbmi, BETA.cortisol)

print(dim(betas))

betas <- as.matrix(betas)
rownames(betas) <- coloc_df$variant
colnames(betas) <-  c("cigday", "whradjbmi", "cortisol")

ses <- coloc_df %>%
  select(SE.cigday, SE.whradjbmi, SE.cortisol)

print(dim(ses))

ses <- as.matrix(ses)
rownames(ses) <- coloc_df$variant
colnames(ses) <-c("cigday", "whradjbmi", "cortisol")

res <- hyprcoloc::hyprcoloc(betas, ses, trait.names=trait.names, snp.id=coloc_df$variant)

#$results
#iteration traits posterior_prob regional_prob candidate_snp posterior_explained_by_snp dropped_trait
#1         1   None             NA        0.0026            NA                         NA        cigday
#2         2   None             NA        0.0001            NA                         NA      cortisol
#
#attr(,"class")
#[1] "hyprcoloc"

#No co-localization according to coloc.
