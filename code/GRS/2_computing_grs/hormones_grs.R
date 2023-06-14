
source("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/code/0_functions/functions_4_mediation.R")

smoking_aligner <- function(query_ss, other_ss){
  
  #We are gonna put the example so that we can know what happens:
  
  #fiadjbmi_ss <- exp_df_found_pos
  #other_ss <- hdl_005
  
  #STEP 0: let's run the matching with TwoSampleMR so we need the data:
  
  #Let's first check if we have the effect_allele_frequency column:
  
  check <- which(colnames(query_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    query_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  exposure <- query_ss %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  
  colnames(exposure) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "rsid.exposure")
  
  exposure$exposure <- "smoking"
  exposure$id.exposure <- "smoking"
  
  #Now with the outcome:
  
  check <- which(colnames(other_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    other_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  outcome <- other_ss %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  
  colnames(outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "rsid.outcome")
  
  outcome$outcome <- "idps"
  outcome$id.outcome <- "idps"
  
  #STEP 1: match the data. This will probably fail with tri-allelic SNPs. Here I think we are gonna be OK...
  
  merged_df <- harmonise_data(exposure, outcome, action=3)
  print(dim(merged_df))
  
  merged_df <- merged_df[which(merged_df$remove == FALSE),] #removing incompatible alleles
  merged_df <- merged_df[which(merged_df$palindromic == FALSE & merged_df$ambiguous == FALSE | merged_df$palindromic == TRUE & merged_df$ambiguous == FALSE),]
  
  #merged_df <- merged_df[which(merged_df$mr_keep),] #removing incompatible alleles
  
  #I checked that all is working great. FIadjBMI betas is positive. The rest is great.
  
  #STEP 3: reorganize the dataframe, cuz we need something clean:
  
  other_ss_aligned <- merged_df %>%
    select("rsid.outcome", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "SNP")
  
  colnames(other_ss_aligned) <- c("variant", "chromosome", "base_pair_location",  "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")
  
  return(other_ss_aligned)
  
}

grs_producer <- function(query_ss, list_of_traits_, type_of_trait_, trait_names_, weighted=FALSE){
  #This code produces GRS for that that has already been matched:
  
  for(index_trait in seq(1, length(list_of_traits_))){
    
    print(trait_names_[index_trait])
    
    aligned_outcome <- as.data.frame(list_of_traits_[index_trait])
    
    #Finally produce the genetic risk scores:
    
    if(weighted == TRUE){
      
      weights <- query_ss$beta
      
      grs_res <- gtx::grs.summary(weights, aligned_outcome$beta, aligned_outcome$standard_error, aligned_outcome$sample_size)
      
      print(grs_res)
      
    } else {
      
      weights <- rep(1, length(aligned_outcome$beta))
      
      grs_res <- gtx::grs.summary(weights, aligned_outcome$beta, aligned_outcome$standard_error, aligned_outcome$sample_size)
      
      print(grs_res)
      
    }
    
    #Let's make the results for the GRS parseable:
    
    #Finally, merge the results:
    
    if(index_trait == 1){
      
      beta_ <- grs_res$ahat
      se_ <- grs_res$aSE
      lower_ci_ <- lower_ci(beta_, se_)
      upper_ci_ <- upper_ci(beta_, se_)
      pval_ <- grs_res$pval
      trait_ <- trait_names_[index_trait]
      type_ <- type_of_trait_[index_trait]
      nsnps_ <- length(aligned_outcome$chr_pos)
      
      grs_vect <- c(beta_, se_, lower_ci_, upper_ci_, pval_, trait_, type_, nsnps_)
      grs_end <- as.data.frame(t(grs_vect))
      colnames(grs_end) <- c("beta", "se", "lower_ci", "upper_ci", "pval", "trait", "type", "nsnps")
      
    } else {
      
      beta_ <- grs_res$ahat
      se_ <- grs_res$aSE
      lower_ci_ <- lower_ci(beta_, se_)
      upper_ci_ <- upper_ci(beta_, se_)
      pval_ <- grs_res$pval
      trait_ <- trait_names_[index_trait]
      type_ <- type_of_trait_[index_trait]
      nsnps_ <- length(aligned_outcome$chr_pos)
      
      grs_vect <- c(beta_, se_, lower_ci_, upper_ci_, pval_, trait_, type_, nsnps_)
      grs_end <- rbind(grs_end, grs_vect)
      
    }
    
  }
  
  #Once everything has been looped we can just return the dataframe
  
  return(grs_end)
  
}

computing_grs <- function(smoking_data, 
                          trait1, trait2, trait3,
                          name1, name2, name3,
                          unweighted_out, weighted_out){
  
  smoking_initiation <- fread(smoking_data)
  
  ######################
  #Loading outcome data#
  ######################
  
  vat <- fread(trait1)
  asat <- fread(trait2)
  gfat <- fread(trait3)
  
  ##############################################
  #STEP 1: matching data for smoking initiation#
  ##############################################
  
  #First volume:
  
  vat_2_initiation <- smoking_aligner(smoking_initiation, vat) #76
  asat_2_initiation <- smoking_aligner(smoking_initiation, asat) #76
  gfat_2_initiation <- smoking_aligner(smoking_initiation, gfat) #76
  
  #Let's match the smoking initiation variants and order all of our variants to that:
  
  smoking_initiation_match <- smoking_initiation[which(smoking_initiation$variant%in%vat_2_initiation$variant),] #76: perfect match
  
  #Let's order all variants:
  
  smoking_initiation_match_ordered <- smoking_initiation_match[order(smoking_initiation_match$variant),]
  
  vat_2_initiation_ordered <- vat_2_initiation[order(vat_2_initiation$variant),]
  asat_2_initiation_ordered <- asat_2_initiation[order(asat_2_initiation$variant),]
  gfat_2_initiation_ordered <- gfat_2_initiation[order(gfat_2_initiation$variant),]
  
  #Let's do some checks on effect alleles:
  
  print(length(which(smoking_initiation_match_ordered$effect_allele == vat_2_initiation_ordered$effect_allele))) #good
  print(length(which(smoking_initiation_match_ordered$effect_allele == asat_2_initiation_ordered$effect_allele))) #good
  print(length(which(smoking_initiation_match_ordered$effect_allele == gfat_2_initiation_ordered$effect_allele))) #good
  
  #######################
  #STEP 2: computing GRS#
  #######################
  
  list_of_traits <- list(vat_2_initiation_ordered, 
                         asat_2_initiation_ordered, 
                         gfat_2_initiation_ordered) 
  
  trait_names <- c(name1, name2, name3)
  type_of_trait <- c(name1, name2, name3)
  
  smoking_initiation_smk_grs_unw <- grs_producer(query_ss = smoking_initiation_match_ordered, list_of_traits_ = list_of_traits, trait_names_ = trait_names,  type_of_trait_ = type_of_trait, weighted = FALSE) 
  
  fwrite(smoking_initiation_smk_grs_unw, unweighted_out)
  
  smoking_initiation_smk_grs_w <- grs_producer(query_ss = smoking_initiation_match_ordered, list_of_traits_ = list_of_traits, trait_names_ = trait_names,  type_of_trait_ = type_of_trait, weighted = TRUE) 
  
  fwrite(smoking_initiation_smk_grs_w, weighted_out)
}


# DATA
smoking_data = smoking_initiation_whr
trait1 = cortisol
trait2 = testosterone
trait3 = shbg

smoking_initiation_whr <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/output/1_curated_data/smoking_initiation_whr_curated.txt"
lifetime_smoking_whr <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/output/1_curated_data/lifetime_smoking_whr_curated.txt"
smoking_initiation_whradjbmi <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/output/1_curated_data/smoking_initiation_whradjbmi_curated.txt"
lifetime_smoking_whradjbmi <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/output/1_curated_data/lifetime_smoking_whradjbmi_curated.txt"

cortisol <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/cortisol.txt"
testosterone <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/testosterone.txt"
shbg <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/shbg.txt"

# CALL
### 3. HORMONES
computing_grs(smoking_initiation_whr, 
              cortisol, testosterone, shbg, 
              "Cortisol", "Testosterone", "SHBG",
              "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_smoking_initiation_whr_grs_res_unweighted.txt",
              "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_smoking_initiation_whr_grs_res_weighted.txt")

computing_grs(smoking_initiation_whradjbmi, 
              cortisol, testosterone, shbg, 
              "Cortisol", "Testosterone", "SHBG",
              "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_smoking_initiation_whradjbmi_grs_res_unweighted.txt",
              "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_smoking_initiation_whradjbmi_grs_res_weighted.txt")

computing_grs(lifetime_smoking_whr, 
              cortisol, testosterone, shbg, 
              "Cortisol", "Testosterone", "SHBG",
              "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_lifetime_whr_grs_res_unweighted.txt",
              "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_lifetime_whr_grs_res_weighted.txt")

computing_grs(lifetime_smoking_whradjbmi, 
              cortisol, testosterone, shbg, 
              "Cortisol", "Testosterone", "SHBG",
              "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_lifetime_whradjbmi_grs_res_unweighted.txt",
              "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_lifetime_whradjbmi_grs_res_weighted.txt")

