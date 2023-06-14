
library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(dplyr)

source("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/code/0_functions/functions_4_mediation.R")

smoking_aligner <- function(query_ss, other_ss){
  
  # query_ss = smoking_initiation
  # other_ss = oestradiol

  query_ss$order_col <- paste(query_ss$chromosome, query_ss$base_pair_location, sep = ":")
  other_ss$order_col <- paste(other_ss$chromosome, other_ss$base_pair_location, sep = ":")
  
  # sort by numeric_part
  query_ss <- query_ss[order(as.numeric(gsub(":.*", "", query_ss$order_col))), ]
  other_ss <- other_ss[order(as.numeric(gsub(":.*", "", other_ss$order_col))), ]
  
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
  
  # missing_rows <- anti_join(exposure, merged_df, by = "SNP")
  # print(dim(missing_rows))
  # fwrite(missing_rows, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/missing_hormones_initiation_whradjbmi_grs_res_weighted.txt")
  
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
                          cortisol, oestradiol, testosterone, shbg,
                          name1, name2, name3, name4,
                          unweighted_out, weighted_out){
  
  smoking_initiation <- fread(smoking_data)
  
  ######################
  #Loading outcome data#
  ######################
  
  cortisol <- fread(cortisol)
  oestradiol <- fread(oestradiol)
  testosterone <- fread(testosterone)
  shbg <- fread(shbg)
  
  ##############################################
  #STEP 1: matching data for smoking initiation#
  ##############################################
  
  #First volume:
  
  cortisol_2_initiation <- smoking_aligner(smoking_initiation, cortisol) #76
  oestradiol_2_initiation <- smoking_aligner(smoking_initiation, oestradiol) #76
  testosterone_2_initiation <- smoking_aligner(smoking_initiation, testosterone) #76
  shbg_2_initiation <- smoking_aligner(smoking_initiation, shbg) #76
  
  # oestradiol_check <- oestradiol[which(oestradiol$variant %in% oestradiol_2_initiation$variant)]
  # smk_check <- smoking_initiation[which(smoking_initiation$variant %in% oestradiol_2_initiation$variant)]
  
  #Let's match the smoking initiation variants and order all of our variants to that:
  
  # smoking_initiation_match <- smoking_initiation[which(smoking_initiation$variant%in%cortisol_2_initiation$variant),] #76: perfect match
  
  #Let's first find the SNPs found in all:
  
  common_snps <- Reduce(intersect, list(smoking_initiation$chr_pos,
                                        cortisol_2_initiation$chr_pos, oestradiol_2_initiation$chr_pos,
                                        testosterone_2_initiation$chr_pos, shbg_2_initiation$chr_pos)) # 110
  
  #We match aligned data to the SNPs that are common for ALL traits:
  smoking_initiation_common <- smoking_initiation[which(smoking_initiation$chr_pos %in% common_snps),]
  cortisol_common <- cortisol_2_initiation[which(cortisol_2_initiation$chr_pos %in% common_snps),]
  oestradiol_common <- oestradiol_2_initiation[which(oestradiol_2_initiation$chr_pos %in% common_snps),]
  testosterone_common <- testosterone_2_initiation[which(testosterone_2_initiation$chr_pos %in% common_snps),]
  shbg_common <- shbg_2_initiation[which(shbg_2_initiation$chr_pos %in% common_snps),]
  
  #Here we order the 110:
  
  smoking_initiation_match_ordered <- smoking_initiation_common[order(as.numeric(smoking_initiation_common$chromosome), as.numeric(smoking_initiation_common$base_pair_location)),]
  cortisol_2_initiation_ordered <- cortisol_common[order(as.numeric(cortisol_common$chromosome), as.numeric(cortisol_common$base_pair_location)),]
  oestradiol_2_initiation_ordered <- oestradiol_common[order(as.numeric(oestradiol_common$chromosome), as.numeric(oestradiol_common$base_pair_location)),]
  testosterone_2_initiation_ordered <- testosterone_common[order(as.numeric(testosterone_common$chromosome), as.numeric(testosterone_common$base_pair_location)),]
  shbg_2_initiation_ordered <- shbg_common[order(as.numeric(shbg_common$chromosome), as.numeric(shbg_common$base_pair_location)),]
  
  #Let's do some checks on effect alleles:
  
  print(length(which(smoking_initiation_match_ordered$effect_allele == cortisol_2_initiation_ordered$effect_allele))) #good
  print(length(which(smoking_initiation_match_ordered$effect_allele == oestradiol_2_initiation_ordered$effect_allele))) #good
  print(length(which(smoking_initiation_match_ordered$effect_allele == testosterone_2_initiation_ordered$effect_allele))) #good
  print(length(which(smoking_initiation_match_ordered$effect_allele == shbg_2_initiation_ordered$effect_allele))) #good
  
  #######################
  #STEP 2: computing GRS#
  #######################
  
  list_of_traits <- list(cortisol_2_initiation_ordered, 
                         oestradiol_2_initiation_ordered, 
                         testosterone_2_initiation_ordered,
                         shbg_2_initiation_ordered) 
  
  trait_names <- c(name1, name2, name3, name4)
  type_of_trait <- c(name1, name2, name3, name4)
  
  smoking_initiation_smk_grs_unw <- grs_producer(query_ss = smoking_initiation_match_ordered, list_of_traits_ = list_of_traits, trait_names_ = trait_names,  type_of_trait_ = type_of_trait, weighted = FALSE) 
  
  fwrite(smoking_initiation_smk_grs_unw, unweighted_out)
  
  smoking_initiation_smk_grs_w <- grs_producer(query_ss = smoking_initiation_match_ordered, list_of_traits_ = list_of_traits, trait_names_ = trait_names,  type_of_trait_ = type_of_trait, weighted = TRUE) 
  
  fwrite(smoking_initiation_smk_grs_w, weighted_out)
}

# DATA
# name1 = "Cortisol"
# name2= "Oestradiol"
# name3="Testosterone"
# name4= "SHBG"

smoking_initiation_whradjbmi <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/output/1_curated_data/smoking_initiation_whradjbmi_curated.txt"
lifetime_smoking_whradjbmi <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/output/1_curated_data/lifetime_smoking_whradjbmi_curated.txt"

cortisol <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/cortisol.txt"
oestradiol <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/curated_estradiol_chrpos_noX.txt"
testosterone <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/testosterone.txt"
shbg <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/shbg.txt"

# CALL
### 3. HORMONES
computing_grs(smoking_initiation_whradjbmi, 
              cortisol, oestradiol, testosterone, shbg, 
              "Cortisol", "Oestradiol", "Testosterone", "SHBG",
              "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_smoking_initiation_whradjbmi_grs_res_unweighted.txt",
              "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_smoking_initiation_whradjbmi_grs_res_weighted.txt")

# unweighted_out = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_smoking_initiation_whradjbmi_grs_res_unweighted.txt"
# weighted_out = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_smoking_initiation_whradjbmi_grs_res_weighted.txt"


computing_grs(lifetime_smoking_whradjbmi, 
              cortisol, oestradiol, testosterone, shbg, 
              "Cortisol", "Oestradiol", "Testosterone", "SHBG",
              "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_lifetime_whradjbmi_grs_res_unweighted.txt",
              "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_lifetime_whradjbmi_grs_res_weighted.txt")

# unweighted_out = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_lifetime_whradjbmi_grs_res_unweighted.txt"
# weighted_out = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_lifetime_whradjbmi_grs_res_weighted.txt"
# 


# rs10000012             C            G                  0.8595  5e-04 -0.0004 
### http://www.phenoscanner.medschl.cam.ac.uk/?query=rs10000012&catalogue=GWAS&p=1e-5&proxies=None&r2=0.8&build=37 

# > df_ordered 
# variant effect_allele other_allele effect_allele_frequency freqse    beta standard_error   p_value direction sample_size chromosome base_pair_location        chr_pos    ---> SEEMS THAT A1 & A2 ARE SWAPPED
# 1:  rs727479             A            C                  0.6491 0.0007  0.1250         0.0094 9.654e-41        ++      311675         15           51534547 chr15:51534547 --> ESTRADIOL (BETA+ ; C (A1) A (A2)) http://www.phenoscanner.medschl.cam.ac.uk/?query=rs727479&catalogue=GWAS&p=1e-5&proxies=None&r2=0.8&build=37
# 2: rs2414095             A            G                  0.3500 0.0007 -0.1246         0.0093 1.304e-40        --      311675         15           51524292 chr15:51524292 --> A G HORMONES (BETA +)
# 3: rs7173595             T            C                  0.6496 0.0007  0.1247         0.0094 1.348e-40        ++      311675         15           51533736 chr15:51533736 --> C T 
# 4: rs7175531             T            C                  0.3503 0.0007 -0.1247         0.0094 1.364e-40        --      311675         15           51534055 chr15:51534055 --> C T 
# 5: rs2414097             A            G                  0.6495 0.0007  0.1245         0.0093 1.593e-40        ++      311675         15           51529835 chr15:51529835 --> A G 
# ---                                                                                                                                                                        
#   7564113: rs9804675             T            C                  0.4048 0.0005  0.0000         0.0090 1.000e+00        +-      311675         11           62819470 chr11:62819470
# 7564114: rs9841832             A            C                  0.2184 0.0015  0.0000         0.0107 1.000e+00        -+      311675          3          184186015 chr3:184186015
# 7564115: rs9843781             T            C                  0.9033 0.0006  0.0000         0.0148 1.000e+00        +-      311675          3            9396766   chr3:9396766
# 7564116: rs9870902             T            C                  0.9850 0.0001  0.0000         0.0367 1.000e+00        -+      311675          3           20994822  chr3:20994822
# 7564117: rs9871963             A            G                  0.6298 0.0004  0.0000         0.0092 1.000e+00        +-      311675          3          141033481 chr3:141033481
