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

#We are going to load the data for all the traits.
#the SNP that we found that the causal effect is significant: rs1317286.

cigday <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/CigDayCurrent/Cig_Day_Current_Curated_FULL.txt")
whradjbmi <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMIsmk/WHRAdjBMI_Smk_Strat_Curated.txt")
shbg <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Colocalization/Curating_data/Curated_data/shbg_curated.txt")
cortisol <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Colocalization/Curating_data/Curated_data/cortisol_curated.txt")

########################################################
#STEP 1: let's define the loci that we are going to use#
########################################################

#First let's get the position and the info of the SNP of interest:

cigday[which(cigday$SNP == "rs1317286"),]

#chr.exposure pos.exposure       SNP other_allele.exposure effect_allele.exposure dot PASS          AF        Headers                                    Headers_info beta.exposure se.exposure    logp eaf.exposure      rsid logp_neg pval.exposure
#1:           15     78896129 rs1317286                     A                      G   . PASS AF=0.343245 ES:SE:LP:AF:ID 0.0766873:0.00573869:40.0044:0.343245:rs1317286     0.0766873  0.00573869 40.0044     0.343245 rs1317286 -40.0044  9.899198e-41

#Let's get the position:

chr <- cigday[which(cigday$SNP == "rs1317286"),]$chr.exposure
pos <- cigday[which(cigday$SNP == "rs1317286"),]$pos.exposure

#And define the borders:

pos_start <- as.numeric(pos)-500000
pos_end <- as.numeric(pos)+500000

#Now let's get all the fellas that are inside this boundary:

cigday_loci <- cigday[which(cigday$chr.exposure == chr & cigday$pos.exposure >= pos_start & cigday$pos.exposure <= pos_end),] #2803 SNPs

#Let's check that all is good:

summary(cigday_loci$chr.exposure) #perfect
summary(cigday_loci$pos.exposure) #perfect

#####################################################
#STEP 2: let's align everything to the positive beta#
#####################################################

#It is easier for us to understand what will be going on. 

cigday_loci$final_beta <- ifelse(as.numeric(cigday_loci$beta.exposure) < 0, as.numeric(cigday_loci$beta.exposure)*(-1), cigday_loci$beta.exposure)
cigday_loci$final_effect_allele <- ifelse(as.numeric(cigday_loci$beta.exposure) < 0, cigday_loci$other_allele.exposure , cigday_loci$effect_allele.exposure)
cigday_loci$final_other_allele <- ifelse(as.numeric(cigday_loci$beta.exposure) < 0, cigday_loci$effect_allele.exposure, cigday_loci$other_allele.exposure)
cigday_loci$final_eaf <- ifelse(as.numeric(cigday_loci$beta.exposure) < 0, 1-cigday_loci$eaf.exposure, cigday_loci$eaf.exposure)

#Let's check if we did everything correctly...

head(cigday_loci) #seems that all is OK!!

cigday_loci_aligned <- cigday_loci %>%
  dplyr::select(chr.exposure, pos.exposure, SNP, final_effect_allele, final_other_allele, final_eaf, final_beta, se.exposure, pval.exposure)

colnames(cigday_loci_aligned) <- c("chr.exposure", "pos.exposure", "rsid", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure")

cigday_loci_aligned$samplesize.exposure <- 33229

cigday_loci_aligned$id.exposure <- "cigday"
cigday_loci_aligned$exposure <- "cigday"

cigday_loci_aligned$SNP <- paste("chr", cigday_loci_aligned$chr.exposure, ":", cigday_loci_aligned$pos.exposure, sep = "")

################################################################
#STEP 2: let's align the data from the other summary statistics#
################################################################

#Before going through the function that will match the data we are going to change colnames to the same in each dataframe.

#First for whradjbmi:

whradjbmi_clean <- whradjbmi %>%
  dplyr::select(chromosome, rs_id, Effect_allele, Other_allele, EAF_HapMapCEU, Effect_SMK, StdErr_SMK, P_value_SMK, N_SMK, chr_pos_37)

colnames(whradjbmi_clean) <- c("chr.outcome", "rsid", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "SNP")

whradjbmi_clean$outcome <- "whradjbmi"
whradjbmi_clean$id.outcome <- "whradjbmi"

#Now for shbg: it is done!! We only need to put the ids

shbg_clean <- shbg

shbg_clean$outcome <- "shbg"
shbg_clean$id.outcome <- "shbg"

colnames(shbg_clean) <- c("rsid", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "SNP", "outcome", "id.outcome")

#Now for cortisol: it is done!! Yay.

cortisol_clean <- cortisol

cortisol_clean$outcome <- "cortisol"
cortisol_clean$id.outcome <- "cortisol"

colnames(cortisol_clean) <- c("rsid", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "SNP", "outcome", "id.outcome")

#And now we can use the function:

whradjbmi_loci_aligned <- loci_aligner(cigday_loci_aligned, whradjbmi_clean) #762
shbg_loci_aligned <- loci_aligner(cigday_loci_aligned, shbg_clean) #2401
cortisol_loci_aligned <- loci_aligner(cigday_loci_aligned, cortisol_clean) #2709

#Now we need to take the loci that match in all of them:
#Let's clean the data for cigday first...

cigday_loci_aligned <- cigday_loci_aligned %>%
  select("rsid", "effect_allele.exposure", "other_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "SNP")

colnames(cigday_loci_aligned) <-  c("variant","effect_allele", "other_allele", "beta", "standard_error", "p_value", "sample_size", "chr_pos")

common_snps <- Reduce(intersect, list(cigday_loci_aligned$chr_pos, whradjbmi_loci_aligned$chr_pos ,shbg_loci_aligned$chr_pos, cortisol_loci_aligned$chr_pos)) #637

#And finally we do the match:

cigday_loci_aligned_match <- cigday_loci_aligned[which(cigday_loci_aligned$chr_pos%in%common_snps),] #we have a duplicate... damn: 639
whradjbmi_loci_aligned_match <- whradjbmi_loci_aligned[which(whradjbmi_loci_aligned$chr_pos%in%common_snps),] #637
shbg_loci_aligned_match <- shbg_loci_aligned[which(shbg_loci_aligned$chr_pos%in%common_snps),] #we have a duplicate... damn: 637
cortisol_loci_aligned_match <- cortisol_loci_aligned[which(cortisol_loci_aligned$chr_pos%in%common_snps),] #we have a duplicate... damn: 637

#Let's try to solve the issue with the duplicate:

dupl <- cigday_loci_aligned_match$chr_pos[which(duplicated(cigday_loci_aligned_match$chr_pos) == TRUE)]

#And now let's compare...

cigday_loci_aligned_match[which(cigday_loci_aligned_match$chr_pos%in%dupl),] #Two triallelic SNPs..., let's check this out.

whradjbmi_loci_aligned_match[which(whradjbmi_loci_aligned_match$chr_pos%in%dupl),]

#Gotcha:

cigday_loci_aligned_match <- cigday_loci_aligned_match[-which(cigday_loci_aligned_match$chr_pos == "chr15:79221662" & cigday_loci_aligned_match$effect_allele == "A" & cigday_loci_aligned_match$other_allele == "G"),]
cigday_loci_aligned_match <- cigday_loci_aligned_match[-which(cigday_loci_aligned_match$chr_pos == "chr15:79309560" & cigday_loci_aligned_match$effect_allele == "G" & cigday_loci_aligned_match$other_allele == "T"),]

#Let's do a small check:

cigday_loci_aligned_match[which(cigday_loci_aligned_match$chr_pos%in%dupl),] #PERFECT
whradjbmi_loci_aligned_match[which(whradjbmi_loci_aligned_match$chr_pos%in%dupl),] #PERFECT

###############################################################################################
#STEP 3: Now we are going to make a dataframe with all the info for the co-localization needed#
###############################################################################################

final_df <- cigday_loci_aligned_match %>%
  dplyr::select( "variant", "chr_pos", "effect_allele",  "other_allele", "beta", "standard_error", "p_value", "sample_size")

colnames(final_df) <- c("variant", "chr_pos", "effect_allele", "other_allele", "BETA.cigday", "SE.cigday", "P.cigday", "N.cigday") 

#And now let's add the information of the other traits:

final_df$BETA.whradjbmi <- whradjbmi_loci_aligned_match$beta
final_df$SE.whradjbmi <- whradjbmi_loci_aligned_match$standard_error
final_df$P.whradjbmi <- whradjbmi_loci_aligned_match$p_value
final_df$N.whradjbmi <- whradjbmi_loci_aligned_match$sample_size

final_df$BETA.shbg <- shbg_loci_aligned_match$beta
final_df$SE.shbg <- shbg_loci_aligned_match$standard_error
final_df$P.shbg <- shbg_loci_aligned_match$p_value
final_df$N.shbg <- shbg_loci_aligned_match$sample_size

final_df$BETA.cortisol <- cortisol_loci_aligned_match$beta
final_df$SE.cortisol <- cortisol_loci_aligned_match$standard_error
final_df$P.cortisol <- cortisol_loci_aligned_match$p_value
final_df$N.cortisol <- cortisol_loci_aligned_match$sample_size

######################
#STEP 4: SAVING DATA:#
######################

fwrite(final_df, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Colocalization/Colocalization/info_4_loci/1_ss_of_all_traits/15_78396129_79396129.txt")
