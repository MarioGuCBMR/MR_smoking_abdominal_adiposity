##############
#INTRODUCTION#
##############

#This code is to get the SNPs for MR-BMA for Smoking Initiation and BMI.

#The number of risk factors will increase for sure.
#And things will change, this is just a test.

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)
library(TwoSampleMR)
library(ieugwasr)

##############
#LOADING DATA#
##############

#For the risk factors we are only going to include SNPs that are in
#all risk factors, taking into account proxies, of course. 
#Things get a bit cumbersome since we need them all.

#For Life-time Smoker

lifetime_smk <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/Lifetime_smk/Lifetime_Smoking_Curated.txt")

lifetime_gw <- lifetime_smk[which(as.numeric(lifetime_smk$P) < 0.00000005),]

lifetime_gw$rsid <- lifetime_gw$SNP
lifetime_gw$pval <- lifetime_gw$P
lifetime_gw_ind <- ieugwasr::ld_clump_local(lifetime_gw, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#And now we get the SNPs that match for the edu one:

edu_out_dat <- extract_outcome_data(
  snps = lifetime_gw_ind$SNP,
  outcomes = 'ukb-b-16489',
  proxies = TRUE
)

#This worked perfectly.
#Let's get the harmonized version:

#We need to change the data so it matches perfectly with the 2SMR version:

colnames(lifetime_gw_ind) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
                                                "other_allele.exposure", "eaf.exposure", "info.exposure", "beta.exposure",
                                                "se.exposure", "pval.exposure", "chr_pos", "rsid", "pval")

lifetime_gw_ind$id.exposure <- "lifetime_smk"
lifetime_gw_ind$exposure <- "lifetime_smk"

dat_1_smk2edu <- harmonise_data(lifetime_gw_ind, edu_out_dat, action = 3)

#Nice, now let's do it the other way around:

dat_1_smk2edu <- dat_1_smk2edu[which(dat_1_smk2edu$palindromic == FALSE),]

##############################################
#Now we are going to do the same, but for edu#
##############################################

#We are gonna extract the data (no udåplicates)

edu_dat <- extract_instruments(outcomes = "ukb-b-16489", p1 = 0.00000005, clump = FALSE) #369

#There is one duplicate: rs2696531.
#Let's see what can we do with this little bastard.

edu_dat[which(edu_dat$SNP == "rs2696531"),] 
lifetime_smk[which(lifetime_smk$SNP == "rs2696531"),] #it is also duplicated here!

#We will probably have to repeat this then..., NOT. One of them is palindromic: C, G. It removes the one we do not like.

edu_dat_2 <- edu_dat[which(edu_dat$pval.exposure < 0.00000005),]
edu_dat_2$rsid <- edu_dat_2$SNP
edu_dat_2$pval <- edu_dat_2$pval.exposure
edu_gw_ind <- ieugwasr::ld_clump_local(edu_dat_2, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#We get 260 SNPs. 
#How many of these are in lifetime?

Lifetime_out_dat <- lifetime_smk[which(lifetime_smk$SNP%in%edu_gw_ind$SNP),] #257/260. Almost!

#We might actually need to do some proxy analysis for this one first.
#But just for three SNPs:

missing_snps <- edu_gw_ind[which(!(edu_gw_ind$SNP%in%Lifetime_out_dat$SNP)),] #Let's get this shit done.

#THEY ARE ALL DEEP IN THE MHC REGION. 
#HENCE WHY THEY WERE NOT OBTAINED.

#That is why the curation beforehand worked perfectly.

#Perfect match.

colnames(Lifetime_out_dat) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome",
                               "other_allele.outcome", "eaf.outcome", "info.outcome", "beta.outcome",
                               "se.outcome", "pval.outcome", "chr_pos")

Lifetime_out_dat$id.outcome <- "lifetime_smk"
Lifetime_out_dat$outcome <- "lifetime_smk"

dat_1_edu2smk <- harmonise_data(edu_gw_ind, Lifetime_out_dat, action = 3)

#Nice, now let's do it the other way around:

dat_1_edu2smk <- dat_1_edu2smk[which(dat_1_edu2smk$palindromic == FALSE),]

###########################################
#Now we are gonna make the final dataframe#
###########################################

dat_1_smk2edu$samplesize.exposure <- 462690

dat_1_smk2edu_clean <- dat_1_smk2edu %>%
  select(chr_pos, chr.exposure, pos.exposure, SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure, samplesize.exposure, eaf.outcome, beta.outcome, se.outcome, pval.outcome, samplesize.outcome)

colnames(dat_1_smk2edu_clean) <- c("chr_pos", "chr", "pos", "SNP", "effect_allele", "other_allele", "eaf.smk", "beta.smk", "se.smk", "pval.smk", "samplesize.smk", "eaf.edu", "beta.edu", "se.edu", "pval.edu", "samplesize.edu")

dat_1_edu2smk$samplesize.outcome <- 462690

dat_1_edu2smk_clean <- dat_1_edu2smk %>%
  select(chr_pos, chr.outcome, pos.outcome, SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure, samplesize.exposure, eaf.outcome, beta.outcome, se.outcome, pval.outcome, samplesize.outcome)

colnames(dat_1_edu2smk_clean) <- c("chr_pos", "chr", "pos", "SNP", "effect_allele", "other_allele", "eaf.edu", "beta.edu", "se.edu", "pval.edu", "samplesize.edu", "eaf.smk", "beta.smk", "se.smk", "pval.smk", "samplesize.smk")

smk_edu_final_df <- rbind(dat_1_smk2edu_clean, dat_1_edu2smk_clean)

#First let's make sure there are no duplicates:

smk_edu_final_df[which(duplicated(smk_edu_final_df$SNP) == TRUE),] #We do have a duplicate.

#Let's investigate it:

smk_edu_final_df[which(smk_edu_final_df$SNP == "rs1050847"),]

#They are exactly the same!
#Where the hell do the come from as duplicated?

edu_dat[which(edu_dat$SNP == "rs1050847"),] #only one here...

lifetime_smk[which(lifetime_smk$SNP == "rs1050847"),] #only one.

#AHÁ. It is because they are genome-wide significant both of them.
#Oh welp, no worries. #We can remove the second one that's it.

smk_edu_final_df <- smk_edu_final_df[which(duplicated(smk_edu_final_df$SNP) == FALSE),] 

######################################################################
#Now we just need to get the independent ones across all the new SNPs#
######################################################################

pval_min <- c()

for(i in seq(1, length(smk_edu_final_df$SNP))){
  
  pval_vect <- c(smk_edu_final_df$pval.smk[i], smk_edu_final_df$pval.edu[i])
  
  index <- which.min(pval_vect)
  
  final_pval <- pval_vect[index]
  
  pval_min <- c(pval_min, final_pval)
  
}

smk_edu_final_df$pval <- pval_min

smk_edu_independent_end <- smk_edu_final_df %>%
  select(SNP, pval)

colnames(smk_edu_independent_end) <- c("rsid", "pval")

#We won't care about pval now. There are mixed so... we ignore it.

independent_end <- ieugwasr::ld_clump_local(smk_edu_independent_end, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.05) #no need to do anything. LD r2 < 0.001 and kb = 10000.

#We end up with 241 SNPs.

smk_edu_final_df <- smk_edu_final_df[which(smk_edu_final_df$SNP%in%independent_end$rsid),]

fwrite(smk_edu_final_df, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/Curated_data/LifeTime_edu_Risk_Factors.txt")
