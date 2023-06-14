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

lifetime_gw <- lifetime_smk[which(lifetime_smk$P < 0.00000005),]

lifetime_gw$rsid <- lifetime_gw$SNP
lifetime_gw$pval <- lifetime_gw$P
lifetime_gw_ind <- ieugwasr::ld_clump_local(lifetime_gw, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#And now we get the SNPs that match for the TDI one:

tdi_out_dat <- extract_outcome_data(
  snps = lifetime_gw_ind$SNP,
  outcomes = 'ukb-b-10011',
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

dat_1_smk2tdi <- harmonise_data(lifetime_gw_ind, tdi_out_dat, action = 3)

#Nice, now let's do it the other way around:

dat_1_smk2tdi <- dat_1_smk2tdi[which(dat_1_smk2tdi$palindromic == FALSE),]

##############################################
#Now we are going to do the same, but for TDI#
##############################################

#We are gonna extract the data (no udåplicates)

tdi_dat <- extract_instruments(outcomes = "ukb-b-10011", p1 = 0.00000005, clump = FALSE) #369
tdi_dat_2 <- tdi_dat[which(tdi_dat$pval.exposure < 0.00000005),]
tdi_dat_2$rsid <- tdi_dat_2$SNP
tdi_dat_2$pval <- tdi_dat_2$pval.exposure
tdi_gw_ind <- ieugwasr::ld_clump_local(tdi_dat_2, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#We get 18 SNPs. 
#How many of these are in lifetime?

Lifetime_out_dat <- lifetime_smk[which(lifetime_smk$SNP%in%tdi_gw_ind$SNP),] #18/18

#Perfect match.

colnames(Lifetime_out_dat) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome",
                               "other_allele.outcome", "eaf.outcome", "info.outcome", "beta.outcome",
                               "se.outcome", "pval.outcome", "chr_pos")

Lifetime_out_dat$id.outcome <- "lifetime_smk"
Lifetime_out_dat$outcome <- "lifetime_smk"

dat_1_tdi2smk <- harmonise_data(tdi_gw_ind, Lifetime_out_dat, action = 3)

#Nice, now let's do it the other way around:

dat_1_tdi2smk <- dat_1_tdi2smk[which(dat_1_tdi2smk$palindromic == FALSE),]

###########################################
#Now we are gonna make the final dataframe#
###########################################

dat_1_smk2tdi$samplesize.exposure <- 462690

dat_1_smk2tdi_clean <- dat_1_smk2tdi %>%
  select(chr_pos, chr.exposure, pos.exposure, SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure, samplesize.exposure, eaf.outcome, beta.outcome, se.outcome, pval.outcome, samplesize.outcome)

colnames(dat_1_smk2tdi_clean) <- c("chr_pos", "chr", "pos", "SNP", "effect_allele", "other_allele", "eaf.smk", "beta.smk", "se.smk", "pval.smk", "samplesize.smk", "eaf.tdi", "beta.tdi", "se.tdi", "pval.tdi", "samplesize.tdi")

dat_1_tdi2smk$samplesize.outcome <- 462690

dat_1_tdi2smk_clean <- dat_1_tdi2smk %>%
  select(chr_pos, chr.outcome, pos.outcome, SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure, samplesize.exposure, eaf.outcome, beta.outcome, se.outcome, pval.outcome, samplesize.outcome)

colnames(dat_1_tdi2smk_clean) <- c("chr_pos", "chr", "pos", "SNP", "effect_allele", "other_allele", "eaf.tdi", "beta.tdi", "se.tdi", "pval.tdi", "samplesize.tdi", "eaf.smk", "beta.smk", "se.smk", "pval.smk", "samplesize.smk")

smk_tdi_final_df <- rbind(dat_1_smk2tdi_clean, dat_1_tdi2smk_clean)

#First let's make sure there are no duplicates:

smk_tdi_final_df[which(duplicated(smk_tdi_final_df$SNP) == TRUE),] #NONE.

######################################################################
#Now we just need to get the independent ones across all the new SNPs#
######################################################################

pval_min <- c()

for(i in seq(1, length(smk_tdi_final_df$SNP))){
  
  pval_vect <- c(smk_tdi_final_df$pval.smk[i], smk_tdi_final_df$pval.tdi[i])
  
  index <- which.min(pval_vect)
  
  final_pval <- pval_vect[index]
  
  pval_min <- c(pval_min, final_pval)
  
}

smk_tdi_final_df$pval <- pval_min

smk_tdi_independent_end <- smk_tdi_final_df %>%
  select(SNP, pval)

colnames(smk_tdi_independent_end) <- c("rsid", "pval")

#We won't care about pval now. There are mixed so... we ignore it.

independent_end <- ieugwasr::ld_clump_local(smk_tdi_independent_end, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.05) #no need to do anything. LD r2 < 0.001 and kb = 10000.

#There were 30 in high LD:

smk_tdi_final_df <- smk_tdi_final_df[which(smk_tdi_final_df$SNP%in%independent_end$rsid),]

fwrite(smk_tdi_final_df, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/Curated_data/LifeTime_TDI_Risk_Factors.txt")
