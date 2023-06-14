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

#For Life-time Smoker:

lifetime_smk <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/Lifetime_smk/Lifetime_Smoking_Curated.txt")

lifetime_gw <- lifetime_smk[which(as.numeric(lifetime_smk$P) < 0.00000005),]

lifetime_gw$rsid <- lifetime_gw$SNP
lifetime_gw$pval <- lifetime_gw$P

pre_independent_set_1 <- lifetime_gw %>%
  select(rsid, pval)

pre_independent_set_1$source <- "Lifetime"

#And now get the genome-wide significant for tdi:

tdi_dat <- extract_instruments(outcomes = "ukb-b-10011", p1 = 0.00000005, clump = FALSE) #369

#In this case there are no duplicates in TDI, so there is no need for the following part of the code.
#Just in case I will leave it here, but it is deprecated.

##########################################
#We continue with the independent setting#
##########################################

pre_independent_set_2 <- tdi_dat %>%
  select(SNP, pval.exposure)

colnames(pre_independent_set_2) <- c("rsid", "pval")

pre_independent_set_2$source <- "tdi"

independent_snps <- rbind(pre_independent_set_1, pre_independent_set_2) #10766 SNPs!

########################################################################
#Careful: we might have duplicates and we didn't take that into account#
########################################################################

length(which(duplicated(independent_snps$rsid) == TRUE)) #40/10766.

#Let's see the examples:

dupl <- independent_snps[which(independent_snps$rsid%in%independent_snps[which(duplicated(independent_snps$rsid) == TRUE),]$rsid),]

dupl <- dupl[order(dupl$rsid),]

#Let's check what is going on with these SNPs.
#First let's order them:

independent_snps <- independent_snps[order(independent_snps$rsid),]

#Let's get the independent SNPs and check how it deals with the duplicates:

final_set_of_snps <- ieugwasr::ld_clump_local(independent_snps, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.00000005) #135...

#####################################
#Let's check how this is dealed with#
#####################################

dupl_in_final_set_of_snps <- final_set_of_snps[which(final_set_of_snps$rsid%in%dupl$rsid),]

#2 of them!!
#Let's check them:

dupl_in_final_set_of_snps <- dupl_in_final_set_of_snps[order(dupl_in_final_set_of_snps$rsid),]

View(dupl_in_final_set_of_snps)

#It does not remove them!!

#Easy, let's check what happens when we order them and remove the duplicates:

#######################
#Repeating the pruning#
#######################

independent_snps <- independent_snps[order(independent_snps$pval),]

independent_snps_no_dupl <- independent_snps[which(duplicated(independent_snps$rsid) == FALSE),]

#Now let's get the independent ones:

final_set_of_snps_no_dupl <- ieugwasr::ld_clump_local(independent_snps_no_dupl, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.00000005) #304...

#134. Very similar. We just removed the duplicate.

##############################################
#SO THIS IS GOING TO BE THE FINAL SET OF SNPs#
##############################################

#Now obtain the SNPs for each of the data.

##############################################################################################
#Before merging the data let's try to take the proxies by retrieving them from tdication data#
##############################################################################################

#And now we get the SNPs that match for the tdi one:

tdi_out_dat <- extract_outcome_data(
  snps = final_set_of_snps_no_dupl$rsid,
  outcomes = 'ukb-b-10011',
  proxies = TRUE
) #No proxies were obtained, btw. All are found.

#PERFECT MATCH.

################################################################
#Now obtain the SNPs for Lifetime from those found in tdication#
################################################################

lifetime_exp_dat <- lifetime_smk[which(lifetime_smk$SNP%in%tdi_out_dat$SNP),] #134!

#We are missing three SNPs, can we get the proxies???

missing_snps <- tdi_out_dat[which(!(tdi_out_dat$SNP%in%lifetime_exp_dat$SNP)),] #Let's get this shit done.

missing_snps #they are deep in the MHC region!!!

#That means that we have finished!!
#We only need to merge the SNPs.

#There is only one issue: the duplicate SNP.
#Let's check:

rsid <- tdi_out_dat$SNP[which(duplicated(tdi_out_dat$SNP) == TRUE)]

#########################################################################
#Checking that duplicate before taking into account the rest of the data#
#########################################################################

lifetime_smk[which(lifetime_smk$SNP == rsid),] #only 1 and it is the palindromic one.

tdi_out_dat[which(tdi_out_dat$SNP == rsid),] #There is gonna be only one, but one of them is palindromic so it will probably be removed, though it is genome-wide...

##################
#Merging the data#
##################

colnames(lifetime_exp_dat) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
                                                "other_allele.exposure", "eaf.exposure", "info.exposure", "beta.exposure",
                                                "se.exposure", "pval.exposure", "chr_pos")

lifetime_exp_dat$id.exposure <- "lifetime_exp_dat"
lifetime_exp_dat$exposure <- "lifetime_exp_dat"

dat_1_smk2tdi <- harmonise_data(lifetime_exp_dat, tdi_out_dat, action = 3)

#Nice, now let's do it the other way around:

dat_1_smk2tdi <- dat_1_smk2tdi[which(dat_1_smk2tdi$remove == FALSE),]
dat_1_smk2tdi <- dat_1_smk2tdi[which(dat_1_smk2tdi$palindromic == FALSE),]

#########################
#We end up with 119 SNPs#
#########################

#Let's obtain where we obtained it from:

final_set_of_snps_no_dupl_clean <- final_set_of_snps_no_dupl[which(final_set_of_snps_no_dupl$rsid%in%dat_1_smk2tdi$SNP),]
final_set_of_snps_no_dupl_clean <- final_set_of_snps_no_dupl_clean[order(match(final_set_of_snps_no_dupl_clean$rsid, dat_1_smk2tdi$SNP)),]

length(which(final_set_of_snps_no_dupl_clean$rsid != dat_1_smk2tdi$SNP)) #perfect
length(which(final_set_of_snps_no_dupl_clean$rsid == dat_1_smk2tdi$SNP)) #perfect

dat_1_smk2tdi$source <- final_set_of_snps_no_dupl_clean$source

######################################
#Let's continue and generate the data#
######################################

dat_1_smk2tdi$samplesize.exposure <- 462690

final_df <- dat_1_smk2tdi %>%
  select(SNP, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, eaf.exposure, pval.exposure, samplesize.exposure, beta.outcome, se.outcome, eaf.outcome, pval.outcome, samplesize.outcome, chr_pos, source)

colnames(final_df) <- c("SNP", "effect_allele", "other_allele", "beta.smk", "se.smk", "eaf.smk", "pval.smk", "samplesize.smk", "beta.tdi", "se.tdi", "eaf.tdi", "pval.tdi", "samplesize.tdi", "chr_pos", "source")

fwrite(final_df, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/Curated_data/LifeTime_TDI_Risk_Factors.txt")

##############################
#Final check about MHC region#
##############################

#Since lifetime smoking has been curated to AVOID MHC region, we shouldn't have any variants after the matching.

#Let's check it out:

chr_6 <- dat_1_smk2tdi[which(dat_1_smk2tdi$chr == 6),]

#Only 6 SNPs! Let's check where they are:

View(chr_6) #none of them in MHC region, as expected.
