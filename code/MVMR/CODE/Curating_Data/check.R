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

#And now get the genome-wide significant for EDU:

edu_dat <- extract_instruments(outcomes = "ukb-b-16489", p1 = 0.00000005, clump = FALSE) #369

#There is only one duplicate. 
#Let's check on lifetime just in case:

lifetime_smk[which(lifetime_smk$SNP == "rs2696531"),] #the duplicate is there!

#It doesn't matter since the LD reference panel does not have triallelic SNPs after all.

#We will do the following:

dupl <- extract_outcome_data(
  snps = "rs2696531",
  outcomes = 'ukb-b-16489',
  proxies = FALSE
)

#Now we got it! Let's remove it from the original:

edu_dat <- edu_dat[-which(edu_dat$SNP == "rs2696531"),]

#And now we merge them, though we have to tweek a bit the data:

colnames(edu_dat)

colnames(dupl) <-  c("SNP", "chr.exposure", "pos.exposure",                 
"beta.exposure",          "se.exposure",            "samplesize.exposure",   
"pval.exposure",          "eaf.exposure",           "effect_allele.exposure",
"other_allele.exposure",  "exposure",               "id.exposure",           
"originalname.exposure",  "exposure.deprecated",    "mr_keep.exposure",      
"data_source.exposure")  

dupl$pval_origin.exposure <- edu_dat$pval_origin.exposure[1]
dupl$id.exposure <- edu_dat$id.exposure[1]

dupl <- dupl %>%
  select(colnames(edu_dat))

edu_dat <- rbind(edu_dat, dupl) #PERFECT.

#We continue with the independent setting:

pre_independent_set_2 <- edu_dat %>%
  select(SNP, pval.exposure)

colnames(pre_independent_set_2) <- c("rsid", "pval")

independent_snps <- rbind(pre_independent_set_1, pre_independent_set_2) #36324 SNPs!

#################################################################################################################
#This a check code to see if we can obtain the same amount of SNPs with the pval_min protocol and Tom's approach#
#################################################################################################################

#To do so, what we need to do is obtain the full dataframe above and obtain the SNPs for both sets of data.

lifetime_check <- lifetime_smk[which(lifetime_smk$SNP%in%independent_snps$rsid),] #34554! Almost all of them.

edu_check <- extract_outcome_data(
  snps = independent_snps$rsid,
  outcomes = 'ukb-b-16489',
  proxies = TRUE
)

#Not all of them are obtained for each dataframe, of course.
#It seems that the bottleneck is actually lifetime smoking.
#Let's try again lifetime_check, but with the proxies above.

lifetime_check <- lifetime_smk[which(lifetime_smk$SNP%in%edu_check$SNP),] #34553! Just one SNP of difference.

###############################
#Let's merge the datasets then#
###############################

colnames(lifetime_check) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
                                "other_allele.exposure", "eaf.exposure", "info.exposure", "beta.exposure",
                                "se.exposure", "pval.exposure", "chr_pos")

lifetime_check$id.exposure <- "lifetime_exp_dat"
lifetime_check$exposure <- "lifetime_exp_dat"

dat_1_smk2edu <- harmonise_data(lifetime_check, edu_check, action = 3)

#And now let's get the independent SNPs...

pval_min <- c()

for(i in seq(1, length(dat_1_smk2edu$SNP))){
  
  pval_vect <- c(dat_1_smk2edu$pval.exposure[i], dat_1_smk2edu$pval.outcome[i])
  
  index <- which.min(pval_vect)
  
  final_pval <- pval_vect[index]
  
  pval_min <- c(pval_min, final_pval)
  
}

independent_snps <- as.data.frame(dat_1_smk2edu$SNP)

independent_snps$pval <- pval_min

colnames(independent_snps) <- c("rsid", "pval")

#Let's get the independent SNPs:

final_set_of_snps <- ieugwasr::ld_clump_local(independent_snps, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.00000005) #300 SNPs.

##############################################
#SO THIS IS GOING TO BE THE FINAL SET OF SNPs#
##############################################

#Now obtain the SNPs for each of the data.

##############################################################################################
#Before merging the data let's try to take the proxies by retrieving them from education data#
##############################################################################################

#And now we get the SNPs that match for the edu one:

edu_out_dat <- extract_outcome_data(
  snps = final_set_of_snps$rsid,
  outcomes = 'ukb-b-16489',
  proxies = TRUE
)

#Even with proxies we only obtained 298 SNPs (a duplicate is retained!).
#Since Education might be the bottleneck (we are already trying proxies here)
#We are gonna, then base the matching with this.

################################################################
#Now obtain the SNPs for Lifetime from those found in education#
################################################################

lifetime_exp_dat <- lifetime_smk[which(lifetime_smk$SNP%in%edu_out_dat$SNP),] #296!

#We are missing three SNPs, can we get the proxies???

missing_snps <- edu_out_dat[which(!(edu_out_dat$SNP%in%lifetime_exp_dat$SNP)),] #Let's get this shit done.

missing_snps 

rsid <- edu_out_dat$SNP[which(duplicated(edu_out_dat$SNP) == TRUE)]

#########################################################################
#Checking that duplicate before taking into account the rest of the data#
#########################################################################

lifetime_smk[which(lifetime_smk$SNP == rsid),] #only 1 and it is the palindromic one.

edu_out_dat[which(edu_out_dat$SNP == rsid),] #There is gonna be only one, but one of them is palindromic so it will probably be removed, though it is genome-wide...

##################
#Merging the data#
##################

colnames(lifetime_exp_dat) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
                                                "other_allele.exposure", "eaf.exposure", "info.exposure", "beta.exposure",
                                                "se.exposure", "pval.exposure", "chr_pos")

lifetime_exp_dat$id.exposure <- "lifetime_exp_dat"
lifetime_exp_dat$exposure <- "lifetime_exp_dat"

dat_1_smk2edu <- harmonise_data(lifetime_exp_dat, edu_out_dat, action = 3)

#Nice, now let's do it the other way around:

dat_1_smk2edu <- dat_1_smk2edu[which(dat_1_smk2edu$remove == FALSE),]
dat_1_smk2edu <- dat_1_smk2edu[which(dat_1_smk2edu$palindromic == FALSE),]

#########################
#We end up with 251 SNPs#
#########################

dat_1_smk2edu$samplesize.exposure <- 462690

final_df <- dat_1_smk2edu %>%
  select(SNP, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, eaf.exposure, pval.exposure, samplesize.exposure, beta.outcome, se.outcome, eaf.outcome, pval.outcome, samplesize.outcome, chr_pos)

colnames(final_df) <- c("SNP", "effect_allele", "other_allele", "beta.smk", "se.smk", "eaf.smk", "pval.smk", "samplesize.smk", "beta.edu", "se.edu", "eaf.edu", "pval.edu", "samplesize.edu", "chr_pos")

############################################################
#Mysteriously we end up with almost the same amount of SNPs#
############################################################

#... I need to think about this.
