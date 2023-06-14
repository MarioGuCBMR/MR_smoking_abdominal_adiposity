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

#######################
#GETTING DATA FOR MVMR#
#######################

###########################################################################
#A. We are going to obtain genome-wide significant SNPs for the two traits#
###########################################################################

#############################
##A.1 For Life-time Smoking##
#############################

#Loading data:

lifetime_smk <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/Lifetime_smk/Lifetime_Smoking_Curated.txt")

#Filtering GW:

lifetime_gw <- lifetime_smk[which(as.numeric(lifetime_smk$P) < 0.00000005),]

#Preparing data to obtain, later the genome-wide significant SNPs.

lifetime_gw$rsid <- lifetime_gw$SNP
lifetime_gw$pval <- lifetime_gw$P

pre_independent_set_1 <- lifetime_gw %>%
  select(rsid, pval)

pre_independent_set_1$source <- "Lifetime"

#####################
##A.2 For Education##
#####################

#And now get the genome-wide significant for EDU:

edu_dat <- extract_instruments(outcomes = "ukb-b-16489", p1 = 0.00000005, clump = FALSE) #25926

#There is only one duplicate, but they are keeping only the first instance. 
#How do we deal with this? 
#We should check both of them and see if it is important.

#1. Is the duplicate genome-wide significant or is it just coincidence?

summary(edu_dat$pval.exposure) #it is genome-wide significant.

#2. Can we remove it since it is not in lifetime smoking?

lifetime_smk[which(lifetime_smk$SNP == "rs2696531"),] #the duplicate is there!

#SNP CHR       BP EFFECT_ALLELE OTHER_ALLELE      EAF     INFO       BETA
#1: rs2696531  17 44355634             C            A 0.787784 0.969450 0.00618641
#2: rs2696531  17 44355634             C            G 0.778039 0.965203 0.00618357
#SE       P        chr_pos
#1: 0.00171913 0.00032 chr17:44355634
#2: 0.00169575 0.00027 chr17:44355634

#3. Let's recover this fella in the original data

dupl <- extract_outcome_data(
  snps = "rs2696531",
  outcomes = 'ukb-b-16489',
  proxies = FALSE
)

dupl

#SNP chr      pos beta.outcome se.outcome samplesize.outcome pval.outcome
#1 rs2696531  17 44355634   -0.0104694 0.00118003             458079  7.19946e-19
#2 rs2696531  17 44355634   -0.0103249 0.00116394             458079  7.29962e-19
#eaf.outcome effect_allele.outcome other_allele.outcome
#1    0.212220                     A                    C
#2    0.221967                     G                    C

#Look at the effects, se and p-values. They are the same, almost. 
#I am afraid of including it, since it is gonna give us issues, for sure.
#Nonetheless, maybe it will be removed since there are other SNPs in high LD with it.
#Let's keep it for the time being.

#Now we got it! Let's do some tricks to include this SNP there.

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

#Finally, let's remove the SNPs in MHC region, just like we did curating lifetime smoking data;

edu_dat_mhc <- edu_dat[which(as.numeric(edu_dat$chr.exposure) == 6 & as.numeric(edu_dat$pos.exposure) >= 26000000 & as.numeric(edu_dat$pos.exposure) <= 34000000),]

summary(as.numeric(edu_dat_mhc$pos.exposure)) #perfect
summary(as.numeric(edu_dat_mhc$chr.exposure)) #perfect

edu_dat <- edu_dat[-which(edu_dat$SNP%in%edu_dat_mhc$SNP),]

#Let's do a quick check that this has worked...

edu_dat_check <- edu_dat[which(edu_dat$chr.exposure == 6),]

edu_dat_check <- edu_dat_check[order(as.numeric(edu_dat_check$pos.exposure)),]

View(edu_dat_check) #PERFECT

#Now let's format this new data to get the independent SNPs:

pre_independent_set_2 <- edu_dat %>%
  select(SNP, pval.exposure)

colnames(pre_independent_set_2) <- c("rsid", "pval")

pre_independent_set_2$source <- "Education"

independent_snps <- rbind(pre_independent_set_1, pre_independent_set_2) #35207 SNPs!

#Now we have the combined set of genome-wide significant SNPs, but we have to make them truly independent.

######################################################
#B investigate the duplicates and perform LD clumping#
######################################################

################################
#B.1 Let's check the duplicates#
################################

no_dupl <- unique(independent_snps$rsid) #we have 34580.

#if we order by p-value and remove them:

independent_snps <- independent_snps[order(independent_snps$pval),]

independent_snps_no_dupl <- independent_snps[-which(duplicated(independent_snps$rsid) == TRUE),] #34580

#Perfect!! We get the same amout of SNPs. That means that we do not have more special cases like the one above.

#There are two cases of duplicates:

#1. The special case in which you found them in the same exposure data.
#2. Those that appear because they are in the two exposures.

#This code above takes care of case 2, since we are taking the p-value of the exposure
#that the SNP is more associated with the SNP. Case 1 is gonna drag us down a bit.

########################################
#B.2 Now let's get the independent ones#
########################################

final_set_of_snps_no_dupl <- ieugwasr::ld_clump_local(independent_snps_no_dupl, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.00000005) #297

#Now obtain the SNPs for each of the data.

###############################################################################################
#C. Merging both datasets so we can have effect sizes of both traits for each independent SNP #
###############################################################################################

##################################################################
#C.1: now we obtain the summary data of all 297 SNPs in Education#
##################################################################

#This one is easy:

edu_out_dat <- extract_outcome_data(
  snps = final_set_of_snps_no_dupl$rsid,
  outcomes = 'ukb-b-16489',
  proxies = FALSE
) #No proxies were obtained, btw. All are found.

#We obtain 298. Is it our friend from case 1?

#####################################################
#C. 2 Now obtain the summary statistics of Lifetime #
#####################################################

lifetime_exp_dat <- lifetime_smk[which(lifetime_smk$SNP%in%edu_out_dat$SNP),] #297!

#There is a missing SNP!

missing_snps <- edu_out_dat[which(!(edu_out_dat$SNP%in%lifetime_exp_dat$SNP)),] #0

#There are no missing SNPs, that means that we are dealing with a duplicate.

dupl <- edu_out_dat$SNP[which(duplicated(edu_out_dat$SNP) == TRUE)]

lifetime_smk[which(lifetime_smk$SNP == dupl),] #it is a different one!

#SNP CHR       BP EFFECT_ALLELE OTHER_ALLELE      EAF     INFO       BETA
#1: rs2606913  16 14916988             T            A 0.177068 0.957933 0.00304362
#SE   P        chr_pos
#1: 0.00185417 0.1 chr16:14916988

#In Lifetime only this one is present.

edu_out_dat[which(edu_out_dat$SNP == dupl),]

#SNP chr      pos beta.outcome se.outcome samplesize.outcome pval.outcome
#171 rs2606913  16 14916988   0.00748826 0.00127209             458079  3.89996e-09
#292 rs2606913  16 14916988   0.00888099 0.00660753             458079  1.80000e-01
#eaf.outcome effect_allele.outcome other_allele.outcome
#171    0.822881                     A                    T
#292    0.008769                     C                    T

#And here we have the reason why. 
#Since we only have one combination, we only need to remove it from edu_dat_out:

edu_out_dat[292,] #This one out!

edu_out_dat <- edu_out_dat[-292,]

#Using data harmonization we won't have any trouble merging the data.

######################
#C.3 Merging the data#
######################

colnames(lifetime_exp_dat) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
                                                "other_allele.exposure", "eaf.exposure", "info.exposure", "beta.exposure",
                                                "se.exposure", "pval.exposure", "chr_pos")

lifetime_exp_dat$id.exposure <- "lifetime_exp_dat"
lifetime_exp_dat$exposure <- "lifetime_exp_dat"

dat_1_smk2edu <- harmonise_data(lifetime_exp_dat, edu_out_dat, action = 3)

#Nice, now let's do it the other way around:

dat_1_smk2edu <- dat_1_smk2edu[which(dat_1_smk2edu$remove == FALSE),]
dat_1_smk2edu <- dat_1_smk2edu[which(dat_1_smk2edu$palindromic == FALSE & dat_1_smk2edu$ambiguous == FALSE | dat_1_smk2edu$palindromic == TRUE & dat_1_smk2edu$ambiguous == FALSE),]

#########################
#We end up with 282 SNPs#
#########################

#Let's obtain where we obtained them from:

final_set_of_snps_no_dupl_clean <- final_set_of_snps_no_dupl[which(final_set_of_snps_no_dupl$rsid%in%dat_1_smk2edu$SNP),]
final_set_of_snps_no_dupl_clean <- final_set_of_snps_no_dupl_clean[order(match(final_set_of_snps_no_dupl_clean$rsid, dat_1_smk2edu$SNP)),]

#There is a mismatch. We have 282 in dat_1_smk, but 281 here. 
#This is probably due to the duplicate fella.

dat_1_smk2edu[which(duplicated(dat_1_smk2edu$SNP) == TRUE),] #none

#That means that all SNPs should be fine.

length(which(final_set_of_snps_no_dupl_clean$rsid != dat_1_smk2edu$SNP)) #perfect
length(which(final_set_of_snps_no_dupl_clean$rsid == dat_1_smk2edu$SNP)) #perfect

dat_1_smk2edu$source <- final_set_of_snps_no_dupl_clean$source #with this we can now where the SNPs come from.

######################################
#Let's continue and generate the data#
######################################

dat_1_smk2edu$samplesize.exposure <- 462690

final_df <- dat_1_smk2edu %>%
  select(SNP, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, eaf.exposure, pval.exposure, samplesize.exposure, beta.outcome, se.outcome, eaf.outcome, pval.outcome, samplesize.outcome, chr_pos, source)

colnames(final_df) <- c("SNP", "effect_allele", "other_allele", "beta.smk", "se.smk", "eaf.smk", "pval.smk", "samplesize.smk", "beta.edu", "se.edu", "eaf.edu", "pval.edu", "samplesize.edu", "chr_pos", "source")

fwrite(final_df, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/Curated_data/LifeTime_edu_Risk_Factors.txt")
