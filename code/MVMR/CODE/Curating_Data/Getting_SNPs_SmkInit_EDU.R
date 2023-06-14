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

SMK <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/RAW_Data/SmokingInitiation.txt.gz")

SMK$INFO <- SMK$EFFECTIVE_N/SMK$N
SMK$chr_pos_37 <- paste("chr", SMK$CHROM, ":", SMK$POS, sep = "")
summary(SMK$INFO) #0.3 is minumum as described by README.

#Filtering GW:

SMK_gw <- SMK[which(as.numeric(SMK$PVALUE) < 0.00000005),] 

#We have to be wary. We need to remove MHC data and those with weird MAFs.Also check INFO if possible.

summary(SMK_gw$PVALUE) #great
summary(SMK_gw$AF) #MAF > 0.01, as w want
summary(SMK_gw$INFO) #>0.70 as we want

#We just need to remove those in the MHC region.

SMK_mhc <- SMK_gw[which(as.numeric(SMK_gw$CHROM) == 6 & as.numeric(SMK_gw$POS) >= 26000000 & as.numeric(SMK_gw$POS) <= 34000000),]

#SMK_mhc == 0. I remember seeing this while working on 2SMR. Let's check if this is true:

SMK_check <- SMK_gw[which(as.numeric(SMK_gw$CHROM) == 6),]

SMK_check <- SMK_check[order(as.numeric(SMK_check$POS)),]

View(SMK_check) #ALL GOOD. We can proceed with SMK_gw.

#Preparing data to obtain, later the genome-wide significant SNPs.

SMK_gw$rsid <- SMK_gw$RSID
SMK_gw$pval <- SMK_gw$PVALUE

#Let's check if we have any weird SNP.

SMK_gw <- SMK_gw[order(SMK_gw$RSID),]

head(SMK_gw) #all good
tail(SMK_gw) #all good

#Finally, let's remove the data weird data

pre_independent_set_1 <- SMK_gw %>%
  select(rsid, pval)

pre_independent_set_1$source <- "Smoking Initiation"

#####################
##A.2 For Education##
#####################

#And now get the genome-wide significant for EDU:

edu_dat <- extract_instruments(outcomes = "ukb-b-16489", p1 = 0.00000005, clump = FALSE) #25926

#There is only one duplicate, but they are keeping only the first instance. 
#How do we deal with this? 
#We should check both of them and see if it is important.

#1. Is the duplicate genome-wide significant or is it just coincidence? Also we should check the data.

summary(edu_dat$pval.exposure) #it is genome-wide significant. 
summary(edu_dat$eaf.exposure) #all good.
colnames(edu_dat) #no info data.

#2. Can we remove it since it is not in Smoking Initiation?

SMK[which(SMK$RSID == "rs2696531"),] #the SNP is not even here, so no worries, yet.

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

independent_snps <- rbind(pre_independent_set_1, pre_independent_set_2) #33146 SNPs!

#Now we have the combined set of genome-wide significant SNPs, but we have to make them truly independent.

######################################################
#B investigate the duplicates and perform LD clumping#
######################################################

################################
#B.1 Let's check the duplicates#
################################

no_dupl <- unique(independent_snps$rsid) #we have 32077.

#if we order by p-value and remove them:

independent_snps <- independent_snps[order(independent_snps$pval),]

independent_snps_no_dupl <- independent_snps[-which(duplicated(independent_snps$rsid) == TRUE),] #32077

#Perfect!! We get the same amout of SNPs. That means that we do not have more special cases like the one above.

#There are two cases of duplicates:

#1. The special case in which you found them in the same exposure data.
#2. Those that appear because they are in the two exposures.

#This code above takes care of case 2, since we are taking the p-value of the exposure
#that the SNP is more associated with the SNP. Case 1 is gonna drag us down a bit.

########################################
#B.2 Now let's get the independent ones#
########################################

final_set_of_snps_no_dupl <- ieugwasr::ld_clump_local(independent_snps_no_dupl, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.00000005) #289

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

#We obtain 290. Is it our friend from case 1?

edu_out_dat$chr_pos_37 <- paste("chr", edu_out_dat$chr, ":", edu_out_dat$pos, sep = "")

#####################################################
#C. 2 Now obtain the summary statistics of Lifetime #
#####################################################

SMK_exp_dat <- SMK[which(SMK$chr_pos_37%in%edu_out_dat$chr_pos_37),] #288!

#There are missing SNPs. But how can that be?

missing_snps <- edu_out_dat[which(!(edu_out_dat$SNP%in%SMK_exp_dat$RSID)),] #2

missing_snps

#SNP chr      pos beta.outcome se.outcome samplesize.outcome pval.outcome eaf.outcome
#173 rs2606913  16 14916988   0.00748826 0.00127209             458079  3.89996e-09    0.822881
#286 rs2606913  16 14916988   0.00888099 0.00660753             458079  1.80000e-01    0.008769

#There is only one missing SNP, really. Our duplicated fella. 
#We will need to get the strongest proxy possible that is available in SMK and in education.

setwd("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/Curating_Data/proxies/SmkInit_EDU/")

getwd()

#Now we get the proxies.

LDlinkR::LDproxy_batch(snp = "rs2606913", pop = "EUR", r2d = "r2", append = TRUE, token = "04cad4ca4374")

proxies_df <- read.table("combined_query_snp_list.txt", header = TRUE)

proxies_df <- proxies_df[which(proxies_df$R2 > 0.8),] #in this way we obtain those that are not in the same loci.

head(proxies_df)

#Now let's get the matches in both, smoking and lifetime.
#Then we will choose the strongest proxy for education (since the SNP comes from there)

SMK_match_proxy <- SMK[which(SMK$chr_pos_37%in%proxies_df$Coord),]  #15/34 are SMK

#Now let's get those that are found here in edu:

edu_proxy <- extract_outcome_data(
  snps = SMK_match_proxy$RSID,
  outcomes = 'ukb-b-16489',
  proxies = FALSE
) #14/15 are in Education.

#Now let's prune them and get the one with the strongest association.
#Let's check the associations though:

summary(edu_proxy$pval.outcome) #not all of them are genome-wide, but most are.

edu_proxy$rsid <- edu_proxy$SNP
edu_proxy$pval <- edu_proxy$pval.outcome

final_proxy <- ieugwasr::ld_clump_local(edu_proxy, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.00000005) #1, as expected.

###################################
#This is the fella that we want!!!#
###################################

#First we remove the original SNP:

edu_out_dat <- edu_out_dat[-which(edu_out_dat$SNP == "rs2606913"),] 

#And now we add the other fella we need to select it.

edu_new_data <- edu_proxy[which(edu_proxy$SNP%in%final_proxy$SNP),]

#And now make two changes to the dataframes:

edu_out_dat$rsid <- edu_out_dat$SNP
edu_out_dat$pval <- edu_out_dat$pval.outcome

#And now add chr_pos_37 to the other one:

edu_new_data$chr_pos_37 <- paste("chr", edu_new_data$chr, ":", edu_new_data$pos, sep = "")

edu_out_dat_end <- rbind(edu_out_dat, edu_new_data) #PERFECT.

#Now let's repeat obtaining the exposure data:

SMK_exp_dat <- SMK[which(SMK$chr_pos_37%in%edu_out_dat_end$chr_pos_37),] #289! WE DID IT.

######################
#C.3 Merging the data#
######################

colnames(SMK_exp_dat) <- c("chr.exposure", "pos.exposure", "SNP", "other_allele.exposure", "effect_allele.exposure", 
                           "eaf.exposure",  "stat.exposure", "pval.exposure", "beta.exposure", "se.exposure", 
                           "samplesize.exposure", "effective_sample.exposure", "studies.exposure", "annotation.exposure", "annotation_full.exposure",
                           "info.exposure", "chr_pos_37")

SMK_exp_dat$id.exposure <- "Smoking Initiation"
SMK_exp_dat$exposure <- "Smoking Initiation"

dat_1_smk2edu <- harmonise_data(SMK_exp_dat, edu_out_dat, action = 3)

#Nice, now let's do it the other way around:

dat_1_smk2edu <- dat_1_smk2edu[which(dat_1_smk2edu$remove == FALSE),]
dat_1_smk2edu <- dat_1_smk2edu[which(dat_1_smk2edu$palindromic == FALSE & dat_1_smk2edu$ambiguous == FALSE | dat_1_smk2edu$palindromic == TRUE & dat_1_smk2edu$ambiguous == FALSE),]

#########################
#We end up with 275 SNPs#
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
  select(SNP, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, eaf.exposure, pval.exposure, samplesize.exposure, beta.outcome, se.outcome, eaf.outcome, pval.outcome, samplesize.outcome, chr_pos_37.x, source)

colnames(final_df) <- c("SNP", "effect_allele", "other_allele", "beta.smk", "se.smk", "eaf.smk", "pval.smk", "samplesize.smk", "beta.edu", "se.edu", "eaf.edu", "pval.edu", "samplesize.edu", "chr_pos", "source")

fwrite(final_df, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/Curated_data/SmkInit_EDU_Risk_Factors.txt")
