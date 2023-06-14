##############
#INTRODUCTION#
##############

#This is a code to run 2SMR on single SNPs that hit the CHRNA5-CHRNA3-CHRNB4 loci.

###################
#Loading libraries#
###################

library(TwoSampleMR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gridExtra)
library(rmarkdown)
library(TwoSampleMR)
library(MRInstruments)
library(jsonlite)
library(httr)
library(data.table)

#########################################
#SECTION A: getting the SNPs of interest#
#########################################

#Note the data is lifetime smoking:

##################
#1. Load the data#
##################

cigday <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/CigDayCurrent/Cig_Day_Current_Curated_FULL.txt")

#########################################
#2. Get the genome-wide significant SNPs#
#########################################

rs1051730_exp <- cigday[which(cigday$SNP == "rs1051730"),]
rs16969968_exp <- cigday[which(cigday$SNP == "rs16969968"),]
rs1317286_exp <- cigday[which(cigday$SNP == "rs1317286"),]

#Checking pvals:

rs1051730_exp

#chr.exposure pos.exposure       SNP other_allele.exposure
#1:           15     78894339 rs1051730                     G
#effect_allele.exposure dot PASS          AF        Headers
#1:                      A   . PASS AF=0.341409 ES:SE:LP:AF:ID
#Headers_info beta.exposure se.exposure
#1: 0.0764471:0.00574399:39.6778:0.341409:rs1051730     0.0764471  0.00574399
#logp eaf.exposure      rsid logp_neg pval.exposure
#1: 39.6778     0.341409 rs1051730 -39.6778  2.099907e-40

#Really strong

rs16969968_exp

#chr.exposure pos.exposure        SNP other_allele.exposure
#1:           15     78882925 rs16969968                     G
#effect_allele.exposure dot PASS          AF        Headers
#1:                      A   . PASS AF=0.340882 ES:SE:LP:AF:ID
#Headers_info beta.exposure se.exposure
#1: 0.0762056:0.00574466:39.4318:0.340882:rs16969968     0.0762056  0.00574466
#logp eaf.exposure       rsid logp_neg pval.exposure
#1: 39.4318     0.340882 rs16969968 -39.4318  3.699985e-40

rs1317286_exp

#chr.exposure pos.exposure       SNP other_allele.exposure
#1:           15     78896129 rs1317286                     A
#effect_allele.exposure dot PASS          AF        Headers
#1:                      G   . PASS AF=0.343245 ES:SE:LP:AF:ID
#Headers_info beta.exposure se.exposure
#1: 0.0766873:0.00573869:40.0044:0.343245:rs1317286     0.0766873  0.00573869
#logp eaf.exposure      rsid logp_neg pval.exposure
#1: 40.0044     0.343245 rs1317286 -40.0044  9.899198e-41

#Really strong!!

#As they said, it would be great to see the effects, so let's continue.

#################################################
#SECTION B: FINDING matching SNPs in the outcome#
#################################################

#First let's load the outcome:

WCAdjBMI <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WCAdjBMIsmk/WCAdjBMI_Smk_Strat_Curated.txt")

##########################################################################################
#1. we are gonna try to match the data with RSID and with chr_pos to get as many variants#
##########################################################################################

WCAdjBMI_2_rs1051730 <- WCAdjBMI[which(WCAdjBMI$rs_id%in%rs1051730_exp$SNP),] #gotcha
WCAdjBMI_2_rs16969968 <- WCAdjBMI[which(WCAdjBMI$rs_id%in%rs16969968_exp$SNP),] #gotcha
WCAdjBMI_2_rs1317286 <- WCAdjBMI[which(WCAdjBMI$rs_id%in%rs1317286_exp$SNP),] #gotcha

##########################################################
#Let's change the dataframes so we can harmonise the data#
##########################################################

#To properly match the data we need to do some tricks here:

rs1051730_exp_clean <- rs1051730_exp
rs16969968_exp_clean <- rs16969968_exp
rs1317286_exp_clean <- rs1317286_exp


#We need to add some stuff:

rs16969968_exp_clean$samplesize.exposure <- 		33229
rs1051730_exp_clean$samplesize.exposure <- 		33229
rs1317286_exp_clean$samplesize.exposure <- 		33229

rs16969968_exp_clean$id.exposure <- "cigday"
rs16969968_exp_clean$exposure <- "cigday"

rs1051730_exp_clean$id.exposure <- "cigday"
rs1051730_exp_clean$exposure <- "cigday"

rs1317286_exp_clean$id.exposure <- "cigday"
rs1317286_exp_clean$exposure <- "cigday"

##############################
#Worked like a freaking charm#
##############################

##############################################
#SECTION C: Merging exposure and outcome data#
##############################################

#We need to format the data for the outcome:

WCAdjBMI_outcome_rs1051730 <- WCAdjBMI_2_rs1051730

colnames(WCAdjBMI_outcome_rs1051730) <- c("chr.outcome", "SNP", "rs_id", "position_hg18",
                                           "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "samplesize.outcome",
                                           "beta.outcome", "se.outcome", "pval.outcome", "N_NonSmk", 
                                           "Effect_NonSMK", "StdErr_NonSMK", "P_value_NonSMK", "chr_pos_37",
                                           "chr_pos_18")

WCAdjBMI_outcome_rs16969968 <- WCAdjBMI_2_rs16969968

colnames(WCAdjBMI_outcome_rs16969968) <- c("chr.outcome", "SNP", "rs_id", "position_hg18",
                                            "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "samplesize.outcome",
                                            "beta.outcome", "se.outcome", "pval.outcome", "N_NonSmk", 
                                            "Effect_NonSMK", "StdErr_NonSMK", "P_value_NonSMK", "chr_pos_37",
                                            "chr_pos_18")


WCAdjBMI_outcome_rs1317286 <- WCAdjBMI_2_rs1317286

colnames(WCAdjBMI_outcome_rs1317286) <- c("chr.outcome", "SNP", "rs_id", "position_hg18",
                                           "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "samplesize.outcome",
                                           "beta.outcome", "se.outcome", "pval.outcome", "N_NonSmk", 
                                           "Effect_NonSMK", "StdErr_NonSMK", "P_value_NonSMK", "chr_pos_37",
                                           "chr_pos_18")


#Getting ids, so troublesome...

WCAdjBMI_outcome_rs16969968$id.outcome <- "WCAdjBMI"
WCAdjBMI_outcome_rs16969968$outcome <- "WCAdjBMI"

WCAdjBMI_outcome_rs1051730$id.outcome <- "WCAdjBMI"
WCAdjBMI_outcome_rs1051730$outcome <- "WCAdjBMI"

WCAdjBMI_outcome_rs1317286$id.outcome <- "WCAdjBMI"
WCAdjBMI_outcome_rs1317286$outcome <- "WCAdjBMI"


dat_1_rs1051730 <- harmonise_data(rs1051730_exp_clean, WCAdjBMI_outcome_rs1051730, action = 3) #76/86: perfect.
dat_1_rs16969968 <- harmonise_data(rs16969968_exp_clean, WCAdjBMI_outcome_rs16969968, action = 3) #76/86: perfect.
dat_1_rs1317286 <- harmonise_data(rs1317286_exp_clean, WCAdjBMI_outcome_rs1317286, action = 3) #76/86: perfect.

#######################
#SECTION D: mF and IVW#
#######################

F_ = ((dat_1_rs1051730$beta.exposure)^2)/((dat_1_rs1051730$se.exposure)^2)
mF  = mean(F_)
print(mF)
#177.131

Isq(dat_1_rs1051730$beta.exposure, dat_1_rs1051730$se.exposure)
#NaN

##NEXT SNP:

F_ = ((dat_1_rs16969968$beta.exposure)^2)/((dat_1_rs16969968$se.exposure)^2)
mF  = mean(F_)
print(mF)
#175.9726

Isq(dat_1_rs1051730$beta.exposure, dat_1_rs1051730$se.exposure)
#NaN

##NEXT SNP:

F_ = ((dat_1_rs1317286$beta.exposure)^2)/((dat_1_rs1317286$se.exposure)^2)
mF  = mean(F_)
print(mF)
#178.5753

Isq(dat_1_rs1317286$beta.exposure, dat_1_rs1317286$se.exposure)

################################
#Let's run the shit out of this#
################################

mr(dat_1_rs1051730)

#id.exposure id.outcome  outcome exposure     method nsnp          b
#1      cigday   WCAdjBMI WCAdjBMI   cigday Wald ratio    1 0.08633421
#se      pval
#1 0.1046475 0.4093716

mr(dat_1_rs16969968)

#id.exposure id.outcome  outcome exposure     method nsnp         b        se
#1      cigday   WCAdjBMI WCAdjBMI   cigday Wald ratio    1 0.1102281 0.1036669
#pval
#1 0.2876499

mr(dat_1_rs1317286)

#id.exposure id.outcome  outcome exposure     method nsnp         b        se
#1      cigday   WCAdjBMI WCAdjBMI   cigday Wald ratio    1 0.2151595 0.1134477
#pval
#1 0.05788713