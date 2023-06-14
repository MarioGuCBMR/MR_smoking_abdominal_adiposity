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

#Note the data is packyears smoking:

##################
#1. Load the data#
##################

packyears_smk <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/PacksPerYear/PackPerYear_Curated_FULL.txt")

#########################################
#2. Get the genome-wide significant SNPs#
#########################################

rs1051730_exp <- packyears_smk[which(packyears_smk$SNP == "rs1051730"),]
rs16969968_exp <- packyears_smk[which(packyears_smk$SNP == "rs16969968"),]
rs1317286_exp <- packyears_smk[which(packyears_smk$SNP == "rs1317286"),]

#Checking pvals:

rs1051730_exp

#chr.exposure pos.exposure       SNP other_allele.exposure effect_allele.exposure dot PASS          AF
#1:           15     78894339 rs1051730                     G                      A   . PASS AF=0.337514
#Headers                                    Headers_info beta.exposure se.exposure    logp
#1: ES:SE:LP:AF:ID 0.0730912:0.00366914:87.5686:0.337514:rs1051730     0.0730912  0.00366914 87.5686
#eaf.exposure      rsid logp_neg pval.exposure
#1:     0.337514 rs1051730 -87.5686  2.700225e-88

#Really strong

rs16969968_exp

#chr.exposure pos.exposure        SNP other_allele.exposure effect_allele.exposure dot PASS          AF
#1:           15     78882925 rs16969968                     G                      A   . PASS AF=0.336998
#Headers                                     Headers_info beta.exposure se.exposure    logp
#1: ES:SE:LP:AF:ID 0.0729731:0.00367071:87.2147:0.336998:rs16969968     0.0729731  0.00367071 87.2147
#eaf.exposure       rsid logp_neg pval.exposure
#1:     0.336998 rs16969968 -87.2147  6.099581e-88

rs1317286_exp

#chr.exposure pos.exposure       SNP other_allele.exposure effect_allele.exposure dot PASS          AF
#1:           15     78896129 rs1317286                     A                      G   . PASS AF=0.339506
#Headers                                   Headers_info beta.exposure se.exposure   logp
#1: ES:SE:LP:AF:ID 0.0734398:0.00366497:88.585:0.339506:rs1317286     0.0734398  0.00366497 88.585
#eaf.exposure      rsid logp_neg pval.exposure
#1:     0.339506 rs1317286  -88.585   2.60016e-89

#Really strong!!

#As they said, it would be great to see the effects, so let's continue.

#################################################
#SECTION B: FINDING matching SNPs in the outcome#
#################################################

#First let's load the outcome:

wc <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WC/WC_combined_Curated.txt")

##########################################################################################
#1. we are gonna try to match the data with RSID and with chr_pos to get as many variants#
##########################################################################################

wc_2_rs1051730 <- wc[which(wc$MarkerName%in%rs1051730_exp$SNP),] #gotcha
wc_2_rs16969968 <- wc[which(wc$MarkerName%in%rs16969968_exp$SNP),] #gotcha
wc_2_rs1317286 <- wc[which(wc$MarkerName%in%rs1317286_exp$SNP),] #gotcha

##########################################################
#Let's change the dataframes so we can harmonise the data#
##########################################################

#To properly match the data we need to do some tricks here:

rs1051730_exp_clean <- rs1051730_exp

rs16969968_exp_clean <- rs16969968_exp

rs1317286_exp_clean <- rs1317286_exp

#We need to add some stuff:

rs16969968_exp_clean$samplesize.exposure <- 142387
rs1051730_exp_clean$samplesize.exposure <- 142387
rs1317286_exp_clean$samplesize.exposure <- 142387

rs16969968_exp_clean$id.exposure <- "packyears_smk"
rs16969968_exp_clean$exposure <- "packyears_smk"

rs1051730_exp_clean$id.exposure <- "packyears_smk"
rs1051730_exp_clean$exposure <- "packyears_smk"

rs1317286_exp_clean$id.exposure <- "packyears_smk"
rs1317286_exp_clean$exposure <- "packyears_smk"

##############################
#Worked like a freaking charm#
##############################

##############################################
#SECTION C: Merging exposure and outcome data#
##############################################

#We need to format the data for the outcome:

wc_outcome_rs1051730 <- wc_2_rs1051730

colnames(wc_outcome_rs1051730) <- c("SNP", "effect_allele.outcome", 
                                     "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                                     "pval.outcome", "samplesize.outcome", "chr_pos_37")

wc_outcome_rs16969968 <- wc_2_rs16969968

colnames(wc_outcome_rs16969968) <- c("SNP","effect_allele.outcome", 
                                      "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                                      "pval.outcome", "samplesize.outcome", "chr_pos_37")

wc_outcome_rs1317286 <- wc_2_rs1317286

colnames(wc_outcome_rs1317286) <- c("SNP","effect_allele.outcome", 
                                     "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                                     "pval.outcome", "samplesize.outcome", "chr_pos_37")


#Getting ids, so troublesome...

wc_outcome_rs16969968$id.outcome <- "wc"
wc_outcome_rs16969968$outcome <- "wc"

wc_outcome_rs1051730$id.outcome <- "wc"
wc_outcome_rs1051730$outcome <- "wc"

wc_outcome_rs1317286$id.outcome <- "wc"
wc_outcome_rs1317286$outcome <- "wc"

dat_1_rs1051730 <- harmonise_data(rs1051730_exp_clean, wc_outcome_rs1051730, action = 3) #76/86: perfect.
dat_1_rs16969968 <- harmonise_data(rs16969968_exp_clean, wc_outcome_rs16969968, action = 3) #76/86: perfect.
dat_1_rs1317286 <- harmonise_data(rs1317286_exp_clean, wc_outcome_rs1317286, action = 3) #76/86: perfect.

#######################
#SECTION D: mF and IVW#
#######################

F_ = ((dat_1_rs1051730$beta.exposure)^2)/((dat_1_rs1051730$se.exposure)^2)
mF  = mean(F_)
print(mF)
#396.8274

Isq(dat_1_rs1051730$beta.exposure, dat_1_rs1051730$se.exposure)
#NaN

##NEXT SNP:

F_ = ((dat_1_rs16969968$beta.exposure)^2)/((dat_1_rs16969968$se.exposure)^2)
mF  = mean(F_)
print(mF)
#395.2077

Isq(dat_1_rs1051730$beta.exposure, dat_1_rs1051730$se.exposure)
#NaN

F_ = ((dat_1_rs1317286$beta.exposure)^2)/((dat_1_rs1317286$se.exposure)^2)
mF  = mean(F_)
print(mF)
#401.5338

Isq(dat_1_rs1317286$beta.exposure, dat_1_rs1317286$se.exposure)
#NaN

################################
#Let's run the shit out of this#
################################

mr(dat_1_rs1051730)

#id.exposure id.outcome outcome      exposure     method nsnp           b         se      pval
#1 packyears_smk         wc      wc packyears_smk Wald ratio    1 -0.07935292 0.04925354 0.1071555

mr(dat_1_rs16969968)

#id.exposure id.outcome outcome      exposure     method nsnp           b         se       pval
#1 packyears_smk         wc      wc packyears_smk Wald ratio    1 -0.08496282 0.04933325 0.08502927

mr(dat_1_rs1317286)

#id.exposure id.outcome outcome      exposure     method nsnp           b         se      pval
#1 packyears_smk         wc      wc packyears_smk Wald ratio    1 -0.07216795 0.06127468 0.2388852
