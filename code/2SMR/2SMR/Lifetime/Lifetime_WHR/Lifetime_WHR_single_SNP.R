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

lifetime_smk <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/Lifetime_smk/Lifetime_Smoking_Curated.txt")

#########################################
#2. Get the genome-wide significant SNPs#
#########################################

rs1051730_exp <- lifetime_smk[which(lifetime_smk$SNP == "rs1051730"),]
rs16969968_exp <- lifetime_smk[which(lifetime_smk$SNP == "rs16969968"),]
rs1317286_exp <- lifetime_smk[which(lifetime_smk$SNP == "rs1317286"),]

#Checking pvals:

rs1051730_exp

#SNP CHR       BP EFFECT_ALLELE OTHER_ALLELE      EAF INFO       BETA         SE       P
#1: rs1051730  15 78894339             G            A 0.668648    1 -0.0183306 0.00147226 1.4e-35
#chr_pos
#1: chr15:78894339

#Really strong

rs16969968_exp

#SNP CHR       BP EFFECT_ALLELE OTHER_ALLELE      EAF INFO       BETA        SE       P
#1: rs16969968  15 78882925             G            A 0.669262    1 -0.0185748 0.0014734 1.9e-36
#chr_pos
#1: chr15:78882925

rs1317286_exp

#SNP CHR       BP EFFECT_ALLELE OTHER_ALLELE      EAF INFO       BETA
#1: rs1317286  15 78896129             A            G 0.666723    1 -0.0184813
#SE       P        chr_pos
#1: 0.00147059 3.2e-36 chr15:78896129

#Interestingly the effects are different!!!

#Really strong!!

#As they said, it would be great to see the effects, so let's continue.

#################################################
#SECTION B: FINDING matching SNPs in the outcome#
#################################################

#First let's load the outcome:

whr <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHR/WHR_combined_Curated.txt")

##########################################################################################
#1. we are gonna try to match the data with RSID and with chr_pos to get as many variants#
##########################################################################################

whr_2_rs1051730 <- whr[which(whr$MarkerName%in%rs1051730_exp$SNP),] #gotcha
whr_2_rs16969968 <- whr[which(whr$MarkerName%in%rs16969968_exp$SNP),] #gotcha
whr_2_rs1317286 <- whr[which(whr$MarkerName%in%rs1317286_exp$SNP),] #gotcha

##########################################################
#Let's change the dataframes so we can harmonise the data#
##########################################################

#To properly match the data we need to do some tricks here:

rs1051730_exp_clean <- rs1051730_exp

colnames(rs1051730_exp_clean) <-  c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
    "other_allele.exposure", "eaf.exposure", "info.exposure", "beta.exposure",
    "se.exposure", "pval.exposure", "chr_pos")

rs16969968_exp_clean <- rs16969968_exp

colnames(rs16969968_exp_clean) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
                                   "other_allele.exposure", "eaf.exposure", "info.exposure", "beta.exposure",
                                   "se.exposure", "pval.exposure", "chr_pos")

#And finally for the SNP we found:

rs1317286_exp_clean <- rs1317286_exp

colnames(rs1317286_exp_clean) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
                                    "other_allele.exposure", "eaf.exposure", "info.exposure", "beta.exposure",
                                    "se.exposure", "pval.exposure", "chr_pos")

#We need to add some stuff:

rs16969968_exp_clean$samplesize.exposure <- 461066
rs1051730_exp_clean$samplesize.exposure <- 461066
rs1317286_exp_clean$samplesize.exposure <- 461066


rs16969968_exp_clean$id.exposure <- "lifetime_smk"
rs16969968_exp_clean$exposure <- "lifetime_smk"

rs1051730_exp_clean$id.exposure <- "lifetime_smk"
rs1051730_exp_clean$exposure <- "lifetime_smk"

rs1317286_exp_clean$id.exposure <- "lifetime_smk"
rs1317286_exp_clean$exposure <- "lifetime_smk"

##############################
#Worked like a freaking charm#
##############################

##############################################
#SECTION C: Merging exposure and outcome data#
##############################################

#We need to format the data for the outcome:

whr_outcome_rs1051730 <- whr_2_rs1051730

colnames(whr_outcome_rs1051730) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", 
                           "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                           "pval.outcome", "samplesize.outcome", "chr_pos_37", "chr_pos_18")

whr_outcome_rs16969968 <- whr_2_rs16969968

colnames(whr_outcome_rs16969968) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", 
                           "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                           "pval.outcome", "samplesize.outcome", "chr_pos_37", "chr_pos_18")

#And finally:

whr_outcome_rs1317286 <- whr_2_rs1317286

colnames(whr_outcome_rs1317286) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", 
                                      "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                                      "pval.outcome", "samplesize.outcome", "chr_pos_37", "chr_pos_18")


#Getting ids, so troublesome...

whr_outcome_rs16969968$id.outcome <- "whr"
whr_outcome_rs16969968$outcome <- "whr"

whr_outcome_rs1051730$id.outcome <- "whr"
whr_outcome_rs1051730$outcome <- "whr"

whr_outcome_rs1317286$id.outcome <- "whr"
whr_outcome_rs1317286$outcome <- "whr"

dat_1_rs1051730 <- harmonise_data(rs1051730_exp_clean, whr_outcome_rs1051730, action = 3) #76/86: perfect.
dat_1_rs16969968 <- harmonise_data(rs16969968_exp_clean, whr_outcome_rs16969968, action = 3) #76/86: perfect.
dat_1_rs1317286 <- harmonise_data(rs1317286_exp_clean, whr_outcome_rs1317286, action = 3) #76/86: perfect.

#######################
#SECTION D: mF and IVW#
#######################

F_ = ((dat_1_rs1051730$beta.exposure)^2)/((dat_1_rs1051730$se.exposure)^2)
mF  = mean(F_)
print(mF)
#155.0188

Isq(dat_1_rs1051730$beta.exposure, dat_1_rs1051730$se.exposure)
#NaN

##NEXT SNP:

F_ = ((dat_1_rs16969968$beta.exposure)^2)/((dat_1_rs16969968$se.exposure)^2)
mF  = mean(F_)
print(mF)
#158.9304

Isq(dat_1_rs16969968$beta.exposure, dat_1_rs16969968$se.exposure)
#NaN

##NEXT SNP:

F_ = ((dat_1_rs1317286$beta.exposure)^2)/((dat_1_rs1317286$se.exposure)^2)
mF  = mean(F_)
print(mF)
#157.9362

Isq(dat_1_rs1317286$beta.exposure, dat_1_rs1317286$se.exposure)
#NaN

################################
#Let's run the shit out of this#
################################

mr(dat_1_rs1051730)

#id.exposure id.outcome outcome     exposure     method nsnp           b        se     pval
#1 lifetime_smk        whr     whr lifetime_smk Wald ratio    1 -0.04364287 0.1909376 0.819202

mr(dat_1_rs16969968)

#id.exposure id.outcome outcome     exposure     method nsnp           b        se     pval
#1 lifetime_smk        whr     whr lifetime_smk Wald ratio    1 -0.05922002 0.1884273 0.753304

mr(dat_1_rs1317286)

#id.exposure id.outcome outcome     exposure     method nsnp        b        se
#1 lifetime_smk        whr     whr lifetime_smk Wald ratio    1 0.129861 0.2380785
#pval
#1 0.5854409

#POSITIVE, but not significant.