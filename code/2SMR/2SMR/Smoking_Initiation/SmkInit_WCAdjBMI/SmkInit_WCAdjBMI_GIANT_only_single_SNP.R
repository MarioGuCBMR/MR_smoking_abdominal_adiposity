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

#We are going to load two types of data and check whether we can find the SNPs with the highest sample size or not.

SMK <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/RAW_Data/SmokingInitiation.txt.gz")
SMK$chr_pos_37 <- paste("chr", SMK$CHROM, ":", SMK$POS, sep = "")
SMK$INFO <- SMK$EFFECTIVE_N/SMK$N

smk_gw <- readxl::read_xlsx("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/RAW_Data/SmkInit_conditionally_independent_snps.xlsx")

colnames(smk_gw) <- smk_gw[1,]

smk_gw <- smk_gw[-1,]

#########################################
#2. Get the genome-wide significant SNPs#
#########################################

rs1051730_exp <- SMK[which(SMK$RSID == "rs1051730"),]
rs16969968_exp <- SMK[which(SMK$RSID == "rs16969968"),]
rs1317286_exp <- SMK[which(SMK$RSID == "rs1317286"),]

#Can we get them with the higher sample size:

check_rs1051730_exp <- smk_gw[which(smk_gw$rsID == "rs1051730"),] #not there
check_rs16969968_exp <- smk_gw[which(smk_gw$rsID == "rs16969968"),] #not there
check_rs1317286_exp <- smk_gw[which(smk_gw$rsID == "rs1317286"),]  #not there

########################################################################
#Perfect, so we just run the analysis with the SmokingInitiation.txt.gz#
########################################################################

#Checking pvals:

rs1051730_exp

#CHROM      POS      RSID REF ALT    AF STAT  PVALUE        BETA          SE      N EFFECTIVE_N
#1:    15 78894339 rs1051730   G   A 0.324 7.37 0.00663 -0.01026409 0.003780826 632802      627497
#Number_of_Studies              ANNO
#1:                34 Synonymous:CHRNA3

#The p-value is not strong at all.

rs16969968_exp

#CHROM      POS       RSID REF ALT    AF STAT PVALUE         BETA          SE      N EFFECTIVE_N
#1:    15 78882925 rs16969968   G   A 0.322 6.45 0.0111 -0.009609521 0.003783745 632802      628725
#Number_of_Studies                 ANNO
#1:                34 Nonsynonymous:CHRNA5
#ANNOFULL
#1: CHRNA5/NM_000745:+:Nonsynonymous(GAT/Asp/D->AAT/Asn/N:Base1193/1407:Codon398/469:Exon5/6):Exon
#chr_pos_37      INFO
#1: chr15:78882925 0.9935572

rs1317286_exp

#CHROM      POS      RSID REF ALT    AF STAT PVALUE         BETA          SE      N EFFECTIVE_N
#1:    15 78896129 rs1317286   A   G 0.332 6.78 0.0092 -0.009814893 0.003769387 632802      627305
#Number_of_Studies          ANNO                                               ANNOFULL
#1:                34 Intron:CHRNA3 CHRNA3/NM_001166694:-:Intron|CHRNA3/NM_000743:-:Intron
#chr_pos_37      INFO
#1: chr15:78896129 0.9913132

#################################################
#SECTION B: FINDING matching SNPs in the outcome#
#################################################

#First let's load the outcome:

WCAdjBMI <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WCAdjBMI/WCAdjBMI_combined_Curated.txt")

##########################################################################################
#1. we are gonna try to match the data with RSID and with chr_pos to get as many variants#
##########################################################################################

WCAdjBMI_2_rs1051730 <- WCAdjBMI[which(WCAdjBMI$MarkerName%in%rs1051730_exp$RSID),] #gotcha
WCAdjBMI_2_rs16969968 <- WCAdjBMI[which(WCAdjBMI$MarkerName%in%rs16969968_exp$RSID),] #gotcha
WCAdjBMI_2_rs1317286 <- WCAdjBMI[which(WCAdjBMI$MarkerName%in%rs1317286_exp$RSID),] #gotcha

##########################################################
#Let's change the dataframes so we can harmonise the data#
##########################################################

#To properly match the data we need to do some tricks here:

rs1051730_exp_clean <- rs1051730_exp

colnames(rs1051730_exp_clean) <-  c("chr.exposure", "pos.exposure", "SNP", "other_allele.exposure",
    "effect_allele.exposure", "eaf.exposure", "stat.exposure", "pval.exposure", "beta.exposure", "se.exposure",
    "samplesize.exposure", "effective_sample.exposure", "studies.exposure", "anno.exposure", "anno_full.exposure",
    "chr_pos_37", "info.exposure")

rs16969968_exp_clean <- rs16969968_exp

colnames(rs16969968_exp_clean) <- c("chr.exposure", "pos.exposure", "SNP", "other_allele.exposure",
                                   "effect_allele.exposure", "eaf.exposure", "stat.exposure", "pval.exposure", "beta.exposure", "se.exposure",
                                   "samplesize.exposure", "effective_sample.exposure", "studies.exposure", "anno.exposure", "anno_full.exposure",
                                   "chr_pos_37", "info.exposure")

#And finally for the SNP we found:

rs1317286_exp_clean <- rs1317286_exp

colnames(rs1317286_exp_clean) <- c("chr.exposure", "pos.exposure", "SNP", "other_allele.exposure",
                                   "effect_allele.exposure", "eaf.exposure", "stat.exposure", "pval.exposure", "beta.exposure", "se.exposure",
                                   "samplesize.exposure", "effective_sample.exposure", "studies.exposure", "anno.exposure", "anno_full.exposure",
                                   "chr_pos_37", "info.exposure")

#We need to add some stuff:

rs16969968_exp_clean$id.exposure <- "SMK"
rs16969968_exp_clean$exposure <- "SMK"

rs1051730_exp_clean$id.exposure <- "SMK"
rs1051730_exp_clean$exposure <- "SMK"

rs1317286_exp_clean$id.exposure <- "SMK"
rs1317286_exp_clean$exposure <- "SMK"

##############################
#Worked like a freaking charm#
##############################

##############################################
#SECTION C: Merging exposure and outcome data#
##############################################

#We need to format the data for the outcome:

WCAdjBMI_outcome_rs1051730 <- WCAdjBMI_2_rs1051730

colnames(WCAdjBMI_outcome_rs1051730) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", 
                           "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                           "pval.outcome", "samplesize.outcome", "chr_pos_37")

WCAdjBMI_outcome_rs16969968 <- WCAdjBMI_2_rs16969968

colnames(WCAdjBMI_outcome_rs16969968) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", 
                                           "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                                           "pval.outcome", "samplesize.outcome", "chr_pos_37")

#And finally:

WCAdjBMI_outcome_rs1317286 <- WCAdjBMI_2_rs1317286

colnames(WCAdjBMI_outcome_rs1317286) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", 
                                          "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                                          "pval.outcome", "samplesize.outcome", "chr_pos_37")


#Getting ids, so troublesome...

WCAdjBMI_outcome_rs16969968$id.outcome <- "WCAdjBMI"
WCAdjBMI_outcome_rs16969968$outcome <- "WCAdjBMI"

WCAdjBMI_outcome_rs1051730$id.outcome <- "WCAdjBMI"
WCAdjBMI_outcome_rs1051730$outcome <- "WCAdjBMI"

WCAdjBMI_outcome_rs1317286$id.outcome <- "WCAdjBMI"
WCAdjBMI_outcome_rs1317286$outcome <- "WCAdjBMI"

dat_1_rs1051730 <- harmonise_data(rs1051730_exp_clean, WCAdjBMI_outcome_rs1051730, action = 3) 
dat_1_rs16969968 <- harmonise_data(rs16969968_exp_clean, WCAdjBMI_outcome_rs16969968, action = 3) 
dat_1_rs1317286 <- harmonise_data(rs1317286_exp_clean, WCAdjBMI_outcome_rs1317286, action = 3) 

#######################
#SECTION D: mF and IVW#
#######################

F_ = ((dat_1_rs1051730$beta.exposure)^2)/((dat_1_rs1051730$se.exposure)^2)
mF  = mean(F_)
print(mF)
#7.37

Isq(dat_1_rs1051730$beta.exposure, dat_1_rs1051730$se.exposure)
#NaN

##NEXT SNP:

F_ = ((dat_1_rs16969968$beta.exposure)^2)/((dat_1_rs16969968$se.exposure)^2)
mF  = mean(F_)
print(mF)
#6.45

Isq(dat_1_rs16969968$beta.exposure, dat_1_rs16969968$se.exposure)
#NaN

##NEXT SNP:

F_ = ((dat_1_rs1317286$beta.exposure)^2)/((dat_1_rs1317286$se.exposure)^2)
mF  = mean(F_)
print(mF)
#6.78

Isq(dat_1_rs1317286$beta.exposure, dat_1_rs1317286$se.exposure)
#NaN

################################
#Let's run the shit out of this#
################################

mr(dat_1_rs1051730)

#id.exposure id.outcome  outcome exposure     method nsnp          b        se      pval
#1         SMK   WCAdjBMI WCAdjBMI      SMK Wald ratio    1 -0.2240822 0.3507374 0.5228952

mr(dat_1_rs16969968)

#id.exposure id.outcome  outcome exposure     method nsnp          b        se      pval
#1         SMK   WCAdjBMI WCAdjBMI      SMK Wald ratio    1 -0.1977206 0.3642221 0.5872282

mr(dat_1_rs1317286)

#id.exposure id.outcome  outcome exposure     method nsnp          b        se     pval
#1         SMK   WCAdjBMI WCAdjBMI      SMK Wald ratio    1 -0.6113159 0.4482983 0.172682


######################################################
#I think that we should not go for this one after all#
######################################################