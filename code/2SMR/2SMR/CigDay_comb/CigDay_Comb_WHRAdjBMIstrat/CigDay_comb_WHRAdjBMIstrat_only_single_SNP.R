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

CIG <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/RAW_Data/CigarettesPerDay.txt.gz")
CIG$chr_pos_37 <- paste("chr", CIG$CHROM, ":", CIG$POS, sep = "")
CIG$INFO <- CIG$EFFECTIVE_N/CIG$N

cig_gw <- readxl::read_xlsx("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/RAW_Data/CigDay_conditionally_independent_snps.xlsx")

colnames(cig_gw) <- cig_gw[1,]

cig_gw <- cig_gw[-1,]

#########################################
#2. Get the genome-wide significant SNPs#
#########################################

rs1051730_exp <- CIG[which(CIG$RSID == "rs1051730"),]
rs16969968_exp <- CIG[which(CIG$RSID == "rs16969968"),]
rs1317286_exp <- CIG[which(CIG$RSID == "rs1317286"),]

#Can we get them with the higher sample size:

check_rs1051730_exp <- cig_gw[which(cig_gw$rsID == "rs1051730"),] #not there
check_rs16969968_exp <- cig_gw[which(cig_gw$rsID == "rs16969968"),] #not there
check_rs1317286_exp <- cig_gw[which(cig_gw$rsID == "rs1317286"),]  #not there

########################################################################
#Perfect, so we just run the analysis with the SmokingInitiation.txt.gz#
########################################################################

#Checking pvals:

rs1051730_exp

#CHROM      POS      RSID REF ALT    AF STAT    PVALUE      BETA          SE      N EFFECTIVE_N
#1:    15 78894339 rs1051730   G   A 0.324  921 2.33e-202 0.1795169 0.005915282 257341      253807
#Number_of_Studies              ANNO
#1:                33 Synonymous:CHRNA3
#ANNOFULL
#1: CHRNA3/NM_001166694:-:Synonymous(TAC/Tyr/Y->TAT/Tyr/Y:Base646/1470:Codon216/490:Exon5/6):Exon|CHRNA3/NM_000743:-:Synonymous(TAC/Tyr/Y->TAT/Tyr/Y:Base646/1518:Codon216/506:Exon5/6):Exon
#chr_pos_37      INFO
#1: chr15:78894339 0.9862672

#The p-value is not strong at all.

rs16969968_exp

#CHROM      POS       RSID REF ALT    AF STAT    PVALUE      BETA          SE      N EFFECTIVE_N
#1:    15 78882925 rs16969968   G   A 0.322  912 2.32e-200 0.1787725 0.005919748 257341      254984
#Number_of_Studies                 ANNO
#1:                33 Nonsynonymous:CHRNA5
#ANNOFULL
#1: CHRNA5/NM_000745:+:Nonsynonymous(GAT/Asp/D->AAT/Asn/N:Base1193/1407:Codon398/469:Exon5/6):Exon
#chr_pos_37      INFO
#1: chr15:78882925 0.9908409

rs1317286_exp

#CHROM      POS      RSID REF ALT    AF STAT   PVALUE      BETA          SE      N EFFECTIVE_N
#1:    15 78896129 rs1317286   A   G 0.332  926 2.6e-203 0.1792113 0.005889249 257341      253212
#Number_of_Studies          ANNO                                               ANNOFULL     chr_pos_37
#1:                33 Intron:CHRNA3 CHRNA3/NM_001166694:-:Intron|CHRNA3/NM_000743:-:Intron chr15:78896129
#INFO
#1: 0.9839551

#################################################
#SECTION B: FINDING matching SNPs in the outcome#
#################################################

#First let's load the outcome:

WHRAdjBMI <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMISmk/WHRAdjBMI_Smk_Strat_Curated.txt")

##########################################################################################
#1. we are gonna try to match the data with RSID and with chr_pos to get as many variants#
##########################################################################################

WHRAdjBMI_2_rs1051730 <- WHRAdjBMI[which(WHRAdjBMI$rs_id%in%rs1051730_exp$RSID),] #gotcha
WHRAdjBMI_2_rs16969968 <- WHRAdjBMI[which(WHRAdjBMI$rs_id%in%rs16969968_exp$RSID),] #gotcha
WHRAdjBMI_2_rs1317286 <- WHRAdjBMI[which(WHRAdjBMI$rs_id%in%rs1317286_exp$RSID),] #gotcha

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

rs16969968_exp_clean$id.exposure <- "CIG"
rs16969968_exp_clean$exposure <- "CIG"

rs1051730_exp_clean$id.exposure <- "CIG"
rs1051730_exp_clean$exposure <- "CIG"

rs1317286_exp_clean$id.exposure <- "CIG"
rs1317286_exp_clean$exposure <- "CIG"

##############################
#Worked like a freaking charm#
##############################

##############################################
#SECTION C: Merging exposure and outcome data#
##############################################

#We need to format the data for the outcome:

WHRAdjBMI_outcome_rs1051730 <- WHRAdjBMI_2_rs1051730

colnames(WHRAdjBMI_outcome_rs1051730) <- c("chr.outcome", "SNP", "marker", "pos_18", "effect_allele.outcome", 
                           "other_allele.outcome", "eaf.outcome", "samplesize.outcome", "beta.outcome", "se.outcome",
                           "pval.outcome", "samplesize_NO.outcome",  "beta_NO.outcome", "se_NO.outcome",
                           "pval_NO.outcome","chr_pos_37", "chr_pos_18")

WHRAdjBMI_outcome_rs16969968 <- WHRAdjBMI_2_rs16969968

colnames(WHRAdjBMI_outcome_rs16969968) <-c("chr.outcome", "SNP", "marker", "pos_18", "effect_allele.outcome", 
                                           "other_allele.outcome", "eaf.outcome", "samplesize.outcome", "beta.outcome", "se.outcome",
                                           "pval.outcome", "samplesize_NO.outcome",  "beta_NO.outcome", "se_NO.outcome",
                                           "pval_NO.outcome","chr_pos_37", "chr_pos_18")

#And finally:

WHRAdjBMI_outcome_rs1317286 <- WHRAdjBMI_2_rs1317286

colnames(WHRAdjBMI_outcome_rs1317286) <- c("chr.outcome", "SNP", "marker", "pos_18", "effect_allele.outcome", 
                                           "other_allele.outcome", "eaf.outcome", "samplesize.outcome", "beta.outcome", "se.outcome",
                                           "pval.outcome", "samplesize_NO.outcome",  "beta_NO.outcome", "se_NO.outcome",
                                           "pval_NO.outcome","chr_pos_37", "chr_pos_18")

#Getting ids, so troublesome...

WHRAdjBMI_outcome_rs16969968$id.outcome <- "WHRAdjBMI"
WHRAdjBMI_outcome_rs16969968$outcome <- "WHRAdjBMI"

WHRAdjBMI_outcome_rs1051730$id.outcome <- "WHRAdjBMI"
WHRAdjBMI_outcome_rs1051730$outcome <- "WHRAdjBMI"

WHRAdjBMI_outcome_rs1317286$id.outcome <- "WHRAdjBMI"
WHRAdjBMI_outcome_rs1317286$outcome <- "WHRAdjBMI"

dat_1_rs1051730 <- harmonise_data(rs1051730_exp_clean, WHRAdjBMI_outcome_rs1051730, action = 3) 
dat_1_rs16969968 <- harmonise_data(rs16969968_exp_clean, WHRAdjBMI_outcome_rs16969968, action = 3) 
dat_1_rs1317286 <- harmonise_data(rs1317286_exp_clean, WHRAdjBMI_outcome_rs1317286, action = 3) 

#######################
#SECTION D: mF and IVW#
#######################

F_ = ((dat_1_rs1051730$beta.exposure)^2)/((dat_1_rs1051730$se.exposure)^2)
mF  = mean(F_)
print(mF)
#921

Isq(dat_1_rs1051730$beta.exposure, dat_1_rs1051730$se.exposure)
#NaN

##NEXT SNP:

F_ = ((dat_1_rs16969968$beta.exposure)^2)/((dat_1_rs16969968$se.exposure)^2)
mF  = mean(F_)
print(mF)
#912

Isq(dat_1_rs16969968$beta.exposure, dat_1_rs16969968$se.exposure)
#NaN

##NEXT SNP:

F_ = ((dat_1_rs1317286$beta.exposure)^2)/((dat_1_rs1317286$se.exposure)^2)
mF  = mean(F_)
print(mF)
#925,99

Isq(dat_1_rs1317286$beta.exposure, dat_1_rs1317286$se.exposure)
#NaN

################################
#Let's run the shit out of this#
################################

mr(dat_1_rs1051730)

#id.exposure id.outcome   outcome exposure     method nsnp           b         se      pval
#1         CIG  WHRAdjBMI WHRAdjBMI      CIG Wald ratio    1 0.007241659 0.02228203 0.7451811

mr(dat_1_rs16969968)
#id.exposure id.outcome   outcome exposure     method nsnp           b         se      pval
#1         CIG  WHRAdjBMI WHRAdjBMI      CIG Wald ratio    1 0.009509294 0.02237481 0.6708367

mr(dat_1_rs1317286)
#id.exposure id.outcome   outcome exposure     method nsnp         b         se       pval
#1         CIG  WHRAdjBMI WHRAdjBMI      CIG Wald ratio    1 0.1009981 0.04966204 0.04198106


######################################################
#I think that we should not go for this one after all#
######################################################