##############
#INTRODUCTION#
##############

#We are going to run CAUSE for Lifetime as exposure and WHR as outcome

memory.limit(size=8000000)

#########
#Library#
#########

library(data.table)
library(tidyverse)
library(cause)
library(TwoSampleMR)

###########
#Functions#
###########

parse_rsid <- function(rsid){
  
  rsid_ <- strsplit(rsid, ":")[[1]][1] #this works for those that do not present the weird formatting :)
  
  return(rsid_)
  
}

##############
#Loading data#
##############

Lifetime <- fread("C:/Users/zlc436/Desktop/SMK_WHRAdjBMI/SmokingCombined_Wootton.txt")
WHRAdjBMI <- fread("C:/Users/zlc436/Downloads/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")

#############################################
#First we are gonna prune some Lifetime data#
#############################################

head(Lifetime)

summary(Lifetime$INFO) #Min here is 0.8. Min should be 0.70, so this is perfect..

#################
#Let's go for AF#
#################

summary(Lifetime$EAF) #already done!!

Lifetime <- Lifetime[which(Lifetime$EAF > 0.01),] #just in case
Lifetime <- Lifetime[which(Lifetime$EAF < 0.99),] #just in case

#Perfect!!

###################################
#Finally let's get the chr and pos#
###################################

Lifetime$chr_pos <- paste(Lifetime$CHR, Lifetime$BP, sep = ":")

###########################################
#Let's filter with the extended MHC region#
###########################################

Lifetime_mhc <- Lifetime[which(Lifetime$CHR == 6 & Lifetime$BP >= 26000000 & Lifetime$BP <= 34000000),]

summary(Lifetime_mhc$CHR) #perfect
summary(Lifetime_mhc$BP) #perfect

Lifetime_ <- Lifetime[which(!(Lifetime$chr_pos%in%Lifetime_mhc$chr_pos)),] #we go 7683352 to 7641395. PERFECT MATCH.

############################
#We check for weird alleles#
############################

yes_vect <- c("A", "G", "T", "C")

Lifetime_ <- Lifetime_[which(Lifetime_$EFFECT_ALLELE%in%yes_vect),]
Lifetime_ <- Lifetime_[which(Lifetime_$OTHER_ALLELE%in%yes_vect),]

#############
#WE ARE DONE#
#############

#####################
#FILTERING WHRAdjBMI#
#####################

WHRAdjBMI$RSID <- as.character(unlist(sapply(WHRAdjBMI$SNP, parse_rsid)))

index <- str_detect(WHRAdjBMI$SNP, ":")

which(index == FALSE)[1]

WHRAdjBMI[2146729,] #PERFECT.

#################
#ONLY GIANT data#
#################

#We are gonna work with GIANT data separately, 
#then merge the data.
#If there are any RSIDs that are duplicated, then I will erase them afterwards.

GIANT <- WHRAdjBMI[which(is.na(WHRAdjBMI$INFO) == TRUE),] #17391 that are erased.

#We are only gonna use the SNPs from GIANT that present allele frequencies in the range that I like.

GIANT <- GIANT[which(GIANT$Freq_Tested_Alle < 0.99),]
GIANT <- GIANT[which(GIANT$Freq_Tested_Alle > 0.01),] #17190

###########################################################
#We are not going to trust the CHR_POS from these bad boys#
###########################################################

#The alleles match with the RSID, but the chr_pos do not. 
#We are going to make one last check:

head(GIANT) #the first one is already a 100% example. I have seen like 10 of them, so we can be sure.

#Thus we are going to get as many chr_pos in build 37 with SNPNexus. 

GIANT$CHR <- NA
GIANT$POS <- NA
GIANT$dbsnp <- "dbsnp"

nexus_query <- GIANT %>%
  select(dbsnp, SNP)

#We got our data!!
#Now let's save it:

write.table(nexus_query, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/snp_nexus/lifetime_WHRAdjBMI_4_nexus.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

#checkio <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/snp_nexus/lifetime_WHRAdjBMI_4_nexus.txt", header = FALSE)

#################
#Retrieving data#
#################

nexus_results <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/snp_nexus/lifetime_WHRAdjBMI_nexus_results.txt")

#Let's do a smol check:

nexus_results <- nexus_results[order(nexus_results$`ALT Allele (IUPAC)`),]

head(nexus_results) #GREAT.
tail(nexus_results) #GREAT.

#Let's get the data from here:

GIANT_nexus_recovered <- GIANT[which(GIANT$SNP%in%nexus_results$dbSNP),]
GIANT_nexus_missed <- GIANT[which(!(GIANT$SNP%in%nexus_results$dbSNP)),] #1395

GIANT_nexus_recovered <- GIANT_nexus_recovered[order(match(GIANT_nexus_recovered$SNP, nexus_results$dbSNP)),]

#Smol check:

length(which(GIANT_nexus_recovered$SNP == nexus_results$dbSNP)) #15795. perfect.

nexus_results$chr_pos <- paste(nexus_results$Chromosome, nexus_results$Position, sep = ":")

GIANT_nexus_recovered$chr_pos <- nexus_results$chr_pos

#########
#PERFECT#
#########

#Let's check the merged allele versions of those missing and that only have RSIDs:

merged_rs <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Romain_Mario/RAW_DATA/MergedRs/RsMergeArch.bcp.gz")

merged_rs$old_rs <- paste("rs", merged_rs$V1, sep = "")
merged_rs$new_rs <- paste("rs", merged_rs$V2, sep = "")

###################
#Now we match them#
###################

index_merged_1 <- which(merged_rs$old_rs%in%GIANT_nexus_missed$RSID) #Only 3! 
index_merged_2 <- which(merged_rs$new_rs%in%GIANT_nexus_missed$RSID) #Most of them. There are a few repeated quite some times. 

#First, which of the old ones are in lifetime?

swapping_snps_1 <- merged_rs[index_merged_1,]
swapping_snps_2 <- merged_rs[index_merged_2,]

which(swapping_snps_1$new_rs%in%Lifetime_$SNP) #0
which(swapping_snps_1$old_rs%in%Lifetime_$SNP) #0

#None, so that means, that only those in GIANT that are in new 
#need to be there.
#Now we need to get the lifetime data that is in either the old or the new column 
#from those that are found in the GIANT data.

swapping_snps_2_in_lifetime <- swapping_snps_2[which(swapping_snps_2$new_rs%in%Lifetime_$SNP),] #109 do not need to be changed. Why? Because GIANT has already those. And lifetime too.

#HENCE, lifetime old are those that need to be changed!!

swapping_snps_2_in_lifetime_old <- swapping_snps_2[which(swapping_snps_2$old_rs%in%Lifetime_$SNP),] #0

#So, let's go:

GIANT_new_2_old <- GIANT_nexus_missed[which(GIANT_nexus_missed$RSID%in%swapping_snps_2_in_lifetime_old$new_rs),]

#Now let's change the RSID:

swapping_snps_2_in_lifetime_old <- swapping_snps_2_in_lifetime_old[order(match(swapping_snps_2_in_lifetime_old$new_rs, GIANT_new_2_old$RSID)),]

head(swapping_snps_2_in_lifetime_old)
head(GIANT_new_2_old)

#PERFECT!!

GIANT_new_2_old$RSID <- swapping_snps_2_in_lifetime_old$old_rs

#####################
#Finally let's merge#
#####################

GIANT_rest <- GIANT_nexus_missed[which(!(GIANT_nexus_missed$RSID%in%swapping_snps_2_in_lifetime_old$new_rs)),]

GIANT_missing_end <- rbind(GIANT_rest, GIANT_new_2_old)

GIANT_new_column <- GIANT_nexus_missed[which(GIANT_nexus_missed$RSID%in%swapping_snps_2$new_rs),] #1010
GIANT_old_column <- GIANT_nexus_missed[which(GIANT_nexus_missed$RSID%in%swapping_snps_2$old_rs),] #0 as expected.

#####################
#Filtering WHRAdjBMI#
#####################

#INFO:

WHRAdjBMI_ <- WHRAdjBMI[which(WHRAdjBMI$INFO > 0.7),]

summary(WHRAdjBMI_$INFO)

#MAF:

summary(WHRAdjBMI_$Freq_Tested_Allele)

WHRAdjBMI_[which.min(WHRAdjBMI_$Freq_Tested_Allele),]

WHRAdjBMI_ <- WHRAdjBMI_[which(as.numeric(WHRAdjBMI_$Freq_Tested_Allele) > 0.01),]
WHRAdjBMI_ <- WHRAdjBMI_[which(as.numeric(WHRAdjBMI_$Freq_Tested_Allele) < 0.99),]

##################
#Finally, alleles#
##################

yes_vect <- c("A", "G", "T", "C")

WHRAdjBMI_ <- WHRAdjBMI_[which(WHRAdjBMI_$Tested_Allele%in%yes_vect),] #all of them, as expected.
WHRAdjBMI_ <- WHRAdjBMI_[which(WHRAdjBMI_$Other_Allele%in%yes_vect),] #all of them, as expected.

##################################################
#We are not gonna do MHC because we did it before#
##################################################

WHRAdjBMI_$chr_pos <- paste(WHRAdjBMI_$CHR, WHRAdjBMI_$POS, sep = ":")

##################
#MERGING DATASETS#
##################

library(cause)
library(tidyverse)

X_UKBB <- gwas_merge(Lifetime_, WHRAdjBMI_, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("BETA", "BETA"), 
                se_cols = c("SE", "SE"), 
                A1_cols = c("EFFECT_ALLELE", "Tested_Allele"), 
                A2_cols = c("OTHER_ALLELE", "Other_Allele"))

X_GIANT_recovered <- gwas_merge(Lifetime_, GIANT_nexus_recovered, snp_name_cols = c("chr_pos", "chr_pos"), 
                                beta_hat_cols = c("BETA", "BETA"), 
                                se_cols = c("SE", "SE"), 
                                A1_cols = c("EFFECT_ALLELE", "Tested_Allele"), 
                                A2_cols = c("OTHER_ALLELE", "Other_Allele"))

X_GIANT_missed <- gwas_merge(Lifetime_, GIANT_missing_end, snp_name_cols = c("SNP", "RSID"), 
                             beta_hat_cols = c("BETA", "BETA"), 
                             se_cols = c("SE", "SE"), 
                             A1_cols = c("EFFECT_ALLELE", "Tested_Allele"), 
                             A2_cols = c("OTHER_ALLELE", "Other_Allele"))

#######################################################################
#We need to understand why when we add more SNPs, we wind up with less#
#######################################################################

#Calculating Nuisance:

###############################
#We need to retrieve the RSIDs#
###############################

head(X_UKBB)

WHRAdjBMI_[which(WHRAdjBMI_$chr_pos == "1:693731"),] #PERFECT.
Lifetime_[which(Lifetime_$chr_pos == "1:693731"),] #PERFECT

WHRAdjBMI_match <- WHRAdjBMI_[which(WHRAdjBMI_$chr_pos%in%X_UKBB$snp),]

WHRAdjBMI_match <- WHRAdjBMI_match[order(match(WHRAdjBMI_match$chr_pos, X_UKBB$snp)),]

length(which(WHRAdjBMI_match$chr_pos != X_UKBB$snp)) #they all match. 
length(which(WHRAdjBMI_match$chr_pos == X_UKBB$snp)) #they all match. 

head(WHRAdjBMI_match$chr_pos)
head(X_UKBB$snp)
tail(WHRAdjBMI_match$chr_pos)
tail(X_UKBB$snp)

##########
#All good#
##########

X_UKBB$snp <- WHRAdjBMI_match$RSID

head(X_UKBB)

length(which(duplicated(X_UKBB$snp) == TRUE))

#######################################
#MERGING GIANT data recovered by dbsnp#
#######################################

###############################
#We need to retrieve the RSIDs#
###############################

head(X_GIANT_recovered)

GIANT_nexus_recovered[which(GIANT_nexus_recovered$chr_pos == "1:2488608"),] #PERFECT.
Lifetime_[which(Lifetime_$chr_pos == "1:2488608"),] #PERFECT

GIANT_nexus_recovered_match <- GIANT_nexus_recovered[which(GIANT_nexus_recovered$chr_pos%in%X_GIANT_recovered$snp),]

GIANT_nexus_recovered_match <- GIANT_nexus_recovered_match[order(match(GIANT_nexus_recovered_match$chr_pos, X_GIANT_recovered$snp)),]

length(which(GIANT_nexus_recovered_match$chr_pos != X_GIANT_recovered$snp)) #they all match. 
length(which(GIANT_nexus_recovered_match$chr_pos == X_GIANT_recovered$snp)) #they all match. 

head(GIANT_nexus_recovered_match$chr_pos)
head(X_GIANT_recovered$snp)
tail(GIANT_nexus_recovered_match$chr_pos)
tail(X_GIANT_recovered$snp)

##########
#All good#
##########

X_GIANT_recovered$snp <- GIANT_nexus_recovered_match$RSID

head(X_GIANT_recovered)

length(which(duplicated(X_GIANT_recovered$snp) == TRUE))

#####################################
#NOW LET'S MERGE THE DATA FROM GIANT#
#####################################

X_GIANT <- rbind(X_GIANT_missed, X_GIANT_recovered)

head(X_GIANT)
tail(X_GIANT)

############################
#Now merging the GIANT data#
############################

#We are gonna divide the data into RSIDs alrady in UKBB and those that are not:

X_GIANT_in_UKBB <- X_GIANT[which(X_GIANT$snp%in%X_UKBB$snp),] #256. SAME.

#We are gonna remove the duplicated in both:

dupl <- X_GIANT_in_UKBB$snp

X_UKBB_check <- X_UKBB[which(X_UKBB$snp%in%X_GIANT_in_UKBB$snp),]

mixed_df <- rbind(X_GIANT_in_UKBB, X_UKBB_check)

mixed_df <- mixed_df[order(mixed_df$snp),]

head(mixed_df)

#They are all the same...
#they present the same alleles, but the effect sizes are different.
#Nonetheless, we should keep those duplicates that have the largest sample size right?
#I think that is fair. If we had triallelic SNPs, that would be another thing
#Since CAUSE removes them.

#However, now that is not the case.
#We are only going to include the SNPs that are not present.

X_GIANT_end <- X_GIANT[which(!(X_GIANT$snp%in%dupl)),]

#And now we merged them_

X_end <- rbind(X_UKBB, X_GIANT_end)

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(params, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/dataframes/params_Lifetime_WHRAdjBMI_JULY.rds")
saveRDS(X_end, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/dataframes/X_Lifetime_WHRAdjBMI_JULY.rds")

params <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/dataframes/params_Lifetime_WHRAdjBMI_JULY.rds")
X_end <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/dataframes/X_Lifetime_WHRAdjBMI_JULY.rds")

#There might be duplicates:

X_end[which(duplicated(X_end$rsid) == TRUE),]

############
#LD pruning#
############

#params has been calculated, so let's go for LD pruning if possible:

library(tidyverse)

variants <- X_end %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))

#1

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr1_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr1_AF0.05_snpdata.RDS")

pruned_1 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#2

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr2_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr2_AF0.05_snpdata.RDS")

pruned_2 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#3

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr3_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr3_AF0.05_snpdata.RDS")

pruned_3 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#4

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr4_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr4_AF0.05_snpdata.RDS")

pruned_4 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#5

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr5_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr5_AF0.05_snpdata.RDS")

pruned_5 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#6

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr6_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr6_AF0.05_snpdata.RDS")

pruned_6 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#7

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr7_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr7_AF0.05_snpdata.RDS")

pruned_7 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#8

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr8_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr8_AF0.05_snpdata.RDS")

pruned_8 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#9

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr9_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr9_AF0.05_snpdata.RDS")

pruned_9 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#10

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr10_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr10_AF0.05_snpdata.RDS")

pruned_10 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#11

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr11_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr11_AF0.05_snpdata.RDS")

pruned_11 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#12

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr12_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr12_AF0.05_snpdata.RDS")

pruned_12 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#13

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr13_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr13_AF0.05_snpdata.RDS")

pruned_13 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#14

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr14_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr14_AF0.05_snpdata.RDS")

pruned_14 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#15

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr15_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr15_AF0.05_snpdata.RDS")

pruned_15 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#16

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr16_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr16_AF0.05_snpdata.RDS")

pruned_16 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))


#17

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr17_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr17_AF0.05_snpdata.RDS")

pruned_17 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#18

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr18_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr18_AF0.05_snpdata.RDS")

pruned_18 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#19

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr19_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr19_AF0.05_snpdata.RDS")

pruned_19 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))


#20

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr20_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr20_AF0.05_snpdata.RDS")

pruned_20 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#21

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr21_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr21_AF0.05_snpdata.RDS")

pruned_21 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))
#22

#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_0.1.RDS?download=1", destfile = "C:/Users/zlc436/Desktop/Leisure_Project/chr22_AF0.05_0.1.RDS")
#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_snpdata.RDS?download=1", destfile="C:/Users/zlc436/Desktop/Leisure_Project/chr22_AF0.05_snpdata.RDS")

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr22_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr22_AF0.05_snpdata.RDS")

#variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
pruned_22 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

pruned_all <- c(pruned_1, pruned_2, pruned_3, pruned_4, pruned_5, pruned_6, pruned_7, pruned_8, pruned_9, pruned_10, pruned_11, pruned_12, pruned_13, pruned_14, pruned_15, pruned_16, pruned_17, pruned_18, pruned_19, pruned_20, pruned_21, pruned_22)

Lifetime_pruned <- Lifetime_[which(Lifetime_$SNP%in%pruned_all),] #they are all here, no need of weird stuff.

write.table(Lifetime_pruned, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/dataframes/Lifetime_WHRAdjBMI_PRUNED_JULY.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

#############
#Small check#
#############

#When the DF are similar, the SNPs should be the same, as we know.
#Indeed, this might be the case between WHR and WHRAdjBMI

just_check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/dataframes/Lifetime_WHRAdjBMI_PRUNED.txt") #exactly the same

length(which(just_check$SNP%in%Lifetime_pruned$SNP)) #all of them.

#Another check:

lifetime_chr_3 <- Lifetime_[which(Lifetime$SNP%in%pruned_3),]

summary(lifetime_chr_3$CHR) #perfect

summary(lifetime_chr_3$P) #perfect.

#Then, adding those SNPs didn't change a thing.
#As expected.

#################
#Let's run CAUSE#
#################

top_Lifetime_pruned_df <- read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/dataframes/Lifetime_WHRAdjBMI_PRUNED_JULY.txt", header = TRUE, stringsAsFactors = FALSE)

top_Lifetime_pruned_snps <- top_Lifetime_pruned_df$SNP

#Full version:

res <- cause(X=X_end, variants = top_Lifetime_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/dataframes/RES_Lifetime_to_WHRAdjBMI_FULL_JULY.rds")

res <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/dataframes/RES_Lifetime_to_WHRAdjBMI_FULL_JULY.rds")

res$elpd

#model1  model2 delta_elpd se_delta_elpd         z
#1    null sharing -85.623210    11.0711000 -7.733939
#2    null  causal -86.883197    11.2659132 -7.712042
#3 sharing  causal  -1.259987     0.2146005 -5.871314

#Sounds really promising:

summary(res, ci_size = 0.95)

#p-value testing that causal model is a better fit:  2.2e-09 
#Posterior medians and  95 % credible intervals:
#  model     gamma               eta                 q                  
#[1,] "Sharing" NA                  "0.19 (0.16, 0.23)" "0.94 (0.81, 0.99)"
#[2,] "Causal"  "0.18 (0.13, 0.24)" "0 (-0.27, 0.24)"   "0.18 (0, 0.86)"   

plot(res)

plot(res, type="data")


