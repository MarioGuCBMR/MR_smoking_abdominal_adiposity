##############
#INTRODUCTION#
##############

#We are going to run CAUSE for SmkInit as exposure and WC as outcome

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

SMK <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/RAW_Data/SmokingInitiation.txt.gz")
WHR <- fread("C:/Users/zlc436/Downloads/whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")

########################################
#First we are gonna prune some SMK data#
########################################

head(SMK)

#We can filter through AF, INFO and STAT... somehow? I dunno how yet.

SMK$INFO <- SMK$EFFECTIVE_N/SMK$N

summary(SMK$INFO) #Min should be 0.30 according to readme.

SMK <- SMK[which(SMK$INFO > 0.70),]

summary(SMK$INFO) #perfect

#################
#Let's go for AF#
#################

summary(SMK$AF) #perfect

SMK <- SMK[which(SMK$AF > 0.01),]
SMK <- SMK[which(SMK$AF < 0.99),]

summary(SMK$AF) #perfect

##########################
#Let's check for the STAT#
##########################

summary(SMK$STAT) #Iit also reflects association. It has nothing to do with our data.

head(SMK)

###################################
#Finally let's get the chr and pos#
###################################

SMK$chr_pos <- paste(SMK$CHROM, SMK$POS, sep = ":")

###########################################
#Let's filter with the extended MHC region#
###########################################

SMK_mhc <- SMK[which(SMK$CHROM == 6 & SMK$POS >= 26000000 & SMK$POS <= 34000000),]

summary(SMK_mhc$CHROM) #perfect
summary(SMK_mhc$POS) #perfect

SMK_ <- SMK[which(!(SMK$chr_pos%in%SMK_mhc$chr_pos)),] #we go 7692565 to 7648667. PERFECT MATCH.

############################
#We check for weird alleles#
############################

yes_vect <- c("A", "G", "T", "C")

SMK_ <- SMK_[which(SMK_$REF%in%yes_vect),]
SMK_ <- SMK_[which(SMK_$ALT%in%yes_vect),]

#############
#WE ARE DONE#
#############

###############
#FILTERING WHR#
###############

WHR$RSID <- as.character(unlist(sapply(WHR$SNP, parse_rsid)))

index <- str_detect(WHR$SNP, ":")

which(index == FALSE)[1]

WHR[2146729,] #PERFECT.

#################
#ONLY GIANT data#
#################

#We are gonna work with GIANT data separately, 
#then merge the data.
#If there are any RSIDs that are duplicated, then I will erase them afterwards.

GIANT <- WHR[which(is.na(WHR$INFO) == TRUE),] #18028 that are erased.

#We are only gonna use the SNPs from GIANT that present allele frequencies in the range that I like.

summary(GIANT$Freq_Tested_Allele)

GIANT <- GIANT[which(GIANT$Freq_Tested_Alle < 0.99),]
GIANT <- GIANT[which(GIANT$Freq_Tested_Alle > 0.01),] #17810

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
#This data is the same as obtained with Lifetime since it only affects WHR.
#Hence, no problemo. We can just obtain it.

#write.table(nexus_query, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/snp_nexus/lifetime_WHR_4_nexus.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

#################
#Retrieving data#
#################

nexus_results <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/snp_nexus/lifetime_WHR_nexus_results.txt")

#Let's do a smol check:

nexus_results <- nexus_results[order(nexus_results$`ALT Allele (IUPAC)`),]

head(nexus_results) #GREAT.
tail(nexus_results) #GREAT.

#Let's get the data from here:

GIANT_nexus_recovered <- GIANT[which(GIANT$SNP%in%nexus_results$dbSNP),]
GIANT_nexus_missed <- GIANT[which(!(GIANT$SNP%in%nexus_results$dbSNP)),] #1446

GIANT_nexus_recovered <- GIANT_nexus_recovered[order(match(GIANT_nexus_recovered$SNP, nexus_results$dbSNP)),]

#Smol check:

length(which(GIANT_nexus_recovered$SNP == nexus_results$dbSNP)) #16364. perfect.

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

which(swapping_snps_1$new_rs%in%SMK_$RSID) #0
which(swapping_snps_1$old_rs%in%SMK_$RSID) #0

#None, so that means, that only those in GIANT that are in new 
#need to be there.
#Now we need to get the lifetime data that is in either the old or the new column 
#from those that are found in the GIANT data.

swapping_snps_2_in_lifetime <- swapping_snps_2[which(swapping_snps_2$new_rs%in%SMK_$RSID),] #135 do not need to be changed. Why? Because GIANT has already those. And lifetime too.

#HENCE, lifetime old are those that need to be changed!!

swapping_snps_2_in_SMK_old <- swapping_snps_2[which(swapping_snps_2$old_rs%in%SMK_$SNP),] #0

#So, let's go:

GIANT_new_2_old <- GIANT_nexus_missed[which(GIANT_nexus_missed$RSID%in%swapping_snps_2_in_SMK_old$new_rs),]

#Now let's change the RSID:

swapping_snps_2_in_SMK_old <- swapping_snps_2_in_SMK_old[order(match(swapping_snps_2_in_SMK_old$new_rs, GIANT_new_2_old$RSID)),]

head(swapping_snps_2_in_SMK_old)
head(GIANT_new_2_old)

#PERFECT!!

GIANT_new_2_old$RSID <- swapping_snps_2_in_SMK_old$old_rs

#####################
#Finally let's merge#
#####################

GIANT_rest <- GIANT_nexus_missed[which(!(GIANT_nexus_missed$RSID%in%swapping_snps_2_in_SMK_old$new_rs)),]

GIANT_missing_end <- rbind(GIANT_rest, GIANT_new_2_old)

GIANT_new_column <- GIANT_nexus_missed[which(GIANT_nexus_missed$RSID%in%swapping_snps_2$new_rs),] #1037
GIANT_old_column <- GIANT_nexus_missed[which(GIANT_nexus_missed$RSID%in%swapping_snps_2$old_rs),] #0 as expected.

###############
#Filtering WHR#
###############

#INFO:

WHR_ <- WHR[which(WHR$INFO > 0.7),]

summary(WHR_$INFO)

#MAF:

summary(WHR_$Freq_Tested_Allele)

WHR_[which.min(WHR_$Freq_Tested_Allele),]

WHR_ <- WHR_[which(as.numeric(WHR_$Freq_Tested_Allele) > 0.01),]
WHR_ <- WHR_[which(as.numeric(WHR_$Freq_Tested_Allele) < 0.99),]

##################
#Finally, alleles#
##################

yes_vect <- c("A", "G", "T", "C")

WHR_ <- WHR_[which(WHR_$Tested_Allele%in%yes_vect),] #all of them, as expected.
WHR_ <- WHR_[which(WHR_$Other_Allele%in%yes_vect),] #all of them, as expected.

##################################################
#We are not gonna do MHC because we did it before#
##################################################

WHR_$chr_pos <- paste(WHR_$CHR, WHR_$POS, sep = ":")

##################
#MERGING DATASETS#
##################

library(cause)
library(tidyverse)

X_UKBB <- gwas_merge(SMK_, WHR_, snp_name_cols = c("chr_pos", "chr_pos"), 
                     beta_hat_cols = c("BETA", "BETA"), 
                     se_cols = c("SE", "SE"), 
                     A1_cols = c("ALT", "Tested_Allele"), 
                     A2_cols = c("REF", "Other_Allele"))

X_GIANT_recovered <- gwas_merge(SMK_, GIANT_nexus_recovered, snp_name_cols = c("chr_pos", "chr_pos"), 
                                beta_hat_cols = c("BETA", "BETA"), 
                                se_cols = c("SE", "SE"), 
                                A1_cols = c("ALT", "Tested_Allele"), 
                                A2_cols = c("REF", "Other_Allele"))

X_GIANT_missed <- gwas_merge(SMK_, GIANT_missing_end, snp_name_cols = c("RSID", "RSID"), 
                             beta_hat_cols = c("BETA", "BETA"), 
                             se_cols = c("SE", "SE"), 
                             A1_cols = c("ALT", "Tested_Allele"), 
                             A2_cols = c("REF", "Other_Allele"))

#######################################################################
#We need to understand why when we add more SNPs, we wind up with less#
#######################################################################

#Calculating Nuisance:

###############################
#We need to retrieve the RSIDs#
###############################

head(X_UKBB)

WHR_[which(WHR_$chr_pos == "1:693731"),] #PERFECT.
SMK_[which(SMK_$chr_pos == "1:693731"),] #PERFECT

WHR_match <- WHR_[which(WHR_$chr_pos%in%X_UKBB$snp),]

WHR_match <- WHR_match[order(match(WHR_match$chr_pos, X_UKBB$snp)),]

which(WHR_match$chr_pos != X_UKBB$snp) #they all match. 

head(WHR_match$chr_pos)
head(X_UKBB$snp)
tail(WHR_match$chr_pos)
tail(X_UKBB$snp)

##########
#All good#
##########

X_UKBB$snp <- WHR_match$RSID

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
SMK_[which(SMK_$chr_pos == "1:2488608"),] #PERFECT

GIANT_nexus_recovered_match <- GIANT_nexus_recovered[which(GIANT_nexus_recovered$chr_pos%in%X_GIANT_recovered$snp),]

GIANT_nexus_recovered_match <- GIANT_nexus_recovered_match[order(match(GIANT_nexus_recovered_match$chr_pos, X_GIANT_recovered$snp)),]

which(GIANT_nexus_recovered_match$chr_pos != X_GIANT_recovered$snp) #they all match. 

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

X_GIANT_in_UKBB <- X_GIANT[which(X_GIANT$snp%in%X_UKBB$snp),] #258. SAME.

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

#Let's check for duplicates, I don't think there should be any left after this.

X_end[which(duplicated(X_end$snp) == TRUE),] #none of them, of course!

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(params, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Smoking_Initiation/dataframes/params_SMK_WHR_JULY.rds")
saveRDS(X_end, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Smoking_Initiation/dataframes/X_SMK_WHR_JULY.rds")

params <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Smoking_Initiation/dataframes/params_SMK_WHR_JULY.rds")
X_end <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Smoking_Initiation/dataframes/X_SMK_WHR_JULY.rds")

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

SMK_pruned <- SMK_[which(SMK_$RSID%in%pruned_all),] #they are all here, no need of weird stuff.

write.table(SMK_pruned, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Smoking_Initiation/dataframes/SMK_WHR_PRUNED_JULY.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

#############
#Small check#
#############

#When the DF are similar, the SNPs should be the same, as we know.
#Indeed, this might be the case between WHR and WHRAdjBMI

just_check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Smoking_Initiation/dataframes/SMK_WHRAdjBMI_PRUNED.txt") #exactly the same

length(which(just_check$SNP%in%SMK_pruned$SNP)) #all of them.

#Then, adding those SNPs didn't change a thing.
#As expected.

#################
#Let's run CAUSE#
#################

top_SMK_pruned_df <- read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Smoking_Initiation/dataframes/SMK_WHR_PRUNED_JULY.txt", header = TRUE, stringsAsFactors = FALSE)

top_SMK_pruned_snps <- top_SMK_pruned_df$RSID

#Full version:

res <- cause(X=X_end, variants = top_SMK_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Smoking_Initiation/dataframes/RES_SMK_to_WHR_FULL_JULY.rds")

res <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Smoking_Initiation/dataframes/RES_SMK_to_WHR_FULL_JULY.rds")

res$elpd

#model1  model2  delta_elpd se_delta_elpd         z
#1    null sharing -111.826854    11.9086681 -9.390374
#2    null  causal -113.198683    12.0692316 -9.379113
#3 sharing  causal   -1.371829     0.1747279 -7.851230

#Sounds really promising:

summary(res)

#p-value testing that causal model is a better fit:  2.1e-15 
#Posterior medians and  95 % credible intervals:
#  model     gamma              eta                 q                  
#[1,] "Sharing" NA                 "0.13 (0.11, 0.15)" "0.96 (0.88, 0.99)"
#[2,] "Causal"  "0.13 (0.1, 0.16)" "0 (-0.15, 0.16)"   "0.18 (0, 0.86)"  

plot(res, type="data")


