##############
#INTRODUCTION#
##############

#We are going to run CAUSE for PacksPerYear as exposure and WHRadjBMI as outcome

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

PacksPerYear <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/PacksPerYear/PackPerYear_Curated_FULL.txt")
WHRadjBMI <- fread("C:/Users/zlc436/Desktop/PULIT_data_reservoir/Combined/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")

#################################################
#First we are gonna prune some PacksPerYear data#
#################################################

#First we need to standardise the data:

head(PacksPerYear)

#There is no INFO data.

#################
#Let's go for AF#
#################

summary(PacksPerYear$eaf.exposure) #already done!!

###################################
#Finally let's get the chr and pos#
###################################

PacksPerYear$chr_pos <- paste(PacksPerYear$chr.exposure, PacksPerYear$pos.exposure, sep = ":")

###########################################
#Let's filter with the extended MHC region#
###########################################

PacksPerYear_mhc <- PacksPerYear[which(PacksPerYear$chr.exposure == 6 & PacksPerYear$pos.exposure >= 26000000 & PacksPerYear$pos.exposure <= 34000000),] #none: they were already remove

PacksPerYear_ <- PacksPerYear

############################
#We check for weird alleles#
############################

yes_vect <- c("A", "G", "T", "C")

PacksPerYear_ <- PacksPerYear_[which(PacksPerYear_$effect_allele.exposure%in%yes_vect),]
PacksPerYear_ <- PacksPerYear_[which(PacksPerYear_$other_allele.exposure%in%yes_vect),]

#############
#WE ARE DONE#
#############

#####################
#FILTERING WHRadjBMI#
#####################

WHRadjBMI$RSID <- as.character(unlist(sapply(WHRadjBMI$SNP, parse_rsid)))

index <- str_detect(WHRadjBMI$SNP, ":")

which(index == FALSE)[1]

WHRadjBMI[2146729,] #PERFECT.

#################
#ONLY GIANT data#
#################

#We are gonna work with GIANT data separately, 
#then merge the data.
#If there are any RSIDs that are duplicated, then I will erase them afterwards.

GIANT <- WHRadjBMI[which(is.na(WHRadjBMI$INFO) == TRUE),] #18028 that are erased.

#We are only gonna use the SNPs from GIANT that present allele frequencies in the range that I like.

GIANT <- GIANT[which(GIANT$Freq_Tested_Alle < 0.99),]
GIANT <- GIANT[which(GIANT$Freq_Tested_Alle > 0.01),] #17810

###########################################################
#We are not going to trust the CHR_POS from these bad boys#
###########################################################

#We already recovered these bad boys in the lifetime-WHRadjBMI set so we are going to use them:

nexus_results <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/snp_nexus/Lifetime_WHRadjBMI_nexus_results.txt")

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

#First, which of the old ones are in PacksPerYear?

swapping_snps_1 <- merged_rs[index_merged_1,]
swapping_snps_2 <- merged_rs[index_merged_2,]

which(swapping_snps_1$new_rs%in%PacksPerYear_$SNP) #0
which(swapping_snps_1$old_rs%in%PacksPerYear_$SNP) #0

#None, so that means, that only those in GIANT that are in new 
#need to be there.
#Now we need to get the PacksPerYear data that is in either the old or the new column 
#from those that are found in the GIANT data.

swapping_snps_2_in_PacksPerYear <- swapping_snps_2[which(swapping_snps_2$new_rs%in%PacksPerYear_$SNP),] #250 do not need to be changed. Why? Because GIANT has already those. And PacksPerYear too.

#HENCE, PacksPerYear old are those that need to be changed!!

swapping_snps_2_in_PacksPerYear_old <- swapping_snps_2[which(swapping_snps_2$old_rs%in%PacksPerYear_$SNP),] #Only 1.

#So, let's go:

GIANT_new_2_old <- GIANT_nexus_missed[which(GIANT_nexus_missed$RSID%in%swapping_snps_2_in_PacksPerYear_old$new_rs),]

#Now let's change the RSID:

swapping_snps_2_in_PacksPerYear_old <- swapping_snps_2_in_PacksPerYear_old[order(match(swapping_snps_2_in_PacksPerYear_old$new_rs, GIANT_new_2_old$RSID)),]

head(swapping_snps_2_in_PacksPerYear_old)
head(GIANT_new_2_old)

#PERFECT!!

GIANT_new_2_old$RSID <- swapping_snps_2_in_PacksPerYear_old$old_rs

#####################
#Finally let's merge#
#####################

GIANT_rest <- GIANT_nexus_missed[which(!(GIANT_nexus_missed$RSID%in%swapping_snps_2_in_PacksPerYear_old$new_rs)),]

GIANT_missing_end <- rbind(GIANT_rest, GIANT_new_2_old)

GIANT_new_column <- GIANT_nexus_missed[which(GIANT_nexus_missed$RSID%in%swapping_snps_2$new_rs),] #1037
GIANT_old_column <- GIANT_nexus_missed[which(GIANT_nexus_missed$RSID%in%swapping_snps_2$old_rs),] #0 as expected.

###############
#Filtering WHRadjBMI#
###############

#INFO:

WHRadjBMI_ <- WHRadjBMI[which(WHRadjBMI$INFO > 0.7),]

summary(WHRadjBMI_$INFO)

#MAF:

summary(WHRadjBMI_$Freq_Tested_Allele)

WHRadjBMI_[which.min(WHRadjBMI_$Freq_Tested_Allele),]

WHRadjBMI_ <- WHRadjBMI_[which(as.numeric(WHRadjBMI_$Freq_Tested_Allele) > 0.01),]
WHRadjBMI_ <- WHRadjBMI_[which(as.numeric(WHRadjBMI_$Freq_Tested_Allele) < 0.99),]

##################
#Finally, alleles#
##################

yes_vect <- c("A", "G", "T", "C")

WHRadjBMI_ <- WHRadjBMI_[which(WHRadjBMI_$Tested_Allele%in%yes_vect),] #all of them, as expected.
WHRadjBMI_ <- WHRadjBMI_[which(WHRadjBMI_$Other_Allele%in%yes_vect),] #all of them, as expected.

##################################################
#We are not gonna do MHC because we did it before#
##################################################

WHRadjBMI_$chr_pos <- paste(WHRadjBMI_$CHR, WHRadjBMI_$POS, sep = ":")

##################
#MERGING DATASETS#
##################

library(cause)
library(tidyverse)

X_UKBB <- gwas_merge(PacksPerYear_, WHRadjBMI_, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("beta.exposure", "BETA"), 
                se_cols = c("se.exposure", "SE"), 
                A1_cols = c("effect_allele.exposure", "Tested_Allele"), 
                A2_cols = c("other_allele.exposure", "Other_Allele"))

X_GIANT_recovered <- gwas_merge(PacksPerYear_, GIANT_nexus_recovered, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("beta.exposure", "BETA"), 
                se_cols = c("se.exposure", "SE"), 
                A1_cols = c("effect_allele.exposure", "Tested_Allele"), 
                A2_cols = c("other_allele.exposure", "Other_Allele"))

X_GIANT_missed <- gwas_merge(PacksPerYear_, GIANT_missing_end, snp_name_cols = c("SNP", "RSID"), 
                                beta_hat_cols = c("beta.exposure", "BETA"), 
                                se_cols = c("se.exposure", "SE"), 
                                A1_cols = c("effect_allele.exposure", "Tested_Allele"), 
                                A2_cols = c("other_allele.exposure", "Other_Allele"))

#######################################################################
#We need to understand why when we add more SNPs, we wind up with less#
#######################################################################

#Calculating Nuisance:

###############################
#We need to retrieve the RSIDs#
###############################

head(X_UKBB)

WHRadjBMI_[which(WHRadjBMI_$chr_pos == "7:92383888"),] #PERFECT.
PacksPerYear_[which(PacksPerYear_$chr_pos == "7:92383888"),] #PERFECT

WHRadjBMI_match <- WHRadjBMI_[which(WHRadjBMI_$chr_pos%in%X_UKBB$snp),]

WHRadjBMI_match <- WHRadjBMI_match[order(match(WHRadjBMI_match$chr_pos, X_UKBB$snp)),]

which(WHRadjBMI_match$chr_pos != X_UKBB$snp) #they all match. 

head(WHRadjBMI_match$chr_pos)
head(X_UKBB$snp)
tail(WHRadjBMI_match$chr_pos)
tail(X_UKBB$snp)

##########
#All good#
##########

X_UKBB$snp <- WHRadjBMI_match$RSID

head(X_UKBB)

length(which(duplicated(X_UKBB$snp) == TRUE))

#######################################
#MERGING GIANT data recovered by dbsnp#
#######################################

###############################
#We need to retrieve the RSIDs#
###############################

head(X_GIANT_recovered)

GIANT_nexus_recovered[which(GIANT_nexus_recovered$chr_pos == "17:21949792"),] #PERFECT.
PacksPerYear_[which(PacksPerYear_$chr_pos == "17:21949792"),] #PERFECT

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

X_GIANT_in_UKBB <- X_GIANT[which(X_GIANT$snp%in%X_UKBB$snp),] #259. SAME.

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

#Checking for any duplicates:

which(duplicated(X_end$snp) == TRUE) #perfect.

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(params, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/PacksPerYear/dataframes/params_PacksPerYear_WHRadjBMI.rds")
saveRDS(X_end, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/PacksPerYear/dataframes/X_PacksPerYear_WHRadjBMI.rds")

X_end <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/PacksPerYear/dataframes/X_PacksPerYear_WHRadjBMI.rds")
params <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/PacksPerYear/dataframes/params_PacksPerYear_WHRadjBMI.rds")

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

PacksPerYear_pruned <- PacksPerYear[which(PacksPerYear$SNP%in%pruned_all),] #they are all here, no need of weird stuff.

write.table(PacksPerYear_pruned, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/PacksPerYear/dataframes/PacksPerYear_WHRadjBMI_PRUNED.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

#################
#Let's run CAUSE#
#################

top_PacksPerYear_pruned_df <- read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/PacksPerYear/dataframes/PacksPerYear_WHRadjBMI_PRUNED.txt", header = TRUE, stringsAsFactors = FALSE)

top_PacksPerYear_pruned_snps <- top_PacksPerYear_pruned_df$SNP

#Full version:

res <- cause(X=X_end, variants = top_PacksPerYear_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/PacksPerYear/dataframes/RES_PacksPerYear_to_WHRadjBMI_FULL.rds")

res <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/PacksPerYear/dataframes/RES_PacksPerYear_to_WHRadjBMI_FULL.rds")

res$elpd

#model1  model2 delta_elpd se_delta_elpd         z
#1    null sharing  -8.805546      3.764534 -2.339080
#2    null  causal  -9.895379      4.441811 -2.227780
#3 sharing  causal  -1.089833      0.979469 -1.112677

#Sounds really promising:

summary(res)

#p-value testing that causal model is a better fit:  0.13 
#Posterior medians and  95 % credible intervals:
#  model     gamma               eta                 q                  
#[1,] "Sharing" NA                  "0.1 (0.05, 0.22)"  "0.61 (0.21, 0.92)"
#[2,] "Causal"  "0.06 (0.01, 0.15)" "0.02 (-0.2, 0.31)" "0.19 (0, 0.85)" 

plot(res, type="data")


