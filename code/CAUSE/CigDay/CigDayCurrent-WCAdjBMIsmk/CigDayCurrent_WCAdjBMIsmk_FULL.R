##############
#INTRODUCTION#
##############

#We are going to run CAUSE for CigDayInit as exposure and WHR as outcome

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

clean_liftover <- function(chr_pos){
  
  tmp_1 <- strsplit(chr_pos, ":")[[1]][1]
  tmp_2 <- strsplit(chr_pos, ":")[[1]][2]
  tmp_3 <- as.numeric(as.character(strsplit(tmp_2, "-")[[1]][1])) +1
  
  final <- paste(tmp_1, tmp_3, sep = ":")
  
  return(final)
  
}

##############
#Loading data#
##############

CigDay <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/CigDayCurrent/Cig_Day_Current_Curated_FULL.txt")
WCAdjBMI <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WCAdjBMIsmk/WCAdjBMI_Smk_Strat_Curated.txt")

###########################################
#First we are gonna prune some CigDay data#
###########################################

head(CigDay)

#We can filter through EAF, INFO and STAT... somehow? I dunno how yet.

##################
#Let's go for EAF#
##################

summary(CigDay$eaf.exposure) #perfect no need to filter

###################################
#Finally let's get the chr and pos#
###################################

chr_ <- paste("chr", CigDay$chr.exposure, sep = "")
CigDay$chr_pos <- paste(chr_, CigDay$pos.exposure, sep = ":")

CigDay <- CigDay[order(CigDay$chr_pos),]

head(CigDay)
tail(CigDay)

#Perfect

###########################################
#Let's filter with the extended MHC region#
###########################################

CigDay_mhc <- CigDay[which(CigDay$chr_pos == 6 & CigDay$pos.exposure >= 26000000 & CigDay$pos.exposure <= 34000000),]

summary(CigDay_mhc$CHROM) #perfect
summary(CigDay_mhc$POS) #perfect

index <- which(CigDay$CHROM == 6 & CigDay$POS >= 26000000 & CigDay$POS <= 34000000)

#There is nothing!! We curated it beforehand :)
#Perfect.

############################
#We check for weird alleles#
############################

yes_vect <- c("A", "G", "T", "C")

CigDay <- CigDay[which(CigDay$effect_allele.exposure%in%yes_vect),]
CigDay <- CigDay[which(CigDay$other_allele.exposure%in%yes_vect),]

#all of them are perfect.

#############
#WE ARE DONE#
#############

#####################
#Filtering WHRAdjBMI#
#####################

head(WCAdjBMI)

#First let's see if we have all the data possible:

WCAdjBMI <- WCAdjBMI[order(WCAdjBMI$chr_pos_37),]

head(WCAdjBMI$chr_pos_37) #we will check later if they match.

#We have 157 SNPs that we could not find the chr_pos at all.

#######################
#Let's make some tests#
#######################

WCAdjBMI_not_liftover_weird <- WCAdjBMI[which(WCAdjBMI$chr_pos_37 == "-"),]

#Let's check whether the RSIDs can be matched with Smoking data:

which(WCAdjBMI_not_liftover_weird$rs_id%in%CigDay$SNP) #none of them. 

#Let's check one more time, but with merged RSIDs
#If this does not work... then they are not here.

merged_rs <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Romain_Mario/RAW_DATA/MergedRs/RsMergeArch.bcp.gz")

merged_rs$old_rs <- paste("rs", merged_rs$V1, sep = "")
merged_rs$new_rs <- paste("rs", merged_rs$V2, sep = "")

###################
#Now we match them#
###################

index_merged_old <- which(merged_rs$old_rs%in%WCAdjBMI_not_liftover_weird$rs_id) #Only 2! 
index_merged_new <- which(merged_rs$new_rs%in%WCAdjBMI_not_liftover_weird$rs_id) #We have 566. That means that we have A LOT of duplicates. All of them were too new it seems!! 

swapping_snps_old_2_new <- merged_rs[index_merged_old,]
swapping_snps_new_2_old <- merged_rs[index_merged_new,]

#Now we are gonna check which of these can be matched.

which(swapping_snps_old_2_new$new_rs%in%CigDay$SNP) #0: that means that we don't need to change these.
which(swapping_snps_new_2_old$old_rs%in%CigDay$SNP) #0: that means that we don't need to change these either.
which(swapping_snps_new_2_old$new_rs%in%CigDay$SNP) #0: that means that we don't need to change these either.

#Nothing that we can do. 
#We will forget about these SNPs.
#It might seem that these SNPs are lost in build 37, from what I could gather from liftover 
#and from dbSNP info. 

#So, let's forget about them, yeah.

################################################
#We already did some filtering, but let's check#
################################################

summary(WCAdjBMI$freq_allele1_hapmapceu) #perfect

#We cannot do anything about INFO and the MHC will disappear once we match the data

##################
#Finally, alleles#
##################

yes_vect <- c("A", "G", "T", "C")

WCAdjBMI <- WCAdjBMI[which(WCAdjBMI$Effect_Allele%in%yes_vect),] #all of them, as expected.
WCAdjBMI <- WCAdjBMI[which(WCAdjBMI$Other_Allele%in%yes_vect),] #all of them, as expected.

#We are not gonna do MHC because we did it before.

##################
#MERGING DATASETS#
##################

library(cause)
library(tidyverse)

X <- gwas_merge(CigDay, WCAdjBMI, snp_name_cols = c("chr_pos", "chr_pos_37"), 
                beta_hat_cols = c("beta.exposure", "Effect_SMK"), 
                se_cols = c("se.exposure", "StdErr_SMK"), 
                A1_cols = c("effect_allele.exposure", "Effect_Allele"), 
                A2_cols = c("other_allele.exposure", "Other_Allele"))

#Calculating Nuisance:

###############################
#We need to retrieve the RSIDs#
###############################

head(X)

WCAdjBMI[which(WCAdjBMI$chr_pos_37 == "chr10:100004906"),] #PERFECT.
CigDay[which(CigDay$chr_pos == "chr10:100004906"),] #PERFECT

####################################################
#We are getting the RSIDs from Cig_Day_ if possible#
####################################################

CigDay_match <- CigDay[which(CigDay$chr_pos%in%X$snp),]

#Smol check:

CigDay_match <- CigDay_match[order(CigDay_match$SNP),]

head(CigDay_match$SNP, 100) #Some of them are the weird SNPs... in this case.
tail(CigDay_match$SNP, 100) #We are fine

#We are going to have to split it and obtain the RSIDs for those in WHRAdjBMI, if they exist.

#.1 retrieve the weird ones from CigDay_match <-

CigDay_match_rare <-  CigDay_match[seq(1,23),]

#Now let's obtain those for WHRAdjBMI:

WHRAdjBMI_rare <- WCAdjBMI[which(WCAdjBMI$chr_pos_37%in%CigDay_match_rare$chr_pos),] #23

#Let's check the RSIDs:

WHRAdjBMI_rare$rs_id #all of them have them. 

#Let's check if they are alright in dbSNP... CORRECT. Also in build 37. ALL GOOD.

X_rare <- X[which(X$snp%in%WHRAdjBMI_rare$chr_pos_37),] #23/23 perfect.

#Now we order them so that we can replace the chromsome and the position by RSID:

WHRAdjBMI_rare <- WHRAdjBMI_rare[order(match(WHRAdjBMI_rare$chr_pos_37, X_rare$snp)),]

#Let's check...

length(which(WHRAdjBMI_rare$chr_pos_37 == X_rare$snp)) #23/23 perfect

head(WHRAdjBMI_rare$chr_pos_37) #awesome.
head(X_rare$snp)  #awesome.

X_rare$snp <- WHRAdjBMI_rare$rs_id

##############################################################################
#Now let's do the same, but for those that are not rare. Then, merge the data#
##############################################################################

#We will do it with CigDay RSIDS, because they are from 2SMR GWAS and they are curated and pass filters of quality.

CigDay_match_normal <- CigDay_match[which(!(CigDay_match$chr_pos%in%WHRAdjBMI_rare$chr_pos_37)),] #perfect. 23 less.

#Let's check that we really did a good job here:

CigDay_match_normal <- CigDay_match_normal[order(CigDay_match_normal$SNP),] #we should have normal RSIDs only.

head(CigDay_match_normal$SNP) #perfect
tail(CigDay_match_normal$SNP) #perfect

#Awesome!! 

X_normal <- X[which(X$snp%in%CigDay_match_normal$chr_pos),]

#And now we order them and we make sure that all is good...

CigDay_match_normal <- CigDay_match_normal[order(match(CigDay_match_normal$chr_pos, X_normal$snp)),]

length(which(CigDay_match_normal$chr_pos == X_normal$snp)) #1948622
length(which(CigDay_match_normal$chr_pos != X_normal$snp)) #0

#Perfect! One more check...

head(CigDay_match_normal$chr_pos)
head(X_normal$snp)

tail(CigDay_match_normal$chr_pos)
tail(X_normal$snp)

#All set to go:

X_normal$snp <- CigDay_match_normal$SNP

#Finally we merge:

X_end <- rbind(X_normal, X_rare) #perfect! It has the same amount of SNPs as the original.

head(X_end)
tail(X_end)

#Checking first and last SNP:

WCAdjBMI[which(WCAdjBMI$rs_id == "rs3750595"),] #PERFECT.
CigDay[which(CigDay$SNP == "rs3750595"),] #PERFECT

#Careful last SNP does not have RSID in CigDay:

WCAdjBMI[which(WCAdjBMI$rs_id == "rs4005780"),] #PERFECT.
CigDay[which(CigDay$SNP == "rs4005780"),] #no RSID as expected..., but
CigDay[which(CigDay$chr_pos == "chr2:90389326"),] #with chr_pos it worked like a freaking charm.

##########
#All good#
##########

###################
#Smol final check #
###################

which(duplicated(X_end$snp)) #PERFECT.

head(X_end)

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE)) #No need to reduce the sample now!!
params <- est_cause_params(X_end, varlist)

saveRDS(params, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/CigDay/CigDayCurrent-WCAdjBMIsmk/dataframes/params_CigDayCurrent_WCAdjBMIsmk_FULL.rds")
saveRDS(X_end, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/CigDay/CigDayCurrent-WCAdjBMIsmk/dataframes/X_CigDayCurrent_WCAdjBMIsmk_FULL.rds")

check <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/CigDay/CigDayCurrent-WCAdjBMIsmk/dataframes/params_CigDayCurrent_WCAdjBMIsmk.rds")

check$rho #they do change a bit...
params$rho #they do change a bit...

check$var #
params$var 

check$prior #-25
params$prior #84

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

CigDay_pruned <- CigDay[which(CigDay$SNP%in%pruned_all),] #they are all here, no need of weird stuff.

write.table(CigDay_pruned, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/CigDay/CigDayCurrent-WCAdjBMIsmk/dataframes/CigDayCurrent_WCAdjBMIsmk_PRUNED_FULL.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")

#################
#Let's run CAUSE#
#################

top_CigDay_pruned_df <- read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/CigDay/CigDayCurrent-WCAdjBMIsmk/dataframes/CigDayCurrent_WCAdjBMIsmk_PRUNED_FULL.txt", header = TRUE, stringsAsFactors = FALSE, sep = ",")

top_CigDay_pruned_snps <- top_CigDay_pruned_df$SNP

#Full version:

res <- cause(X=X_end, variants = top_CigDay_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/CigDay/CigDayCurrent-WCAdjBMIsmk/dataframes/RES_CigDayCurrent_to_WCAdjBMIsmk_FULL_VCF.rds")

res <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/CigDay/CigDayCurrent-WCAdjBMIsmk/dataframes/RES_CigDayCurrent_to_WCAdjBMIsmk_FULL_VCF.rds")

res$elpd

#model1  model2 delta_elpd se_delta_elpd        z
#1    null sharing  0.4390632     0.1173100 3.742761
#2    null  causal  1.0057515     0.5366074 1.874278
#3 sharing  causal  0.5666882     0.4291566 1.320470

summary(res)

#p-value testing that causal model is a better fit:  0.91 
#Posterior medians and  95 % credible intervals:
#  model     gamma                eta                   q                  
#[1,] "Sharing" NA                   "0.06 (-0.87, 0.85)"  "0.16 (0, 0.75)"   
#[2,] "Causal"  "0.06 (-0.27, 0.29)" "-0.03 (-0.89, 0.74)" "0.22 (0.01, 0.85)"

plot(res, type="data")


