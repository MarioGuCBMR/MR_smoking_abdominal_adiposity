##############
#INTRODUCTION#
##############

#We are going to run CAUSE for SmkInit as exposure and HIPadjBMI as outcome

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
HIPadjBMI <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/HIPAdjBMI_comb/hipadjbmi_combined_Curated.txt")

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

#####################
#Filtering HIPadjBMI#
#####################

#We don't need to filter much because we have already done it before!

#We do not have info and I checked the frequencies. All good:

summary(HIPadjBMI$FreqAllele1HapMapCEU) #perfect.

#The main issue is gonna be... the matching.

##################
#Finally, alleles#
##################

yes_vect <- c("A", "G", "T", "C")

HIPadjBMI <- HIPadjBMI[which(HIPadjBMI$Allele1%in%yes_vect),] #all of them, as expected.
HIPadjBMI <- HIPadjBMI[which(HIPadjBMI$Allele2%in%yes_vect),] #all of them, as expected.

###############################################################################
#Let's do first a match of those SNPs that do not have chromosome and position#
###############################################################################

HIPadjBMI_missing <- HIPadjBMI[which(HIPadjBMI$chr_pos_37 == "-"),]

which(HIPadjBMI_missing$MarkerName%in%SMK_$RSID) #none of them. 

#Let's check the merged allele versions of those RSIDs:

merged_rs <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Romain_Mario/RAW_DATA/MergedRs/RsMergeArch.bcp.gz")

merged_rs$old_rs <- paste("rs", merged_rs$V1, sep = "")
merged_rs$new_rs <- paste("rs", merged_rs$V2, sep = "")

###################
#Now we match them#
###################

index_merged_1 <- which(merged_rs$old_rs%in%HIPadjBMI_missing$MarkerName) #Only 2! 
index_merged_2 <- which(merged_rs$new_rs%in%HIPadjBMI_missing$MarkerName) #Only 2! 

swapping_snps_1 <- merged_rs[index_merged_1,]
swapping_snps_2 <- merged_rs[index_merged_2,]

which(swapping_snps_1$new_rs%in%SMK_$RSID) #0
which(swapping_snps_1$old_rs%in%SMK_$RSID) #0

which(swapping_snps_2$new_rs%in%SMK_$RSID) #0
which(swapping_snps_2$old_rs%in%SMK_$RSID) #0

#It doesn't matter which combination we do, this still gets out of luck.
#Hence, we are just gonna work with what we have: the chromosome and positions.

##################
#MERGING DATASETS#
##################

library(cause)
library(tidyverse)

#Careful:

SMK_$chr_pos <- paste("chr", SMK_$chr_pos, sep = "")

X <- gwas_merge(SMK_, HIPadjBMI, snp_name_cols = c("chr_pos", "chr_pos_37"), 
                beta_hat_cols = c("BETA", "b"), 
                se_cols = c("SE", "se"), 
                A1_cols = c("ALT", "Allele1"), 
                A2_cols = c("REF", "Allele2"))

#Calculating Nuisance:

###############################
#We need to retrieve the RSIDs#
###############################

head(X)

HIPadjBMI[which(HIPadjBMI$chr_pos_37 == "chr1:752566"),] #PERFECT. Also checked in PS.
SMK_[which(SMK_$chr_pos == "chr1:752566"),] #PERFECT

SMK_match <- SMK_[which(SMK_$chr_pos%in%X$snp),]

SMK_match <- SMK_match[order(match(SMK_match$chr_pos, X$snp)),]

which(SMK_match$chr_pos != X$snp) #they all match. 

head(SMK_match)
head(X)
tail(SMK_match)
tail(X)

#Careful. Smk_Init might have some chr_pos that do not present RSID.
#Let's see if that is the case:

check <- SMK_match[order(SMK_match$RSID),]

head(check$RSID, 100) #exactly. for the first 23 SNPs.

#We need RSIDs, so we are gonna do it the following way:

#######################################
#First we are gonna get the clean ones#
#######################################

index_clean <- which(SMK_match$RSID != ".")

X_clean <- X[index_clean,]
SMK_clean <- SMK_match[index_clean,]

which(X_clean$snp != SMK_clean$chr_pos)

head(X_clean$snp)
head(SMK_clean$chr_pos)

X_clean$snp <- SMK_clean$RSID #now we can return the RSID.

##########
#All good#
##########

index_not_clean <- which(SMK_match$RSID == ".")

X_not_clean <- X[index_not_clean,]

#And now we are gonna retrieve the chromosome and position from 
#our dear HIPadjBMI data:

HIPadjBMI_not_clean <- HIPadjBMI[which(HIPadjBMI$chr_pos_37%in%X_not_clean$snp),]

HIPadjBMI_not_clean <- HIPadjBMI_not_clean[order(match(HIPadjBMI_not_clean$chr_pos_37, X_not_clean$snp)),]

#Let's do some checks:

which(HIPadjBMI_not_clean$chr_pos_37 != X_not_clean$snp) #cool.

head(HIPadjBMI_not_clean$chr_pos_37)
head(X_not_clean$snp) 

#Perfect, now we retrieve the RSIDs:

X_not_clean$snp <- HIPadjBMI_not_clean$MarkerName

#Finally, let's check in PhenoScanner if these exist:

check_ps <- phenoscanner::phenoscanner(X_not_clean$snp)

results <- check_ps$snps #they are not in PS.

#Then, the chances that those RSIDs are found in the CAUSE df is low...
#Let's try with the chromosome and the positions:

check_ps <- phenoscanner::phenoscanner(HIPadjBMI_not_clean$chr_pos_37)

results <- check_ps$snps #still 0.

################################################
#I checked in dbSNP and the SNPs do exist there#
################################################

#We will have to make do with this:

X_end <- rbind(X_clean, X_not_clean)

head(X_end)
tail(X_end)

#PERFECT.

##########################
#RUNNING CAUSE PARAMETERS#
##########################

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(params, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Smoking_Initiation/dataframes/params_SMKInit_HIPadjBMI.rds")
saveRDS(X_end, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Smoking_Initiation/dataframes/X_SMKInit_HIPadjBMI.rds")

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

write.table(SMK_pruned, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Smoking_Initiation/dataframes/SMKInit_HIPadjBMI_PRUNED.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

#################
#Let's run CAUSE#
#################

top_smk_pruned_df <- read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Smoking_Initiation/dataframes/SMKInit_HIPadjBMI_PRUNED.txt", header = TRUE, stringsAsFactors = FALSE)

top_smk_pruned_snps <- top_smk_pruned_df$RSID

#Full version:

res <- cause(X=X_end, variants = top_smk_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Smoking_Initiation/dataframes/RES_SMKinit_to_HIPadjBMI_FULL.rds")

res <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Smoking_Initiation/dataframes/RES_SMKinit_to_HIPadjBMI_FULL.rds")

res$elpd

#model1  model2 delta_elpd se_delta_elpd          z
#1    null sharing -0.3435779     1.2480973 -0.2752813
#2    null  causal -0.1968079     1.7720938 -0.1110595
#3 sharing  causal  0.1467700     0.5310463  0.2763789

#This is really bad:

summary(res)

#p-value testing that causal model is a better fit:  0.61 
#Posterior medians and  95 % credible intervals:
#  model     gamma                eta                   q               
#[1,] "Sharing" NA                   "-0.07 (-0.58, 0.08)" "0.17 (0, 0.75)"
#[2,] "Causal"  "-0.02 (-0.1, 0.05)" "0 (-0.39, 0.36)"     "0.18 (0, 0.84)"

#OH MY. The causal effects are negative, but non-causal. Also the shared pathways are really small.
#However, the results are not that surprising Mario. They are what we expected from 2SMR!

plot(res, type="data")


