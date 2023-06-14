##############
#INTRODUCTION#
##############

#We are going to run CAUSE for CigDayInit as exposure and WHR as outcome

memory.limit(size=800000000)

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

CigDay <- fread("C:/Users/zlc436/Desktop/SMK_WHRAdjBMI/CigarettesPerDay.txt.gz")
WCAdjBMI <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WCAdjBMIsmk/WCAdjBMI_Smk_Strat_Curated.txt")

###########################################
#First we are gonna prune some CigDay data#
###########################################

head(CigDay)

#We can filter through AF, INFO and STAT... somehow? I dunno how yet.

CigDay$INFO <- CigDay$EFFECTIVE_N/CigDay$N

summary(CigDay$INFO) #Min should be 0.30 according to readme.

CigDay <- CigDay[which(CigDay$INFO > 0.70),]

summary(CigDay$INFO) #perfect

#################
#Let's go for AF#
#################

CigDay <- CigDay[which(CigDay$AF > 0.01),]
CigDay <- CigDay[which(CigDay$AF < 0.99),]

summary(CigDay$AF) #perfect

##########################
#Let's check for the STAT#
##########################

summary(CigDay$STAT) #Iit also reflects association. It has nothing to do with our data.

head(CigDay)

###################################
#Finally let's get the chr and pos#
###################################

chr_ <- paste("chr", CigDay$CHROM, sep = "")
CigDay$chr_pos <- paste(chr_, CigDay$POS, sep = ":")

CigDay <- CigDay[order(CigDay$chr_pos),]

head(CigDay)
tail(CigDay)

###########################################
#Let's filter with the extended MHC region#
###########################################

CigDay_mhc <- CigDay[which(CigDay$CHROM == 6 & CigDay$POS >= 26000000 & CigDay$POS <= 34000000),]

summary(CigDay_mhc$CHROM) #perfect
summary(CigDay_mhc$POS) #perfect

index <- which(CigDay$CHROM == 6 & CigDay$POS >= 26000000 & CigDay$POS <= 34000000)

CigDay_ <- CigDay[-index,] #we go 7748746 to 7704844. PERFECT MATCH.

which(CigDay_$CHROM == 6 & CigDay_$POS >= 26000000 & CigDay_$POS <= 34000000)

#Perfect.

############################
#We check for weird alleles#
############################

yes_vect <- c("A", "G", "T", "C")

CigDay_ <- CigDay_[which(CigDay_$REF%in%yes_vect),]
CigDay_ <- CigDay_[which(CigDay_$ALT%in%yes_vect),]

#############
#WE ARE DONE#
#############

##############
#Filtering WC#
##############

head(WCAdjBMI)

#First let's see if we have all the data possible:

WCAdjBMI <- WCAdjBMI[order(WCAdjBMI$chr_pos_37),]

head(WCAdjBMI$chr_pos_37) #we will check later if they match.

#We have 159 SNPs that we could not find the chr_pos at all.

#######################
#Let's make some tests#
#######################

WCAdjBMI_not_liftover_weird <- WCAdjBMI[which(WCAdjBMI$chr_pos_37 == "-"),]

#Let's check whether the RSIDs can be matched with Smoking data:

which(WCAdjBMI_not_liftover_weird$rs_id%in%CigDay_$RSID) #none of them. 

#Let's check one more time, but with merged RSIDs
#If this does not work... then they are not here.

merged_rs <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Romain_Mario/RAW_DATA/MergedRs/RsMergeArch.bcp.gz")

merged_rs$old_rs <- paste("rs", merged_rs$V1, sep = "")
merged_rs$new_rs <- paste("rs", merged_rs$V2, sep = "")

###################
#Now we match them#
###################

index_merged <- which(merged_rs$old_rs%in%WCAdjBMI_not_liftover_weird$rs_id) #Only 2! 

swapping_snps <- merged_rs[index_merged,]

###################
#Now we match them#
###################

index_merged_old <- which(merged_rs$old_rs%in%WCAdjBMI_not_liftover_weird$rs_id) #Only 2! 
index_merged_new <- which(merged_rs$new_rs%in%WCAdjBMI_not_liftover_weird$rs_id) #We have 556. That means that we have A LOT of duplicates. All of them were too new it seems!! 

swapping_snps_old_2_new <- merged_rs[index_merged_old,] #2
swapping_snps_new_2_old <- merged_rs[index_merged_new,] #556

#Now we are gonna check which of these can be matched.

#First we match which of the old RSIDs can be match to CigDay with their updated RSIDs:

which(swapping_snps_old_2_new$new_rs%in%CigDay_$RSID) #none.

#Let's check again if those old are in CigDay, but we saw before that there are 0:

which(swapping_snps_new_2_old$old_rs%in%CigDay_$RSID) #Exactly.

#Finally, 553 of them present new RSIDs, meaning that maybe we cannot match them because they are too new:

which(swapping_snps_new_2_old$old_rs%in%CigDay_$RSID) #Nope, they are not here.

#And lastly, let's do the same and check whether the new ones that are coolio.

which(swapping_snps_new_2_old$new_rs%in%CigDay_$RSID) #Nope, they are not here.

#Nothing that we can do. 
#We will forget about these SNPs.

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
WCAdjBMI <- WCAdjBMI[which(WCAdjBMI$Effect_Allele%in%yes_vect),] #all of them, as expected.

#We are not gonna do MHC because we did it before.

##################
#MERGING DATASETS#
##################

library(cause)
library(tidyverse)

X <- gwas_merge(CigDay_, WCAdjBMI, snp_name_cols = c("chr_pos", "chr_pos_37"), 
                beta_hat_cols = c("BETA", "Effect_SMK"), 
                se_cols = c("SE", "StdErr_SMK"), 
                A1_cols = c("ALT", "Effect_Allele"), 
                A2_cols = c("REF", "Other_Allele"))

#Calculating Nuisance:

###############################
#We need to retrieve the RSIDs#
###############################

head(X)

WCAdjBMI[which(WCAdjBMI$chr_pos_37 == "chr10:100004906"),] #PERFECT.
CigDay_[which(CigDay_$chr_pos == "chr10:100004906"),] #PERFECT

####################################################
#We are getting the RSIDs from Cig_Day_ if possible#
####################################################

CigDay_match <- CigDay_[which(CigDay_$chr_pos%in%X$snp),]

#Smol check:

CigDay_match <- CigDay_match[order(CigDay_match$RSID),]

head(CigDay_match$RSID, 100) #We are only missing around 20SNPs.

#Let's check whether they can be removed due to pval:

CigDay_match_test <- CigDay_match[which(CigDay_match$RSID == "."),]

summary(CigDay_match_test$PVALUE) #they have P > 0.001, the threshold above. They can totally be removed.

##############################################
#Making sure that we remove the SNPs properly#
##############################################

snps_2_remove <- CigDay_match_test$chr_pos

X_end <- X[which(!(X$snp%in%snps_2_remove)),] #Perfect.

#And now we get the RSIDs from the exposure data:

CigDay_match_end <- CigDay_[which(CigDay_$chr_pos%in%X_end$snp),]

CigDay_match_end <- CigDay_match_end[order(match(CigDay_match_end$chr_pos, X_end$snp)),]

which(CigDay_match_end$chr_pos != X_end$snp) #they all match. 

head(CigDay_match_end$chr_pos)
head(X_end$snp)
tail(CigDay_match_end$chr_pos)
tail(X_end$snp)

##########
#All good#
##########

X_end$snp <- CigDay_match_end$RSID

head(X_end)

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(params, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/CigDay/CigDay-WCAdjBMIsmk/dataframes/params_CigDayInit_WCAdjBMIsmk.rds")
saveRDS(X_end, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/CigDay/CigDay-WCAdjBMIsmk/dataframes/X_CigDayInit_WCAdiBMIsmk.rds")

check <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/CigDay/CigDay-WCAdjBMIsmk/dataframes/X_CigDayInit_WCAdiBMIsmk.rds") 

length(which(check$snp%in%X_end$snp)) #all of them.

#Let's check that the data is the same

head(X_end, 1)
head(check, 1)

#0 problems found here!!

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

CigDay_pruned <- CigDay_[which(CigDay_$RSID%in%pruned_all),] #they are all here, no need of weird stuff.

write.table(CigDay_pruned, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/CigDay/CigDay-WCAdjBMIsmk/dataframes/CigDayInit_WCAdjBMIsmk_PRUNED.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

#################
#Let's run CAUSE#
#################

top_CigDay_pruned_df <- read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/CigDay/CigDay-WCAdjBMIsmk/dataframes/CigDayInit_WCAdjBMIsmk_PRUNED.txt", header = TRUE, stringsAsFactors = FALSE)

top_CigDay_pruned_snps <- top_CigDay_pruned_df$RSID

#Full version:

res <- cause(X=X_end, variants = top_CigDay_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/CigDay/CigDay-WCAdjBMIsmk/dataframes/RES_CigDay_to_WCAdjBMIsmk_FULL.rds")

res <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/CigDay/CigDay-WCAdjBMIsmk/dataframes/RES_CigDay_to_WCAdjBMIsmk_FULL.rds")

res$elpd

#model1  model2 delta_elpd se_delta_elpd         z
#1    null sharing  0.4204150     0.3295496 1.2757260
#2    null  causal  0.7581764     0.9386353 0.8077433
#3 sharing  causal  0.3377614     0.6142236 0.5498998

#Sounds really promising:

summary(res)

#p-value testing that causal model is a better fit:  0.71 
#Posterior medians and  95 % credible intervals:
#  model     gamma                eta                  q               
#[1,] "Sharing" NA                   "0.04 (-0.46, 0.53)" "0.14 (0, 0.75)"
#[2,] "Causal"  "0.03 (-0.11, 0.11)" "-0.02 (-0.48, 0.4)" "0.19 (0, 0.86)"

plot(res, type="data")


