##############
#INTRODUCTION#
##############

#We are going to run CAUSE for Lifetime smoking as exposure and WC as outcome

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

Lifetime <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/RAW_Data/SmokingCombined_Wootton.txt")
WC <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WC/WC_UKBB_Combined_Curated_FULL.txt")

#############################################
#First we are gonna prune some Lifetime data#
#############################################

head(Lifetime)

#We can filter through EAF and INFO:

summary(Lifetime$INFO) #Min is 0.8. No need to filter a thing.

#################
#Let's go for AF#
#################

summary(Lifetime$EAF) #Min is 0.01 and max 0.99. I think these are approximations, but let's check:

Lifetime <- Lifetime[which(Lifetime$EAF > 0.01),]
Lifetime <- Lifetime[which(Lifetime$EAF < 0.99),]

summary(Lifetime$EAF) #They were approximations!

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

Lifetime_ <- Lifetime[which(!(Lifetime$chr_pos%in%Lifetime_mhc$chr_pos)),] #we go from 7683352 to 7641395. PERFECT MATCH.

############
#Smol check#
############

check <- Lifetime_[which(Lifetime_$CHR == 6),]

check <- check[order(check$BP),]

View(check) #worked like a freaking charm.

############################
#We check for weird alleles#
############################

yes_vect <- c("A", "G", "T", "C")

Lifetime_ <- Lifetime_[which(Lifetime_$EFFECT_ALLELE%in%yes_vect),]
Lifetime_ <- Lifetime_[which(Lifetime_$OTHER_ALLELE%in%yes_vect),]

#############
#WE ARE DONE#
#############

###############
#Filtering WC#
###############

#We don't need to filter much because we have already done it before!

#We do not have info and I checked the frequencies. All good:

summary(WC$eaf.exposure) #perfect.

#For the matching we are going to use chromosome and position, instead of RSID.

##################
#Finally, alleles#
##################

yes_vect <- c("A", "G", "T", "C")

WC <- WC[which(WC$effect_allele.exposure%in%yes_vect),] #all of them, as expected.
WC <- WC[which(WC$other_allele.exposure%in%yes_vect),] #all of them, as expected.

################################################################
#Let's first check whether we have chr_pos for all the variants#
################################################################

WC$chr_pos <- paste(WC$chr.exposure, ":", WC$pos.exposure, sep = "")

WC <- WC[order(WC$chr_pos),]

head(WC) #seems all is OK
tail(WC) #seems all is OK

##################
#MERGING DATASETS#
##################

library(cause)
library(tidyverse)

#NOTE: DUE TO 2SMR FORMATTING THE COLUMNS ARE NAMED AS EXPOSURE, BUT THEY ARE THE OUTCOME
#IN THIS CASE.

X <- gwas_merge(Lifetime_, WC, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("BETA", "beta.exposure"), 
                se_cols = c("SE", "se.exposure"), 
                A1_cols = c("EFFECT_ALLELE", "effect_allele.exposure"), 
                A2_cols = c("OTHER_ALLELE", "other_allele.exposure"))

#Calculating Nuisance:

###############################
#We need to retrieve the RSIDs#
###############################

head(X)

WC[which(WC$chr_pos == "1:693731"),] #PERFECT. Also checked in PS.
Lifetime_[which(Lifetime_$chr_pos == "1:693731"),] #PERFECT

Lifetime_match <- Lifetime_[which(Lifetime_$chr_pos%in%X$snp),]

Lifetime_match <- Lifetime_match[order(match(Lifetime_match$chr_pos, X$snp)),]

which(Lifetime_match$chr_pos != X$snp) #they all match. 

head(Lifetime_match)
head(X)
tail(Lifetime_match)
tail(X)

#Careful. Lifetime_ might have some chr_pos that do not present RSID.
#Let's see if that is the case:

check <- Lifetime_match[order(Lifetime_match$SNP),]

head(check$SNP, 1000) #exactly. for the first 174 SNPs.

#For the first 174 we have them in another format!
#They do not present an RSID, but present the chromsome and position format.

#We need RSIDs, so we are gonna do it the following way:

#######################################
#First we are gonna get the clean ones#
#######################################

index_clean <- which(str_detect(Lifetime_match$SNP, "rs") != FALSE)

X_clean <- X[index_clean,]
Lifetime_clean <- Lifetime_match[index_clean,]

which(X_clean$snp != Lifetime_clean$chr_pos)

head(X_clean$snp)
head(Lifetime_clean$chr_pos)

X_clean$snp <- Lifetime_clean$SNP #now we can return the RSID.

##########
#All good#
##########

index_not_clean <- which(str_detect(Lifetime_match$SNP, "rs") == FALSE)

X_not_clean <- X[index_not_clean,]

#And now we are gonna retrieve the chromosome and position from 
#our dear WC data:

WC_not_clean <- WC[which(WC$chr_pos%in%X_not_clean$snp),]

WC_not_clean <- WC_not_clean[order(match(WC_not_clean$chr_pos, X_not_clean$snp)),]

#Let's do some checks:

which(WC_not_clean$chr_pos != X_not_clean$snp) #cool.

head(WC_not_clean$chr_pos)
head(X_not_clean$snp) 

#Unfortunately, these fellas are also in the same format in WC data!
#Let's see if we can obtain them with phenoscanner.
#Those that we cannot obtain, we will end up removing them.

X_not_clean_copy <- X_not_clean

X_not_clean_copy$snp_ <- paste("chr", X_not_clean_copy$snp, sep = "")

#Finally, let's check in PhenoScanner if these exist:

check_ps_1 <- phenoscanner::phenoscanner(X_not_clean_copy$snp_[seq(1,100)])
check_ps_2 <- phenoscanner::phenoscanner(X_not_clean_copy$snp_[seq(101,174)])

results <- rbind(check_ps_1$snps, check_ps_2$snps) #we recovered one! 

#But the alleles do not match. It is considered an insertion in PS, while we only have standard alleles for CAUSE.

################################################
#I checked in dbSNP and the SNPs do exist there#
################################################

#We will have to make do with this. 
#We will still introduce the non-clean, with their original IDS.
#Maybe CAUSE presents these IDs too, you never now.

X_not_clean$snp <- WC_not_clean$SNP

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

saveRDS(params, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/dataframes/params_Lifetime_WC_UKBB.rds")
saveRDS(X_end, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/dataframes/X_LifetimeInit_WC_UKBB.rds")

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

write.table(Lifetime_pruned, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/dataframes/Lifetime_WC_UKBB_PRUNED.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

#################
#Let's run CAUSE#
#################

top_Lifetime_pruned_df <- read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/dataframes/Lifetime_WC_UKBB_PRUNED.txt", header = TRUE, stringsAsFactors = FALSE)

top_Lifetime_pruned_snps <- top_Lifetime_pruned_df$SNP

#Full version:

res <- cause(X=X_end, variants = top_Lifetime_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/dataframes/RES_Lifetime_to_WC_UKBB_FULL.rds")

res <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/CAUSE/Lifetime_Smoking/dataframes/RES_Lifetime_to_WC_UKBB_FULL.rds")

res$elpd

#model1  model2  delta_elpd se_delta_elpd          z
#1    null sharing -184.214639    15.2227920 -12.101239
#2    null  causal -185.653218    15.4236607 -12.036910
#3 sharing  causal   -1.438579     0.2133242  -6.743628

#Sounds really promising:

summary(res)

#p-value testing that causal model is a better fit:  7.7e-12 
#Posterior medians and  95 % credible intervals:
#  model     gamma               eta                q               
#[1,] "Sharing" NA                  "0.36 (0.32, 0.4)" "0.98 (0.92, 1)"
#[2,] "Causal"  "0.35 (0.29, 0.41)" "0 (-0.29, 0.28)"  "0.18 (0, 0.86)"

plot(res, type="data")


