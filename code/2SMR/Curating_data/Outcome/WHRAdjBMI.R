##############
#INTRODUCTION#
##############

#This is a quick code to check how much do we need to clean the WHRAdjBMI stratified data.

###########
#Libraries#
###########

library(data.table)
library(tidyverse)
library(TwoSampleMR)

##############
#Loading data#
##############

whradjbmi <- fread("C:/Users/zlc436/Desktop/SMK_WHRAdjBMI/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt") 

#Let's check what columns do we have:

whradjbmi <- whradjbmi[order(whradjbmi$MarkerName),]

head(whradjbmi, 1000)

#Here we just have the SNPs. 

##############################
#Filtering data and saving it#
##############################

#What we are going to do is the following:

#a) Filter for MAF/INFO if possible.
#b) remove MHC region SNPs here, if possible. If not, in the outcome.
#c) Save the data so we do not have to wait for the long query in the database.

whradjbmi_maf <- whradjbmi[which(whradjbmi$FreqAllele1HapMapCEU > 0.01),]
whradjbmi_maf <- whradjbmi_maf[which(whradjbmi_maf$FreqAllele1HapMapCEU < 0.99),]

summary(whradjbmi_maf$FreqAllele1HapMapCEU) #perfect.

#Now the INFO, do we have that data?

colnames(whradjbmi_maf) #we do not. 

##########################
#REMOVING MHC region SNPs#
##########################

#WAIT. We cannot do this because we use the filters from build 37 from Warrington et al.
#Hence, what we are gonna do is... just include them.
#Since we removed them in the exposure, they will be removed automatically.

################################################
#Getting as many chr_pos possible from BMI data#
################################################

BMI <- fread("C:/Users/zlc436/Desktop/Leisure_Project/BMI_treated.txt")

BMI <- BMI[which(is.na(BMI$INFO) == FALSE),]
BMI <- BMI[which(is.na(BMI$CHR) == FALSE),]

summary(BMI$CHR) #perfect
summary(BMI$INFO) #perfect

#We removed the weird SNPs.

chr_ <- paste("chr", BMI$CHR, sep = "")
BMI$chr_pos <- paste(chr_, BMI$POS, sep = ":")

#Now we are gonna match by SNP, 
#but we might miss some of them because they were merged.
#But we can do that later.

BMI_match <- BMI[which(BMI$SNP%in%whradjbmi_maf$MarkerName),] #most of them:

whradjbmi_Pulit <- whradjbmi_maf[which(whradjbmi_maf$MarkerName%in%BMI_match$SNP),]

#We have duplicates, but they are bound to have the same chr_pos, so we do not care
#about them. Let's check just in case:

dupl <- BMI_match$SNP[which(duplicated(BMI_match$SNP) == TRUE)]

BMI_dupl <- BMI_match[which(BMI_match$SNP%in%dupl),]

#Let's order them and get ready the dataframe:

BMI_dupl <- BMI_dupl[order(BMI_dupl$SNP),]

head(BMI_dupl) #same chr_pos
tail(BMI_dupl) #same chr_pos

#PERFECT. We can remove duplicates:

BMI_match <- BMI_match[which(duplicated(BMI_match$SNP) == FALSE),] #now they match.

BMI_match <- BMI_match[order(match(BMI_match$SNP, whradjbmi_Pulit$MarkerName)),]

#Let's check if we did it properly:

which(BMI_match$SNP != whradjbmi_Pulit$MarkerName) #perfect.

head(BMI_match$SNP)
head(whradjbmi_Pulit$MarkerName)

#We can match the chr_pos:

whradjbmi_Pulit$chr_pos_37 <- BMI_match$chr_pos

#########
#PERFECT#
#########

#Let's go and get those missing:

whradjbmi_missing <- whradjbmi_maf[which(!(whradjbmi_maf$MarkerName%in%BMI_match$SNP)),]

#Let's see how many of them can be recovered from SNP Nexus:

for_nexus <- whradjbmi_missing

for_nexus$dbsnp <- "dbsnp"

for_nexus <- for_nexus %>%
  select(dbsnp, MarkerName)

for_nexus <- for_nexus[order(for_nexus$MarkerName),]

head(for_nexus)
tail(for_nexus)

########################################################
#Important we have three SNP that are in another format#
########################################################

#lifted_over manually:

#chr1:161981977-161981979
#chr6:13159524-13159526
#chr6:26020935-26020937

gotcha_vect <- c("chr1:161981978", "chr6:13159525", "chr6:26020936")

which(BMI$chr_pos%in%gotcha_vect) #can we recover them with BMI? YES.

BMI_gotcha <- BMI[which(BMI$chr_pos%in%gotcha_vect),]

###############################################
#Let's prepare the data for getting these snps#
###############################################

whradjbmi_missing <- whradjbmi_missing[order(whradjbmi_missing$MarkerName),]

whradjbmi_missing_gotcha <- whradjbmi_missing[seq(1,3),]

#After checking, the SNP match. 
#Hence:

whradjbmi_missing_gotcha$MarkerName <- BMI_gotcha$SNP

whradjbmi_missing_gotcha$chr_pos_37 <- BMI_gotcha$chr_pos

#########
#PERFECT#
#########

#With that being settled..., we only need to run SNPNexus:

head(for_nexus)

for_nexus <- for_nexus[-seq(1,3),]

head(for_nexus) #PERFECT:

write.table(for_nexus, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMI_comb/WHRAdjBMI_combined_4_Nexus.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

####################
#Let's recover them#
####################

recovered_snps <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMI_comb/WHRAdjBMI_combined_in_Nexus.txt")

recovered_snps <- recovered_snps[order(recovered_snps$`REF Allele`),]

head(recovered_snps) #one deletion
tail(recovered_snps) #nothing.

#It seems we only need to remove the deletion. 
#To do so, since we have IUPAC alleles we will just take all the alleles
#that are letters, to avoid insertions and deletions.

good_alleles <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K",
                  "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U",
                  "V", "W", "X", "Y", "Z")

recovered_snps <- recovered_snps[which(recovered_snps$`REF Allele`%in%good_alleles),] #19596->19595. PERFECT.
recovered_snps <- recovered_snps[which(recovered_snps$`ALT Allele (IUPAC)`%in%good_alleles),] #19596->19595. PERFECT.

whradjbmi_recovered <- whradjbmi_missing[which(whradjbmi_missing$MarkerName%in%recovered_snps$dbSNP),]

recovered_snps <- recovered_snps[order(match(recovered_snps$dbSNP, whradjbmi_recovered$MarkerName)),]

#Doing some checkity checks:

which(recovered_snps$dbSNP != whradjbmi_recovered$MarkerName)

head(recovered_snps$dbSNP)
head(whradjbmi_recovered$MarkerName)

#PERFECT:

chr_ <- paste("chr", recovered_snps$Chromosome, sep = "")

recovered_snps$chr_pos <- paste(chr_, recovered_snps$Position, sep = ":")

whradjbmi_recovered$chr_pos_37 <- recovered_snps$chr_pos

###############################
#Let's check the missing ones:#
###############################

whradjbmi_missing_2 <- whradjbmi_missing[which(!(whradjbmi_missing$MarkerName%in%recovered_snps$dbSNP)),]

View(whradjbmi_missing_2)

#The first three were already recovered, so let's remove them:

whradjbmi_missing_2 <- whradjbmi_missing_2[-seq(1,3),]

head(whradjbmi_missing_2)

################################################################################################################
#For the missing 2000 we are gonna try to match them with the curated WHR, since most of them were solved there# 
################################################################################################################

WHR_curated <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHR/WHR_combined_Curated.txt")

whradjbmi_missing_2_new <- whradjbmi_missing_2[which(whradjbmi_missing_2$MarkerName%in%WHR_curated$MarkerName),] #most of them!!

WHR_curated_match <- WHR_curated[which(WHR_curated$MarkerName%in%whradjbmi_missing_2_new$MarkerName),] #perfect match.

######################################################
#But we are only gonna use those that have chr_pos_37#
######################################################

WHR_curated_match <- WHR_curated_match[order(WHR_curated_match$chr_pos_37),]

head(WHR_curated_match) #as expected

WHR_curated_match <- WHR_curated_match[which(WHR_curated_match$chr_pos_37 != "-"),] #we removed those that could not be found in WHR. The number makes sense.

whradjbmi_missing_2_new_match <- whradjbmi_missing_2_new[which(whradjbmi_missing_2_new$MarkerName%in%WHR_curated_match$MarkerName),]

#############################################################
#Now let's order them and retrieve the proper amount of data#
#############################################################

whradjbmi_missing_2_new_match <- whradjbmi_missing_2_new_match[order(match(whradjbmi_missing_2_new_match$MarkerName, WHR_curated_match$MarkerName)),]

head(whradjbmi_missing_2_new_match$MarkerName)
head(WHR_curated_match$MarkerName)

which(whradjbmi_missing_2_new_match$MarkerName != WHR_curated_match$MarkerName) #PERFECT.

whradjbmi_missing_2_new_match$chr_pos_37 <- WHR_curated_match$chr_pos_37

#Done!

##########################################################
#Now let's get those that we could not retrieve before...#
##########################################################

whradjbmi_missing_3 <- whradjbmi_missing_2[which(!(whradjbmi_missing_2$MarkerName%in%whradjbmi_missing_2_new_match$MarkerName)),] #the number makes sense?

#2471 - 2286 #it does

#Let's take a quick look at them:

whradjbmi_missing_3[order(whradjbmi_missing_3$MarkerName),] #they present RSIDs... 

#Let's try with phenoscanner. Our last resort:

results_1 <- phenoscanner::phenoscanner(whradjbmi_missing_3$MarkerName[seq(1,100)])
results_2 <- phenoscanner::phenoscanner(whradjbmi_missing_3$MarkerName[seq(100,184)])

results <- rbind(results_1$snps, results_2$snps) #We recovered 7 SNPs!!!

whradjbmi_missing_3_match <- whradjbmi_missing_3[which(whradjbmi_missing_3$MarkerName%in%results$snp),]

results <- results[order(match(results$snp, whradjbmi_missing_3_match$MarkerName)),]

head(results$snp)
head(whradjbmi_missing_3_match$MarkerName)

whradjbmi_missing_3_match$chr_pos_37 <- results$hg19_coordinates

#PERFECT.

################################
#Are some of the SNPs merged???#
################################

#Let's try our final luck:

whradjbmi_missing_4 <- whradjbmi_missing_3[which(!(whradjbmi_missing_3$MarkerName%in%whradjbmi_missing_3_match$MarkerName)),]

#7 less SNPs, as expected:
#Let's see if the 177 SNPs are merged, or not.

merged_rs <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Romain_Mario/RAW_DATA/MergedRs/RsMergeArch.bcp.gz")

merged_rs$old_rs <- paste("rs", merged_rs$V1, sep = "")
merged_rs$new_rs <- paste("rs", merged_rs$V2, sep = "")

index_merged <- which(merged_rs$old_rs%in%whradjbmi_missing_4$MarkerName) #2 of them:

swapping_snps <- merged_rs[index_merged,]

swapping_snps

#V1        V2  V3 V4                      V5                      V6        V7 V8            V9
#1:      668 281865545 136  1 2014-08-25 17:55:49.607 2014-08-25 17:55:49.607 281865545  1 JIRA SNP-6845
#2: 12020931   9604473 131  0   2009-12-02 15:53:00.0   2010-03-15 13:21:00.0   9604473  0           rsm
#old_rs      new_rs
#1:      rs668 rs281865545
#2: rs12020931   rs9604473

#I checked them manually.
#They all have weird positions in build 37 that do not allow to locate them properly.
#So far, we are just gonna give them the "-" treatment. If the RSIDs are found while matching...
#Then better.

whradjbmi_missing_4$chr_pos_37 <- "-"

#############
#CONCLUSIONS#
#############

#We solved most of the SNPs, except for 2000, let's see how can we match them afterwards playing with the merged SNPs.

whradjbmi_end <- rbind(whradjbmi_Pulit, whradjbmi_recovered, whradjbmi_missing_gotcha, whradjbmi_missing_2_new_match, whradjbmi_missing_3_match, whradjbmi_missing_4) #PERFECT number. It matches the data perfectily.

dim(whradjbmi_end)
dim(whradjbmi_maf)

fwrite(whradjbmi_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMI_comb/WHRAdjBMI_combined_Curated.txt")

############################################
#Finally, let's see if we can load the data#
############################################

check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMI_comb/WHRAdjBMI_combined_Curated.txt")
