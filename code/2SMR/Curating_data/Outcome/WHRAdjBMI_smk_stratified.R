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

###################
#Loading functions#
###################

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

whradjbmi <- fread("C:/Users/zlc436/Desktop/SMK_WHRAdjBMI/WHRadjBMI.Stratified/WHRadjBMI.StratifiedBySMK.CombinedSexes.EuropeanOnly.txt") 

#Let's check what columns do we have:

head(whradjbmi)

#We have data from the hg18 chromsome and position which is annoyting.
#First we are gonna remove the SNPs that we do not care about.
#Then, we are gonna check if we can recover that data somehow.

##############################
#Filtering data and saving it#
##############################

#What we are going to do is the following:

#a) Filter for MAF/INFO if possible.
#b) remove MHC region SNPs here, if possible. If not, in the outcome.
#c) Save the data so we do not have to wait for the long query in the database.

summary(whradjbmi$EAF_HapMapCEU) #perfect.

whradjbmi_maf <- whradjbmi[which(whradjbmi$EAF_HapMapCEU > 0.01),]
whradjbmi_maf <- whradjbmi_maf[which(whradjbmi_maf$EAF_HapMapCEU < 0.99),]

summary(whradjbmi_maf$EAF_HapMapCEU) #perfect.

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

BMI_match <- BMI[which(BMI$SNP%in%whradjbmi_maf$rs_id),] #most of them:

whradjbmi_Pulit <- whradjbmi_maf[which(whradjbmi_maf$rs_id%in%BMI_match$SNP),]

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

BMI_match <- BMI_match[order(match(BMI_match$SNP, whradjbmi_Pulit$rs_id)),]

#Let's check if we did it properly:

which(BMI_match$SNP != whradjbmi_Pulit$rs_id) #perfect.

head(BMI_match$SNP)
head(whradjbmi_Pulit$rs_id)

#We can match the chr_pos:

whradjbmi_Pulit$chr_pos_37 <- BMI_match$chr_pos

#########
#PERFECT#
#########

#Let's go and get those missing:

whradjbmi_missing <- whradjbmi_maf[which(!(whradjbmi_maf$rs_id%in%BMI_match$SNP)),]

#Let's see how many of them can be recovered from SNP Nexus:

for_nexus <- whradjbmi_missing

for_nexus$dbsnp <- "dbsnp"

for_nexus <- for_nexus %>%
  select(dbsnp, rs_id)

for_nexus <- for_nexus[order(for_nexus$rs_id),]

head(for_nexus)

write.table(for_nexus, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMISmk/WHRAdjBMI_Smk_Strat_4_Nexus.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

####################
#Let's recover them#
####################

recovered_snps <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMISmk/WHRAdjBMI_Smk_Strat_in_Nexus.txt")

recovered_snps <- recovered_snps[order(recovered_snps$`REF Allele`),]

head(recovered_snps) #one deletion
tail(recovered_snps) #nothing.

#It seems we only need to remove the deletion. 
#To do so, since we have IUPAC alleles we will just take all the alleles
#that are letters, to avoid insertions and deletions.

good_alleles <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K",
                  "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U",
                  "V", "W", "X", "Y", "Z")

recovered_snps <- recovered_snps[which(recovered_snps$`REF Allele`%in%good_alleles),] #18082->18081. PERFECT.
recovered_snps <- recovered_snps[which(recovered_snps$`ALT Allele (IUPAC)`%in%good_alleles),] #18082->18081. PERFECT.

whradjbmi_recovered <- whradjbmi_missing[which(whradjbmi_missing$rs_id%in%recovered_snps$dbSNP),]

recovered_snps <- recovered_snps[order(match(recovered_snps$dbSNP, whradjbmi_recovered$rs_id)),]

#Doing some checkity checks:

which(recovered_snps$dbSNP != whradjbmi_recovered$rs_id)

head(recovered_snps$dbSNP)
head(whradjbmi_recovered$rs_id)

#PERFECT:

chr_ <- paste("chr", recovered_snps$Chromosome, sep = "")

recovered_snps$chr_pos <- paste(chr_, recovered_snps$Position, sep = ":")

whradjbmi_recovered$chr_pos_37 <- recovered_snps$chr_pos

###############################
#Let's check the missing ones:#
###############################

whradjbmi_missing_2 <- whradjbmi_missing[which(!(whradjbmi_missing$rs_id%in%recovered_snps$dbSNP)),]

#############################################################
#For the missing 2000 we are gonna lift over the coordinates#
#############################################################

chr_ <- paste("chr", whradjbmi_missing_2$chromosome, sep = "")

pos <- paste(whradjbmi_missing_2$position_hg18 -1, "-", whradjbmi_missing_2$position_hg18+1, sep = "")

chr_pos_18 <- paste(chr_, pos, sep = ":")

whradjbmi_missing_2$chr_pos_18 <- chr_pos_18

#############################
#Okay, let's do the liftover#
#############################

write.table(chr_pos_18, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMISmk/whradjbmi_Smk_Strat_4_liftover.txt", quote = FALSE, row.names = FALSE, col.names =  FALSE)

#########################################################################################
#Now we are gonna redo Liftover without the failing conversions so we can match the data#
#########################################################################################

failed_conversions <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMISmk/Failed_Conversions.txt")

whradjbmi_missing_2_new <- whradjbmi_missing_2[which(!(whradjbmi_missing_2$chr_pos_18%in%failed_conversions$`#Deleted in new`)),]

#The numbers match: we had 157 failed conversions.

write.table(whradjbmi_missing_2_new$chr_pos_18, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMISmk/whradjbmi_Smk_Strat_4_liftover_curated.txt", quote = FALSE, row.names = FALSE, col.names =  FALSE)

converted_data <- as.data.frame(read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMISmk/whradjbmi_Smk_Strat_liftover_correct.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

##########################
#Let's get the chr_pos_37#
##########################

chr_pos_37 <- as.character(unlist(sapply(converted_data$V1, clean_liftover)))

whradjbmi_missing_2_new$chr_pos_37 <- chr_pos_37

head(whradjbmi_missing_2_new)
tail(whradjbmi_missing_2_new)

#Perfect!!!

######################################################
#Now... what are the ones that we could not retrieve?#
######################################################

whradjbmi_missing_3 <- whradjbmi_missing_2[which(!(whradjbmi_missing_2$rs_id%in%whradjbmi_missing_2_new$rs_id)),]

#Let's take a quick look at them:

whradjbmi_missing_3[order(whradjbmi_missing_3$rs_id),] #they present RSIDs... 

#Let's try with phenoscanner. Our last resort:

results_1 <- phenoscanner::phenoscanner(whradjbmi_missing_3$rs_id[seq(1,100)])
results_2 <- phenoscanner::phenoscanner(whradjbmi_missing_3$rs_id[seq(100,157)])

results <- rbind(results_1$snps, results_2$snps) #0 SNPs found in phenoscanner. Maybe because it is down.

################################
#Are some of the SNPs merged???#
################################

merged_rs <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Romain_Mario/RAW_DATA/MergedRs/RsMergeArch.bcp.gz")

merged_rs$old_rs <- paste("rs", merged_rs$V1, sep = "")
merged_rs$new_rs <- paste("rs", merged_rs$V2, sep = "")

index_merged <- which(merged_rs$old_rs%in%whradjbmi_missing_3$rs_id) #2 of them:

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

whradjbmi_missing_3$chr_pos_37 <- "-"

#############
#CONCLUSIONS#
#############

#We solved most of the SNPs, except for 159, let's see how can we match them afterwards playing with the merged SNPs.
#We need to add a new column, since for missing_2 and missing_3 we used liftover to recover SNPs.

whradjbmi_Pulit$chr_pos_18 <- "-"
whradjbmi_recovered$chr_pos_18 <- "-"

whradjbmi_end <- rbind(whradjbmi_Pulit, whradjbmi_recovered, whradjbmi_missing_2_new, whradjbmi_missing_3) #PERFECT number. It matches the data perfectily.

dim(whradjbmi_end) #PERFECT
dim(whradjbmi_maf) #PERFECT

fwrite(whradjbmi_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMISmk/WHRAdjBMI_Smk_Strat_Curated.txt", na = "-")

############################################
#Finally, let's see if we can load the data#
############################################

check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMISmk/WHRAdjBMI_Smk_Strat_Curated.txt")
