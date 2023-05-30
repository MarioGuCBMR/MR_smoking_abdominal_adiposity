##############
#INTRODUCTION#
##############

#This is a quick code to check how much do we need to clean the wcadjbmi stratified data.

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

wcadjbmi <- fread("C:/Users/zlc436/Downloads/WCadjBMI.StratifiedBySMK.CombinedSexes.EuropeanOnly.txt") 

#Let's check what columns do we have:

head(wcadjbmi)

#We have data from the hg18 chromsome and position which is annoyting.
#First we are gonna remove the SNPs that we do not care about.
#Then, we are gonna check if we can recover that data somehow.

##############################
#Filtering data and saving it#
##############################

#What we are going to do is the following:

#Special post-hoc case: we can NAs in some SMK cases. 
#We are gonna remove those.

wcadjbmi <- wcadjbmi[which(is.na(wcadjbmi$N_SMK) == FALSE),] #Did this do it?

summary(wcadjbmi$Effect_SMK) #perfect we have effect everywhere.

#a) Filter for MAF/INFO if possible.
#b) remove MHC region SNPs here, if possible. If not, in the outcome.
#c) Save the data so we do not have to wait for the long query in the database.

wcadjbmi_maf <- wcadjbmi[which(wcadjbmi$freq_allele1_hapmapceu > 0.01),]
wcadjbmi_maf <- wcadjbmi_maf[which(wcadjbmi_maf$freq_allele1_hapmapceu  < 0.99),]

summary(wcadjbmi_maf$freq_allele1_hapmapceu) #perfect.

#Now the INFO, do we have that data?

colnames(wcadjbmi_maf) #we do not. 

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

BMI_match <- BMI[which(BMI$SNP%in%wcadjbmi_maf$rs_id),] #most of them:

wcadjbmi_Pulit <- wcadjbmi_maf[which(wcadjbmi_maf$rs_id%in%BMI_match$SNP),]

#We have duplicates, but they are bound to have the same chr_pos, so we do not care
#about them. Let's check just in case:

dupl <- BMI_match$SNP[which(duplicated(BMI_match$SNP) == TRUE)]

BMI_dupl <- BMI_match[which(BMI_match$SNP%in%dupl),]

#Let's order them and get ready the dataframe:

BMI_dupl <- BMI_dupl[order(BMI_dupl$SNP),]

head(BMI_dupl) #same chr_pos
tail(BMI_dupl) #same chr_pos

#PERFECT. We can remove duplicates:

BMI_match <- BMI_match[which(duplicated(BMI_match$SNP) == FALSE),] #now match!! 

###################################
#Let's check if we did it properly#
###################################

BMI_match <- BMI_match[order(match(as.character(BMI_match$SNP), as.character(wcadjbmi_Pulit$rs_id))),]

which(BMI_match$SNP != wcadjbmi_Pulit$rs_id) #perfect.

head(BMI_match$SNP)
head(wcadjbmi_Pulit$rs_id)

#We can match the chr_pos:

wcadjbmi_Pulit$chr_pos_37 <- BMI_match$chr_pos

#########
#PERFECT#
#########

#Let's go and get those missing:

wcadjbmi_missing <- wcadjbmi_maf[which(!(wcadjbmi_maf$rs_id%in%BMI_match$SNP)),]

#Let's see how many of them can be recovered from SNP Nexus:

for_nexus <- wcadjbmi_missing

for_nexus$dbsnp <- "dbsnp"

for_nexus <- for_nexus %>%
  select(dbsnp, rs_id)

for_nexus <- for_nexus[order(for_nexus$rs_id),]

head(for_nexus)

write.table(for_nexus, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WCAdjBMIsmk/wccadjbmi_Smk_Strat_4_Nexus.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

####################
#Let's recover them#
####################

recovered_snps <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WCAdjBMIsmk/wcadjbmi_Smk_Strat_in_Nexus.txt")

recovered_snps <- recovered_snps[order(recovered_snps$`REF Allele`),]

head(recovered_snps) #one deletion
tail(recovered_snps) #nothing.

#It seems we only need to remove the deletion. 
#To do so, since we have IUPAC alleles we will just take all the alleles
#that are letters, to avoid insertions and deletions.

good_alleles <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K",
                  "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U",
                  "V", "W", "X", "Y", "Z")

recovered_snps <- recovered_snps[which(recovered_snps$`REF Allele`%in%good_alleles),] #17679->17679. PERFECT.
recovered_snps <- recovered_snps[which(recovered_snps$`ALT Allele (IUPAC)`%in%good_alleles),] #17679->17679. PERFECT.

wcadjbmi_recovered <- wcadjbmi_missing[which(wcadjbmi_missing$rs_id%in%recovered_snps$dbSNP),]

recovered_snps <- recovered_snps[order(match(recovered_snps$dbSNP, wcadjbmi_recovered$rs_id)),]

#Doing some checkity checks:

which(recovered_snps$dbSNP != wcadjbmi_recovered$rs_id)

head(recovered_snps$dbSNP)
head(wcadjbmi_recovered$rs_id)

#PERFECT:

chr_ <- paste("chr", recovered_snps$Chromosome, sep = "")

recovered_snps$chr_pos <- paste(chr_, recovered_snps$Position, sep = ":")

wcadjbmi_recovered$chr_pos_37 <- recovered_snps$chr_pos

###############################
#Let's check the missing ones:#
###############################

wcadjbmi_missing_2 <- wcadjbmi_missing[which(!(wcadjbmi_missing$rs_id%in%recovered_snps$dbSNP)),]

View(wcadjbmi_missing_2)

#############################################################
#For the missing 2000 we are gonna lift over the coordinates#
#############################################################

chr_ <- paste("chr", wcadjbmi_missing_2$chromosome, sep = "")

pos <- paste(wcadjbmi_missing_2$position_hg18 -1, "-", wcadjbmi_missing_2$position_hg18+1, sep = "")

chr_pos_18 <- paste(chr_, pos, sep = ":")

wcadjbmi_missing_2$chr_pos_18 <- chr_pos_18

#############################
#Okay, let's do the liftover#
#############################

write.table(chr_pos_18, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WCAdjBMIsmk/WCAdjBMI_Smk_Strat_4_liftover.txt", quote = FALSE, row.names = FALSE, col.names =  FALSE)

#########################################################################################
#Now we are gonna redo Liftover without the failing conversions so we can match the data#
#########################################################################################

failed_conversions <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WCAdjBMIsmk/Failed_Conversions.txt")

wcadjbmi_missing_2_new <- wcadjbmi_missing_2[which(!(wcadjbmi_missing_2$chr_pos_18%in%failed_conversions$`#Partially deleted in new`)),]

#The numbers match: we had 159 failed conversions.

write.table(wcadjbmi_missing_2_new$chr_pos_18, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WCAdjBMIsmk/WCAdjBMI_Smk_Strat_4_liftover_curated.txt", quote = FALSE, row.names = FALSE, col.names =  FALSE)

converted_data <- as.data.frame(read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WCAdjBMIsmk/WCAdjBMI_Conversion_4_liftover_correct.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

##########################
#Let's get the chr_pos_37#
##########################

chr_pos_37 <- as.character(unlist(sapply(converted_data$V1, clean_liftover)))

wcadjbmi_missing_2_new$chr_pos_37 <- chr_pos_37

head(wcadjbmi_missing_2_new)
tail(wcadjbmi_missing_2_new)

#Perfect!!!

######################################################
#Now... what are the ones that we could not retrieve?#
######################################################

wcadjbmi_missing_3 <- wcadjbmi_missing_2[which(!(wcadjbmi_missing_2$rs_id%in%wcadjbmi_missing_2_new$rs_id)),]

#Let's take a quick look at them:

wcadjbmi_missing_3[order(wcadjbmi_missing_3$rs_id),] #they present RSIDs... 

#Let's try with phenoscanner. Our last resort:

results_1 <- phenoscanner::phenoscanner(wcadjbmi_missing_3$rs_id[seq(1,100)])
results_2 <- phenoscanner::phenoscanner(wcadjbmi_missing_3$rs_id[seq(100,159)])

results <- rbind(results_1$snps, results_2$snps) #0 SNPs found in phenoscanner. Maybe because it is down.

################################
#Are some of the SNPs merged???#
################################

merged_rs <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Romain_Mario/RAW_DATA/MergedRs/RsMergeArch.bcp.gz")

merged_rs$old_rs <- paste("rs", merged_rs$V1, sep = "")
merged_rs$new_rs <- paste("rs", merged_rs$V2, sep = "")

index_merged <- which(merged_rs$old_rs%in%wcadjbmi_missing_3$rs_id) #2 of them:

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

wcadjbmi_missing_3$chr_pos_37 <- "-"

#############
#CONCLUSIONS#
#############

#We solved most of the SNPs, except for 159, let's see how can we match them afterwards playing with the merged SNPs.
#We need to add a new column, since for missing_2 and missing_3 we used liftover to recover SNPs.

wcadjbmi_Pulit$chr_pos_18 <- "-"
wcadjbmi_recovered$chr_pos_18 <- "-"

wcadjbmi_end <- rbind(wcadjbmi_Pulit, wcadjbmi_recovered, wcadjbmi_missing_2_new, wcadjbmi_missing_3) #PERFECT number. It matches the data perfectily.

dim(wcadjbmi_end) #PERFECT
dim(wcadjbmi_maf) #PERFECT

fwrite(wcadjbmi_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WCAdjBMIsmk/WCAdjBMI_Smk_Strat_Curated.txt")

############################################
#Finally, let's see if we can load the data#
############################################

check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WCAdjBMIsmk/WCAdjBMI_Smk_Strat_Curated.txt")
