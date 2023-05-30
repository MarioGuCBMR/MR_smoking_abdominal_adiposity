##############
#INTRODUCTION#
##############

#This is a quick code to check how much do we need to clean smoking initiation data.
#Since this data is gonna be exposure data, we only need to take into account
#that the data is well-prepared before the ld-pruning to ensure that all is good.

###########
#Libraries#
###########

library(data.table)
library(tidyverse)
library(TwoSampleMR)

##############
#Loading data#
##############

#To avoid sample overlap we are using UKB-b-20261. I didn't remember this!!
#Then, the curation is gonna be really simple. 
#We are gonna do it here, just once to avoid waiting too long when we load the data.

#Nonetheless, we will need to the whole dataset if we are looking for proxies later...
#Hence, I don't know if it is the best idea actually.

smoking <- fread("C:/Users/zlc436/Desktop/SMK_WHRAdjBMI/SmokingCombined_Wootton.txt") #Smoking Combined Wootton. 

#Okay, let's check the RSIDs:

smoking <- smoking[order(smoking$SNP),]

head(smoking, 1500) #we have some that are weird.

##############################
#Filtering data and saving it#
##############################

#What we are going to do is the following:

#a) Filter for MAF/INFO if possible.
#b) remove MHC region SNPs here, if possible. If not, in the outcome.
#c) Save the data so we do not have to wait for the long query in the database.

smoking_maf <- smoking[which(smoking$EAF > 0.01),]
smoking_maf <- smoking_maf[which(smoking_maf$EAF < 0.99),]

summary(smoking_maf$EAF) #perfect.

#Now the INFO, do we have that data?

summary(smoking_maf$INFO) #Freaking perfect. All of them > 0.8

##########################
#REMOVING MHC region SNPs#
##########################

smoking_maf_mhc <- smoking_maf[which(as.numeric(smoking_maf$CHR) == 6 & as.numeric(smoking_maf$BP) >= 26000000 & as.numeric(smoking_maf$BP) <= 34000000),]

summary(as.numeric(smoking_maf_mhc$CHR)) #perfect.
summary(as.numeric(smoking_maf_mhc$BP)) #perfect.

#Now let's check if we had any interesting variants there:

summary(as.numeric(smoking_maf_mhc$P)) #we actually do. But we said we should remove it so...

smoking_maf_no_mhc <- smoking_maf[which(!(smoking_maf$SNP%in%smoking_maf_mhc$SNP)),]

dim(smoking_maf)[1] - dim(smoking_maf_no_mhc)[1] #PERFECT. The numbers match.

######################
#Let's check the RSID#
######################

smoking_maf_no_mhc <- smoking_maf_no_mhc[order(smoking_maf_no_mhc$SNP),]

head(smoking_maf_no_mhc)
tail(smoking_maf_no_mhc)

##########
#CAREFUL!#
##########

#We do have weird SNPs. All the freaking exposure are the same, man.
#Alright, let's check which ones have RSIDs and which ones do not.

rsid_index <- which(str_detect(smoking_maf_no_mhc$SNP, "rs") == TRUE)
no_rsid_index <- which(str_detect(smoking_maf_no_mhc$SNP, "rs") == FALSE)

#Only 205 that are not found!
#This one is gonna be easy, then.

smoking_RSID <- smoking_maf_no_mhc[rsid_index,]

#Let's check that this has worked:

head(smoking_RSID) #perfect.

#Now let' do the other one:

smoking_NO_RSID <- smoking_maf_no_mhc[no_rsid_index,]

######################
#Let's check that out#
######################

View(smoking_NO_RSID) #The code detected them properly, as expected.

#Let's get those RSIDs!

chr <- paste("chr", smoking_NO_RSID$CHR, sep = "")
smoking_NO_RSID$chr_pos <- paste(chr, smoking_NO_RSID$BP, sep = ":")

head(smoking_NO_RSID) #perfect.

#Let's match them in phenoscanner

results_1 <- phenoscanner::phenoscanner(smoking_NO_RSID$chr_pos[seq(1,100)])
results_2 <- phenoscanner::phenoscanner(smoking_NO_RSID$chr_pos[seq(101,200)])
results_3 <- phenoscanner::phenoscanner(smoking_NO_RSID$chr_pos[seq(201,205)])

results_data <- rbind(results_1$snps, results_2$snps, results_3$snps) #0 variables!!!

#Phenoscanner is DEAD.

#########################################
#I am afraid they are not in build 37...#
#########################################

#Let's check the ones with RSID:

head(smoking_RSID)

#I checked them in dbSNP. THEY are in build 37!
#It is just that phenoscanner is NOT working today.

#We do not have a choice. 
#Let's try with SNPNexus

#############################
#Matching data with SNPNexus#
#############################

smoking_NO_RSID_4_nexus <- smoking_NO_RSID

smoking_NO_RSID_4_nexus$Chromosome <- "Chromosome"

smoking_NO_RSID_4_nexus_clean <- smoking_NO_RSID_4_nexus %>%
  select(Chromosome, CHR, BP, EFFECT_ALLELE, OTHER_ALLELE)

smoking_NO_RSID_4_nexus_clean$value <- 1

write.table(smoking_NO_RSID_4_nexus_clean, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/Lifetime_smk/Lifetime_Smk_4_Nexus.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

################################################
#Let's take those SNPs that have been recovered#
################################################

recovered_snps <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/Lifetime_smk/Lifetime_Smk_in_Nexus.txt")

#Not all of them have RSIDs...

recovered_snps <- recovered_snps[which(recovered_snps$dbSNP != "None"),]

#We need to match properly the data:

recovered_snps$chr_pos <- paste(recovered_snps$Chromosome, recovered_snps$Position, sep = ":")

smoking_maf_no_mhc$chr_pos <- paste(smoking_maf_no_mhc$CHR, smoking_maf_no_mhc$BP, sep = ":")

smoking_recovered <- smoking_maf_no_mhc[which(smoking_maf_no_mhc$chr_pos%in%recovered_snps$chr_pos),]

#We have duplicates in recovered_snps. 
#As long as the RSIDs are the same, we do not care of removing them.
#Let's check them out:

recovered_snps[which(duplicated(recovered_snps$dbSNP) == TRUE),] #none...
recovered_snps[which(duplicated(recovered_snps$chr_pos) == TRUE),] #1. Let's check differences.

recovered_snps[which(recovered_snps$chr_pos == "10:48698435"),] #none...

#Both SNPs seem pretty unusual. 

#Variation ID        dbSNP Chromosome Position REF Allele ALT Allele (IUPAC)
#1: chr10:48698435:A/G:1  rs879951753         10 48698435          A                  G
#2: chr10:48698435:A/G:1 rs1041404501         10 48698435          A                  G
#Minor Allele Minor Allele Global Frequency     Contig Contig Position   Band
#1:         None                          None GL000098.1          592728 q11.22
#2:         None                          None GL000098.1          592728 q11.22
#chr_pos
#1: 10:48698435
#2: 10:48698435

#Let's check them in dbSNP.
#One of them is merged hence why we can erase it:

recovered_snps <- recovered_snps[which(recovered_snps$dbSNP != "rs1041404501"),]

###########################################
#Finally, let's match the hell out of this#
###########################################

recovered_snps <- recovered_snps[order(match(recovered_snps$chr_pos, smoking_recovered$chr_pos)),]

which(recovered_snps$chr_pos != smoking_recovered$chr_pos)

head(recovered_snps)
head(smoking_recovered)

#So far, so good.

smoking_recovered$SNP <- recovered_snps$dbSNP

#PERFECT.

##############################################################
#Let's try getting some of the SNPs with pulit reference data#
##############################################################

smoking_recovered$chr_pos <- paste("chr", smoking_recovered$chr_pos, sep = "")

smoking_missing <- smoking_NO_RSID[which(!(smoking_NO_RSID$chr_pos%in%smoking_recovered$chr_pos)),] #88, make sense.

##################
#Getting BMI data#
##################

BMI <- fread("C:/Users/zlc436/Desktop/Leisure_Project/BMI_treated.txt")

BMI <- BMI[which(is.na(BMI$INFO) == FALSE),]
BMI <- BMI[which(is.na(BMI$CHR) == FALSE),]

summary(BMI$CHR) #perfect
summary(BMI$INFO) #perfect

#We removed the weird SNPs.

chr_ <- paste("chr", BMI$CHR, sep = "")
BMI$chr_pos <- paste(chr_, BMI$POS, sep = ":")

BMI_match <- BMI[which(BMI$chr_pos%in%smoking_missing$chr_pos),]

#There are none.

#############
#CONCLUSIONS#
#############

#We are gonna save even those that do not have RSID, just in case we can match them somehow.
#In any case, this is done.

chr_ <- paste("chr", smoking_RSID$CHR, sep = "")
smoking_RSID$chr_pos <- paste(chr_, smoking_RSID$BP, sep = ":")

head(smoking_RSID$chr_pos)

smk_end <- rbind(smoking_missing, smoking_recovered, smoking_RSID)

dim(smk_end)
dim(smoking_maf_no_mhc)

fwrite(smk_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/Lifetime_smk/Lifetime_Smoking_Curated.txt")

############################################
#Finally, let's see if we can load the data#
############################################

check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Lifetime_Smoking_Curated.txt")

check <- check[order(check$chr_pos),]

head(check)
tail(check)
