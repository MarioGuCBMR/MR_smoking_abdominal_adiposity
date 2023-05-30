##############
#INTRODUCTION#
##############

#This is a quick code to check how much do we need to clean the hip data.

###########
#Libraries#
###########

library(data.table)
library(tidyverse)
library(TwoSampleMR)

##################
#Loading function#
##################

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

hip <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/RAW_Data/GIANT_2015_HIP_COMBINED_EUR.txt") 

#Let's check what columns do we have:

hip <- hip[order(hip$MarkerName),]

head(hip, 1000)

#Just checked FTO on PhenoScanner. 
#Everything seems to match up!

#Here we just have the chromsomes and the positions...
#But we do not have the MAF so we will remove them.

##############################
#Filtering data and saving it#
##############################

#What we are going to do is the following:

#a) Filter for MAF/INFO if possible.
#b) remove MHC region SNPs here, if possible. If not, in the outcome.
#c) Save the data so we do not have to wait for the long query in the database.

#The hip data has some of the SNPs without allele frequency.
#We cannot know if the imputation was done correctly...
#Hence, we are gonna remove them.

hip_maf <- hip[which(hip$FreqAllele1HapMapCEU > 0.01),]
hip_maf <- hip_maf[which(hip_maf$FreqAllele1HapMapCEU < 0.99),]

summary(hip_maf$FreqAllele1HapMapCEU) #perfect.

#Now the INFO, do we have that data?

colnames(hip_maf) #we do not. 

##########################
#REMOVING MHC region SNPs#
##########################

#WAIT. We cannot do this because we use the filters from build 37 from Warrington et al and we do not have CHR and POS.
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

BMI_match <- BMI[which(BMI$SNP%in%hip_maf$MarkerName),] #most of them:

hip_Pulit <- hip_maf[which(hip_maf$MarkerName%in%BMI_match$SNP),]

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

BMI_match <- BMI_match[order(match(BMI_match$SNP, hip_Pulit$MarkerName)),]

#Let's check if we did it properly:

which(BMI_match$SNP != hip_Pulit$MarkerName) #perfect.

head(BMI_match$SNP)
head(hip_Pulit$MarkerName)

#We can match the chr_pos:

hip_Pulit$chr_pos_37 <- BMI_match$chr_pos

#########
#PERFECT#
#########

#Let's go and get those missing:

hip_missing <- hip_maf[which(!(hip_maf$MarkerName%in%BMI_match$SNP)),]

#Let's see how many of them can be recovered from SNP Nexus:

for_nexus <- hip_missing

for_nexus$dbsnp <- "dbsnp"

for_nexus <- for_nexus %>%
  select(dbsnp, MarkerName)

for_nexus <- for_nexus[order(for_nexus$MarkerName),]

head(for_nexus)
tail(for_nexus)

########################################################
#Important we have three SNP that are in another format#
########################################################

#I manually transformed them in liftover since there are so few.
#Let's check them:

check_3_snps <- as.data.frame(read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WC/check_3_SNPs_in_liftover.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

#We are gonna use liftover to transform them.
#Though technically there is no need to,
#They are the same as we had in hipAdjBMI combined.

#chr1:161981977-161981979
#chr6:13159524-13159526
#chr6:26020935-26020937

gotcha_vect <- c("chr1:161981978", "chr6:13159525", "chr6:26020936")

which(BMI$chr_pos%in%gotcha_vect) #can we recover them with BMI? YES.

BMI_gotcha <- BMI[which(BMI$chr_pos%in%gotcha_vect),]

###############################################
#Let's prepare the data for getting these snps#
###############################################

hip_missing <- hip_missing[order(hip_missing$MarkerName),]

hip_missing_gotcha <- hip_missing[seq(1,3),]

#After checking, the SNP match. 
#Hence:

hip_missing_gotcha$MarkerName <- BMI_gotcha$SNP

hip_missing_gotcha$chr_pos_37 <- BMI_gotcha$chr_pos

#########
#PERFECT#
#########

#With that being settled..., we only need to run SNPNexus:

head(for_nexus)

for_nexus <- for_nexus[-seq(1,3),]

head(for_nexus) #PERFECT:

write.table(for_nexus, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/HIP/hip_combined_4_Nexus.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

####################
#Let's recover them#
####################

recovered_snps <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/hip/hip_combined_in_Nexus.txt")

recovered_snps <- recovered_snps[order(recovered_snps$`REF Allele`),]

head(recovered_snps) #one deletion
tail(recovered_snps) #nothing.

#It seems we only need to remove the deletion. 
#To do so, since we have IUPAC alleles we will just take all the alleles
#that are letters, to avoid insertions and deletions.

good_alleles <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K",
                  "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U",
                  "V", "W", "X", "Y", "Z")

recovered_snps <- recovered_snps[which(recovered_snps$`REF Allele`%in%good_alleles),] #19573->19572. PERFECT.
recovered_snps <- recovered_snps[which(recovered_snps$`ALT Allele (IUPAC)`%in%good_alleles),] #19572->19572. PERFECT.

##########

hip_recovered <- hip_missing[which(hip_missing$MarkerName%in%recovered_snps$dbSNP),]

recovered_snps <- recovered_snps[order(match(recovered_snps$dbSNP, hip_recovered$MarkerName)),]

#Doing some checkity checks:

which(recovered_snps$dbSNP != hip_recovered$MarkerName)

head(recovered_snps$dbSNP)
head(hip_recovered$MarkerName)

#PERFECT:

chr_ <- paste("chr", recovered_snps$Chromosome, sep = "")

recovered_snps$chr_pos <- paste(chr_, recovered_snps$Position, sep = ":")

hip_recovered$chr_pos_37 <- recovered_snps$chr_pos

###############################
#Let's check the missing ones:#
###############################

hip_missing_2 <- hip_missing[which(!(hip_missing$MarkerName%in%recovered_snps$dbSNP)),]

View(hip_missing_2)

#The first three were already recovered, so let's remove them:

hip_missing_2 <- hip_missing_2[-seq(1,3),]

head(hip_missing_2)

##############################################################################################################
#We did liftover of many of these with WHR data that had CHR and POS so maybe we can retrieve them from there#
##############################################################################################################

whr_curated <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHR/WHR_combined_Curated.txt")

whr_liftovered <- whr_curated[which(whr_curated$MarkerName%in%hip_missing_2$MarkerName),] #2463/2469. SO CLOSE!

head(whr_liftovered)
tail(whr_liftovered) #not all of them could be recovered from liftover...

whr_liftovered <- whr_liftovered[which(whr_liftovered$chr_pos_37 != "-"),]

hip_missing_2_matched <- hip_missing_2[which(hip_missing_2$MarkerName%in%whr_liftovered$MarkerName),]

#Now we reorder them and get the chr and position:

hip_missing_2_matched <- hip_missing_2_matched[order(match(hip_missing_2_matched$MarkerName, whr_liftovered$MarkerName)),]

which(hip_missing_2_matched$MarkerName != whr_liftovered$MarkerName) #perfect!!

head(hip_missing_2_matched$MarkerName) #perfect
head(whr_liftovered$MarkerName) #perfect

hip_missing_2_matched$chr_pos_37 <- whr_liftovered$chr_pos_37

#COOL! Now let's try to solve the ones that we don't have:

################################################
#Let's try phenoscanner and the merged SNP list#
################################################

hip_missing_3 <- hip_missing_2[which(!(hip_missing_2$MarkerName%in%hip_missing_2_matched$MarkerName)),]

#175 SNPs! Not bad

results_1 <- phenoscanner::phenoscanner(hip_missing_3$MarkerName[seq(1,100)])
results_2 <- phenoscanner::phenoscanner(hip_missing_3$MarkerName[seq(101,175)])

results <- rbind(results_1$snps, results_2$snps) #we got only 0 SNPs back. Usually we recovered like 5, but it is not working now. 

hip_missing_3_matched <- hip_missing_3[which(hip_missing_3$MarkerName%in%results$snp),]

#Let's try to match them:

hip_missing_3_matched <- hip_missing_3_matched[order(match(hip_missing_3_matched$MarkerName, results$snp)),]

which(hip_missing_3_matched$MarkerName != results$snp)

head(hip_missing_3_matched$MarkerName)
head(results$snp)

hip_missing_3_matched$chr_pos_37 <- results$hg19_coordinates

###################################################
#We will do our best in the CAUSE pipeline then...#
###################################################

hip_missing_4 <- hip_missing_3[which(!(hip_missing_3$MarkerName%in%results$snp)),]

hip_missing_4$chr_pos_37 <- "-" #to include all variants just in case.

#############
#CONCLUSIONS#
#############

#We solved most of the SNPs, except for 175.
#The hip_missing_3_matched did not work for hip.

hip_end <- rbind(hip_Pulit, hip_recovered, hip_missing_gotcha, hip_missing_2_matched, hip_missing_4) #PERFECT number. It matches the data perfectily.

dim(hip_end)
dim(hip_maf)

fwrite(hip_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/HIP/Hip_combined_Curated.txt")

############################################
#Finally, let's see if we can load the data#
############################################

check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/hip/hip_combined_Curated.txt")
