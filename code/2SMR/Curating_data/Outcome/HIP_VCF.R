##############
#INTRODUCTION#
##############

#This code is to curated Ever-smoker data from the downloaded vcf file.

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)
library(TwoSampleMR)

###################
#Loading functions#
###################

beta_retriever <- function(header_info){
  
  beta <- strsplit(header_info, ":")[[1]][1]
  
  return(beta)
  
}

se_retriever <- function(header_info){
  
  se <- strsplit(header_info, ":")[[1]][2]
  
  return(se)
  
}

logp_retriever <- function(header_info){
  
  logp <- strsplit(header_info, ":")[[1]][3]
  
  return(logp)
  
}

EAF_retriever <- function(header_info){
  
  EAF <- strsplit(header_info, ":")[[1]][4]
  
  return(EAF)
  
}


rsid_retriever <- function(header_info){
  
  rsid <- strsplit(header_info, ":")[[1]][5]
  
  return(rsid)
  
}

##############
#Loading data#
##############

hip_data<-read.table("C:/Users/zlc436/Downloads/ukb-b-15590.vcf.gz", stringsAsFactors = FALSE)

#We are also gonna load the genome-wide significant SNPs to be able to see if the data
#matches and check whether we have weird results or not with this format.

hip_check <- extract_instruments(outcomes = "ukb-b-15590", p1 = 0.00000005, clump = FALSE) #40892

#############
#Small check#
#############

hip_check[1,]

#beta.exposure pval.exposure samplesize.exposure pos.exposure chr.exposure se.exposure
#1     0.0182154   5.99998e-10              462117    209490527            1  0.00294241
#id.exposure       SNP effect_allele.exposure other_allele.exposure eaf.exposure
#1 ukb-b-15590 rs1395753                      T                     A     0.131647
#exposure mr_keep.exposure pval_origin.exposure
#1 Hip circumference || id:ukb-b-15590             TRUE             reported
#data_source.exposure
#1                  igd

head(hip_data) #the RSID is in V3

hip_data[which(hip_data$V3 == "rs1395753"),]

#V1        V2        V3 V4 V5 V6   V7          V8             V9
#611271  1 209490527 rs1395753  A  T  . PASS AF=0.131647 ES:SE:LP:AF:ID
#V10
#611271 0.0182154:0.00294241:9.22185:0.131647:rs1395753

#It seems that the effect allele is the second one!!
#I just checked in phenoscanner. And exactly. 

#Let's see if I can transform and de transform the p-value:

abs(log10(hip_check$pval.exposure[which(hip_check$SNP == "rs1395753")]))  #perfect

#Let's play now how to detransform this:

log10(0.00000005) #-7.3

10^-7.30 #AHA.

#The important thing was the minus.
#That is why...
#We might have an issue here. 

log10(0.05) #it is always gonna be negative.

log10(0.85) #always gonna be negative.

#Hence why people often do the -log10.

#######################
#REFORMATTING THE DATA#
#######################

#Now that we know the basics of the data.
#Let's start by making some functions to make it work:

hip_data_copy <- hip_data

colnames(hip_data_copy) <- c("chr.exposure", "pos.exposure", "SNP", "other_allele.exposure", "effect_allele.exposure", "dot", "PASS", "AF", "Headers", "Headers_info")

#Now we are gonna make functions to retrieve the info on headers info:

hip_data_copy$beta.exposure <- as.numeric(as.character(unlist(sapply(hip_data_copy$Headers_info, beta_retriever))))
hip_data_copy$se.exposure <- as.numeric(as.character(unlist(sapply(hip_data_copy$Headers_info, se_retriever))))
hip_data_copy$logp <- as.numeric(as.character(unlist(sapply(hip_data_copy$Headers_info, logp_retriever))))
hip_data_copy$eaf.exposure <- as.numeric(as.character(unlist(sapply(hip_data_copy$Headers_info, EAF_retriever))))
hip_data_copy$rsid <- as.character(unlist(sapply(hip_data_copy$Headers_info, rsid_retriever)))

#To convert to p-value I need to transform them again:

hip_data_copy$logp_neg <- hip_data_copy$logp*(-1)

#And thus...

hip_data_copy$pval.exposure <- 10^(hip_data_copy$logp_neg)

############
#FINISHED!!#
############

####################################################################
#Now we start the curation process, let's do some checks beforehand#
####################################################################

#1. Are they OR or logODDs??

summary(hip_data_copy$beta.exposure) #No 1s, so no Odd Ratios. Moreover we have SE so...

#2. There is a filter with a dot. What is that?

head(hip_data_copy$dot)
tail(hip_data_copy$dot)

length(which(hip_data_copy$dot == ".")) #all of them.

#something from the vcf file I guess.

#.3 We have a column that says: PASS.

length(which(hip_data_copy$PASS == "PASS")) #all of them.

#That means that it was a way to retrieve the good SNPs for the authors.

#4. Let's check again if the info is correct:

hip_check[1,]

#beta.exposure pval.exposure samplesize.exposure pos.exposure chr.exposure se.exposure
#1     0.0182154   5.99998e-10              462117    209490527            1  0.00294241
#id.exposure       SNP effect_allele.exposure other_allele.exposure eaf.exposure
#1 ukb-b-15590 rs1395753                      T                     A     0.131647
#exposure mr_keep.exposure pval_origin.exposure
#1 Hip circumference || id:ukb-b-15590             TRUE             reported
#data_source.exposure
#1                  igd

head(hip_data) #the RSID is in V3

hip_data_copy[which(hip_data_copy$SNP == "rs1395753"),]

#chr.exposure pos.exposure       SNP other_allele.exposure effect_allele.exposure dot
#611271            1    209490527 rs1395753                     A                      T   .
#PASS          AF        Headers                                    Headers_info
#611271 PASS AF=0.131647 ES:SE:LP:AF:ID 0.0182154:0.00294241:9.22185:0.131647:rs1395753
#beta.exposure se.exposure    logp eaf.exposure      rsid logp_neg pval.exposure
#611271     0.0182154  0.00294241 9.22185     0.131647 rs1395753 -9.22185  5.999983e-10

###############
#Perfect match#
###############

#Let's try another one:

hip_check$SNP[2]

#"rs6681466"

hip_check$beta.exposure[2]

#0.0175975

hip_check$effect_allele.exposure[2]

#"T"

hip_check$pval.exposure[2]

#1.6e-09

#And now we check:

hip_data_copy[which(hip_data_copy$SNP == "rs6681466"),]

##################
#FREAKING PERFECT#
##################

##########
#CURATION#
##########

#We don't have INFO, so we cannot rely on that. 
#1. But we can remove those variants that are MAF < 0.01

hip_data_maf <- hip_data_copy[which(hip_data_copy$eaf.exposure > 0.01),]
hip_data_maf <- hip_data_maf[which(hip_data_maf$eaf.exposure < 0.99),]

#We go from 9.2M to 7.7M. Good enough!!

#2. Now we remove the MHC region:

hip_data_maf_mhc <- hip_data_maf[which(as.numeric(hip_data_maf$chr.exposure) == 6 & as.numeric(hip_data_maf$pos.exposure) >= 26000000 & as.numeric(hip_data_maf$pos.exposure) <= 34000000),]

summary(as.numeric(hip_data_maf_mhc$chr.exposure)) #perfect.
summary(as.numeric(hip_data_maf_mhc$pos.exposure)) #perfect.

#Now let's check if we had any interesting variants there:

summary(as.numeric(hip_data_maf_mhc$pval.exposure)) #we actually do. But we said we should remove it so...

hip_data_maf_no_mhc <- hip_data_maf[which(!(hip_data_maf$SNP%in%hip_data_maf_mhc$SNP)),]

#3. We should have done that before. But let's check whether we are in build 37. I already know that it is, hence why I just check it here.

head(hip_data_maf_no_mhc)

#Perfect.

#4. Let's check the RSIDs:

hip_data_maf_no_mhc <- hip_data_maf_no_mhc[order(hip_data_maf_no_mhc$SNP),]

head(hip_data_maf_no_mhc)

#We know that some of them are not in the correct format. 
#So we will try to obtain them.

##########
#CAREFUL!#
##########

#We do have weird SNPs. All the freaking exposure are the same, man.
#Alright, let's check which ones have RSIDs and which ones do not.

rsid_index <- which(str_detect(hip_data_maf_no_mhc$SNP, "rs") == TRUE)
no_rsid_index <- which(str_detect(hip_data_maf_no_mhc$SNP, "rs") == FALSE) #we have 531

#We have 531, as with the other...
#This one is gonna be easy, then.

#In any case, we do not have RSID for these ones.
#But that is not as important.
#We will match with chr_pos. We are gonna check a couple of more no RSID SNPs in dbSNP.
#Just to see what the hell is going on and why there are no RSIDs on them.

hip_data_RSID <- hip_data_maf_no_mhc[rsid_index,]

#Let's check that this has worked:

head(hip_data_RSID) #perfect.

#Now let's do the other one:

hip_data_NO_RSID <- hip_data_maf_no_mhc[no_rsid_index,]

######################
#Let's check that out#
######################

View(hip_data_NO_RSID) #The code detected them properly, as expected.

#Let's get those RSIDs!

chr <- paste("chr", hip_data_NO_RSID$chr.exposure, sep = "")
hip_data_NO_RSID$chr_pos <- paste(chr, hip_data_NO_RSID$pos.exposure, sep = ":")

head(hip_data_NO_RSID, 150) #the first one is found in build 37!
tail(hip_data_NO_RSID) #the first one is found in build 37!

#Let's see in dbsnp... 
#It seems that we can trust the chr_pos!

hip_data_NO_RSID[150,] #ONLY FOUND IN DBSNP FOR BUILD 37. FINAL EVIDENCE MAN.

#We can use the chr_pos to match the data, so no issues here.

######################################################################################
#Thus, let's save the data and match as much as possible with chromosome and position#
######################################################################################

#In the whole process we will be very careful if one of the SNPs RSID is not there.
#Because if that is the case we will have a problem when clumping.

fwrite(hip_data_maf_no_mhc, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/HIP/HIP_UKBB_Combined_Curated_FULL.txt")

############################################
#Finally, let's see if we can load the data#
############################################

check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/HIP/HIP_UKBB_Combined_Curated_FULL.txt")

#And let's check if the data holds up with SNPs at other p-value thresholds:

hip_check_2 <- extract_instruments(outcomes = "ukb-b-15590", p1 = 0.00005, clump = FALSE) #369

head(hip_check_2)

hip_check_2[1,]

#pval.exposure beta.exposure se.exposure samplesize.exposure chr.exposure pos.exposure
#1   2.59998e-09      0.010713  0.00179964              462166            1      7727854
#id.exposure       SNP effect_allele.exposure other_allele.exposure eaf.exposure
#1  ukb-b-9405 rs1891215                      C                     T      0.45662
#exposure mr_keep.exposure pval_origin.exposure
#1 Waist circumference || id:ukb-b-9405             TRUE             reported
#data_source.exposure
#1                  igd

#And now let's get the data:

check[which(check$SNP == "rs1891215"),]

#chr.exposure pos.exposure       SNP other_allele.exposure effect_allele.exposure dot PASS
#1:            1      7727854 rs1891215                     T                      C   . PASS
#AF        Headers                                  Headers_info beta.exposure
#1: AF=0.45662 ES:SE:LP:AF:ID 0.010713:0.00179964:8.58503:0.45662:rs1891215      0.010713
#se.exposure    logp eaf.exposure      rsid logp_neg pval.exposure
#1:  0.00179964 8.58503      0.45662 rs1891215 -8.58503   2.59998e-09

#It's a perfect match.