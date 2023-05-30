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

smk_data<-read.table("C:/Users/zlc436/Downloads/ukb-b-469.vcf.gz", stringsAsFactors = FALSE)

#We are also gonna load the genome-wide significant SNPs to be able to see if the data
#matches and check whether we have weird results or not with this format.

smk_check <- extract_instruments(outcomes = "ukb-b-469", p1 = 0.00000005, clump = FALSE) #369

#############
#Small check#
#############

smk_check[1,]

#beta.exposure se.exposure pval.exposure pos.exposure samplesize.exposure chr.exposure
#1     0.0359918  0.00650565       3.2e-08     42545827               33229            8
#id.exposure       SNP effect_allele.exposure other_allele.exposure eaf.exposure
#1   ukb-b-469 rs9792277                      A                     G     0.776236
#exposure
#1 Number of cigarettes currently smoked daily (current cigarette smokers) || id:ukb-b-469
#mr_keep.exposure pval_origin.exposure data_source.exposure
#1             TRUE             reported                  igd

head(smk_data) #the RSID is in V3

smk_data[which(smk_data$V3 == "rs9792277"),]

#V1       V2        V3 V4 V5 V6   V7          V8             V9
#3949864  8 42545827 rs9792277  G  A  . PASS AF=0.776236 ES:SE:LP:AF:ID
#V10
#3949864 0.0359918:0.00650565:7.49485:0.776236:rs9792277

#It seems that the effect allele is the second one!!
#I just checked in phenoscanner. And exactly. 

#Let's see if I can transform and de transform the p-value:

abs(log10(smk_check$pval.exposure[which(smk_check$SNP == "rs11265238")])) 

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

smk_data_copy <- smk_data

colnames(smk_data_copy) <- c("chr.exposure", "pos.exposure", "SNP", "other_allele.exposure", "effect_allele.exposure", "dot", "PASS", "AF", "Headers", "Headers_info")

#Now we are gonna make functions to retrieve the info on headers info:

smk_data_copy$beta.exposure <- as.numeric(as.character(unlist(sapply(smk_data_copy$Headers_info, beta_retriever))))
smk_data_copy$se.exposure <- as.numeric(as.character(unlist(sapply(smk_data_copy$Headers_info, se_retriever))))
smk_data_copy$logp <- as.numeric(as.character(unlist(sapply(smk_data_copy$Headers_info, logp_retriever))))
smk_data_copy$eaf.exposure <- as.numeric(as.character(unlist(sapply(smk_data_copy$Headers_info, EAF_retriever))))
smk_data_copy$rsid <- as.character(unlist(sapply(smk_data_copy$Headers_info, rsid_retriever)))

#To convert to p-value I need to transform them again:

smk_data_copy$logp_neg <- smk_data_copy$logp*(-1)

#And thus...

smk_data_copy$pval.exposure <- 10^(smk_data_copy$logp_neg)

############
#FINISHED!!#
############

####################################################################
#Now we start the curation process, let's do some checks beforehand#
####################################################################

#1. Are they OR or logODDs??

summary(smk_data_copy$beta.exposure) #No 1s, so no Odd Ratios. Moreover we have SE so...

#2. There is a filter with a dot. What is that?

head(smk_data_copy$dot)
tail(smk_data_copy$dot)

length(which(smk_data_copy$dot == ".")) #all of them.

#something from the vcf file I guess.

#.3 We have a column that says: PASS.

length(which(smk_data_copy$PASS == "PASS")) #all of them.

#That means that it was a way to retrieve the good SNPs for the authors.

#4. Let's check again if the info is correct:

smk_check[1,]

#beta.exposure se.exposure pval.exposure pos.exposure samplesize.exposure chr.exposure
#1     0.0359918  0.00650565       3.2e-08     42545827               33229            8
#id.exposure       SNP effect_allele.exposure other_allele.exposure eaf.exposure
#1   ukb-b-469 rs9792277                      A                     G     0.776236
#exposure
#1 Number of cigarettes currently smoked daily (current cigarette smokers) || id:ukb-b-469
#mr_keep.exposure pval_origin.exposure data_source.exposure
#1             TRUE             reported                  igd

head(smk_data) #the RSID is in V3

smk_data_copy[which(smk_data_copy$SNP == "rs9792277"),]

#chr.exposure pos.exposure       SNP other_allele.exposure effect_allele.exposure
#3949864            8     42545827 rs9792277                     G                      A
#dot PASS          AF        Headers
#3949864   . PASS AF=0.776236 ES:SE:LP:AF:ID
#Headers_info beta.exposure se.exposure    logp
#3949864 0.0359918:0.00650565:7.49485:0.776236:rs9792277     0.0359918  0.00650565 7.49485
#eaf.exposure      rsid logp_neg pval.exposure
#3949864     0.776236 rs9792277 -7.49485       3.2e-08

###############
#Perfect match#
###############

#Let's try another one:

smk_check$SNP[2]

#"rs4736835"

smk_check$beta.exposure[2]

#0.0359934

smk_check$effect_allele.exposure[2]

#"C"

smk_check$pval.exposure[2]

#3.09999e-08

#And now we check:

smk_data_copy[which(smk_data_copy$SNP == "rs4736835"),]

##################
#FREAKING PERFECT#
##################

##########
#CURATION#
##########

#We don't have INFO, so we cannot rely on that. 
#1. But we can remove those variants that are MAF < 0.01

smk_data_maf <- smk_data_copy[which(smk_data_copy$eaf.exposure > 0.01),]
smk_data_maf <- smk_data_maf[which(smk_data_maf$eaf.exposure < 0.99),]

#We go from 9.7M to 7.7M. Good enough!!

#2. Now we remove the MHC region:

smk_data_maf_mhc <- smk_data_maf[which(as.numeric(smk_data_maf$chr.exposure) == 6 & as.numeric(smk_data_maf$pos.exposure) >= 26000000 & as.numeric(smk_data_maf$pos.exposure) <= 34000000),]

summary(as.numeric(smk_data_maf_mhc$chr.exposure)) #perfect.
summary(as.numeric(smk_data_maf_mhc$pos.exposure)) #perfect.

#Now let's check if we had any interesting variants there:

summary(as.numeric(smk_data_maf_mhc$pval.exposure)) #we actually do. But we said we should remove it so...

smk_data_maf_no_mhc <- smk_data_maf[which(!(smk_data_maf$SNP%in%smk_data_maf_mhc$SNP)),]

#3. We should have done that before. But let's check whether we are in build 37. I already know that it is, hence why I just check it here.

head(smk_data_maf_no_mhc)

#Perfect.

#4. Let's check the RSIDs:

smk_data_maf_no_mhc <- smk_data_maf_no_mhc[order(smk_data_maf_no_mhc$SNP),]

#We know that some of them are not in the correct format. 
#So we will try to obtain them.

##########
#CAREFUL!#
##########

#We do have weird SNPs. All the freaking exposure are the same, man.
#Alright, let's check which ones have RSIDs and which ones do not.

rsid_index <- which(str_detect(smk_data_maf_no_mhc$SNP, "rs") == TRUE)
no_rsid_index <- which(str_detect(smk_data_maf_no_mhc$SNP, "rs") == FALSE) #we have 531

#Way more than the 73 found with the p< 0.05 that are not found!
#This one is gonna be easy, then.

#In any case, we do not have RSID for these ones.
#But that is not as important.
#We will match with chr_pos. We are gonna check a couple of more no RSID SNPs in dbSNP.
#Just to see what the hell is going on and why there are no RSIDs on them.

smk_data_RSID <- smk_data_maf_no_mhc[rsid_index,]

#Let's check that this has worked:

head(smk_data_RSID) #perfect.

#Now let' do the other one:

smk_data_NO_RSID <- smk_data_maf_no_mhc[no_rsid_index,]

######################
#Let's check that out#
######################

View(smk_data_NO_RSID) #The code detected them properly, as expected.

#Let's get those RSIDs!

chr <- paste("chr", smk_data_NO_RSID$chr.exposure, sep = "")
smk_data_NO_RSID$chr_pos <- paste(chr, smk_data_NO_RSID$pos.exposure, sep = ":")

head(smk_data_NO_RSID, 150) #the first one is found in build 37!
tail(smk_data_NO_RSID) #the first one is found in build 37!

#Let's see in dbsnp... 
#It seems that we can trust the chr_pos!

smk_data_NO_RSID[150,] #ONLY FOUND IN DBSNP FOR BUILD 37. FINAL EVIDENCE MAN.

######################################################################################
#Thus, let's save the data and match as much as possible with chromosome and position#
######################################################################################

#In the whole process we will be very careful if one of the SNPs RSID is not there.
#Because if that is the case we will have a problem when clumping.

fwrite(smk_data_maf_no_mhc, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/CigDayCurrent/Cig_Day_Current_Curated_FULL.txt")

############################################
#Finally, let's see if we can load the data#
############################################

check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/CigDayCurrent/Cig_Day_Current_Curated_FULL.txt")

#And let's check if the data holds up with SNPs at other p-value thresholds:

smk_check_2 <- extract_instruments(outcomes = "ukb-b-469", p1 = 0.00005, clump = FALSE) #369

head(smk_check_2)

smk_check_2[1,]

#beta.exposure se.exposure pval.exposure pos.exposure samplesize.exposure chr.exposure
#1      0.023754  0.00584976   4.90004e-05    110989638               33229            1
#id.exposure        SNP effect_allele.exposure other_allele.exposure eaf.exposure
#1   ukb-b-469 rs12405338                      A                     G     0.321633
#exposure
#1 Number of cigarettes currently smoked daily (current cigarette smokers) || id:ukb-b-469
#mr_keep.exposure pval_origin.exposure data_source.exposure
#1             TRUE             reported                  igd

#And now let's get the data:

check[which(check$SNP == "rs12405338"),]

#chr.exposure pos.exposure        SNP other_allele.exposure effect_allele.exposure dot
#1:            1    110989638 rs12405338                     G                      A   .
#PASS          AF        Headers                                   Headers_info
#1: PASS AF=0.321633 ES:SE:LP:AF:ID 0.023754:0.00584976:4.3098:0.321633:rs12405338
#beta.exposure se.exposure   logp eaf.exposure       rsid logp_neg pval.exposure
#1:      0.023754  0.00584976 4.3098     0.321633 rs12405338  -4.3098  4.900044e-05

#It's a perfect match.