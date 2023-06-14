##############
#INTRODUCTION#
##############

#This is a quick code to check how much do we need to clean education initiation data.
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

#Here we are gonna use: college or university degree as proxy for education.

education <- extract_instruments(outcomes = "ukb-b-16489", p1 = 0.005, clump = FALSE) #Education

#It removes duplicates...
#That is so frustrating... can we not turn it off??

#We cannot. Oh, welp. What can we do.

summary(education$pval.exposure) #the filtering has worked just fine.

##############################
#Filtering data and saving it#
##############################

#What we are going to do is the following:

#a) Filter for MAF/INFO if possible.
#b) remove MHC region SNPs here, if possible. If not, in the outcome.
#c) Save the data so we do not have to wait for the long query in the database.

education_maf <- education[which(education$eaf.exposure > 0.01),]
education_maf <- education_maf[which(education_maf$eaf.exposure < 0.99),]

summary(education_maf$eaf.exposure) #perfect.

#Now the INFO, do we have that data?

colnames(education_maf) #we do not. 

##########################
#REMOVING MHC region SNPs#
##########################

education_maf_mhc <- education_maf[which(as.numeric(education_maf$chr.exposure) == 6 & as.numeric(education_maf$pos.exposure) >= 26000000 & as.numeric(education_maf$pos.exposure) <= 34000000),]

summary(as.numeric(education_maf_mhc$chr.exposure)) #perfect.
summary(as.numeric(education_maf_mhc$pos.exposure)) #perfect.

#Now let's check if we had any interesting variants there:

summary(as.numeric(education_maf_mhc$pval.exposure)) #we actually do. But we said we should remove it so...

education_maf_no_mhc <- education_maf[which(!(education_maf$SNP%in%education_maf_mhc$SNP)),]

#Smol check:

dim(education_maf)[1] - dim(education_maf_no_mhc)[1] #PERFECT.

######################
#Let's check the RSID#
######################

education_maf_no_mhc <- education_maf_no_mhc[order(education_maf_no_mhc$SNP),]

head(education_maf_no_mhc)
tail(education_maf_no_mhc)

##########
#CAREFUL!#
##########

#We do have weird SNPs. All the freaking exposure are the same, man.
#Alright, let's check which ones have RSIDs and which ones do not.

rsid_index <- which(str_detect(education_maf_no_mhc$SNP, "rs") == TRUE)
no_rsid_index <- which(str_detect(education_maf_no_mhc$SNP, "rs") == FALSE)

#Only 73 that are not found!
#This one is gonna be easy, then.

education_RSID <- education_maf_no_mhc[rsid_index,]

#Let's check that this has worked:

head(education_RSID) #perfect.

#Now let' do the other one:

education_NO_RSID <- education_maf_no_mhc[no_rsid_index,]

######################
#Let's check that out#
######################

View(education_NO_RSID) #The code detected them properly, as expected.

#Let's get those RSIDs!

chr <- paste("chr", education_NO_RSID$chr.exposure, sep = "")
education_NO_RSID$chr_pos <- paste(chr, education_NO_RSID$pos.exposure, sep = ":")

head(education_NO_RSID) #perfect.

#Let's match them in phenoscanner

results <- phenoscanner::phenoscanner(education_NO_RSID$chr_pos)

results_data <- results$snps #0 variables!!!

#########################################
#I am afraid they are not in build 37...#
#########################################

#Let's check the ones with RSID:

head(education_RSID)

#I checked them in dbSNP. THEY are in build 37!
#It is just that phenoscanner does not have them.

#We do not have a choice. 
#Let's try with SNPNexus

#############################
#Matching data with SNPNexus#
#############################

education_NO_RSID_4_nexus <- education_NO_RSID

education_NO_RSID_4_nexus$Chromosome <- "Chromosome"

education_NO_RSID_4_nexus_clean <- education_NO_RSID_4_nexus %>%
  select(Chromosome, chr.exposure, pos.exposure, effect_allele.exposure, other_allele.exposure)

education_NO_RSID_4_nexus_clean$value <- 1

write.table(education_NO_RSID_4_nexus_clean, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/education_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/SmkInit/Ever_Never_4_Nexus.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

#######
#JESUS#
#######

#None of them have RSIDs!!!! Not even in SNPNexus.
#Then, there is nothing that we can do...
#Unless they are matched by chromsome and position in the outcome.
#And we can retrieve their data there.

#So far, so good.

#############
#CONCLUSIONS#
#############

#We can save the data as it is.

fwrite(education_maf_no_mhc, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/education_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/SmkInit/Ever_Never_Smoker_Curated.txt")

############################################
#Finally, let's see if we can load the data#
############################################

check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/education_WHR/CODE_4_Replication/2SMR/Curating_data/Ever_Never_Smoker_Curated.txt")
