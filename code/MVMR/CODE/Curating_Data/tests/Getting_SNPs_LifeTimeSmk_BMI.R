##############
#INTRODUCTION#
##############

#This code is to get the SNPs for MR-BMA for Smoking Initiation and BMI.

#The number of risk factors will increase for sure.
#And things will change, this is just a test.

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)
library(TwoSampleMR)
library(ieugwasr)
library(jsonlite)
library(httr)

###########
#FUNCTIONS#
###########

ld_proxies_MGU <- function(snp){
  #Enter vector of SNPs it will output all its LD friends from 1KG European LD friends with r2 > 0.8
  #in 500kb.
  
  fake_df <- t(as.data.frame(c("d_prime", "variation2", "population_name", "r2", "variation1")))
  
  colnames(fake_df) <- fake_df[1,]
  rownames(fake_df) <- c(1)
  
  #Setting the server:
  
  server <- "http://grch37.rest.ensembl.org"
  
  for(i in snp){
    
    ext_1 <- paste("/ld/human/", i, sep = "")
    ext_2 <- paste(ext_1, "/1000GENOMES:phase_3:EUR", sep = "")
    
    r <- GET(paste(server, ext_2, sep = ""), content_type("application/json"))
    new <- fromJSON(toJSON(content(r)))
    
    fake_df <- rbind(fake_df, new)
    
  }
  
  #Now filtering for those that are in high LD:
  
  final_df <- fake_df[which(as.numeric(fake_df$r2) > 0.8),] #The NAs by coercion are the rows from the fake_df, ignore them!
  
  return(final_df)
  
}

##############
#LOADING DATA#
##############

#For the risk factors we are only going to include SNPs that are in
#all risk factors, taking into account proxies, of course. 
#Things get a bit cumbersome since we need them all.

#For Life-time Smoker

exposure_dat <- read_exposure_data(
  filename = "C:/Users/zlc436/Desktop/SMK_WHRAdjBMI/SmokingCombined_Wootton.txt",
  sep = " ",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  eaf_col = "EAF",
  pval_col = "P",
  #samplesize_col = "N"
)

exposure_dat_2 <- exposure_dat[which(exposure_dat$pval.exposure <= 0.00000005),]
exposure_dat_2$rsid <- exposure_dat_2$SNP
exposure_dat_2$pval <- exposure_dat_2$pval.exposure
life_time_smk_rf <- ieugwasr::ld_clump(exposure_dat_2)

#And now we get the SNPs that match for the educational one:

edu_out_dat <- extract_outcome_data(
  snps = life_time_smk_rf$SNP,
  outcomes = 'ukb-b-16489',
  proxies = TRUE
)

#This worked perfectly.
#Let's get the harmonized version:

dat_1_smk2edu <- harmonise_data(life_time_smk_rf, edu_out_dat, action = 3)

#Nice, now let's do it the other way around:

dat_1_smk2edu <- dat_1_smk2edu[which(dat_1_smk2edu$palindromic == FALSE),]

#Same for educational attainment.
#Which we don't have, but we will use Age of finished education as a proxy.

edu_dat <- extract_instruments(outcomes = "ukb-b-16489", p1 = 0.05, clump = FALSE) #Age completed full time education 
edu_dat_2 <- edu_dat[which(edu_dat$pval.exposure <= 0.00000005),]
edu_dat_2$rsid <- edu_dat_2$SNP
edu_dat_2$pval <- edu_dat_2$pval.exposure
ea_rf <- ieugwasr::ld_clump(edu_dat_2)

life_time_smk_out_dat <- read_outcome_data(
  snps = ea_rf$SNP,
  filename = "C:/Users/zlc436/Desktop/SMK_WHRAdjBMI/SmokingCombined_Wootton.txt",
  sep = " ",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  eaf_col = "EAF",
  pval_col = "P",
  #samplesize_col = "N"
)

dat_1_edu2smk <- harmonise_data(ea_rf, life_time_smk_out_dat, action = 3)

#Nice, now let's do it the other way around:

dat_1_edu2smk <- dat_1_edu2smk[which(dat_1_edu2smk$palindromic == FALSE),]

#We also want for physical activity and tdi intake.
#So let's go!

#The first thing we need to do is updating this...
#We only need beta, se and pval for each. 
#We need to be wary.
#Let's go:

dat_1_smk2edu_clean <- dat_1_smk2edu %>%
  select(SNP, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, pval.exposure, beta.outcome, se.outcome, pval.outcome)

colnames(dat_1_smk2edu_clean) <- c("SNP", "effect_allele", "other_allele", "beta.smk", "se.smk", "pval.smk", "beta.ea", "se.ea", "pval.ea")

dat_1_edu2smk_clean <- dat_1_edu2smk %>%
  select(SNP, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, pval.exposure, beta.outcome, se.outcome, pval.outcome)

colnames(dat_1_edu2smk_clean) <- c("SNP", "effect_allele", "other_allele", "beta.ea", "se.ea", "pval.ea", "beta.smk", "se.smk", "pval.smk")

smk_edu_final_df <- rbind(dat_1_smk2edu_clean, dat_1_edu2smk_clean)

#Now we just need to get the independent ones:

smk_edu_independent_end <- smk_edu_final_df %>%
  select(SNP, pval.smk)

colnames(smk_edu_independent_end) <- c("rsid", "pval")

#We won't care about pval now. There are mixed so... we ignore it.

independent_end <- ieugwasr::ld_clump(smk_edu_independent_end, clump_p = 0.05)

#There were 30 in high LD:

smk_edu_final_df <- smk_edu_final_df[which(smk_edu_final_df$SNP%in%independent_end$rsid),]

fwrite(smk_edu_final_df, "N:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/Multivariate_Draft/CODE/Multivariate_Mario/OUTPUT/LifeTime_Edu_Risk_Factors.txt")

cool <- fread("N:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/Multivariate_Draft/CODE/Multivariate_Mario/OUTPUT/LifeTime_Edu_Risk_Factors.txt")

#This is perfect.
#Let's go for TDI

tdi_dat <- extract_instruments(outcomes = "ukb-b-10011", p1 = 0.05, clump = FALSE) #tdi intake 
tdi_dat_2 <- tdi_dat[which(tdi_dat$pval.exposure <= 0.00000005),]
tdi_dat_2$rsid <- tdi_dat_2$SNP
tdi_dat_2$pval <- tdi_dat_2$pval.exposure
tdi_rf <- ieugwasr::ld_clump(tdi_dat_2)

#Lots of them, actually.
#Let's check how many of these are independent from our main set:

#Perfect, now we need those in the other exposures:

life_time_smk_out_dat <- read_outcome_data(
  snps = tdi_rf$SNP,
  filename = "C:/Users/zlc436/Desktop/SMK_WHRAdjBMI/SmokingCombined_Wootton.txt",
  sep = " ",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  eaf_col = "EAF",
  pval_col = "P",
  #samplesize_col = "N"
)

ea_out_dat <- extract_outcome_data(
  snps = tdi_rf$SNP,
  outcomes = 'ukb-b-16489',
  proxies = TRUE
)

#This is perfect. Let's merge them in the same df. 
#To do so, we will align all data to tdi data.

smk_tdi <- harmonise_data(tdi_rf, life_time_smk_out_dat, action = 3)

ea_tdi <- harmonise_data(tdi_rf, ea_out_dat, action = 3)

#Perfect, now we can take the data from smk_tdi:

final_tdi_gwi_df <- smk_tdi %>%
  select(SNP, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, pval.exposure, beta.outcome, se.outcome, pval.outcome)

colnames(final_tdi_gwi_df) <- c("SNP", "effect_allele", "other_allele", "beta.tdi", "se.tdi", "pval.tdi", "beta.smk", "se.smk", "pval.smk")

#Let's do one small check:

which(final_tdi_gwi_df$effect_allele != ea_tdi$effect_allele.exposure) #As expected, they all match. Awesome.

final_tdi_gwi_df$beta.ea <- ea_tdi$beta.outcome
final_tdi_gwi_df$se.ea <- ea_tdi$pval.outcome
final_tdi_gwi_df$pval.ea <- ea_tdi$pval.outcome

#Now we can merge it?

#We cannot yet. Because we need to add to cool the effects of those SNP.

tdi_out_dat <- extract_outcome_data(
  snps = cool$SNP,
  outcomes = 'ukb-b-10011',
  proxies = TRUE
)

#In this case, even with proxies, there is one SNP that is not found.
#Hence we are removing it from the final list:

cool <- cool[which(cool$SNP%in%tdi_out_dat$SNP),]

#Oh, that means that we probably have a duplicated in the original df:

cool <- cool[which(duplicated(cool$SNP) == FALSE),] #exactly!

tdi_out_dat <- tdi_out_dat[order(match(tdi_out_dat$SNP, cool$SNP)),]

which(tdi_out_dat$SNP != cool$SNP)

#Once you get the gist of it, this is easy peasy lemon squeazy.

new_beta <- ifelse(tdi_out_dat$effect_allele.outcome != cool$effect_allele, tdi_out_dat$beta.outcome*(-1), tdi_out_dat$beta.outcome)
new_a1 <- ifelse(tdi_out_dat$effect_allele.outcome != cool$effect_allele, tdi_out_dat$other_allele.outcome, tdi_out_dat$effect_allele.outcome)

#Small check:

which(new_a1 != cool$effect_allele) #Perfect.

cool$beta.tdi <- new_beta
cool$se.tdi <- tdi_out_dat$se.outcome
cool$pval.tdi <- tdi_out_dat$pval.outcome

#And now, we can merge both dataframes.

final_smk_ea_tdi_df <- rbind(cool, final_tdi_gwi_df)

#And now we have to be wary because we have to make them independent.

independent_check <- final_smk_ea_tdi_df %>%
  select(SNP, pval.smk)

colnames(independent_check) <- c("rsid", "pval")

independent_snps <- ieugwasr::ld_clump(independent_check, clump_p = 0.05)

#Lots of them were removed. As expected!
#These are the ones that we should keep:

final_smk_ea_tdi_df_independent  <- final_smk_ea_tdi_df[which(final_smk_ea_tdi_df$SNP%in%independent_snps$rsid),]

fwrite(final_smk_ea_tdi_df_independent, "N:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/Multivariate_Draft/CODE/Multivariate_Mario/OUTPUT/LifeTime_Edu_TDI_Risk_Factors.txt")

