##############
#INTRODUCTION#
##############

#This is the 2SMR code to run Lifetime smoking with WHR (GIANT)

###################
#Loading libraries#
###################

library(TwoSampleMR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gridExtra)
library(rmarkdown)
library(TwoSampleMR)
library(MRInstruments)
library(jsonlite)
library(httr)
library(data.table)
library(readxl)
library(tidyverse)

###################
#Loading functions#
###################

mr_plots <- function(dat)
{
  require(TwoSampleMR)
  require(ggplot2)
  require(ggrepel)
  require(dplyr)
  
  temp <- subset(dat, outcome == outcome[1] & exposure == exposure[1])
  exposure_name <- temp$exposure[1]
  outcome_name <- temp$outcome[1]
  
  if(! "labels" %in% names(dat)) dat$labels <- NA
  
  exposure_units <- temp$units.exposure[1]
  outcome_units <- temp$units.outcome[1]
  
  mrs <- mr_singlesnp(temp, all_method=c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))
  mrl <- mr_leaveoneout(temp)
  
  mrs$index <- 1:nrow(mrs)
  mrl$index <- 1:nrow(mrl)
  
  mrs <- dplyr::arrange(merge(mrs, select(temp, SNP, labels), all.x=TRUE), index)
  mrl <- dplyr::arrange(merge(mrl, select(temp, SNP, labels), all.x=TRUE), index)
  
  mrres <- mr(temp, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))
  
  gridExtra::grid.arrange(
    mr_forest_plot(mrs)[[1]] +
      ggplot2::labs(
        title="a)",
        x=paste0("MR estimate (", outcome_units, " per ", exposure_units, ")")
      ),
    mr_scatter_plot(mrres, temp)[[1]] +
      ggplot2::labs(
        title="b)",
        x=paste0("SNP effect on ", exposure_name),
        y=paste0("SNP effect on ", outcome_name)
      ) +
      geom_label_repel(ggplot2::aes(label=labels), point.padding = unit(0.5, "lines")),
    mr_leaveoneout_plot(mrl)[[1]] +
      ggplot2::labs(
        title="c)",
        x=paste0("MR estimate (", outcome_units, " per ", exposure_units, ")"),
        y="Excluded variant"
      ),
    mr_funnel_plot(mrs)[[1]] +
      ggplot2::labs(title="d)") +
      ggplot2::theme(legend.position="none") +
      ggrepel::geom_label_repel(ggplot2::aes(label=labels), point.padding = unit(0.5, "lines")),
    ncol=2
  )
}

parse_rsid <- function(rsid){
  
  rsid_ <- strsplit(rsid, ":")[[1]][1] #this works for those that do not present the weird formatting :)
  
  return(rsid_)
  
}

######################################################
#SECTION A: FINDING INDEPENDENT, GW SMK SNPs#
######################################################

#Data from Wotton et al.

##################
#1. Load the data#
##################

SMK <- read_xlsx("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/RAW_Data/SmkInit_conditionally_independent_snps.xlsx")

#Let's get the data correctly:

colnames(SMK) <- SMK[1,]

SMK <- SMK[-1,]

View(SMK)

SMK$chr_pos_37 <- paste("chr", SMK$Chr, ":", SMK$Position, sep = "")
SMK$INFO <- SMK$`Effective N / max N`

summary(as.numeric(SMK$INFO)) #>0.30. 3 NAs, but still reported by them.

#########################################
#2. Get the genome-wide significant SNPs#
#########################################

summary(as.numeric(SMK$Pvalue))

SMK_gw <- SMK[which(as.numeric(SMK$Pvalue) <= 0.00000005),] #378: correct. The total number reported is this. The others are weird rows.

SMK_gw$`Alternate Allele Frequency` <- as.numeric(SMK_gw$`Alternate Allele Frequency`)
SMK_gw$Beta <- as.numeric(SMK_gw$Beta)
SMK_gw$SE <- as.numeric(SMK_gw$SE)
SMK_gw$Pvalue <- as.numeric(SMK_gw$Pvalue)
SMK_gw$N <- as.numeric(SMK_gw$N)

SMK_gw <- SMK_gw[-which(SMK_gw$Chr == 6 & SMK_gw$Position >= 26000000 & SMK_gw$Position <= 34000000),] #374

View(SMK_gw)

#################################################
#SECTION B: FINDING matching SNPs in the outcome#
#################################################

#First let's load the outcome:

whr <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHR/WHR_combined_Curated.txt")

#During the matching we are going to be EXTREMELY WARY.
#In other versions we curated the data beforehand because we didn't have enough RSIDs.
#Here we have all and we have other parameters that we need to be wary of.

##########################################################################################
#1. we are gonna try to match the data with RSID and with chr_pos to get as many variants#
##########################################################################################

whr_match_by_rsid <- whr[which(whr$MarkerName%in%SMK_gw$rsID),] #173/378

#Do we have a triallelic SNP that is chosen twice?

dupl <- whr_match_by_rsid$MarkerName[which(duplicated(whr_match_by_rsid$MarkerName) == TRUE)]

whr_match_by_rsid[which(whr_match_by_rsid$MarkerName == dupl),]

#Empty data.table (0 rows and 12 cols): MarkerName,Chr,Pos,Allele1,Allele2,FreqAllele1HapMapCEU...

summary(whr_match_by_rsid$FreqAllele1HapMapCEU) #all good. No need to do anything else

#################################
#Doing some checks with the data#
#################################

#1: do the chromosome and positions match? They should.

SMK_gw$chr_pos_37 <- paste("chr", SMK_gw$Chr, ":", SMK_gw$Position, sep = "")

length(which(whr_match_by_rsid$chr_pos_37%in%SMK_gw$chr_pos_37)) #173/378. Exactly the same as with RSID.

#We lose many SNPs that we need to recover.
#No problem. We will do so after some quick analysis.

#2. Check sample size, INFO and allele frequency

summary(whr_match_by_rsid$N) #no NAs, all good.
#summary(whr_match_by_rsid$INFO) #no info in GIANT data here.
summary(whr_match_by_rsid$FreqAllele1HapMapCEU) #perfect.

#We can merge the data without any worries.

######################################################
#SECTION D.2: obtainin proxies for these missing SNPs#
######################################################

#First we are gonna get the SNPs that are missing.

SMK_gw_missing <- SMK_gw[which(!(SMK_gw$chr_pos_37%in%whr_match_by_rsid$chr_pos_37)),]

#We are gonna save these SNPs for checking the rest of the traits.

fwrite(SMK_gw_missing, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Smoking_Initiation/SmkInit_WHR/proxies/Conditionally_independent_snps/conditionally_missing_snps.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

#And now we are going to get the proxies:

#####################################################################################
#THIS SECTION IS DEPRECATED SINCE WE CANNOT RECOVER THE SNPS WITH THEIR MAXIMUM LOCI#
#####################################################################################

#First we set the directory where we save them.

#setwd("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Smoking_Initiation/SmkInit_WHR/proxies/Conditionally_independent_snps/")

#getwd()

#Now we run the analysis.

#LDlinkR::LDproxy_batch(snp = SMK_gw_missing$chr_pos_37, pop = "EUR", r2d = "r2", append = TRUE, token = "04cad4ca4374")

#We retrieve them: I manually add a column for the rownames since they cannot be read properly with fread.
#And I set up the pathway.

#proxies_df <- read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Smoking_Initiation/SmkInit_WHR/proxies/combined_query_snp_list.txt", header = TRUE)

#######################################################
#The proxies dataframe needs to be cleaned in two ways#
#######################################################

#a) we need to remove the same snp that we are querying by taking all those SNPs with Distance = 0 & R2 == 1.

#weird_snps <-  proxies_df[which(proxies_df$query_snp == proxies_df$Coord & proxies_df$Distance == 0 & proxies_df$R2 == 1),]

#length(which(weird_snps$query_snp%in%SMK_gw_ind_missing$chr_pos_37)) #ALL IS GOOD. They are the query SNPs.

#Let's clean the data of them, then...

#proxies_df_clean <- proxies_df[which(!(proxies_df$query_snp == proxies_df$Coord & proxies_df$Distance == 0 & proxies_df$R2 == 1)),]

#head(proxies_df_clean)

#Now we are only interested in those with a r2 > 0.8

#proxies_df_clean_curated <- proxies_df_clean[which(proxies_df_clean$R2 > 0.8),]

#Just in case, let's check the MAF:

#summary(proxies_df_clean_curated$MAF) #perfect. No need to remove any.

#Cool, we are going to take the proxies with chr_pos, only.
#Just in case, since we have RSIDs that are a bit weird.

#proxies_chr_pos <- proxies_df_clean_curated$Coord

#Let's check in PS that they are correct.

#head(proxies_chr_pos)

#Mostly they are in buld 37, but we have many that are merged.
#We need to be wary. Actually some of them are not even correct.
#By a couple of positions, possible due to the changes...

#SMK_proxies <- SMK[which(SMK$chr_pos_37%in%proxies_chr_pos),] #1473/1781! Makes sense
#whr_proxies <- whr[which(whr$chr_pos_37%in%proxies_chr_pos),] #695/1781. Makes sense

#Now let's match the data:

#SMK_proxies_match <- SMK_proxies[which(SMK_proxies$chr_pos_37%in%whr_proxies$chr_pos_37),] #689/695. We miss some of them
#whr_proxies_match <- whr_proxies[which(whr_proxies$chr_pos_37%in%SMK_proxies$chr_pos_37),] #689/695.

#There are 5 that are missing in SMK.
#Let's check them just in case:

#whr_proxies[which(!(whr_proxies$chr_pos_37%in%SMK_proxies$chr_pos_37)),] #689/695.

#They are just missing. Nothing weird here.

#689 is the best that we can do!

#SMK_proxies_match <- SMK_proxies_match[order(SMK_proxies_match$RSID),]

#head(SMK_proxies_match) #ALL good we have no weird SNPs.

#Let's check for any MHC:

#SMK_proxies_match_chr6 <- SMK_proxies_match[which(SMK_proxies_match$CHROM == 6),]

#SMK_proxies_match_chr6 <- SMK_proxies_match_chr6[order(as.numeric(SMK_proxies_match_chr6$POS)),]

#View(SMK_proxies_match_chr6) #ALL GOOD.

###################################################################
#SECTION D: combined pre and post proxy data and get the lead SNPs# #DEPRECATED.
###################################################################

#SMK_proxies_match$rsid <- SMK_proxies_match$RSID
#SMK_proxies_match$pval <- as.numeric(SMK_proxies_match$PVALUE)

#Now, be wary..., we need to get those that are not missing to avoid removing the proxies we just got.

#SMK_gw_ind_not_missing <- SMK_gw_ind[which(SMK_gw_ind$chr_pos_37%in%whr_match_by_rsid$chr_pos_37),]

#SMK_gw_ind_plus_proxies <- rbind(SMK_gw_ind_not_missing, SMK_proxies_match) #we get 53 + 689 SNPs = 742. Makes total sense.

#When we run the ld_clump we will get only the lead SNPs in each loci!
#Which is exactly what we want. 

#SMK_gw_ind_post_proxies <- ieugwasr::ld_clump_local(SMK_gw_ind_plus_proxies, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#We end up with 86/93. This is actually really good.

##############################################
#SECTION E: Merging exposure and outcome data#
##############################################

#We need to format the data for the outcome.
#If we repeat the code we should erase chr_ column since I redid the code
#after changing all.
#However it took SO MUCH to get the matching that I won't repeat it today.

colnames(SMK_gw) <- c("locus.exposure", "chr.exposure", "pos.exposure", "SNP", "other_allele.exposure", "effect_allele.exposure",
                          "eaf.exposure", "STAT.exposure", "pval.exposure", "beta.exposure", "se.exposure",
                          "samplesize.exposure", "directions.exposure", "Cochran_Q.exposure", "Q_pval.exposure", "Q_df.exposure", "I2",  "effect_sample.exposure", "maf.exposure", "info.exposure",
                          "studies.exposure", "gene.exposure", "annotation.exposure", "genes_with_variants_r2_03", "annotation_full.exposure",
                          "Conditional_Variant_1",                                   
                          "Conditional_Variant_2",                                   
                          "Conditional_Variant_3",                                   
                          "Conditional_Variant_4",                                   
                          "Conditional_Variant_5",                                   
                          "Conditional_Variant_6",                                   
                          "Conditional_Variant_7",                                   
                          "Conditionally independent or sentinel (1=sentinel)",  
                          "chr_pos_37", "info.exposure_2")


whr_outcome <- whr 

colnames(whr_outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", 
                           "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                           "pval.outcome", "samplesize.outcome", "chr_pos_37", "chr_pos_18")

#Getting ids, so troublesome...

SMK_gw$id.exposure <- "smoking initiation"
whr_outcome$id.outcome <- "WHR"
SMK_gw$exposure <- "smoking initiation"
whr_outcome$outcome <- "WHR"

dat_1_pre <- harmonise_data(SMK_gw, whr_outcome, action = 3) #173: perfect.

########################################
#SECTION F: filtering harmonise results#
########################################

dat_1_filt <- dat_1_pre[which(dat_1_pre$remove == FALSE),] #173/173
dat_1_filt <- dat_1_filt[which(dat_1_filt$palindromic == FALSE & dat_1_filt$ambiguous == FALSE | dat_1_filt$palindromic == TRUE & dat_1_filt$ambiguous == FALSE),] #173/173

#So we are gonna use the steiger filtering like that:

dat_1_filt$units.exposure <- "log odds"
dat_1_filt$units.outcome <- "SD"
dat_1_filt$ncase.exposure <- 557337
dat_1_filt$ncontrol.exposure <- 674753

#We are not gonna set a prevalence because I have no clue what it is. 
#They say steiger is not too sensitive about it, so let's keep it this way.

dat_1_filt <- steiger_filtering(dat_1_filt) 
dat_1_filt_1 <- dat_1_filt[which(dat_1_filt$steiger_dir == TRUE),] #122/175. Not bad!!

dat_1_filt <- dat_1_filt_1

dat_1_filt$mr_keep <- TRUE #this allows us to keep the unambiguous.

#We have 121/175

fwrite(dat_1_filt, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/output/1_curated_data/smoking_initiation_whr_curated.txt")

##########################
#Final checks on SMK data#
##########################

summary(dat_1_filt$eaf.exposure) #perfect
summary(as.numeric(dat_1_filt$info.exposure_2)) #perfect.

##############################################
#SECTION G: Checking validity of the variants#
##############################################

#Are the instruments valid??
#Let's calculate mF and Isq to check it.

F_ = ((dat_1_filt$beta.exposure)^2)/((dat_1_filt$se.exposure)^2)
mF  = mean(F_)

print(mF)
#44.50083

Isq(dat_1_filt$beta.exposure, dat_1_filt$se.exposure)
#0.9770934

#########################################
#SECTION H: First round of 2SMR analysis#
#########################################

mr(dat_1_filt)

#        id.exposure id.outcome outcome           exposure                    method nsnp          b         se         pval
#1 smoking initiation        WHR     WHR smoking initiation                  MR Egger  121 0.18985367 0.09956550 0.0589546966
#2 smoking initiation        WHR     WHR smoking initiation           Weighted median  121 0.11621749 0.03014001 0.0001152946
#3 smoking initiation        WHR     WHR smoking initiation Inverse variance weighted  121 0.07307571 0.02200239 0.0008961109
#4 smoking initiation        WHR     WHR smoking initiation               Simple mode  121 0.20126647 0.08471833 0.0190986355
#5 smoking initiation        WHR     WHR smoking initiation             Weighted mode  121 0.19161415 0.08649330 0.0286221445

#####################################################
#SECTION I: ANALYSIS ON PLEIOTROPY AND HETEROGENEITY#
#####################################################

library(meta)

dat_1_filt$mr <- dat_1_filt$beta.outcome/dat_1_filt$beta.exposure
dat_1_filt$mr_se <- ((dat_1_filt$mr*((dat_1_filt$se.exposure/dat_1_filt$beta.exposure)^2+(dat_1_filt$se.outcome/dat_1_filt$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_filt$mr, dat_1_filt$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0 [<0.0000; <0.0000]; tau = 0 [<0.0000; <0.0000];
#I^2 = 0.0% [0.0%; 0.0%]; H = 1.00 [1.00; 1.00]

#Test of heterogeneity:
#  Q d.f. p-value
#70.77  119  0.9999

#Lots of heterogeneity.

mr_rucker(dat_1_filt)

#[[1]]$intercept
#Method     Estimate          SE       CI_low       CI_upp         P
#1  Egger fixed effects -0.002222733 0.001432823 -0.005031015 0.0005855484 0.1208307
#2 Egger random effects -0.002222733 0.001432823 -0.005031015 0.0005855484 0.9395847

#[[1]]$Q
#Method         Q  df         P
#1   Q_ivw 72.962581 120 0.9997797
#2 Q_egger 71.516316 119 0.9998256
#3  Q_diff  1.446264   1 0.2291283

#[[1]]$res
#[1] "A"

#[[1]]$selected
#Method nsnp   Estimate         SE     CI_low    CI_upp            P
#1 Rucker  121 0.07307571 0.01715652 0.03944954 0.1067019 2.050181e-05

#We need to remove variants with Egger Radial.

#######################################
#SECTION J: GENERATE SENSITIVITY PLOTS#
#######################################

dat_1_filt$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Smoking_Initiation/SmkInit_WHR/SMK_conditionally_independent_WHR_GIANT_all.tiff ", units="in", width=10, height=10, res=300)
mr_plots(dat_1_filt)
dev.off()

############################
#NO NEED TO REMOVE ANYTHING#
############################

##########################################
#SECTION K: CHECKING BY REMOVING OUTLIERS#
##########################################

library(RadialMR)

radial_input <- format_radial(dat_1_filt$beta.exposure, dat_1_filt$beta.outcome, dat_1_filt$se.exposure, dat_1_filt$se.outcome, dat_1_filt$SNP)

radial_output <- ivw_radial(radial_input, alpha = 0.05, weights = 3) 

outliers <- radial_output$outliers$SNP #NO OUTLIERS DETECTED.
