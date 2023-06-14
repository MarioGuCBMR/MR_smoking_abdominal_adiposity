##############
#INTRODUCTION#
##############

#This is the 2SMR code to run Lifetime smoking with WHR (Pulit)

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

#############################################
#SECTION A: FINDING INDEPENDENT, GW CIG SNPs#
#############################################

#Data from Wotton et al.

##################
#1. Load the data#
##################

CIG <- read_xlsx("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/RAW_Data/CigDay_conditionally_independent_snps.xlsx")

#Let's get the data correctly:

colnames(CIG) <- CIG[1,]

CIG <- CIG[-1,]

View(CIG)

CIG$chr_pos_37 <- paste("chr", CIG$Chr, ":", CIG$Position, sep = "")
CIG$INFO <- CIG$`Effective N / max N`

summary(as.numeric(CIG$INFO)) #>0.30 as expected.

#########################################
#2. Get the genome-wide significant SNPs#
#########################################

summary(as.numeric(CIG$Pvalue))

CIG_gw <- CIG[which(as.numeric(CIG$Pvalue) <= 0.00000005),] #55: correct. The total number reported is this. The others are weird rows.

CIG_gw$`Alternate Allele Frequency` <- as.numeric(CIG_gw$`Alternate Allele Frequency`)
CIG_gw$Beta <- as.numeric(CIG_gw$Beta)
CIG_gw$SE <- as.numeric(CIG_gw$SE)
CIG_gw$Pvalue <- as.numeric(CIG_gw$Pvalue)
CIG_gw$N <- as.numeric(CIG_gw$N)

View(CIG_gw)

#################################################
#SECTION B: FINDING matching SNPs in the outcome#
#################################################

#First let's load the outcome:

wcadjbmi <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WCAdjBMISmk/WCAdjBMI_Smk_Strat_Curated.txt")

#During the matching we are going to be EXTREMELY WARY.
#In other versions we curated the data beforehand because we didn't have enough RSIDs.
#Here we have all and we have other parameters that we need to be wary of.

##########################################################################################
#1. we are gonna try to match the data with RSID and with chr_pos to get as many variants#
##########################################################################################

wcadjbmi_match_by_rsid <- wcadjbmi[which(wcadjbmi$rs_id%in%CIG_gw$rsID),] #15/55

#Do we have a triallelic SNP that is chosen twice?

dupl <- wcadjbmi_match_by_rsid$rs_id[which(duplicated(wcadjbmi_match_by_rsid$rs_id) == TRUE)]

wcadjbmi_match_by_rsid[which(wcadjbmi_match_by_rsid$rs_id == dupl),]

#Empty data.table (0 rows and 12 cols): MarkerName,Chr,Pos,Allele1,Allele2,FreqAllele1HapMapCEU...

summary(wcadjbmi_match_by_rsid$freq_allele1_hapmapceu) #all good. No need to do anything else

#################################
#Doing some checks with the data#
#################################

#1: do the chromosome and positions match? They should.

CIG_gw$chr_pos_37 <- paste("chr", CIG_gw$Chr, ":", CIG_gw$Position, sep = "")

length(which(wcadjbmi_match_by_rsid$chr_pos_37%in%CIG_gw$chr_pos_37)) #15/55. Exactly the same as with RSID.

#We lose many SNPs that we need to recover.
#No problem. We will do so after some quick analysis.

#2. Check sample size, INFO and allele frequency

summary(wcadjbmi_match_by_rsid$N_SMK) #no NAs, all good.
#summary(whr_match_by_rsid$INFO) #no info in GIANT data here.
summary(wcadjbmi_match_by_rsid$freq_allele1_hapmapceu) #perfect.

#We can merge the data without any worries.

######################################################
#SECTION D.2: obtainin proxies for these missing SNPs#
######################################################

#First we are gonna get the SNPs that are missing.

#CIG_gw_missing <- CIG_gw[which(!(CIG_gw$chr_pos_37%in%whr_match_by_rsid$chr_pos_37)),]

#We are gonna save these SNPs for checking the rest of the traits.

#fwrite(CIG_gw_missing, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Smoking_Initiation/CIGInit_WHR/proxies/Conditionally_independent_snps/conditionally_missing_snps.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

#And now we are going to get the proxies:

#####################################################################################
#THIS SECTION IS DEPRECATED SINCE WE CANNOT RECOVER THE SNPS WITH THEIR MAXIMUM LOCI#
#####################################################################################

#First we set the directory where we save them.

#setwd("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Smoking_Initiation/CIGInit_WHR/proxies/Conditionally_independent_snps/")

#getwd()

#Now we run the analysis.

#LDlinkR::LDproxy_batch(snp = CIG_gw_missing$chr_pos_37, pop = "EUR", r2d = "r2", append = TRUE, token = "04cad4ca4374")

#We retrieve them: I manually add a column for the rownames since they cannot be read properly with fread.
#And I set up the pathway.

#proxies_df <- read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Smoking_Initiation/CIGInit_WHR/proxies/combined_query_snp_list.txt", header = TRUE)

#######################################################
#The proxies dataframe needs to be cleaned in two ways#
#######################################################

#a) we need to remove the same snp that we are querying by taking all those SNPs with Distance = 0 & R2 == 1.

#weird_snps <-  proxies_df[which(proxies_df$query_snp == proxies_df$Coord & proxies_df$Distance == 0 & proxies_df$R2 == 1),]

#length(which(weird_snps$query_snp%in%CIG_gw_ind_missing$chr_pos_37)) #ALL IS GOOD. They are the query SNPs.

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

#CIG_proxies <- CIG[which(CIG$chr_pos_37%in%proxies_chr_pos),] #1473/1781! Makes sense
#whr_proxies <- whr[which(whr$chr_pos_37%in%proxies_chr_pos),] #695/1781. Makes sense

#Now let's match the data:

#CIG_proxies_match <- CIG_proxies[which(CIG_proxies$chr_pos_37%in%whr_proxies$chr_pos_37),] #689/695. We miss some of them
#whr_proxies_match <- whr_proxies[which(whr_proxies$chr_pos_37%in%CIG_proxies$chr_pos_37),] #689/695.

#There are 5 that are missing in CIG.
#Let's check them just in case:

#whr_proxies[which(!(whr_proxies$chr_pos_37%in%CIG_proxies$chr_pos_37)),] #689/695.

#They are just missing. Nothing weird here.

#689 is the best that we can do!

#CIG_proxies_match <- CIG_proxies_match[order(CIG_proxies_match$RSID),]

#head(CIG_proxies_match) #ALL good we have no weird SNPs.

#Let's check for any MHC:

#CIG_proxies_match_chr6 <- CIG_proxies_match[which(CIG_proxies_match$CHROM == 6),]

#CIG_proxies_match_chr6 <- CIG_proxies_match_chr6[order(as.numeric(CIG_proxies_match_chr6$POS)),]

#View(CIG_proxies_match_chr6) #ALL GOOD.

###################################################################
#SECTION D: combined pre and post proxy data and get the lead SNPs# #DEPRECATED.
###################################################################

#CIG_proxies_match$rsid <- CIG_proxies_match$RSID
#CIG_proxies_match$pval <- as.numeric(CIG_proxies_match$PVALUE)

#Now, be wary..., we need to get those that are not missing to avoid removing the proxies we just got.

#CIG_gw_ind_not_missing <- CIG_gw_ind[which(CIG_gw_ind$chr_pos_37%in%whr_match_by_rsid$chr_pos_37),]

#CIG_gw_ind_plus_proxies <- rbind(CIG_gw_ind_not_missing, CIG_proxies_match) #we get 53 + 689 SNPs = 742. Makes total sense.

#When we run the ld_clump we will get only the lead SNPs in each loci!
#Which is exactly what we want. 

#CIG_gw_ind_post_proxies <- ieugwasr::ld_clump_local(CIG_gw_ind_plus_proxies, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#We end up with 86/93. This is actually really good.

##############################################
#SECTION E: Merging exposure and outcome data#
##############################################

#We need to format the data for the outcome.
#If we repeat the code we should erase chr_ column since I redid the code
#after changing all.
#However it took SO MUCH to get the matching that I won't repeat it today.

colnames(CIG_gw) <- c("locus.exposure", "chr.exposure", "pos.exposure", "SNP", "other_allele.exposure", "effect_allele.exposure",
                          "eaf.exposure", "STAT.exposure", "pval.exposure", "beta.exposure", "se.exposure",
                          "samplesize.exposure", "directions.exposure", "Cochran_Q.exposure", "Q_pval.exposure", "Cochran_Q_df", "I2",  "effect_sample.exposure", "maf.exposure", "info.exposure",
                          "studies.exposure", "gene.exposure", "annotation.exposure", "genes_with_variants_r2_03", "sign_pascal",
                          "Conditional_Variant_1",                                   
                          "Conditional_Variant_2",                                   
                          "Conditional_Variant_3",                                   
                          "Conditional_Variant_4",                                   
                          "Conditional_Variant_5",                                   
                          "Conditional_Variant_6",                                   
                          "Conditional_Variant_7",    
                          "Conditional_Variant_8",                                   
                          "Conditionally independent or sentinel (1=sentinel)",  
                          "chr_pos_37", "info.exposure_2")


wcadjbmi_outcome <- wcadjbmi 

colnames(wcadjbmi_outcome) <- c("chr.outcome", "SNP", "markername", "position_hg18",
                                 "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "samplesize.outcome",
                                 "beta.outcome", "se.outcome", "pval.outcome", "N_NonSmk", 
                                 "Effect_NonSMK", "StdErr_NonSMK", "P_value_NonSMK", "chr_pos_37",
                                 "chr_pos_18")


#Getting ids, so troublesome...

CIG_gw$id.exposure <- "cigarettes per day (current/past)"
wcadjbmi_outcome$id.outcome <- "WCadjBMI"
CIG_gw$exposure <- "cigarettes per day (current/past)"
wcadjbmi_outcome$outcome <- "WCadjBMI"

dat_1_pre <- harmonise_data(CIG_gw, wcadjbmi_outcome, action = 3) #15: perfect.

########################################
#SECTION F: filtering harmonise results#
########################################

dat_1_filt <- dat_1_pre[which(dat_1_pre$remove == FALSE),] #15/15
dat_1_filt <- dat_1_filt[which(dat_1_filt$palindromic == FALSE & dat_1_filt$ambiguous == FALSE | dat_1_filt$palindromic == TRUE & dat_1_filt$ambiguous == FALSE),] #15/15

#So we are gonna use the steiger filtering like that:

dat_1_filt$units.exposure <- "SD"
dat_1_filt$units.outcome <- "SD"

#We are not gonna set a prevalence because I have no clue what it is. 
#They say steiger is not too sensitive about it, so let's keep it this way.

dat_1_filt <- steiger_filtering(dat_1_filt) 
dat_1_filt_1 <- dat_1_filt[which(dat_1_filt$steiger_dir == TRUE),] #13/15. Not bad!!

dat_1_filt <- dat_1_filt_1

dat_1_filt$mr_keep <- TRUE #this allows us to keep the unambiguous.

#We have 13/15

######################
#Let's check the SNPs#
######################

which(dat_1_filt$chr.exposure == 6) #perfect 
summary(as.numeric(dat_1_filt$info.exposure_2)) #perfect.
summary(as.numeric(dat_1_filt$eaf.exposure)) #perfect.

##############################################
#SECTION G: Checking validity of the variants#
##############################################

#Are the instruments valid??
#Let's calculate mF and Isq to check it.

F_ = ((dat_1_filt$beta.exposure)^2)/((dat_1_filt$se.exposure)^2)
mF  = mean(F_)

print(mF)
#149.4615

Isq(dat_1_filt$beta.exposure, dat_1_filt$se.exposure)
#0.9934471

#########################################
#SECTION H: First round of 2SMR analysis#
#########################################

mr(dat_1_filt)

#id.exposure id.outcome  outcome           exposure                    method nsnp          b         se      pval
#1 Cigarettes per day   WCadjBMI WCadjBMI Smoking initiation                  MR Egger   13 0.12297765 0.10422272 0.2629085
#2 Cigarettes per day   WCadjBMI WCadjBMI Smoking initiation           Weighted median   13 0.11472025 0.08029220 0.1530662
#3 Cigarettes per day   WCadjBMI WCadjBMI Smoking initiation Inverse variance weighted   13 0.10438079 0.07031942 0.1377080
#4 Cigarettes per day   WCadjBMI WCadjBMI Smoking initiation               Simple mode   13 0.07352216 0.22788375 0.7525273
#5 Cigarettes per day   WCadjBMI WCadjBMI Smoking initiation             Weighted mode   13 0.11292625 0.07734675 0.1699701

#####################################################
#SECTION I: ANALYSIS ON PLEIOTROPY AND HETEROGENEITY#
#####################################################

library(meta)

dat_1_filt$mr <- dat_1_filt$beta.outcome/dat_1_filt$beta.exposure
dat_1_filt$mr_se <- ((dat_1_filt$mr*((dat_1_filt$se.exposure/dat_1_filt$beta.exposure)^2+(dat_1_filt$se.outcome/dat_1_filt$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_filt$mr, dat_1_filt$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0 [0.0000; 0.2932]; tau = 0 [0.0000; 0.5414];
#I^2 = 0.0% [0.0%; 44.3%]; H = 1.00 [1.00; 1.34]

#Test of heterogeneity:
#  Q d.f. p-value
#9.35   12  0.6729

#No of heterogeneity.

mr_rucker(dat_1_filt)

#[[1]]$intercept
#Method      Estimate          SE       CI_low      CI_upp         P
#1  Egger fixed effects -0.0009323811 0.003622687 -0.008032718 0.006167956 0.7968910
#2 Egger random effects -0.0009323811 0.003622687 -0.008032718 0.006167956 0.6015545

#[[1]]$Q
#Method          Q df         P
#1   Q_ivw 9.76363550 12 0.6366880
#2 Q_egger 9.70519194 11 0.5570841
#3  Q_diff 0.05844356  1 0.8089731

#[[1]]$res
#[1] "A"

#[[1]]$selected
#Method nsnp  Estimate         SE     CI_low    CI_upp          P
#1 Rucker   13 0.1043808 0.06342938 -0.0199385 0.2287001 0.09984156

#We need to remove variants with Egger Radial.

#######################################
#SECTION J: GENERATE SENSITIVITY PLOTS#
#######################################

dat_1_filt$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/CigDay_comb/CigDay_Comb_WCadjBMIstrat/CigDay_comb_conditionally_independent_WCAdjBMIstrat_all.tiff ", units="in", width=10, height=10, res=300)
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
