##############
#INTRODUCTION#
##############

#This is the 2SMR code to run Pack Years > whradjbmi.

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

############
#LD PROXIES#
############

#We need to make a function that gets the info I need for LD proxies.

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

######################################################"
#SECTION A: FINDING INDEPENDENT, GW packyears_smk SNPs#
######################################################"

#Data from GWAS MRCIEU

##################
#1. Load the data#
##################

packyears_smk <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/PacksPerYear/PackPerYear_Curated_FULL.txt")

#########################################
#2. Get the genome-wide significant SNPs#
#########################################

packyears_smk_gw <- packyears_smk[which(packyears_smk$pval.exposure <= 0.00000005),] #1068.

#Let's check that the filtering is done correctly...

summary(packyears_smk_gw$pval.exposure) #perfect.

#And that all the SNPs present RSIDs that we can use:

packyears_smk_gw <- packyears_smk_gw[order(packyears_smk_gw$SNP),]

head(packyears_smk_gw) #all of them with RSIDs.
tail(packyears_smk_gw) #all of them with RSIDs.

#Let's check the eaf for a sec, just in case:

summary(packyears_smk_gw$eaf.exposure) #perfect.

#And we do a final check to avoid SNPs in the MHC region:

which(packyears_smk_gw$chr.exposure == 6 & packyears_smk_gw$pos.exposure >= 26000000 & packyears_smk_gw$pos.exposure <= 34000000) #as expected.

#We removed them properly in the curation stage.
#Perfect!!

#########################################################
#3. calculate which SNPs are independent from each other#
#########################################################

#we get the data for the pruning

packyears_smk_gw$rsid <- packyears_smk_gw$SNP
packyears_smk_gw$pval <- packyears_smk_gw$pval.exposure

#And we perform the pruning:

packyears_smk_gw_ind <- ieugwasr::ld_clump_local(packyears_smk_gw, bfile = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/2SMR/1000G/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/R-4.1.2/library/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) 

#################################################
#SECTION B: FINDING matching SNPs in the outcome#
#################################################

#First let's load the outcome:

whradjbmi <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMI_comb/WHRAdjBMI_Combined_Curated.txt")

##########################################################################################
#1. we are gonna try to match the data with RSID and with chr_pos to get as many variants#
##########################################################################################

whradjbmi_match_by_rsid <- whradjbmi[which(whradjbmi$MarkerName%in%packyears_smk_gw_ind$SNP),] #7/13

############################################################################
#Small code to get the chr_pos for packyears_smk (full) and the gw_ind SNPs#
############################################################################

packyears_smk$chr_ <- paste("chr", packyears_smk$chr.exposure, sep = "")
packyears_smk$chr_pos_37 <- paste(packyears_smk$chr_, packyears_smk$pos.exposure, sep = ":")

head(packyears_smk)

packyears_smk_gw_ind$chr_ <- paste("chr", packyears_smk_gw_ind$chr.exposure, sep = "")
packyears_smk_gw_ind$chr_pos_37 <- paste(packyears_smk_gw_ind$chr_, packyears_smk_gw_ind$pos.exposure, sep = ":")

packyears_smk_gw_ind

##############################
#Now let's match with chr_pos#
##############################

whradjbmi_match_by_chr_pos <- whradjbmi[which(whradjbmi$chr_pos_37%in%packyears_smk_gw_ind$chr_pos_37),] #7/13

######################################################
#SECTION C: get proxies for the SNPs that are missing#
######################################################

#First we are gonna get the SNPs that are missing.

packyears_smk_gw_ind_missing <- packyears_smk_gw_ind[which(!(packyears_smk_gw_ind$chr_pos_37%in%whradjbmi_match_by_chr_pos$chr_pos_37)),]

setwd("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/PackYears/PackYears_WHRAdjBMI/proxies/")

getwd()

LDlinkR::LDproxy_batch(snp = packyears_smk_gw_ind_missing$chr_pos_37, pop = "EUR", r2d = "r2", append = TRUE, token = "04cad4ca4374")

#And now we are going to get the proxies from WHR since they are bound to be the same as WHRAjdBMI

proxies_df <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/PackYears/PackYears_WHRAdjBMI/proxies/combined_query_snp_list.txt")

#Let's check all the missing SNPs are in the query_snps:

query_snps <- proxies_df[which(proxies_df$query_snp == proxies_df$Coord & proxies_df$Distance == 0 & proxies_df$R2 == 1),]

length(which(query_snps$query_snp%in%packyears_smk_gw_ind_missing$chr_pos_37)) #73/73 we can proceed.

#LET'S CLEAN THE DATA THEN:

#a) we need to remove the same snp that we are querying by taking all those SNPs with Distance = 0 & R2 == 1.

proxies_df_clean <- proxies_df[which(!(proxies_df$query_snp == proxies_df$Coord & proxies_df$Distance == 0 & proxies_df$R2 == 1)),]

head(proxies_df_clean)

#Now we are only interested in those with a r2 > 0.8

proxies_df_clean_curated <- proxies_df_clean[which(proxies_df_clean$R2 > 0.8),]

#Just in case, let's check the MAF:

summary(proxies_df_clean_curated$MAF) #perfect. No need to remove any.

#Mostly they are in buld 37, but we have many that are merged.
#We need to be wary. Actually some of them are not even correct.
#By a couple of positions, possible due to the changes...

#Let's find the proxies in both df:

packyears_smk_proxies <- packyears_smk[which(packyears_smk$chr_pos_37%in%proxies_df_clean_curated$Coord),] #158
whradjbmi_proxies <- whradjbmi[which(whradjbmi$chr_pos_37%in%proxies_df_clean_curated$Coord),] #38

#Now let's match the data:

packyears_smk_proxies_match <- packyears_smk_proxies[which(packyears_smk_proxies$chr_pos_37%in%whradjbmi_proxies$chr_pos_37),] #38/38
whradjbmi_proxies_match <- whradjbmi_proxies[which(whradjbmi_proxies$chr_pos_37%in%packyears_smk_proxies$chr_pos_37),] #38/38

#whr is a bottleneck and this is the best we can do right now!! However, there are 4 SNPs that are not found in Smk.

packyears_smk_proxies_match_2 <- packyears_smk_proxies[which(packyears_smk_proxies$chr_pos_37%in%whradjbmi_proxies$chr_pos_37),] #691/691

#It is the same! The mismatch is NOT due to RSIDs, but just because the variant is no there.

##################
#CHECKING ALLELES#
##################

proxies_check <- proxies_df_clean_curated[which(proxies_df_clean_curated$Coord%in%packyears_smk_proxies_match_2$chr_pos_37),]

proxies_check <- proxies_check[order(proxies_check$Alleles),]

head(proxies_check)
tail(proxies_check)

###################################################################
#SECTION D: combined pre and post proxy data and get the lead SNPs#
###################################################################

#We need to add a column to merge the 7 not missing variants with the proxies and perform the clumping.

packyears_smk_proxies_match$pval <- packyears_smk_proxies_match$pval.exposure

#Now, be wary..., we need to get those that are not missing to avoid removing the proxies we just got.

packyears_smk_gw_ind_not_missing <- packyears_smk_gw_ind[which(packyears_smk_gw_ind$chr_pos_37%in%whradjbmi_match_by_chr_pos$chr_pos_37),]

packyears_smk_gw_ind_plus_proxies <- rbind(packyears_smk_gw_ind_not_missing, packyears_smk_proxies_match) #we get 7 + 38 = 45 SNPs. Makes total sense.

#When we run the ld_clump we will get only the lead SNPs in each loci!
#Which is exactly what we want. 

packyears_smk_gw_ind_plus_proxies$rsid <- packyears_smk_gw_ind_plus_proxies$SNP
packyears_smk_gw_ind_plus_proxies$pval <- packyears_smk_gw_ind_plus_proxies$pval.exposure

packyears_smk_gw_ind_post_proxies <- ieugwasr::ld_clump_local(packyears_smk_gw_ind_plus_proxies, bfile = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/2SMR/1000G/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/R-4.1.2/library/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) 

#We do not recover all, but we obtain 12/13. That is super!

##############################################
#SECTION E: Merging exposure and outcome data#
##############################################

#We need to format the data for the outcome:

whradjbmi_outcome <- whradjbmi

colnames(whradjbmi_outcome) <- c("SNP", "effect_allele.outcome", 
                           "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                           "pval.outcome", "samplesize.outcome", "chr_pos_37")


#Getting ids, so troublesome...

packyears_smk_gw_ind_post_proxies$id.exposure <- "pack-years"
whradjbmi_outcome$id.outcome <- "WHRadjBMI"
packyears_smk_gw_ind_post_proxies$exposure <- "pack-years"
whradjbmi_outcome$outcome <- "WHRadjBMI"

dat_1_pre <- harmonise_data(packyears_smk_gw_ind_post_proxies, whradjbmi_outcome, action = 3) #12/12: perfect.

########################################
#SECTION F: filtering harmonise results#
########################################

dat_1_filt <- dat_1_pre[which(dat_1_pre$remove == FALSE),] #12/12
dat_1_filt <- dat_1_filt[which(dat_1_filt$palindromic == FALSE & dat_1_filt$ambiguous == FALSE | dat_1_filt$palindromic == TRUE & dat_1_filt$ambiguous == FALSE),] #10/12

#So we are gonna use the steiger filtering like that:

dat_1_filt$samplesize.exposure <- 	142387
dat_1_filt$units.exposure <- "SD"
dat_1_filt$units.outcome <- "SD"

dat_1_filt <- steiger_filtering(dat_1_filt) 
dat_1_filt_1 <- dat_1_filt[which(dat_1_filt$steiger_dir == TRUE),] #10/10

dat_1_filt <- dat_1_filt_1

dat_1_filt$mr_keep = TRUE #this is to allow the palindromic unambiguous to be kept.

##############################################
#SECTION G: Checking validity of the variants#
##############################################

#Are the instruments valid??
#Let's calculate mF and Isq to check it.

F_ = ((dat_1_filt$beta.exposure)^2)/((dat_1_filt$se.exposure)^2)
mF  = mean(F_)

print(mF)
#77.73074

Isq(dat_1_filt$beta.exposure, dat_1_filt$se.exposure)
#0.9883244

#########################################
#SECTION H: First round of 2SMR analysis#
#########################################

mr(dat_1_filt)

#id.exposure id.outcome   outcome                 exposure                    method nsnp           b         se       pval
#1 Pack years adult smoking  WHRadjBMI WHRadjBMI Pack years adult smoking                  MR Egger   10 0.007076867 0.07579261 0.92790453
#2 Pack years adult smoking  WHRadjBMI WHRadjBMI Pack years adult smoking           Weighted median   10 0.054698724 0.04306501 0.20403371
#3 Pack years adult smoking  WHRadjBMI WHRadjBMI Pack years adult smoking Inverse variance weighted   10 0.081513091 0.03879782 0.03564352
#4 Pack years adult smoking  WHRadjBMI WHRadjBMI Pack years adult smoking               Simple mode   10 0.028491411 0.10127857 0.78482784
#5 Pack years adult smoking  WHRadjBMI WHRadjBMI Pack years adult smoking             Weighted mode   10 0.047185285 0.04721821 0.34375603

#####################################################
#SECTION I: ANALYSIS ON PLEIOTROPY AND HETEROGENEITY#
#####################################################

library(meta)

dat_1_filt$mr <- dat_1_filt$beta.outcome/dat_1_filt$beta.exposure
dat_1_filt$mr_se <- ((dat_1_filt$mr*((dat_1_filt$se.exposure/dat_1_filt$beta.exposure)^2+(dat_1_filt$se.outcome/dat_1_filt$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_filt$mr, dat_1_filt$mr_se)

#Quantifying heterogeneity:
#  tau^2 < 0.0001 [0.0000; 0.0835]; tau = 0.0017 [0.0000; 0.2890]
#I^2 = 0.0% [0.0%; 62.4%]; H = 1.00 [1.00; 1.63]

#Test of heterogeneity:
#  Q d.f. p-value
#8.82    9  0.4540

#No heterogeneity. 
#Let's check the rucker framework:

mr_rucker(dat_1_filt)

#[[1]]$intercept
#Method    Estimate          SE       CI_low      CI_upp         P
#1  Egger fixed effects 0.003370218 0.002925263 -0.002363191 0.009103627 0.2492768
#2 Egger random effects 0.003370218 0.002964368 -0.002439837 0.009180273 0.1277881

#[[1]]$Q
#Method        Q df         P
#1   Q_ivw 9.542675  9 0.3887591
#2 Q_egger 8.215323  8 0.4127246
#3  Q_diff 1.327352  1 0.2492768

#[[1]]$res
#[1] "A"

#[[1]]$selected
#Method nsnp   Estimate         SE      CI_low    CI_upp          P
#1 Rucker   10 0.08151309 0.03767849 0.007664605 0.1553616 0.03051156

#######################################
#SECTION J: GENERATE SENSITIVITY PLOTS#
#######################################

dat_1_filt$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/PackYears/PackYears_WHRAdjBMI/packyears_smk_whradjbmi_original_20082021.tiff", units="in", width=10, height=10, res=300)
mr_plots(dat_1_filt)
dev.off()

#There seems to be a SNP that is too powerful and might be inducing pleiotropy.
#Interesting. Here it seems the outliers seem to be coming from the positive side.

##############################
#SECTION K: Removing outliers#
##############################

library(RadialMR)

radial_input <- format_radial(dat_1_filt$beta.exposure, dat_1_filt$beta.outcome, dat_1_filt$se.exposure, dat_1_filt$se.outcome, dat_1_filt$SNP)

radial_output <- ivw_radial(radial_input, alpha = 0.05, weights = 3, tol = 0.0001) 

#NO OUTLIERS.