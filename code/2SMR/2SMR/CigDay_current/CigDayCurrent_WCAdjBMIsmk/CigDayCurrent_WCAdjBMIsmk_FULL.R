##############
#INTRODUCTION#
##############

#This is the 2SMR code to run cigday (current) > WCAdjBMI.

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

#################################################
#SECTION A: FINDING INDEPENDENT, GW cigday SNPs#
#################################################

#Note the data is Packs per Year (previously + currently), we are using it as a proxy.

##################
#1. Load the data#
##################

cigday <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/CigDayCurrent/Cig_Day_Current_Curated_FULL.txt")

#########################################
#2. Get the genome-wide significant SNPs#
#########################################

#Only the genome-wide gets me 1 SNP... consid

cigday_gw <- extract_instruments(outcomes = "ukb-b-469", p1 = 0.0000005, clump = FALSE) #657

#Only those in the huge dataset we curated are gonna go through the pipeline:

cigday_gw <- cigday_gw[which(cigday_gw$SNP%in%cigday$SNP),] #657 all of them.

#Checking for duplicates:

cigday_gw[which(duplicated(cigday_gw$SNP) == TRUE),] #NONE.

#Let's check that the filtering is done correctly...

summary(cigday_gw$pval.exposure) #perfect.

#And that all the SNPs present RSIDs that we can use:

cigday_gw <- cigday_gw[order(cigday_gw$SNP),]

head(cigday_gw) #all of them with RSIDs.
tail(cigday_gw) #all of them with RSIDs.

#Let's check the eaf for a sec, just in case:

summary(cigday_gw$eaf.exposure) #perfect.

#And we do a final check to avoid SNPs in the MHC region:

#In this case we do have them, so we need to remomve them:

which(cigday_gw$chr.exposure == 6 & cigday_gw$pos.exposure >= 26000000 & cigday_gw$pos.exposure <= 34000000)

#########################################################
#3. calculate which SNPs are independent from each other#
#########################################################

#we get the data for the pruning

cigday_gw$rsid <- cigday_gw$SNP
cigday_gw$pval <- cigday_gw$pval.exposure

#And we perform the pruning:

cigday_gw_ind <- ieugwasr::ld_clump_local(cigday_gw, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#We end up with 4! It is true that many of the variants are removed from the panel.
#But hey, we are using the same panel for everything, so nothing we can do here, pal.

View(cigday_gw_ind) #none of those SNPs are present...

#################################################
#SECTION B: FINDING matching SNPs in the outcome#
#################################################

#First let's load the outcome:

WCAdjBMI <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WCAdjBMIsmk/WCAdjBMI_Smk_Strat_Curated.txt")

##########################################################################################
#1. we are gonna try to match the data with RSID and with chr_pos to get as many variants#
##########################################################################################

WCAdjBMI_match_by_rsid <- WCAdjBMI[which(WCAdjBMI$rs_id%in%cigday_gw_ind$SNP),] #1/4

########################################################################
#Small code to get the chr_pos for cigday (current) and the gw_ind SNPs#
########################################################################

cigday$chr_ <- paste("chr", cigday$chr.exposure, sep = "")
cigday$chr_pos_37 <- paste(cigday$chr_, cigday$pos.exposure, sep = ":")

head(cigday)

cigday_gw_ind$chr_ <- paste("chr", cigday_gw_ind$chr.exposure, sep = "")
cigday_gw_ind$chr_pos_37 <- paste(cigday_gw_ind$chr_, cigday_gw_ind$pos.exposure, sep = ":")

cigday_gw_ind

##############################
#Now let's match with chr_pos#
##############################

WCAdjBMI_match_by_chr_pos <- WCAdjBMI[which(WCAdjBMI$chr_pos_37%in%cigday_gw_ind$chr_pos_3 ),] #1/4

######################################################
#SECTION C: get proxies for the SNPs that are missing#
######################################################

#First we are gonna get the SNPs that are missing.

cigday_gw_ind_missing <- cigday_gw_ind[which(!(cigday_gw_ind$chr_pos_37%in%WCAdjBMI_match_by_chr_pos$chr_pos_37)),]

getwd()

#And now we are going to get the proxies:

LDlinkR::LDproxy_batch(snp = cigday_gw_ind_missing$chr_pos_37, pop = "EUR", r2d = "r2", append = TRUE, token = "04cad4ca4374")

proxies_df <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/CigDay_current/CigDayCurrent_WCAdjBMIsmk/proxies/combined_query_snp_list_full.txt", fill = TRUE)

#######################################################
#The proxies dataframe needs to be cleaned in two ways#
#######################################################

#a) we need to remove the same snp that we are querying by taking all those SNPs with Distance = 0 & R2 == 1.

proxies_df_clean <- proxies_df[which(!(proxies_df$query_snp == proxies_df$Coord & proxies_df$Distance == 0 & proxies_df$R2 == 1)),]

head(proxies_df_clean)

#Now we are only interested in those with a r2 > 0.8

proxies_df_clean_curated <- proxies_df_clean[which(proxies_df_clean$R2 > 0.8),]

summary(proxies_df_clean_curated$R2)

#Just in case, let's check the MAF:

summary(proxies_df_clean_curated$MAF) #perfect. No need to remove any.

cigday_proxies <- cigday[which(cigday$SNP%in%proxies_df_clean_curated$RS_Number),] #84/102
WCAdjBMI_proxies <- WCAdjBMI[which(WCAdjBMI$rs_id%in%proxies_df_clean_curated$RS_Number),] #33/102

#Now let's match the data:

cigday_proxies_match <- cigday_proxies[which(cigday_proxies$SNP%in%WCAdjBMI_proxies$rs_id),] #33/33
WCAdjBMI_proxies_match <- WCAdjBMI_proxies[which(WCAdjBMI_proxies$rs_id%in%cigday_proxies$SNP),] #33/33

#WCAdjBMI is a bottleneck and this is the best we can do right now!! However, there are 5 SNPs that are not found in Smk.
#And that makes me wonder if it is due to chr_pos...

cigday_proxies_match_2 <- cigday_proxies[which(cigday_proxies$chr_pos_37%in%WCAdjBMI_proxies$chr_pos_37),] #33/33

###############################################
#We are gonna run this again, for chr_pos only#
###############################################

cigday_proxies_chr_pos <- cigday[which(cigday$chr_pos_37%in%proxies_df_clean_curated$Coord),] #86/102
WCAdjBMI_proxies_chr_pos <- WCAdjBMI[which(WCAdjBMI$chr_pos_37%in%proxies_df_clean_curated$Coord),] #33/110

#Now let's match the data:

cigday_proxies_match_chr_pos <- cigday_proxies_chr_pos[which(cigday_proxies_chr_pos$chr_pos_37%in%WCAdjBMI_proxies_chr_pos$chr_pos_37),] #33/33
WCAdjBMI_proxies_match_chr_pos <- WCAdjBMI_proxies_chr_pos[which(WCAdjBMI_proxies_chr_pos$chr_pos_37%in%cigday_proxies_chr_pos$chr_pos_37),] #33/33

#We have no doubt about it. 33 it is.

##################
#CHECKING ALLELES#
##################

proxies_check <- proxies_df_clean_curated[which(proxies_df_clean_curated$Coord%in%WCAdjBMI_proxies_match_chr_pos$chr_pos_37),]

proxies_check <- proxies_check[order(proxies_check$Alleles),]

head(proxies_check)
tail(proxies_check)

#It is the same! The mismatch is NOT due to RSIDs, but just because the variant is no there.

###################################################################
#SECTION D: combined pre and post proxy data and get the lead SNPs#
###################################################################

cigday_proxies_match_chr_pos$rsid <- cigday_proxies_match_chr_pos$SNP
cigday_proxies_match_chr_pos$pval <- cigday_proxies_match_chr_pos$pval.exposure

#Now, be wary..., we need to get those that are not missing to avoid removing the proxies we just got.

cigday_gw_ind_not_missing <- cigday_gw_ind[which(cigday_gw_ind$chr_pos_37%in%WCAdjBMI_match_by_chr_pos$chr_pos_37),]

cigday_proxies_match_chr_pos_clean <- cigday_proxies_match_chr_pos %>%
  select(pval.exposure, chr.exposure, se.exposure, pos.exposure, beta.exposure,
         SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure,
         rsid, pval, chr_, chr_pos_37)

#We need to add some stuff:

cigday_proxies_match_chr_pos_clean$samplesize.exposure <- 	33229
cigday_proxies_match_chr_pos_clean$id.exposure <- cigday_gw_ind_not_missing$id.exposure[1]
cigday_proxies_match_chr_pos_clean$exposure <- cigday_gw_ind_not_missing$exposure[1]
cigday_proxies_match_chr_pos_clean$pval_origin.exposure <- cigday_gw_ind_not_missing$pval_origin.exposure[1]
cigday_proxies_match_chr_pos_clean$data_source.exposure <- cigday_gw_ind_not_missing$data_source.exposure[1]
cigday_proxies_match_chr_pos_clean$mr_keep.exposure <- cigday_gw_ind_not_missing$mr_keep.exposure[1]

cigday_gw_ind_plus_proxies <- rbind(cigday_gw_ind_not_missing, cigday_proxies_match_chr_pos_clean) #we get 1 + 33 SNPs = 34. Makes total sense.

#When we run the ld_clump we will get only the lead SNPs in each loci!
#Which is exactly what we want. 

cigday_gw_ind_post_proxies <- ieugwasr::ld_clump_local(cigday_gw_ind_plus_proxies, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#We end up with 3/4 SNPs! Not even sad.

################################
#Checking that this is the case#
################################

test_original <- cigday_gw_ind[order(cigday_gw_ind$chr_pos_37),]
test_post_proxy <- cigday_gw_ind_post_proxies[order(cigday_gw_ind_post_proxies$chr_pos_37),]

test_original$chr_pos_37
test_post_proxy$chr_pos_37

##############################
#Worked like a freaking charm#
##############################

##############################################
#SECTION E: Merging exposure and outcome data#
##############################################

#We need to format the data for the outcome:

WCAdjBMI_outcome <- WCAdjBMI

colnames(WCAdjBMI_outcome) <- c("chr.outcome", "SNP", "rs_id", "position_hg18",
                                 "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "samplesize.outcome",
                                 "beta.outcome", "se.outcome", "pval.outcome", "N_NonSmk", 
                                 "Effect_NonSMK", "StdErr_NonSMK", "P_value_NonSMK", "chr_pos_37",
                                 "chr_pos_18")
#Getting ids, so troublesome...

cigday_gw_ind_post_proxies$id.exposure <- "cigarettes per day (current)"
WCAdjBMI_outcome$id.outcome <- "WCadjBMI"
cigday_gw_ind_post_proxies$exposure <- "cigarettes per day (current)"
WCAdjBMI_outcome$outcome <- "WCadjBMI"

dat_1_pre <- harmonise_data(cigday_gw_ind_post_proxies, WCAdjBMI_outcome, action = 3) #3/3: perfect.

########################################
#SECTION F: filtering harmonise results#
########################################

dat_1_filt <- dat_1_pre[which(dat_1_pre$remove == FALSE),] #3/3
dat_1_filt <- dat_1_filt[which(dat_1_filt$palindromic == FALSE & dat_1_filt$ambiguous == FALSE | dat_1_filt$palindromic == TRUE & dat_1_filt$ambiguous == FALSE),] #3/3

#Despite being a binary trait... the units that we have are SD. 
#So we are gonna use the steiger filtering like that.

dat_1_filt$units.exposure <- "SD"
dat_1_filt$units.outcome <- "SD"
dat_1_filt <- steiger_filtering(dat_1_filt) 
dat_1_filt <- dat_1_filt[which(dat_1_filt$steiger_dir == TRUE),] #3/3

#We might have to repeat this.

##############################################
#SECTION G: Checking validity of the variants#
##############################################

#Are the instruments valid??
#Let's calculate mF and Isq to check it.

F_ = ((dat_1_filt$beta.exposure)^2)/((dat_1_filt$se.exposure)^2)
mF  = mean(F_)

print(mF)
#79.60415

Isq(dat_1_filt$beta.exposure, dat_1_filt$se.exposure)
#0.9875545

#########################################
#SECTION H: First round of 2SMR analysis#
#########################################

mr(dat_1_filt)

#id.exposure id.outcome  outcome exposure                    method nsnp         b        se       pval
#1      cigday   WCAdjBMI WCAdjBMI   cigday                  MR Egger    3 0.6040342 0.4954444 0.43732760
#2      cigday   WCAdjBMI WCAdjBMI   cigday           Weighted median    3 0.1956761 0.1095432 0.07405212
#3      cigday   WCAdjBMI WCAdjBMI   cigday Inverse variance weighted    3 0.1114326 0.1771586 0.52934960
#4      cigday   WCAdjBMI WCAdjBMI   cigday               Simple mode    3 0.1791616 0.1722066 0.40741563
#5      cigday   WCAdjBMI WCAdjBMI   cigday             Weighted mode    3 0.2082756 0.1100587 0.19896702

#No significant association, but they all agree on the causal direction.

#####################################################
#SECTION I: ANALYSIS ON PLEIOTROPY AND HETEROGENEITY#
#####################################################

library(meta)

dat_1_filt$mr <- dat_1_filt$beta.outcome/dat_1_filt$beta.exposure
dat_1_filt$mr_se <- ((dat_1_filt$mr*((dat_1_filt$se.exposure/dat_1_filt$beta.exposure)^2+(dat_1_filt$se.outcome/dat_1_filt$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_filt$mr, dat_1_filt$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0.0971 [0.0000; 6.5144]; tau = 0.3116 [0.0000; 2.5523];
#I^2 = 65.3% [0.0%; 90.0%]; H = 1.70 [1.00; 3.17]

#Test of heterogeneity:
#  Q d.f. p-value
#5.76    2  0.0562

#No heterogeneity. 
#Let's check the rucker framework:

mr_rucker(dat_1_filt)

#[[1]]$intercept
#Method  Estimate         SE      CI_low      CI_upp          P
#1  Egger fixed effects -0.029828 0.01613989 -0.06146159 0.001805598 0.06458901
#2 Egger random effects -0.029828 0.02813618 -0.08497390 0.025317904 0.85545717

#[[1]]$Q
#Method        Q df          P
#1   Q_ivw 6.454436  2 0.03966771
#2 Q_egger 3.038991  1 0.08128644
#3  Q_diff 3.415445  1 0.06458901

#[[1]]$res
#[1] "B"

#[[1]]$selected
#Method nsnp  Estimate        SE     CI_low    CI_upp         P
#2 Rucker    3 0.1114326 0.1771586 -0.2357919 0.4586571 0.5293496

#######################################
#SECTION J: GENERATE SENSITIVITY PLOTS#
#######################################

dat_1_filt$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/CigDay_current/CigDayCurrent_WCAdjBMIsmk/CigDayCurrent_WCAdjBMIsmk_original.tiff", units="in", width=10, height=10, res=300)
mr_plots(dat_1_filt)
dev.off()

#There seems to be a SNP that is too powerful and might be inducing pleiotropy.
#Let's check with RadialMR:

##############################
#SECTION K: Removing outliers#
##############################

library(RadialMR)

radial_input <- format_radial(dat_1_filt$beta.exposure, dat_1_filt$beta.outcome, dat_1_filt$se.exposure, dat_1_filt$se.outcome, dat_1_filt$SNP)

radial_output <- ivw_radial(radial_input, alpha = 0.05, weights = 3, tol = 0.0001) 

outliers <- radial_output$outliers$SNP #NO OUTLIER.

#There is an extreme allele, but nothing that we can do about it. 
#It is not significant, either way. 
