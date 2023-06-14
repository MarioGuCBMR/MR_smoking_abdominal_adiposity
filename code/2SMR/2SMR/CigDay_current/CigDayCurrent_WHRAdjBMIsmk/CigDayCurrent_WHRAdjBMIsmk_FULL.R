##############
#INTRODUCTION#
##############

#This is the 2SMR code to run cigday (Ever Never Smoker) > WHRAdjBMI.

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

cigday[which(cigday$SNP == "rs11438004"),] #not here!

#########################################
#2. Get the genome-wide significant SNPs#
#########################################

#With genome-wide significant SNPs we get very few SNPs to get all the analysis going so we are using 5E-07.

cigday_gw <- extract_instruments(outcomes = "ukb-b-469", p1 = 0.0000005, clump = FALSE) #657

cigday_replication <- cigday[which(as.numeric(cigday$pval.exposure) <= 0.0000005),] #PERFECT:

#Only those in the huge dataset we curated are gonna go through the pipeline:

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

WHRAdjBMI <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMIsmk/WHRAdjBMI_Smk_Strat_Curated.txt")

##########################################################################################
#1. we are gonna try to match the data with RSID and with chr_pos to get as many variants#
##########################################################################################

WHRAdjBMI_match_by_rsid <- WHRAdjBMI[which(WHRAdjBMI$rs_id%in%cigday_gw_ind$SNP),] #1/4

#####################################################################
#Small code to get the chr_pos for cigday (full) and the gw_ind SNPs#
#####################################################################

cigday$chr_ <- paste("chr", cigday$chr.exposure, sep = "")
cigday$chr_pos_37 <- paste(cigday$chr_, cigday$pos.exposure, sep = ":")

head(cigday)

cigday_gw_ind$chr_ <- paste("chr", cigday_gw_ind$chr.exposure, sep = "")
cigday_gw_ind$chr_pos_37 <- paste(cigday_gw_ind$chr_, cigday_gw_ind$pos.exposure, sep = ":")

cigday_gw_ind

##############################
#Now let's match with chr_pos#
##############################

WHRAdjBMI_match_by_chr_pos <- WHRAdjBMI[which(WHRAdjBMI$chr_pos_37%in%cigday_gw_ind$chr_pos_3 ),] #1/4

######################################################
#SECTION C: get proxies for the SNPs that are missing#
######################################################

#First we are gonna get the SNPs that are missing.

cigday_gw_ind_missing <- cigday_gw_ind[which(!(cigday_gw_ind$chr_pos_37%in%WHRAdjBMI_match_by_chr_pos$chr_pos_37)),]

getwd()

#And now we are going to get the proxies:

LDlinkR::LDproxy_batch(snp = cigday_gw_ind_missing$chr_pos_37, pop = "EUR", r2d = "r2", append = TRUE, token = "04cad4ca4374")

proxies_df <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/CigDay_current/CigDayCurrent_WHRAdjBMIsmk/proxies/combined_query_snp_list_full.txt", fill = TRUE)

#Let's check that the data is OK:

query_snps <- proxies_df[which(proxies_df$query_snp == proxies_df$Coord & proxies_df$Distance == 0 & proxies_df$R2 == 1),]
  
length(which(query_snps$query_snp%in%cigday_gw_ind_missing$chr_pos_37)) #3/3

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

#Mostly they are in buld 37, but we have many that are merged.
#We need to be wary. Actually some of them are not even correct.
#By a couple of positions, possible due to the changes...

#Let's find the proxies in both df:

cigday_proxies <- cigday[which(cigday$SNP%in%proxies_df_clean_curated$RS_Number),] #84/102
WHRAdjBMI_proxies <- WHRAdjBMI[which(WHRAdjBMI$rs_id%in%proxies_df_clean_curated$RS_Number),] #33/102

#Now let's match the data:

cigday_proxies_match <- cigday_proxies[which(cigday_proxies$SNP%in%WHRAdjBMI_proxies$rs_id),] #33/33
WHRAdjBMI_proxies_match <- WHRAdjBMI_proxies[which(WHRAdjBMI_proxies$rs_id%in%cigday_proxies$SNP),] #33/33

#WHRAdjBMI is a bottleneck and this is the best we can do right now!! However, there are 5 SNPs that are not found in Smk.
#And that makes me wonder if it is due to chr_pos...

cigday_proxies_match_2 <- cigday_proxies[which(cigday_proxies$chr_pos_37%in%WHRAdjBMI_proxies$chr_pos_37),] #33/33

###############################################
#We are gonna run this again, for chr_pos only#
###############################################

cigday_proxies_chr_pos <- cigday[which(cigday$chr_pos_37%in%proxies_df_clean_curated$Coord),] #86/102
WHRAdjBMI_proxies_chr_pos <- WHRAdjBMI[which(WHRAdjBMI$chr_pos_37%in%proxies_df_clean_curated$Coord),] #33/110

#Now let's match the data:

cigday_proxies_match_chr_pos <- cigday_proxies_chr_pos[which(cigday_proxies_chr_pos$chr_pos_37%in%WHRAdjBMI_proxies_chr_pos$chr_pos_37),] #33/33
WHRAdjBMI_proxies_match_chr_pos <- WHRAdjBMI_proxies_chr_pos[which(WHRAdjBMI_proxies_chr_pos$chr_pos_37%in%cigday_proxies_chr_pos$chr_pos_37),] #33/33

#We have no doubt about it. 33 it is.

##################
#CHECKING ALLELES#
##################

proxies_check <- proxies_df_clean_curated[which(proxies_df_clean_curated$Coord%in%WHRAdjBMI_proxies_match_chr_pos$chr_pos_37),]

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

cigday_gw_ind_not_missing <- cigday_gw_ind[which(cigday_gw_ind$chr_pos_37%in%WHRAdjBMI_match_by_chr_pos$chr_pos_37),]

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

#We end up with 3/4 SNPs! Not even mad.

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

WHRAdjBMI_outcome <- WHRAdjBMI

colnames(WHRAdjBMI_outcome) <- c("chr.outcome", "SNP", "rs_id", "position_hg18",
                                 "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "samplesize.outcome",
                                 "beta.outcome", "se.outcome", "pval.outcome", "N_NonSmk", 
                                 "Effect_NonSMK", "StdErr_NonSMK", "P_value_NonSMK", "chr_pos_37",
                                 "chr_pos_18")
#Getting ids, so troublesome...

cigday_gw_ind_post_proxies$id.exposure <- "cigarettes per day (current)"
WHRAdjBMI_outcome$id.outcome <- "WHRadjBMI"
cigday_gw_ind_post_proxies$exposure <- "cigarettes per day (current)"
WHRAdjBMI_outcome$outcome <- "WHRadjBMI"

dat_1_pre <- harmonise_data(cigday_gw_ind_post_proxies, WHRAdjBMI_outcome, action = 3) #3/3: perfect.

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

#id.exposure id.outcome   outcome exposure                    method nsnp         b        se       pval
#1      cigday  WHRAdjBMI WHRAdjBMI   cigday                  MR Egger    3 0.2224246 0.2920199 0.58560378
#2      cigday  WHRAdjBMI WHRAdjBMI   cigday           Weighted median    3 0.2357308 0.1047012 0.02435621
#3      cigday  WHRAdjBMI WHRAdjBMI   cigday Inverse variance weighted    3 0.2395651 0.1011420 0.01785545
#4      cigday  WHRAdjBMI WHRAdjBMI   cigday               Simple mode    3 0.2176637 0.1456677 0.27371012
#5      cigday  WHRAdjBMI WHRAdjBMI   cigday             Weighted mode    3 0.2334495 0.1155911 0.18086020

#No significant association, but they all agree on the causal direction.

#####################################################
#SECTION I: ANALYSIS ON PLEIOTROPY AND HETEROGENEITY#
#####################################################

library(meta)

dat_1_filt$mr <- dat_1_filt$beta.outcome/dat_1_filt$beta.exposure
dat_1_filt$mr_se <- ((dat_1_filt$mr*((dat_1_filt$se.exposure/dat_1_filt$beta.exposure)^2+(dat_1_filt$se.outcome/dat_1_filt$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_filt$mr, dat_1_filt$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0 [0.0000; 0.1605]; tau = 0 [0.0000; 0.4007];
#I^2 = 0.0% [0.0%; 0.0%]; H = 1.00 [1.00; 1.00]

#Test of heterogeneity:
#  Q d.f. p-value
#0.17    2  0.9199

#No heterogeneity. 
#Let's check the rucker framework:

mr_rucker(dat_1_filt)

#[[1]]$intercept
#Method    Estimate          SE      CI_low     CI_upp         P
#1  Egger fixed effects 0.001041373 0.006826048 -0.01233743 0.01442018 0.8787463
#2 Egger random effects 0.001041373 0.006826048 -0.01233743 0.01442018 0.4393731

#[[1]]$Q
#Method          Q df         P
#1   Q_ivw 0.17212370  2 0.9175375
#2 Q_egger 0.16820878  1 0.6817086
#3  Q_diff 0.00391492  1 0.9501095

#[[1]]$res
#[1] "A"

#[[1]]$selected
#Method nsnp  Estimate         SE    CI_low    CI_upp            P
#1 Rucker    3 0.2395651 0.02967132 0.1814104 0.2977198 6.805236e-16

#######################################
#SECTION J: GENERATE SENSITIVITY PLOTS#
#######################################

dat_1_filt$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/CigDay_current/CigDayCurrent_WHRAdjBMIsmk/CigDayCurrent_WHRAdjBMIsmk_original.tiff", units="in", width=10, height=10, res=300)
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
