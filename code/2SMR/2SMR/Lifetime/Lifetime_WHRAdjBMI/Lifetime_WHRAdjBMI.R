##############
#INTRODUCTION#
##############

#This is the 2SMR code to run lifetime_smk (Ever Never Smoker) > whr.

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

######################################################
#SECTION A: FINDING INDEPENDENT, GW lifetime_smk SNPs#
######################################################

#Data from Wotton et al.

##################
#1. Load the data#
##################

lifetime_smk <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/Lifetime_smk/Lifetime_Smoking_Curated.txt")

#########################################
#2. Get the genome-wide significant SNPs#
#########################################

lifetime_smk_gw <- lifetime_smk[which(lifetime_smk$P <= 0.00000005),] #10397 if <, 10441 if <=. This is an expectional case since we lose the 126th SNP reported.

#Let's check that the filtering is done correctly...

summary(lifetime_smk_gw$P) #perfect.

#And that all the SNPs present RSIDs that we can use:

lifetime_smk_gw <- lifetime_smk_gw[order(lifetime_smk_gw$SNP),]

head(lifetime_smk_gw) #all of them with RSIDs.
tail(lifetime_smk_gw) #all of them with RSIDs.

#Let's check the eaf for a sec, just in case:

summary(lifetime_smk_gw$EAF) #perfect.

#And we do a final check to avoid SNPs in the MHC region:

which(lifetime_smk_gw$CHR == 6 & lifetime_smk_gw$BP >= 26000000 & lifetime_smk_gw$BP <= 34000000) #as expected.

#We removed them properly in the curation stage.
#Perfect!!

#########################################################
#3. calculate which SNPs are independent from each other#
#########################################################

#we get the data for the pruning

lifetime_smk_gw$rsid <- lifetime_smk_gw$SNP
lifetime_smk_gw$pval <- lifetime_smk_gw$P

#And we perform the pruning:

lifetime_smk_gw_ind <- ieugwasr::ld_clump_local(lifetime_smk_gw, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/R-4.1.3/library/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#We have some variants that are not present in the reference panel, but it is a necessary evil, I guess.

#We end up with 125 variants. However we are supposed to get 126.
#Let's check the data reported by Wottoon:

lifetime_smk_gw_og <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/RAW_Data/Lifetime_gw.txt")

#Let's do some checks:

length(which(lifetime_smk_gw_og$SNP%in%lifetime_smk_gw_ind$SNP)) #116/126. This can happen depending on the 1000G reference panel used. So no worries.
length(which(lifetime_smk_gw_og$SNP%in%lifetime_smk_gw$SNP)) #125/126. This can happen depending on the 1000G reference panel used. So no worries.

#There is one that cannot be found there.
#Most probably due to MAF or INFO. Let's check.

lifetime_smk_gw_og[which(!(lifetime_smk_gw_og$SNP%in%lifetime_smk_gw$SNP)),] #125/126. This can happen depending on the 1000G reference panel used. So no worries.

#SNP CHR       BP EFFECT_ALLELE OTHER_ALLELE      EAF INFO      BETA         SE       P
#1: rs6935954   6 26255451             A            G 0.421108    1 0.0095817 0.00140192 8.2e-12

#IT IS IN THE MHC REGION.

########################################################
#We are gonna use the 126 SNPs, but removing this fella#
########################################################

lifetime_smk_gw_ind <- lifetime_smk_gw_og[which(lifetime_smk_gw_og$SNP%in%lifetime_smk_gw$SNP),] 

summary(lifetime_smk_gw_ind$EAF) #perfect
summary(lifetime_smk_gw_ind$INFO) #perfect.

#################################################
#SECTION B: FINDING matching SNPs in the outcome#
#################################################

#First let's load the outcome:

whradjbmi <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHRAdjBMI_comb/WHRAdjBMI_Combined_Curated.txt")

##########################################################################################
#1. we are gonna try to match the data with RSID and with chr_pos to get as many variants#
##########################################################################################

whradjbmi_match_by_rsid <- whradjbmi[which(whradjbmi$MarkerName%in%lifetime_smk_gw_ind$SNP),] #52/125

###########################################################################
#Small code to get the chr_pos for lifetime_smk (full) and the gw_ind SNPs#
###########################################################################

lifetime_smk$chr_ <- paste("chr", lifetime_smk$CHR, sep = "")
lifetime_smk$chr_pos_37 <- paste(lifetime_smk$chr_, lifetime_smk$BP, sep = ":")

head(lifetime_smk)

lifetime_smk_gw_ind$chr_ <- paste("chr", lifetime_smk_gw_ind$CHR, sep = "")
lifetime_smk_gw_ind$chr_pos_37 <- paste(lifetime_smk_gw_ind$chr_, lifetime_smk_gw_ind$BP, sep = ":")

lifetime_smk_gw_ind

##############################
#Now let's match with chr_pos#
##############################

whradjbmi_match_by_chr_pos <- whradjbmi[which(whradjbmi$chr_pos_37%in%lifetime_smk_gw_ind$chr_pos_37),] #52/125

######################################################
#SECTION C: get proxies for the SNPs that are missing#
######################################################

#First we are gonna get the SNPs that are missing.

lifetime_smk_gw_ind_missing <- lifetime_smk_gw_ind[which(!(lifetime_smk_gw_ind$chr_pos_37%in%whradjbmi_match_by_chr_pos$chr_pos_37)),]

#And now we are going to get the proxies from WHR since they are bound to be the same as WHRAjdBMI

proxies_df <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Lifetime/Lifetime_WHR/proxies/combined_query_snp_list.txt", fill = TRUE)

#Let's check all the missing SNPs are in the query_snps:

query_snps <- proxies_df[which(proxies_df$query_snp == proxies_df$Coord & proxies_df$Distance == 0 & proxies_df$R2 == 1),]

length(which(query_snps$query_snp%in%lifetime_smk_gw_ind_missing$chr_pos_37)) #73/73 we can proceed.

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

lifetime_smk_proxies <- lifetime_smk[which(lifetime_smk$chr_pos_37%in%proxies_df_clean_curated$Coord),] #1917/2321
whradjbmi_proxies <- whradjbmi[which(whradjbmi$chr_pos_37%in%proxies_df_clean_curated$Coord),] #696/2321

#Now let's match the data:

lifetime_smk_proxies_match <- lifetime_smk_proxies[which(lifetime_smk_proxies$chr_pos_37%in%whradjbmi_proxies$chr_pos_37),] #691/696
whradjbmi_proxies_match <- whradjbmi_proxies[which(whradjbmi_proxies$chr_pos_37%in%lifetime_smk_proxies$chr_pos_37),] #691/691

#whr is a bottleneck and this is the best we can do right now!! However, there are 4 SNPs that are not found in Smk.

lifetime_smk_proxies_match_2 <- lifetime_smk_proxies[which(lifetime_smk_proxies$chr_pos_37%in%whradjbmi_proxies$chr_pos_37),] #691/691

#It is the same! The mismatch is NOT due to RSIDs, but just because the variant is no there.

##################
#CHECKING ALLELES#
##################

proxies_check <- proxies_df_clean_curated[which(proxies_df_clean_curated$Coord%in%lifetime_smk_proxies_match_2$chr_pos_37),]

proxies_check <- proxies_check[order(proxies_check$Alleles),]

head(proxies_check)
tail(proxies_check)

###################################################################
#SECTION D: combined pre and post proxy data and get the lead SNPs#
###################################################################

lifetime_smk_proxies_match <- lifetime_smk_proxies_match %>%
  select("SNP",           "CHR",           "BP",            "EFFECT_ALLELE", "OTHER_ALLELE",  "EAF",          
         "INFO",          "BETA",          "SE",            "P",             "chr_",          "chr_pos_37")  

#Now, be wary..., we need to get those that are not missing to avoid removing the proxies we just got.

lifetime_smk_gw_ind_not_missing <- lifetime_smk_gw_ind[which(lifetime_smk_gw_ind$chr_pos_37%in%whradjbmi_match_by_chr_pos$chr_pos_37),]

lifetime_smk_gw_ind_plus_proxies <- rbind(lifetime_smk_gw_ind_not_missing, lifetime_smk_proxies_match) #we get 52 + 699 = 751 SNPs. Makes total sense.

#When we run the ld_clump we will get only the lead SNPs in each loci!
#Which is exactly what we want. 

lifetime_smk_gw_ind_plus_proxies$rsid <- lifetime_smk_gw_ind_plus_proxies$SNP
lifetime_smk_gw_ind_plus_proxies$pval <- lifetime_smk_gw_ind_plus_proxies$P

lifetime_smk_gw_ind_post_proxies <- ieugwasr::ld_clump_local(lifetime_smk_gw_ind_plus_proxies, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/R-4.1.3/library/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#We end up with 104 SNPs! We could not retrieve 125, but I think this is rather good.

##############################################
#SECTION E: Merging exposure and outcome data#
##############################################

#We need to format the data for the outcome:

colnames(lifetime_smk_gw_ind_post_proxies) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
                                                "other_allele.exposure", "eaf.exposure", "info.exposure", "beta.exposure",
                                                "se.exposure", "pval.exposure", "chr_", "chr_pos", "rsid", "pval")


whradjbmi_outcome <- whradjbmi

colnames(whradjbmi_outcome) <- c("SNP", "effect_allele.outcome", 
                           "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                           "pval.outcome", "samplesize.outcome", "chr_pos_37")


#Getting ids, so troublesome...

lifetime_smk_gw_ind_post_proxies$id.exposure <- "lifetime smoking"
whradjbmi_outcome$id.outcome <- "WHRadjBMI"
lifetime_smk_gw_ind_post_proxies$exposure <- "lifetime smoking"
whradjbmi_outcome$outcome <- "WHRadjBMI"

dat_1_pre <- harmonise_data(lifetime_smk_gw_ind_post_proxies, whradjbmi_outcome, action = 3) #105/105: perfect.

########################################
#SECTION F: filtering harmonise results#
########################################

dat_1_filt <- dat_1_pre[which(dat_1_pre$remove == FALSE),] #104/104
dat_1_filt <- dat_1_filt[which(dat_1_filt$palindromic == FALSE & dat_1_filt$ambiguous == FALSE | dat_1_filt$palindromic == TRUE & dat_1_filt$ambiguous == FALSE),] #100/104

#So we are gonna use the steiger filtering like that:

dat_1_filt$samplesize.exposure <- 462690
dat_1_filt$units.exposure <- "SD"
dat_1_filt$units.outcome <- "SD"

dat_1_filt <- steiger_filtering(dat_1_filt) 
dat_1_filt_1 <- dat_1_filt[which(dat_1_filt$steiger_dir == TRUE),] #88/102

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
#45.72774

Isq(dat_1_filt$beta.exposure, dat_1_filt$se.exposure)
#0.9782794

#########################################
#SECTION H: First round of 2SMR analysis#
#########################################

mr(dat_1_filt)

#id.exposure id.outcome   outcome         exposure                    method nsnp           b         se       pval
#1 Lifetime smoking  WHRadjBMI WHRadjBMI Lifetime smoking                  MR Egger   88  0.25327887 0.21828603 0.24913469
#2 Lifetime smoking  WHRadjBMI WHRadjBMI Lifetime smoking           Weighted median   88  0.07532565 0.07633708 0.32376494
#3 Lifetime smoking  WHRadjBMI WHRadjBMI Lifetime smoking Inverse variance weighted   88  0.13543279 0.05322250 0.01093863
#4 Lifetime smoking  WHRadjBMI WHRadjBMI Lifetime smoking               Simple mode   88 -0.18804095 0.22756785 0.41089023
#5 Lifetime smoking  WHRadjBMI WHRadjBMI Lifetime smoking             Weighted mode   88 -0.06731967 0.17492961 0.70129615

#####################################################
#SECTION I: ANALYSIS ON PLEIOTROPY AND HETEROGENEITY#
#####################################################

library(meta)

dat_1_filt$mr <- dat_1_filt$beta.outcome/dat_1_filt$beta.exposure
dat_1_filt$mr_se <- ((dat_1_filt$mr*((dat_1_filt$se.exposure/dat_1_filt$beta.exposure)^2+(dat_1_filt$se.outcome/dat_1_filt$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_filt$mr, dat_1_filt$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0.0429 [0.0000; 0.1421]; tau = 0.2072 [0.0000; 0.3770];
#I^2 = 17.6% [0.0%; 37.5%]; H = 1.10 [1.00; 1.26]

#Test of heterogeneity:
#  Q d.f. p-value
#104.32   86  0.0872

#No heterogeneity. 
#Let's check the rucker framework:

mr_rucker(dat_1_filt)

#[[1]]$intercept
#Method     Estimate          SE       CI_low      CI_upp         P
#1  Egger fixed effects -0.001265675 0.002005578 -0.005196535 0.002665185 0.5279899
#2 Egger random effects -0.001265675 0.002273080 -0.005720829 0.003189479 0.7111716

#[[1]]$Q
#Method           Q df          P
#1   Q_ivw 110.8693663 87 0.04309531
#2 Q_egger 110.4711075 86 0.03892776
#3  Q_diff   0.3982588  1 0.52798985

#[[1]]$res
#[1] "B"

#[[1]]$selected
#Method nsnp  Estimate        SE     CI_low   CI_upp          P
#2 Rucker   88 0.1354328 0.0532225 0.03111862 0.239747 0.01093863

#Rucker's way of calculating the IVW here seems funny. MR and meta both agree that the IVW is 0.15~, but
#not significant. Also Rucker suggest that there is pleiotropy, but what is going on is that we have three
#SNPs only. 
#Very difficult to say.

#######################################
#SECTION J: GENERATE SENSITIVITY PLOTS#
#######################################

dat_1_filt$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Lifetime/Lifetime_WHRAdjBMI/Lifetime_smk_whradjbmi_original_20082021.tiff", units="in", width=10, height=10, res=300)
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

outliers <- radial_output$outliers$SNP

#rs2867112: outlier. BMI
#rs4543592: not an outlier. Current smoking (CARM1P1)
#rs6778080: outlier. Impedance.
#rs8078228: outlier? Related to asthma.
#rs889398: arm fat. outlier.

dat_1_post <- dat_1_filt[which(!(dat_1_filt$SNP%in%outliers)),]

fwrite(dat_1_post, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/output/1_curated_data/lifetime_smoking_whradjbmi_curated.txt")

################################################
#SECTION G.2: Checking validity of the variants#
################################################

#Are the instruments valid??
#Let's calculate mF and Isq to check it.

F_ = ((dat_1_post$beta.exposure)^2)/((dat_1_post$se.exposure)^2)
mF  = mean(F_)

print(mF)
#45.77518

Isq(dat_1_post$beta.exposure, dat_1_post$se.exposure)
#0.9781689

#########################################
#SECTION H.2: 2nd round of 2SMR analysis#
#########################################

mr(dat_1_post)

#id.exposure id.outcome   outcome         exposure                    method nsnp           b         se         pval
#1 Lifetime smoking  WHRadjBMI WHRadjBMI Lifetime smoking                  MR Egger   83  0.32163569 0.20223641 0.1156415390
#2 Lifetime smoking  WHRadjBMI WHRadjBMI Lifetime smoking           Weighted median   83  0.12080979 0.07629541 0.1133194359
#3 Lifetime smoking  WHRadjBMI WHRadjBMI Lifetime smoking Inverse variance weighted   83  0.17628354 0.04978219 0.0003984681
#4 Lifetime smoking  WHRadjBMI WHRadjBMI Lifetime smoking               Simple mode   83 -0.17743587 0.21675010 0.4153761951
#5 Lifetime smoking  WHRadjBMI WHRadjBMI Lifetime smoking             Weighted mode   83 -0.05890863 0.16311136 0.7189117311

#####################################################
#SECTION I: ANALYSIS ON PLEIOTROPY AND HETEROGENEITY#
#####################################################

library(meta)

dat_1_post$mr <- dat_1_post$beta.outcome/dat_1_post$beta.exposure
dat_1_post$mr_se <- ((dat_1_post$mr*((dat_1_post$se.exposure/dat_1_post$beta.exposure)^2+(dat_1_post$se.outcome/dat_1_post$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_post$mr, dat_1_post$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0 [0.0000; 0.0890]; tau = 0 [0.0000; 0.2983];
#I^2 = 0.0% [0.0%; 26.4%]; H = 1.00 [1.00; 1.17]

#Test of heterogeneity:
#  Q d.f. p-value
#80.69   81  0.4887

#No heterogeneity. 
#Let's check the rucker framework:

mr_rucker(dat_1_post)

#[[1]]$intercept
#Method     Estimate          SE       CI_low      CI_upp         P
#1  Egger fixed effects -0.001563381 0.002068367 -0.005617306 0.002490545 0.4497377
#2 Egger random effects -0.001563381 0.002107909 -0.005694807 0.002568046 0.7708574

#[[1]]$Q
#Method          Q df         P
#1   Q_ivw 84.6979617 82 0.3972948
#2 Q_egger 84.1266487 81 0.3840717
#3  Q_diff  0.5713129  1 0.4497377

#[[1]]$res
#[1] "A"

#[[1]]$selected
#Method nsnp  Estimate        SE     CI_low    CI_upp            P
#1 Rucker   83 0.1762835 0.0489829 0.08027882 0.2722883 0.0003195913

#Rucker's way of calculating the IVW here seems funny. MR and meta both agree that the IVW is 0.15~, but
#not significant. Also Rucker suggest that there is pleiotropy, but what is going on is that we have three
#SNPs only. 
#Very difficult to say.

#######################################
#SECTION J: GENERATE SENSITIVITY PLOTS#
#######################################

dat_1_post$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Lifetime/Lifetime_WHRAdjBMI/Lifetime_smk_whradjbmi_post_20082021.tiff", units="in", width=10, height=10, res=300)
mr_plots(dat_1_post)
dev.off()

###############################
#IMPORTANT TO INVESTIGATE THEM#
###############################


