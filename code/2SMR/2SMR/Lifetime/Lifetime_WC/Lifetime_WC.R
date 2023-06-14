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

wc <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WC/WC_combined_Curated.txt")

##########################################################################################
#1. we are gonna try to match the data with RSID and with chr_pos to get as many variants#
##########################################################################################

wc_match_by_rsid <- wc[which(wc$MarkerName%in%lifetime_smk_gw_ind$SNP),] #52/125

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

wc_match_by_chr_pos <- wc[which(wc$chr_pos_37%in%lifetime_smk_gw_ind$chr_pos_37),] #52/125

######################################################
#SECTION C: get proxies for the SNPs that are missing#
######################################################

#First we are gonna get the SNPs that are missing.

lifetime_smk_gw_ind_missing <- lifetime_smk_gw_ind[which(!(lifetime_smk_gw_ind$chr_pos_37%in%wc_match_by_chr_pos$chr_pos_37)),]

#We are gonna obtain the data from WHR since the SNPs are gonna be the same.

proxies_df <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Lifetime/Lifetime_WHR/proxies/combined_query_snp_list.txt", fill = TRUE)

#Let's check if the query is cool:

query_snps <- proxies_df[which(proxies_df$query_snp == proxies_df$Coord & proxies_df$Distance == 0 & proxies_df$R2 == 1),]

length(which(query_snps$query_snp%in%lifetime_smk_gw_ind_missing$chr_pos_37)) #73/73. ALL GOOD.

#The proxies dataframe needs to be cleaned in two ways:

#a) we need to remove the same snp that we are querying by taking all those SNPs with Distance = 0 & R2 == 1.

proxies_df_clean <- proxies_df[which(!(proxies_df$query_snp == proxies_df$Coord & proxies_df$Distance == 0 & proxies_df$R2 == 1)),]

head(proxies_df_clean)

#Now we are only interested in those with a r2 > 0.8

proxies_df_clean_curated <- proxies_df_clean[which(proxies_df_clean$R2 > 0.8),]

#Just in case, let's check the MAF:

summary(proxies_df_clean_curated$MAF) #perfect. No need to remove any.

################################################
#Let's check if there are any in the MHC region#
################################################

proxies_chr6 <- proxies_df_clean_curated[order(proxies_df_clean_curated$Coord),]

View(proxies_chr6) #seems like all is good.

#Mostly they are in buld 37, but we have many that are merged.
#We need to be wary. Actually some of them are not even correct.
#By a couple of positions, possible due to the changes...

#Let's find the proxies in both df:

lifetime_smk_proxies <- lifetime_smk[which(lifetime_smk$chr_pos_37%in%proxies_df_clean_curated$Coord),] #1893/2308
wc_proxies <- wc[which(wc$chr_pos_37%in%proxies_df_clean_curated$Coord),] #703/2308

#Now let's match the data:

lifetime_smk_proxies_match <- lifetime_smk_proxies[which(lifetime_smk_proxies$chr_pos_37%in%wc_proxies$chr_pos_37),] #699/703
wc_proxies_match <- wc_proxies[which(wc_proxies$chr_pos_37%in%lifetime_smk_proxies$chr_pos_37),] #699/703

#There are five that cannot be found.
#Why?

#c is a bottleneck and this is the best we can do right now!! However, there are 5 SNPs that are not found in Smk.

lifetime_smk_proxies_match_2 <- lifetime_smk_proxies[which(lifetime_smk_proxies$chr_pos_37%in%wc_proxies$chr_pos_37),] #699/703

#It is the same! The mismatch is NOT due to RSIDs, but just because the variant is no there.

##################
#Checking alleles#
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

lifetime_smk_gw_ind_not_missing <- lifetime_smk_gw_ind[which(lifetime_smk_gw_ind$chr_pos_37%in%wc_match_by_chr_pos$chr_pos_37),]

lifetime_smk_gw_ind_plus_proxies <- rbind(lifetime_smk_gw_ind_not_missing, lifetime_smk_proxies_match) #we get 52 + 699 = 751 SNPs. Makes total sense.

#When we run the ld_clump we will get only the lead SNPs in each loci!
#Which is exactly what we want. 

lifetime_smk_gw_ind_plus_proxies$rsid <- lifetime_smk_gw_ind_plus_proxies$SNP
lifetime_smk_gw_ind_plus_proxies$pval <- lifetime_smk_gw_ind_plus_proxies$P

#When we run the ld_clump we will get only the lead SNPs in each loci!
#Which is exactly what we want. 

lifetime_smk_gw_ind_post_proxies <- ieugwasr::ld_clump_local(lifetime_smk_gw_ind_plus_proxies, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/R-4.1.3/library/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#We end up with 104 SNPs! We could not retrieve 125, but I think this is rather good.

##############################################
#SECTION E: Merging exposure and outcome data#
##############################################

#We need to format the data for the outcome:

colnames(lifetime_smk_gw_ind_post_proxies) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
                                                "other_allele.exposure", "eaf.exposure", "info.exposure", "beta.exposure",
                                                "se.exposure", "pval.exposure", "chr_", "chr_pos", "rsid", "pval")


wc_outcome <- wc

colnames(wc_outcome) <- c("SNP", "effect_allele.outcome", 
                           "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                           "pval.outcome", "samplesize.outcome", "chr_pos_37")


#Getting ids, so troublesome...

lifetime_smk_gw_ind_post_proxies$id.exposure <- "lifetime smoking"
wc_outcome$id.outcome <- "WC"
lifetime_smk_gw_ind_post_proxies$exposure <- "lifetime smoking"
wc_outcome$outcome <- "WC"

dat_1_pre <- harmonise_data(lifetime_smk_gw_ind_post_proxies, wc_outcome, action = 3) #104/104: perfect.

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
dat_1_filt_1 <- dat_1_filt[which(dat_1_filt$steiger_dir == TRUE),] #80/100

dat_1_filt <- dat_1_filt_1

dat_1_filt$mr_keep = TRUE #to keep unambiguous palindromes

##############################################
#SECTION G: Checking validity of the variants#
##############################################

#Are the instruments valid??
#Let's calculate mF and Isq to check it.

F_ = ((dat_1_filt$beta.exposure)^2)/((dat_1_filt$se.exposure)^2)
mF  = mean(F_)

print(mF)
#44.97587

Isq(dat_1_filt$beta.exposure, dat_1_filt$se.exposure)
#0.977956

#########################################
#SECTION H: First round of 2SMR analysis#
#########################################

mr(dat_1_filt)

#id.exposure id.outcome outcome         exposure                    method nsnp           b         se      pval
#1 Lifetime smoking         WC      WC Lifetime smoking                  MR Egger   80 -0.18173681 0.21619094 0.4031242
#2 Lifetime smoking         WC      WC Lifetime smoking           Weighted median   80  0.02023217 0.07925368 0.7985040
#3 Lifetime smoking         WC      WC Lifetime smoking Inverse variance weighted   80  0.05509914 0.05608067 0.3258546
#4 Lifetime smoking         WC      WC Lifetime smoking               Simple mode   80  0.16733247 0.22844418 0.4660389
#5 Lifetime smoking         WC      WC Lifetime smoking             Weighted mode   80 -0.12309700 0.17123629 0.4743407

#####################################################
#SECTION I: ANALYSIS ON PLEIOTROPY AND HETEROGENEITY#
#####################################################

library(meta)

dat_1_filt$mr <- dat_1_filt$beta.outcome/dat_1_filt$beta.exposure
dat_1_filt$mr_se <- ((dat_1_filt$mr*((dat_1_filt$se.exposure/dat_1_filt$beta.exposure)^2+(dat_1_filt$se.outcome/dat_1_filt$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_filt$mr, dat_1_filt$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0.0319 [0.0000; 0.1177]; tau = 0.1785 [0.0000; 0.3430];
#I^2 = 13.0% [0.0%; 34.8%]; H = 1.07 [1.00; 1.24]

#Test of heterogeneity:
#  Q d.f. p-value
#90.85   79  0.1706

#Very few heterogeneity.

mr_rucker(dat_1_filt)

#[[1]]$intercept
#Method    Estimate          SE       CI_low      CI_upp         P
#1  Egger fixed effects 0.002524146 0.002017186 -0.001429466 0.006477757 0.2108176
#2 Egger random effects 0.002524146 0.002225534 -0.001837820 0.006886112 0.1283605

#[[1]]$Q
#Method         Q df          P
#1   Q_ivw 96.510589 79 0.08789279
#2 Q_egger 94.944786 78 0.09312316
#3  Q_diff  1.565803  1 0.21081758

#[[1]]$res
#[1] "A"

#[[1]]$selected
#Method nsnp   Estimate         SE      CI_low    CI_upp         P
#1 Rucker   80 0.05509914 0.05073869 -0.04434687 0.1545452 0.2775059

#######################################
#SECTION J: GENERATE SENSITIVITY PLOTS#
#######################################

dat_1_filt$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Lifetime/Lifetime_WC/Lifetime_smk_wc_original_20082021.tiff", units="in", width=10, height=10, res=300)
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

outliers <- radial_output$outliers$SNP

#rs963354: sleep duration/past smoking. Maybe not it: PQLC2L.
#rs8078228: outlier! Asthma and (PGAP3)
#rs9919670: not much of an outlier: (NCAM1).

dat_1_post <- dat_1_filt[which(!(dat_1_filt$SNP%in%outliers)),]

fwrite(dat_1_post, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/output/1_curated_data/lifetime_smoking_wc_curated.txt")

################################################
#SECTION G.2: Checking validity of the variants#
################################################

#Are the instruments valid??
#Let's calculate mF and Isq to check it.

F_ = ((dat_1_post$beta.exposure)^2)/((dat_1_post$se.exposure)^2)
mF  = mean(F_)

print(mF)
#44.37687

Isq(dat_1_post$beta.exposure, dat_1_post$se.exposure)
#0.9777174

#########################################
#SECTION H.2: 2nd round of 2SMR analysis#
#########################################

mr(dat_1_post)

#id.exposure id.outcome outcome         exposure                    method nsnp           b         se       pval
#1 Lifetime smoking         WC      WC Lifetime smoking                  MR Egger   77 -0.09866058 0.20303940 0.62844172
#2 Lifetime smoking         WC      WC Lifetime smoking           Weighted median   77  0.05137629 0.07902031 0.51558526
#3 Lifetime smoking         WC      WC Lifetime smoking Inverse variance weighted   77  0.11369601 0.05232016 0.02977413
#4 Lifetime smoking         WC      WC Lifetime smoking               Simple mode   77  0.17160297 0.21421505 0.42558384
#5 Lifetime smoking         WC      WC Lifetime smoking             Weighted mode   77 -0.13855911 0.16984730 0.41717178

#####################################################
#SECTION I: ANALYSIS ON PLEIOTROPY AND HETEROGENEITY#
#####################################################

library(meta)

dat_1_post$mr <- dat_1_post$beta.outcome/dat_1_post$beta.exposure
dat_1_post$mr_se <- ((dat_1_post$mr*((dat_1_post$se.exposure/dat_1_post$beta.exposure)^2+(dat_1_post$se.outcome/dat_1_post$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_post$mr, dat_1_post$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0 [0.0000; 0.0686]; tau = 0 [0.0000; 0.2618];
#I^2 = 0.0% [0.0%; 23.1%]; H = 1.00 [1.00; 1.14]

#Test of heterogeneity:
#  Q d.f. p-value
#71.77   76  0.6161

#No heterogeneity. 
#Let's check the rucker framework:

mr_rucker(dat_1_post)

#[[1]]$intercept
#Method    Estimate          SE       CI_low      CI_upp         P
#1  Egger fixed effects 0.002258453 0.002075482 -0.001809418 0.006326323 0.2765253
#2 Egger random effects 0.002258453 0.002075482 -0.001809418 0.006326323 0.1382627

#[[1]]$Q
#Method         Q df         P
#1   Q_ivw 75.386038 76 0.4983203
#2 Q_egger 74.214354 75 0.5039327
#3  Q_diff  1.171684  1 0.2790554

#[[1]]$res
#[1] "A"

#[[1]]$selected
#Method nsnp Estimate        SE     CI_low    CI_upp          P
#1 Rucker   77 0.113696 0.0521084 0.01156542 0.2158266 0.02911593

#######################################
#SECTION J: GENERATE SENSITIVITY PLOTS#
#######################################

dat_1_post$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Lifetime/Lifetime_WC/Lifetime_smk_wc_post_20082021.tiff", units="in", width=10, height=10, res=300)
mr_plots(dat_1_post)
dev.off()
