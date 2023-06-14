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

wcadjbmi <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WCAdjBMI/WCAdjBMI_Combined_Curated.txt")

##########################################################################################
#1. we are gonna try to match the data with RSID and with chr_pos to get as many variants#
##########################################################################################

wcadjbmi_match_by_rsid <- wcadjbmi[which(wcadjbmi$MarkerName%in%lifetime_smk_gw_ind$SNP),] #52/125

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

wcadjbmi_match_by_chr_pos <- wcadjbmi[which(wcadjbmi$chr_pos_37%in%lifetime_smk_gw_ind$chr_pos_37),] #52/125

######################################################
#SECTION C: get proxies for the SNPs that are missing#
######################################################

#First we are gonna get the SNPs that are missing.

lifetime_smk_gw_ind_missing <- lifetime_smk_gw_ind[which(!(lifetime_smk_gw_ind$chr_pos_37%in%wcadjbmi_match_by_chr_pos$chr_pos_37)),]

#Getting proxies from WHR, because they are the same.

proxies_df <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Lifetime/Lifetime_WHR/proxies/combined_query_snp_list.txt", fill = TRUE)

#We are gonna make sure the query_snps are the same:

query_snps <- proxies_df[which(proxies_df$query_snp == proxies_df$Coord & proxies_df$Distance == 0 & proxies_df$R2 == 1),]

length(which(query_snps$query_snp%in%lifetime_smk_gw_ind_missing$chr_pos_37)) #73/73

#The proxies dataframe needs to be cleaned in two ways:

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
wcadjbmi_proxies <- wcadjbmi[which(wcadjbmi$chr_pos_37%in%proxies_df_clean_curated$Coord),] #696/2321

#Now let's match the data:

lifetime_smk_proxies_match <- lifetime_smk_proxies[which(lifetime_smk_proxies$chr_pos_37%in%wcadjbmi_proxies$chr_pos_37),] #691/695
wcadjbmi_proxies_match <- wcadjbmi_proxies[which(wcadjbmi_proxies$chr_pos_37%in%lifetime_smk_proxies$chr_pos_37),] #691/695

#whr is a bottleneck and this is the best we can do right now!! However, there are 4 SNPs that are not found in Smk.

lifetime_smk_proxies_match_2 <- lifetime_smk_proxies[which(lifetime_smk_proxies$chr_pos_37%in%wcadjbmi_proxies$chr_pos_37),] #691/695

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

lifetime_smk_gw_ind_not_missing <- lifetime_smk_gw_ind[which(lifetime_smk_gw_ind$chr_pos_37%in%wcadjbmi_match_by_chr_pos$chr_pos_37),]

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


wcadjbmi_outcome <- wcadjbmi

colnames(wcadjbmi_outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", 
                           "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                           "pval.outcome", "samplesize.outcome", "chr_pos_37")


#Getting ids, so troublesome...

lifetime_smk_gw_ind_post_proxies$id.exposure <- "lifetime smoking"
wcadjbmi_outcome$id.outcome <- "WCadjBMI"
lifetime_smk_gw_ind_post_proxies$exposure <- "lifetime smoking"
wcadjbmi_outcome$outcome <- "WCadjBMI"

dat_1_pre <- harmonise_data(lifetime_smk_gw_ind_post_proxies, wcadjbmi_outcome, action = 3) #104/104: perfect.

########################################
#SECTION F: filtering harmonise results#
########################################

dat_1_filt <- dat_1_pre[which(dat_1_pre$remove == FALSE),] #104/104
dat_1_filt <- dat_1_filt[which(dat_1_filt$palindromic == FALSE & dat_1_filt$ambiguous == FALSE | dat_1_filt$palindromic == TRUE & dat_1_filt$ambiguous == FALSE),] #100/100

#So we are gonna use the steiger filtering like that:

dat_1_filt$samplesize.exposure <- 462690
dat_1_filt$units.exposure <- "SD"
dat_1_filt$units.outcome <- "SD"

dat_1_filt <- steiger_filtering(dat_1_filt) 
dat_1_filt_1 <- dat_1_filt[which(dat_1_filt$steiger_dir == TRUE),] #89/100

dat_1_filt <- dat_1_filt_1

dat_1_filt$mr_keep = TRUE #this keeps the palindromes out.

##############################################
#SECTION G: Checking validity of the variants#
##############################################

#Are the instruments valid??
#Let's calculate mF and Isq to check it.

F_ = ((dat_1_filt$beta.exposure)^2)/((dat_1_filt$se.exposure)^2)
mF  = mean(F_)

print(mF)
#45.43091

Isq(dat_1_filt$beta.exposure, dat_1_filt$se.exposure)
#0.9781544

#########################################
#SECTION H: First round of 2SMR analysis#
#########################################

mr(dat_1_filt)

#id.exposure id.outcome  outcome         exposure                    method nsnp         b         se        pval
#1 Lifetime smoking   WCadjBMI WCadjBMI Lifetime smoking                  MR Egger   89 0.3691878 0.19263339 0.058580521
#2 Lifetime smoking   WCadjBMI WCadjBMI Lifetime smoking           Weighted median   89 0.1727334 0.07126656 0.015360508
#3 Lifetime smoking   WCadjBMI WCadjBMI Lifetime smoking Inverse variance weighted   89 0.1372780 0.04791408 0.004168939
#4 Lifetime smoking   WCadjBMI WCadjBMI Lifetime smoking               Simple mode   89 0.1684762 0.19295187 0.384956882
#5 Lifetime smoking   WCadjBMI WCadjBMI Lifetime smoking             Weighted mode   89 0.2282976 0.14406701 0.116630305

#####################################################
#SECTION I: ANALYSIS ON PLEIOTROPY AND HETEROGENEITY#
#####################################################

library(meta)

dat_1_filt$mr <- dat_1_filt$beta.outcome/dat_1_filt$beta.exposure
dat_1_filt$mr_se <- ((dat_1_filt$mr*((dat_1_filt$se.exposure/dat_1_filt$beta.exposure)^2+(dat_1_filt$se.outcome/dat_1_filt$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_filt$mr, dat_1_filt$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0 [0.0000; 0.0807]; tau = 0 [0.0000; 0.2841];
#I^2 = 0.0% [0.0%; 24.7%]; H = 1.00 [1.00; 1.15]

#Test of heterogeneity:
#  Q d.f. p-value
#86.79   88  0.5165

#No heterogeneity. 
#Let's check the rucker framework:

mr_rucker(dat_1_filt)

#[[1]]$intercept
#Method     Estimate          SE       CI_low      CI_upp         P
#1  Egger fixed effects -0.002517562 0.001992832 -0.006423440 0.001388317 0.2064783
#2 Egger random effects -0.002517562 0.002025875 -0.006488204 0.001453080 0.8930114

#[[1]]$Q
#Method         Q df         P
#1   Q_ivw 91.504974 88 0.3779445
#2 Q_egger 89.909026 87 0.3941380
#3  Q_diff  1.595949  1 0.2064783

#[[1]]$res
#[1] "A"

#[[1]]$selected
#Method nsnp Estimate         SE     CI_low    CI_upp           P
#1 Rucker   89 0.137278 0.04698748 0.04518427 0.2293718 0.003482523

#######################################
#SECTION J: GENERATE SENSITIVITY PLOTS#
#######################################

dat_1_filt$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Lifetime/Lifetime_WCAdjBMI/Lifetime_smk_wcadjbmi_original_20082021.tiff", units="in", width=10, height=10, res=300)
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

#rs1050847: (outlier: height.) -ZCCHC14-
#rs17309874: weight/BMI: BDNF-AS 
#rs6778080: impedance: USP4.

dat_1_post <- dat_1_filt[which(!(dat_1_filt$SNP%in%outliers)),]

fwrite(dat_1_post, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/output/1_curated_data/lifetime_smoking_wcadjbmi_curated.txt")

################################################
#SECTION G.2: Checking validity of the variants#
################################################

#Are the instruments valid??
#Let's calculate mF and Isq to check it.

F_ = ((dat_1_post$beta.exposure)^2)/((dat_1_post$se.exposure)^2)
mF  = mean(F_)

print(mF)
#45.46403

Isq(dat_1_post$beta.exposure, dat_1_post$se.exposure)
#0.9781403

#########################################
#SECTION H.2: 2nd round of 2SMR analysis#
#########################################

mr(dat_1_post)

#id.exposure id.outcome  outcome         exposure                    method nsnp         b         se        pval
#1 Lifetime smoking   WCadjBMI WCadjBMI Lifetime smoking                  MR Egger   86 0.3257999 0.19082837 0.091461597
#2 Lifetime smoking   WCadjBMI WCadjBMI Lifetime smoking           Weighted median   86 0.1728547 0.07163334 0.015819740
#3 Lifetime smoking   WCadjBMI WCadjBMI Lifetime smoking Inverse variance weighted   86 0.1436221 0.04794191 0.002737675
#4 Lifetime smoking   WCadjBMI WCadjBMI Lifetime smoking               Simple mode   86 0.1672467 0.20021081 0.405860830
#5 Lifetime smoking   WCadjBMI WCadjBMI Lifetime smoking             Weighted mode   86 0.2317830 0.13197272 0.082639879

#####################################################
#SECTION I: ANALYSIS ON PLEIOTROPY AND HETEROGENEITY#
#####################################################

library(meta)

dat_1_post$mr <- dat_1_post$beta.outcome/dat_1_post$beta.exposure
dat_1_post$mr_se <- ((dat_1_post$mr*((dat_1_post$se.exposure/dat_1_post$beta.exposure)^2+(dat_1_post$se.outcome/dat_1_post$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_post$mr, dat_1_post$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0 [0.0000; 0.0482]; tau = 0 [0.0000; 0.2196];
#I^2 = 0.0% [0.0%; 14.5%]; H = 1.00 [1.00; 1.08]

#Test of heterogeneity:
#  Q d.f. p-value
#73.48   85  0.8093

#No heterogeneity. 
#Let's check the rucker framework:

mr_rucker(dat_1_post)

#[[1]]$intercept
#Method     Estimate          SE       CI_low      CI_upp         P
#1  Egger fixed effects -0.001981207 0.001908055 -0.005720926 0.001758511 0.2991123
#2 Egger random effects -0.001981207 0.001908055 -0.005720926 0.001758511 0.8504439

#[[1]]$Q
#Method          Q df         P
#1   Q_ivw 76.7642247 85 0.7263159
#2 Q_egger 75.7914345 84 0.7268459
#3  Q_diff  0.9727902  1 0.3239853

#[[1]]$res
#[1] "A"

#[[1]]$selected
#Method nsnp  Estimate         SE     CI_low    CI_upp           P
#1 Rucker   86 0.1436221 0.04556017 0.05432586 0.2329184 0.001619552

#######################################
#SECTION J: GENERATE SENSITIVITY PLOTS#
#######################################

dat_1_post$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/Lifetime/Lifetime_WCAdjBMI/Lifetime_smk_wcadjbmi_post_20082021.tiff", units="in", width=10, height=10, res=300)
mr_plots(dat_1_post)
dev.off()
