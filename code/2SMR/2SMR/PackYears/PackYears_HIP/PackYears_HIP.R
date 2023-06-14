##############
#INTRODUCTION#
##############

#This is the 2SMR code to run packyears_smk > hip.

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
#SECTION A: FINDING INDEPENDENT, GW packyears_smk SNPs#
######################################################

#Data from Wotton et al.

##################
#1. Load the data#
##################

packyears_smk <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/PacksPerYear/PackPerYear_Curated_FULL.txt")

#########################################
#2. Get the genome-wide significant SNPs#
#########################################

packyears_smk_gw <- packyears_smk[which(packyears_smk$pval.exposure <= 0.00000005),] 

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

hip <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/HIP/hip_combined_Curated.txt")

##########################################################################################
#1. we are gonna try to match the data with RSID and with chr_pos to get as many variants#
##########################################################################################

hip_match_by_rsid <- hip[which(hip$MarkerName%in%packyears_smk_gw_ind$SNP),] #7/13

###########################################################################
#Small code to get the chr_pos for packyears_smk (full) and the gw_ind SNPs#
###########################################################################

packyears_smk$chr_ <- paste("chr", packyears_smk$chr.exposure, sep = "")
packyears_smk$chr_pos_37 <- paste(packyears_smk$chr_, packyears_smk$pos.exposure, sep = ":")

head(packyears_smk)

packyears_smk_gw_ind$chr_ <- paste("chr", packyears_smk_gw_ind$chr.exposure, sep = "")
packyears_smk_gw_ind$chr_pos_37 <- paste(packyears_smk_gw_ind$chr_, packyears_smk_gw_ind$pos.exposure, sep = ":")

packyears_smk_gw_ind

##############################
#Now let's match with chr_pos#
##############################

hip_match_by_chr_pos <- hip[which(hip$chr_pos_37%in%packyears_smk_gw_ind$chr_pos_37),] #52/125

######################################################
#SECTION C: get proxies for the SNPs that are missing#
######################################################

#First we are gonna get the SNPs that are missing.

packyears_smk_gw_ind_missing <- packyears_smk_gw_ind[which(!(packyears_smk_gw_ind$chr_pos_37%in%hip_match_by_chr_pos$chr_pos_37)),]

#We are gonna obtain the data from WHR since the SNPs are gonna be the same.

proxies_df <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/PackYears/PackYears_WHRAdjBMI/proxies/combined_query_snp_list.txt")

#Let's check if the query is cool:

query_snps <- proxies_df[which(proxies_df$query_snp == proxies_df$Coord & proxies_df$Distance == 0 & proxies_df$R2 == 1),]

length(which(query_snps$query_snp%in%packyears_smk_gw_ind_missing$chr_pos_37)) #73/73. ALL GOOD.

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

packyears_smk_proxies <- packyears_smk[which(packyears_smk$chr_pos_37%in%proxies_df_clean_curated$Coord),] 
hip_proxies <- hip[which(hip$chr_pos_37%in%proxies_df_clean_curated$Coord),] 

#Now let's match the data:

packyears_smk_proxies_match <- packyears_smk_proxies[which(packyears_smk_proxies$chr_pos_37%in%hip_proxies$chr_pos_37),] #38
hip_proxies_match <- hip_proxies[which(hip_proxies$chr_pos_37%in%packyears_smk_proxies$chr_pos_37),] #38

#There are five that cannot be found.
#Why?

#c is a bottleneck and this is the best we can do right now!! However, there are 5 SNPs that are not found in Smk.

packyears_smk_proxies_match_2 <- packyears_smk_proxies[which(packyears_smk_proxies$chr_pos_37%in%hip_proxies$chr_pos_37),] #38

#It is the same! The mismatch is NOT due to RSIDs, but just because the variant is no there.

##################
#Checking alleles#
##################

proxies_check <- proxies_df_clean_curated[which(proxies_df_clean_curated$Coord%in%packyears_smk_proxies_match_2$chr_pos_37),]

proxies_check <- proxies_check[order(proxies_check$Alleles),]

head(proxies_check)
tail(proxies_check)

###################################################################
#SECTION D: combined pre and post proxy data and get the lead SNPs#
###################################################################

packyears_smk_proxies_match$pval <- packyears_smk_proxies_match$pval.exposure 

#Now, be wary..., we need to get those that are not missing to avoid removing the proxies we just got.

packyears_smk_gw_ind_not_missing <- packyears_smk_gw_ind[which(packyears_smk_gw_ind$chr_pos_37%in%hip_match_by_chr_pos$chr_pos_37),]

packyears_smk_gw_ind_plus_proxies <- rbind(packyears_smk_gw_ind_not_missing, packyears_smk_proxies_match) #we get 7 + 38 = 45 SNPs. Makes total sense.

#When we run the ld_clump we will get only the lead SNPs in each loci!
#Which is exactly what we want. 

packyears_smk_gw_ind_plus_proxies$rsid <- packyears_smk_gw_ind_plus_proxies$SNP
packyears_smk_gw_ind_plus_proxies$pval <- packyears_smk_gw_ind_plus_proxies$pval.exposure

#When we run the ld_clump we will get only the lead SNPs in each loci!
#Which is exactly what we want. 

packyears_smk_gw_ind_post_proxies <- ieugwasr::ld_clump_local(packyears_smk_gw_ind_plus_proxies, bfile = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/2SMR/1000G/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/R-4.1.2/library/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) 

#We end up with 12/13.

##############################################
#SECTION E: Merging exposure and outcome data#
##############################################

#We need to format the data for the outcome:

hip_outcome <- hip

colnames(hip_outcome) <- c("SNP", "effect_allele.outcome", 
                           "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                           "pval.outcome", "samplesize.outcome", "chr_pos_37")


#Getting ids, so troublesome...

packyears_smk_gw_ind_post_proxies$id.exposure <- "pack-years"
hip_outcome$id.outcome <- "HC"
packyears_smk_gw_ind_post_proxies$exposure <- "pack-years"
hip_outcome$outcome <- "HC"

dat_1_pre <- harmonise_data(packyears_smk_gw_ind_post_proxies, hip_outcome, action = 3) #104/104: perfect.

########################################
#SECTION F: filtering harmonise results#
########################################

dat_1_filt <- dat_1_pre[which(dat_1_pre$remove == FALSE),] #12/12
dat_1_filt <- dat_1_filt[which(dat_1_filt$palindromic == FALSE & dat_1_filt$ambiguous == FALSE | dat_1_filt$palindromic == TRUE & dat_1_filt$ambiguous == FALSE),] #10/10

#So we are gonna use the steiger filtering like that:

dat_1_filt$samplesize.exposure <- 142387
dat_1_filt$units.exposure <- "SD"
dat_1_filt$units.outcome <- "SD"

dat_1_filt <- steiger_filtering(dat_1_filt) 
dat_1_filt_1 <- dat_1_filt[which(dat_1_filt$steiger_dir == TRUE),] #81/100

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
#77.73074

Isq(dat_1_filt$beta.exposure, dat_1_filt$se.exposure)
#0.9883244

#########################################
#SECTION H: First round of 2SMR analysis#
#########################################

mr(dat_1_filt)

#id.exposure id.outcome outcome                 exposure                    method nsnp           b         se       pval
#1 Pack years adult smoking         HC      HC Pack years adult smoking                  MR Egger   10 -0.20914486 0.12689208 0.13792199
#2 Pack years adult smoking         HC      HC Pack years adult smoking           Weighted median   10 -0.10987663 0.04925751 0.02570382
#3 Pack years adult smoking         HC      HC Pack years adult smoking Inverse variance weighted   10 -0.04342099 0.06821665 0.52444014
#4 Pack years adult smoking         HC      HC Pack years adult smoking               Simple mode   10 -0.12182380 0.12533710 0.35646215
#5 Pack years adult smoking         HC      HC Pack years adult smoking             Weighted mode   10 -0.11328411 0.04821129 0.04332335

#####################################################
#SECTION I: ANALYSIS ON PLEIOTROPY AND HETEROGENEITY#
#####################################################

library(meta)

dat_1_filt$mr <- dat_1_filt$beta.outcome/dat_1_filt$beta.exposure
dat_1_filt$mr_se <- ((dat_1_filt$mr*((dat_1_filt$se.exposure/dat_1_filt$beta.exposure)^2+(dat_1_filt$se.outcome/dat_1_filt$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_filt$mr, dat_1_filt$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0.0487 [0.0064; 0.3081]; tau = 0.2208 [0.0799; 0.5551]
#I^2 = 61.7% [23.8%; 80.8%]; H = 1.62 [1.15; 2.28]

#Test of heterogeneity:
#  Q d.f. p-value
#23.50    9  0.0052

#Very few heterogeneity.

mr_rucker(dat_1_filt)

#[[1]]$intercept
#Method    Estimate          SE       CI_low     CI_upp          P
#1  Egger fixed effects 0.007443959 0.003056764  0.001452812 0.01343511 0.01488185
#2 Egger random effects 0.007443959 0.004926403 -0.002211613 0.01709953 0.06538996

#[[1]]$Q
#Method         Q df           P
#1   Q_ivw 26.709463  9 0.001562173
#2 Q_egger 20.779061  8 0.007758075
#3  Q_diff  5.930402  1 0.014881850

#[[1]]$res
#[1] "D"

#[[1]]$selected
#Method nsnp   Estimate        SE     CI_low     CI_upp        P
#4 Rucker   10 -0.2091449 0.1268921 -0.4578488 0.03955904 0.137922

#######################################
#SECTION J: GENERATE SENSITIVITY PLOTS#
#######################################

dat_1_filt$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/PackYears/PackYears_HIP/packyears_smk_hip_original_13092021.tiff", units="in", width=10, height=10, res=300)
mr_plots(dat_1_filt)
dev.off()

#There seems to be a SNP that is too powerful and might be inducing pleiotropy.
#Let's check with RadialMR:

##############################
#SECTION K: Removing outliers#
##############################

library(RadialMR)

radial_input <- format_radial(dat_1_filt$beta.exposure, dat_1_filt$beta.outcome, dat_1_filt$se.exposure, dat_1_filt$se.outcome, dat_1_filt$SNP)

radial_output <- egger_radial(radial_input, alpha = 0.05, weights = 3) 

outliers <- radial_output$outliers$SNP

#rs1364063: #arm fat percentage. It is really an outlier.
#rs150353:  Arm fat-free mass right

#The only real outlier might be rs2160316 and it might actually just 
#a smoking variants that shares pathways with hip.

dat_1_post <- dat_1_filt[which(!(dat_1_filt$SNP%in%outliers)),]

################################################
#SECTION G.2: Checking validity of the variants#
################################################

#Are the instruments valid??
#Let's calculate mF and Isq to check it.

F_ = ((dat_1_post$beta.exposure)^2)/((dat_1_post$se.exposure)^2)
mF  = mean(F_)

print(mF)
#89.58663

Isq(dat_1_post$beta.exposure, dat_1_post$se.exposure)
#0.9901934

#########################################
#SECTION H.2: 2nd round of 2SMR analysis#
#########################################

mr(dat_1_post)

#id.exposure id.outcome outcome                 exposure                    method nsnp           b         se       pval
#1 Pack years adult smoking         HC      HC Pack years adult smoking                  MR Egger    8 -0.21496659 0.11154721 0.10224949
#2 Pack years adult smoking         HC      HC Pack years adult smoking           Weighted median    8 -0.10981388 0.04902805 0.02510279
#3 Pack years adult smoking         HC      HC Pack years adult smoking Inverse variance weighted    8 -0.05565872 0.05871784 0.34317972
#4 Pack years adult smoking         HC      HC Pack years adult smoking               Simple mode    8 -0.12269632 0.12017356 0.34124092
#5 Pack years adult smoking         HC      HC Pack years adult smoking             Weighted mode    8 -0.11447199 0.04925477 0.05307462

#####################################################
#SECTION I: ANALYSIS ON PLEIOTROPY AND HETEROGENEITY#
#####################################################

library(meta)

dat_1_post$mr <- dat_1_post$beta.outcome/dat_1_post$beta.exposure
dat_1_post$mr_se <- ((dat_1_post$mr*((dat_1_post$se.exposure/dat_1_post$beta.exposure)^2+(dat_1_post$se.outcome/dat_1_post$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_post$mr, dat_1_post$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0.0234 [0.0000; 0.1983]; tau = 0.1529 [0.0000; 0.4453]
#I^2 = 48.8% [0.0%; 77.2%]; H = 1.40 [1.00; 2.09]

#Test of heterogeneity:
#  Q d.f. p-value
#13.68    7  0.0572

#No heterogeneity. 
#Let's check the rucker framework:

mr_rucker(dat_1_post)

#[[1]]$intercept
#Method   Estimate          SE        CI_low     CI_upp          P
#1  Egger fixed effects 0.00780705 0.003716988  0.0005218877 0.01509221 0.03569630
#2 Egger random effects 0.00780705 0.004813205 -0.0016266595 0.01724076 0.05240099

#[[1]]$Q
#Method         Q df          P
#1   Q_ivw 14.472477  7 0.04338839
#2 Q_egger 10.060922  6 0.12210978
#3  Q_diff  4.411554  1 0.03569630

#[[1]]$res
#[1] "C"

#[[1]]$selected
#Method nsnp   Estimate         SE    CI_low      CI_upp          P
#3 Rucker    8 -0.2149666 0.08614209 -0.383802 -0.04613119 0.04681315

#######################################
#SECTION J: GENERATE SENSITIVITY PLOTS#
#######################################

dat_1_post$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/PackYears/PackYears_HIP/packyears_smk_hip_post_13092021.tiff", units="in", width=10, height=10, res=300)
mr_plots(dat_1_post)
dev.off()
