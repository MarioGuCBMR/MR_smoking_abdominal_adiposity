##############
#INTRODUCTION#
##############

#This is the 2SMR code to run packyears_smk  > whr.

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

match_ss_data <- function(query_chr_pos, query_alleles_1, query_alleles_2, bim_data){
  #This function takes a chr_pos and its alleles and checks whether they are present in bim data.
  #It returns the index of the bim data where the allele that you are supposed to use is located.
  
  #STEP 1: getting the data to test (it will be commented later)
  
  #query_chr_pos = "chr7:130465054"
  #query_alleles_1 = "T_G"
  #query_alleles_2 = "G_T"
  #bim_data = dictionary
  
  #STEP 2: Let's create the ref_alleles in bim_data:
  
  #bim_data$ref_alleles <- paste(bim_data$V6, "_", bim_data$V5, sep = "")
  
  #STEP 3: Now let's perform the query:
  
  bim_data_match <- bim_data[which(bim_data$chr_pos%in%query_chr_pos & bim_data$ref_alleles == query_alleles_1 | bim_data$chr_pos%in%query_chr_pos & bim_data$ref_alleles == query_alleles_2),]
  
  #Finally, let's get the index:  
  
  if(rlang::is_empty(bim_data_match) == TRUE){
    
    next()
    
  }
  
  index_end <- which(bim_data$chr_pos == bim_data_match$chr_pos & bim_data$ref_alleles == bim_data_match$ref_alleles)
  
  #And we return the index: 
  
  return(index_end)
  
}

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

whr <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WHR/WHR_Combined_Curated.txt")

##########################################################################################
#1. we are gonna try to match the data with RSID and with chr_pos to get as many variants#
##########################################################################################

whr_match_by_rsid <- whr[which(whr$MarkerName%in%packyears_smk_gw_ind$SNP),] #7/13

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

whr_match_by_chr_pos <- whr[which(whr$chr_pos_37%in%packyears_smk_gw_ind$chr_pos_37),] #7/7

######################################################
#SECTION C: get proxies for the SNPs that are missing#
######################################################

#First we are gonna get the SNPs that are missing.

packyears_smk_gw_ind_missing <- packyears_smk_gw_ind[which(!(packyears_smk_gw_ind$chr_pos_37%in%whr_match_by_chr_pos$chr_pos_37)),]

#We already have the proxies for these SNPs. 
#Let's check whether they are the same missing SNPs as with WHRadjBMI, just in case.
#To do so, let's read the WHRadjBMI proxies and look for the query SNPs. 

proxies_df <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/PackYears/PackYears_WHRAdjBMI/proxies/combined_query_snp_list.txt")

#Let's get the query SNPs:

query_snps <- proxies_df[which(proxies_df$query_snp == proxies_df$Coord & proxies_df$Distance == 0 & proxies_df$R2 == 1),]

#Let's check if they are the same:

length(which(packyears_smk_gw_ind_missing$chr_pos_37%in%proxies_df$query_snp)) #6/6: all of them. This is perfect

#LET'S GET THOSE PROXIES.

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

packyears_smk_proxies <- packyears_smk[which(packyears_smk$chr_pos_37%in%proxies_df_clean_curated$Coord),] #1893/2308
whr_proxies <- whr[which(whr$chr_pos_37%in%proxies_df_clean_curated$Coord),] #703/2308

#Now let's match the data:

packyears_smk_proxies_match <- packyears_smk_proxies[which(packyears_smk_proxies$chr_pos_37%in%whr_proxies$chr_pos_37),] #699/703
whr_proxies_match <- whr_proxies[which(whr_proxies$chr_pos_37%in%packyears_smk_proxies$chr_pos_37),] #699/703

#whr is a bottleneck and this is the best we can do right now!! However, there are 4 SNPs that are not found in Smk.

packyears_smk_proxies_match_2 <- packyears_smk_proxies[which(packyears_smk_proxies$chr_pos_37%in%whr_proxies$chr_pos_37),] #699/699

#It is the same! The mismatch is NOT due to RSIDs, but just because the variant is no there.

##################
#CHECKING ALLELES#
##################

proxies_check <- proxies_df_clean_curated[which(proxies_df_clean_curated$Coord%in%packyears_smk_proxies_match_2$chr_pos_37),]

proxies_check <- proxies_check[order(proxies_check$Alleles),]

head(proxies_check)
tail(proxies_check)

check <- packyears_smk_proxies_match_2

check <- check[order(check$chr_pos_37),]
proxies_check <- proxies_check[order(proxies_check$Coord),]

head(check$SNP)
head(proxies_check$RS_Number)

head(proxies_check$Alleles)

head(check$effect_allele.exposure)
head(check$other_allele.exposure)

tail(proxies_check$Alleles)

tail(check$effect_allele.exposure)
tail(check$other_allele.exposure)

#It seems that all is good.
#We just need the final comprobation. These are all bi-allelic. So they should be fine!!

#We need to do some checks on the data:

check$alleles_1 <- paste("(", check$effect_allele.exposure, "/", check$other_allele.exposure, ")", sep = "")
check$alleles_2 <- paste("(", check$other_allele.exposure, "/", check$effect_allele.exposure, ")", sep = "")

proxies_check$chr_pos <- proxies_check$Coord
proxies_check$ref_alleles <- proxies_check$Alleles

#And to the proxies_df

index_vect <- c()

for(i in seq(1,length(check$chr_pos))){
  
  print(i)
  
  index_tmp <- as.numeric(as.character(unlist(sapply(check$chr_pos[i], match_ss_data, query_alleles_1 = check$alleles_1[i], query_alleles_2 = check$alleles_2[i], bim_data = proxies_check))))
  
  index_vect <- c(index_vect, index_tmp)
  
}

length(index_vect) #ALL OF THEM!!

#That means that they all match!! No worries about this :)
#We do not need to do this for the rest of packyears.

###################################################################
#SECTION D: combined pre and post proxy data and get the lead SNPs#
###################################################################

packyears_smk_proxies_match$pval <- packyears_smk_proxies_match$pval.exposure

#Now, be wary..., we need to get those that are not missing to avoid removing the proxies we just got.

packyears_smk_gw_ind_not_missing <- packyears_smk_gw_ind[which(packyears_smk_gw_ind$chr_pos_37%in%whr_match_by_chr_pos$chr_pos_37),]

packyears_smk_gw_ind_plus_proxies <- rbind(packyears_smk_gw_ind_not_missing, packyears_smk_proxies_match) #we get 7 + 38 = 45 SNPs. Makes total sense.

#When we run the ld_clump we will get only the lead SNPs in each loci!
#Which is exactly what we want. 

packyears_smk_gw_ind_plus_proxies$rsid <- packyears_smk_gw_ind_plus_proxies$SNP
packyears_smk_gw_ind_plus_proxies$pval <- packyears_smk_gw_ind_plus_proxies$pval.exposure

packyears_smk_gw_ind_post_proxies <- ieugwasr::ld_clump_local(packyears_smk_gw_ind_plus_proxies, bfile = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/2SMR/1000G/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/R-4.1.2/library/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) 

#We end up with 12/13.

##############################################
#SECTION E: Merging exposure and outcome data#
##############################################

#We need to format the data for the outcome:

whr_outcome <- whr

colnames(whr_outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", 
                           "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                           "pval.outcome", "samplesize.outcome", "chr_pos_37", "chr_pos_18")


#Getting ids, so troublesome...

packyears_smk_gw_ind_post_proxies$id.exposure <- "pack-years"
whr_outcome$id.outcome <- "WHR"
packyears_smk_gw_ind_post_proxies$exposure <- "pack-years"
whr_outcome$outcome <- "WHR"

dat_1_pre <- harmonise_data(packyears_smk_gw_ind_post_proxies, whr_outcome, action = 3) #12/12: perfect.

########################################
#SECTION F: filtering harmonise results#
########################################

dat_1_filt <- dat_1_pre[which(dat_1_pre$remove == FALSE),] #12/12
dat_1_filt <- dat_1_filt[which(dat_1_filt$palindromic == FALSE & dat_1_filt$ambiguous == FALSE | dat_1_filt$palindromic == TRUE & dat_1_filt$ambiguous == FALSE),] #10/12

#So we are gonna use the steiger filtering like that:

dat_1_filt$samplesize.exposure <- 142387
dat_1_filt$units.exposure <- "SD"
dat_1_filt$units.outcome <- "SD"

dat_1_filt$mr_keep = TRUE #this is to avoid issues with the palindromes afterwards.

dat_1_filt <- steiger_filtering(dat_1_filt) 
dat_1_filt_1 <- dat_1_filt[which(dat_1_filt$steiger_dir == TRUE),] #83/100

dat_1_filt <- dat_1_filt_1

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

#id.exposure id.outcome outcome                 exposure                    method nsnp            b         se      pval
#1 Pack years adult smoking        WHR     WHR Pack years adult smoking                  MR Egger   10 -0.155674042 0.09939761 0.1559411
#2 Pack years adult smoking        WHR     WHR Pack years adult smoking           Weighted median   10  0.001253659 0.04453187 0.9775410
#3 Pack years adult smoking        WHR     WHR Pack years adult smoking Inverse variance weighted   10  0.050934778 0.06198703 0.4112473
#4 Pack years adult smoking        WHR     WHR Pack years adult smoking               Simple mode   10  0.151024861 0.11516154 0.2221844
#5 Pack years adult smoking        WHR     WHR Pack years adult smoking             Weighted mode   10 -0.014923483 0.04363463 0.7401941

#####################################################
#SECTION I: ANALYSIS ON PLEIOTROPY AND HETEROGENEITY#
#####################################################

library(meta)

dat_1_filt$mr <- dat_1_filt$beta.outcome/dat_1_filt$beta.exposure
dat_1_filt$mr_se <- ((dat_1_filt$mr*((dat_1_filt$se.exposure/dat_1_filt$beta.exposure)^2+(dat_1_filt$se.outcome/dat_1_filt$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_filt$mr, dat_1_filt$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0.0300 [0.0014; 0.2254]; tau = 0.1731 [0.0379; 0.4748]
#I^2 = 56.2% [11.3%; 78.4%]; H = 1.51 [1.06; 2.15]

#Test of heterogeneity:
#  Q d.f. p-value
#20.57    9  0.0147

#No heterogeneity. 
#Let's check the rucker framework:

mr_rucker(dat_1_filt)

#[[1]]$intercept
#Method    Estimate          SE      CI_low     CI_upp           P
#1  Egger fixed effects 0.009448428 0.002884976 0.003793979 0.01510288 0.001056449
#2 Egger random effects 0.009448428 0.003926536 0.001752558 0.01714430 0.008057490

#[[1]]$Q
#Method        Q df           P
#1   Q_ivw 25.54513  9 0.002423791
#2 Q_egger 14.81920  8 0.062757774
#3  Q_diff 10.72592  1 0.001056449

#[[1]]$res
#[1] "C"

#[[1]]$selected
#Method nsnp  Estimate         SE     CI_low      CI_upp          P
#3 Rucker   10 -0.155674 0.07303121 -0.2988126 -0.01253549 0.06562828

#######################################
#SECTION J: GENERATE SENSITIVITY PLOTS#
#######################################

dat_1_filt$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/PackYears/PackYears_WHR/packyears_smk_whr_original_20082021.tiff", units="in", width=10, height=10, res=300)
mr_plots(dat_1_filt)
dev.off()

##############################
#SECTION K: Removing outliers#
##############################

library(RadialMR)

radial_input <- format_radial(dat_1_filt$beta.exposure, dat_1_filt$beta.outcome, dat_1_filt$se.exposure, dat_1_filt$se.outcome, dat_1_filt$SNP)

radial_output <- egger_radial(radial_input, alpha = 0.05, weights = 3) 

outliers <- radial_output$outliers$SNP

#rs1364063: #arm fat percentage. It is really an outlier.

dat_1_post <- dat_1_filt[which(!(dat_1_filt$SNP%in%outliers)),]

################################################
#SECTION G.2: Checking validity of the variants#
################################################

#Are the instruments valid??
#Let's calculate mF and Isq to check it.

F_ = ((dat_1_post$beta.exposure)^2)/((dat_1_post$se.exposure)^2)
mF  = mean(F_)

print(mF)
#83.05064

Isq(dat_1_post$beta.exposure, dat_1_post$se.exposure)
#0.9892917

#########################################
#SECTION H.2: 2nd round of 2SMR analysis#
#########################################

mr(dat_1_post)

#id.exposure id.outcome outcome                 exposure                    method nsnp            b         se      pval
#1 Pack years adult smoking        WHR     WHR Pack years adult smoking                  MR Egger    9 -0.103215016 0.08343033 0.2559241
#2 Pack years adult smoking        WHR     WHR Pack years adult smoking           Weighted median    9 -0.002438654 0.04303646 0.9548122
#3 Pack years adult smoking        WHR     WHR Pack years adult smoking Inverse variance weighted    9  0.027976349 0.04645668 0.5470391
#4 Pack years adult smoking        WHR     WHR Pack years adult smoking               Simple mode    9  0.182808196 0.12792676 0.1908633
#5 Pack years adult smoking        WHR     WHR Pack years adult smoking             Weighted mode    9 -0.021594982 0.04234643 0.6238475

#After removing the outliers we get better results.
#I am sure that it is due, mostly, to those variants associated with smoking.

#####################################################
#SECTION I: ANALYSIS ON PLEIOTROPY AND HETEROGENEITY#
#####################################################

library(meta)

dat_1_post$mr <- dat_1_post$beta.outcome/dat_1_post$beta.exposure
dat_1_post$mr_se <- ((dat_1_post$mr*((dat_1_post$se.exposure/dat_1_post$beta.exposure)^2+(dat_1_post$se.outcome/dat_1_post$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_post$mr, dat_1_post$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0.0106 [0.0000; 0.1124]; tau = 0.1029 [0.0000; 0.3353]
#I^2 = 30.8% [0.0%; 68.0%]; H = 1.20 [1.00; 1.77]

#Test of heterogeneity:
#  Q d.f. p-value
#11.56    8  0.1717

#No heterogeneity. 
#Let's check the rucker framework:

mr_rucker(dat_1_post)

#[[1]]$intercept
#Method    Estimate          SE        CI_low     CI_upp          P
#1  Egger fixed effects 0.006251678 0.003150833  7.615878e-05 0.01242720 0.04724079
#2 Egger random effects 0.006251678 0.003461976 -5.336705e-04 0.01303703 0.03547386

#[[1]]$Q
#Method         Q df          P
#1   Q_ivw 12.387545  8 0.13473208
#2 Q_egger  8.450753  7 0.29452815
#3  Q_diff  3.936792  1 0.04724079

#[[1]]$res
#[1] "A"

#[[1]]$selected
#Method nsnp   Estimate         SE      CI_low   CI_upp         P
#1 Rucker    9 0.02797635 0.03733366 -0.04519628 0.101149 0.4536403

#######################################
#SECTION J: GENERATE SENSITIVITY PLOTS#
#######################################

dat_1_post$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/2SMR/PackYears/PackYears_WHR/packyears_smk_whr_post_20082021.tiff", units="in", width=10, height=10, res=300)
mr_plots(dat_1_post)
dev.off()

###############################
#IMPORTANT TO INVESTIGATE THEM#
###############################


