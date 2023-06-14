##############
#INTRODUCTION#
##############

#This is a small test to check how to deal with MVMR
#results despite having the issues of sample overlap
#between the risk factors.

#In this test we are gonna make it simple and we are 
#just gonna use edu as an additional risk factor.

#Well, then let's go.

###########
#libraries#
###########

library(TwoSampleMR)
library(tidyverse)
library(data.table)

###################
#Loading functions#
###################

mvmr_Mario <- function(dat, outcome, path_){
  #This functions runs MVMR for different ranges of correlation between phenotypes.
  #It runs three times: 1) without taking into account outliers, 2) taking into account only extreme
  #outliers and 3) taking into account all potential outliers. 
  
  #The steps to performed are:
  
  #1) Generate covariance matrices with fake phenotypic correlation matrices and format input data.
  #2) Run conditional mF and conditional Cochran's Q.
  #3) Run MVMR IVW
  #4) Repeat 1-3 after detecting outliers with mvmr_cochran_Q
  
  #########################################################
  #STEP 1: construct the covariance matrix and format data#
  #########################################################
  
  #We first need a covariance matrix for each type of phenotypic correlation.
  
  Pcov_0 <- matrix(data = c(1, 0, 0 , 1), nrow = 2, ncol = 2) #we make the fake phenotypic correlations
  Pcov_0.05 <- matrix(data = c(1, -0.05, -0.05 , 1), nrow = 2, ncol = 2) #we make the fake phenotypic correlations
  Pcov_0.10 <- matrix(data = c(1, -0.10, -0.10 , 1), nrow = 2, ncol = 2) #we make the fake phenotypic correlations
  Pcov_0.15 <- matrix(data = c(1, -0.15, -0.15 , 1), nrow = 2, ncol = 2) #we make the fake phenotypic correlations
  Pcov_0.20 <- matrix(data = c(1, -0.20, -0.20 , 1), nrow = 2, ncol = 2) #we make the fake phenotypic correlations
  Pcov_0.25 <- matrix(data = c(1, -0.25, -0.25 , 1), nrow = 2, ncol = 2) #we make the fake phenotypic correlations
  
  #Let's change the names of the matrices:
  
  colnames(Pcov_0) <- c("1", "2")
  rownames(Pcov_0) <- c("1", "2")
  
  colnames(Pcov_0.05) <- c("1", "2")
  rownames(Pcov_0.05) <- c("1", "2")
  
  colnames(Pcov_0.10) <- c("1", "2")
  rownames(Pcov_0.10) <- c("1", "2")
  
  colnames(Pcov_0.15) <- c("1", "2")
  rownames(Pcov_0.15) <- c("1", "2")
  
  colnames(Pcov_0.20) <- c("1", "2")
  rownames(Pcov_0.20) <- c("1", "2")
  
  colnames(Pcov_0.25) <- c("1", "2")
  rownames(Pcov_0.25) <- c("1", "2")  
  
  #Now we need a dataframe with the SEs.
  #It needs to have the double of rows because MVMR is only made for <1 SNP.
  
  seBXGs = cbind(dat$se.smk, dat$se.edu)
  
  colnames(seBXGs) <- c("1", "2")
  
  #And now we run the covariance matrix:
  
  covariance_matrix_0 <- phenocov_mvmr(Pcov = Pcov_0, seBXGs = seBXGs)
  covariance_matrix_0.05 <- phenocov_mvmr(Pcov = Pcov_0.05, seBXGs = seBXGs)
  covariance_matrix_0.10 <- phenocov_mvmr(Pcov = Pcov_0.10, seBXGs = seBXGs)
  covariance_matrix_0.15 <- phenocov_mvmr(Pcov = Pcov_0.15, seBXGs = seBXGs)
  covariance_matrix_0.20 <- phenocov_mvmr(Pcov = Pcov_0.20, seBXGs = seBXGs)
  covariance_matrix_0.25 <- phenocov_mvmr(Pcov = Pcov_0.25, seBXGs = seBXGs)
  
  #Maybe there is something wrong with my matrix. 
  #I can't be 100% sure SINCE THERE ARE NO EXAMPLES ANYWHERE
  
  ##########################################
  #Let's format the summary statistics data#
  ##########################################
  
  BXGs <- cbind(dat$beta.smk, dat$beta.edu)
  
  colnames(BXGs) <- c("1", "2")
  
  BYG <-outcome$final_beta
  
  seBYG <- outcome$se
  
  final_df <- format_mvmr(BXGs = BXGs, seBXGs = seBXGs, BYG = BYG, seBYG = seBYG, RSID = dat$SNP)
  
  ############################################
  #STEP 2: Run conditional mF and Cochran's Q#
  ############################################
  
  IV_0 <- strength_mvmr(final_df, covariance_matrix_0) 
  IV_0.05 <- strength_mvmr(final_df, covariance_matrix_0.05)
  IV_0.10 <- strength_mvmr(final_df, covariance_matrix_0.10)
  IV_0.15 <- strength_mvmr(final_df, covariance_matrix_0.15)
  IV_0.20 <- strength_mvmr(final_df, covariance_matrix_0.20)
  IV_0.25 <- strength_mvmr(final_df, covariance_matrix_0.25)
  
  #Now let's get the df with the info:
  
  IV_df <- rbind(IV_0, IV_0.05, IV_0.10, IV_0.15, IV_0.20, IV_0.25)
  
  IV_df$SNP <- c(length(covariance_matrix_0), length(covariance_matrix_0.05), length(covariance_matrix_0.10), length(covariance_matrix_0.15), length(covariance_matrix_0.20), length(covariance_matrix_0.25))
  
  rownames(IV_df) <- c("Correlation_0", "Correlation_0.05", "Correlation_0.10", "Correlation_0.15", "Correlation_0.20", "Correlation_0.25")
  colnames(IV_df) <- c("Smoking_mF", "Edu_mF", "Number of SNPs")
  
  #Now we do the same for pleiotropy:
  
  pleiotropy_0 <- pleiotropy_mvmr(final_df, covariance_matrix_0) 
  pleiotropy_0.05 <- pleiotropy_mvmr(final_df, covariance_matrix_0.05) 
  pleiotropy_0.10 <- pleiotropy_mvmr(final_df, covariance_matrix_0.10) 
  pleiotropy_0.15 <- pleiotropy_mvmr(final_df, covariance_matrix_0.15) 
  pleiotropy_0.20 <- pleiotropy_mvmr(final_df, covariance_matrix_0.20)
  pleiotropy_0.25 <- pleiotropy_mvmr(final_df, covariance_matrix_0.25)
  
  pleiotropy_df <- rbind(pleiotropy_0, pleiotropy_0.05, pleiotropy_0.10, pleiotropy_0.15, pleiotropy_0.20, pleiotropy_0.25)
  
  rownames(pleiotropy_df) <- c("Correlation_0", "Correlation_0.05", "Correlation_0.10", "Correlation_0.15", "Correlation_0.20", "Correlation_0.25")
  colnames(pleiotropy_df) <- c("Qstat_original", "Qpval_original")
  
  ######################
  #STEP 3: Run IVW MVMR#
  ######################
  
  results_ivw_original_0 <- ivw_mvmr(final_df,gencov = covariance_matrix_0)
  results_ivw_original_0.05 <- ivw_mvmr(final_df,gencov = covariance_matrix_0.05)
  results_ivw_original_0.10 <- ivw_mvmr(final_df,gencov = covariance_matrix_0.10)
  results_ivw_original_0.15 <- ivw_mvmr(final_df,gencov = covariance_matrix_0.15)
  results_ivw_original_0.20 <- ivw_mvmr(final_df,gencov = covariance_matrix_0.20)
  results_ivw_original_0.25 <- ivw_mvmr(final_df,gencov = covariance_matrix_0.25)
  
  results_ivw_original_all <- rbind(results_ivw_original_0, results_ivw_original_0.05, results_ivw_original_0.10 , results_ivw_original_0.15, results_ivw_original_0.20, results_ivw_original_0.25)
  
  colnames(results_ivw_original_all) <-  c("Estimate_original", "Std. Error_original", "t value_original", "Pr(>|t|)_original")
  rownames(results_ivw_original_all) <-  c("SMK_0", "EDU_0", "SMK_0.05", "EDU_0.05", "SMK_0.10", "EDU_0.10", "SMK_0.15", "EDU_0.15", "SMK_20", "EDU_20", "SMK_25", "EDU_25")
  
  #########################
  #STEP 4: detect outliers#
  #########################
  
  outliers_0 <- as.data.frame(unlist(mvmr_cochran_Q(final_df, covariance_matrix_0)))
  outliers_0.05 <- as.data.frame(unlist(mvmr_cochran_Q(final_df, covariance_matrix_0.05)))
  outliers_0.10 <- as.data.frame(unlist(mvmr_cochran_Q(final_df, covariance_matrix_0.10)))
  outliers_0.15 <- as.data.frame(unlist(mvmr_cochran_Q(final_df, covariance_matrix_0.15)))
  outliers_0.20 <- as.data.frame(unlist(mvmr_cochran_Q(final_df, covariance_matrix_0.20)))
  outliers_0.25 <- as.data.frame(unlist(mvmr_cochran_Q(final_df, covariance_matrix_0.25)))
  
  outliers_0$SNP <- final_df$SNP
  outliers_0.05$SNP <- final_df$SNP
  outliers_0.10$SNP <- final_df$SNP
  outliers_0.15$SNP <- final_df$SNP
  outliers_0.20$SNP <- final_df$SNP
  outliers_0.25$SNP <- final_df$SNP
  
  colnames(outliers_0) <- c("Qstat", "SNP")
  colnames(outliers_0.05) <- c("Qstat", "SNP")
  colnames(outliers_0.10) <- c("Qstat", "SNP")
  colnames(outliers_0.15) <- c("Qstat", "SNP")
  colnames(outliers_0.20) <- c("Qstat", "SNP")
  colnames(outliers_0.25) <- c("Qstat", "SNP")
  
  #Now we are ready to redo the analysis taking into account different types of outliers!!
  
  ########################################################################
  # REDOING ANALYSIS WITHOUT EXTREME OUTLIERS THAT HAVE COCHRAN'S Q > 10 #
  ########################################################################
  
  ##################################################################
  #STEP 1: get the outliers for each type of phenotypic correlation#
  ##################################################################
  
  outlier_Q_10_0 <- outliers_0[which(outliers_0$Qstat > 10),] 
  outlier_Q_10_0.05 <- outliers_0.05[which(outliers_0.05$Qstat > 10),] 
  outlier_Q_10_0.10 <- outliers_0.10[which(outliers_0.10$Qstat > 10),] 
  outlier_Q_10_0.15 <- outliers_0.15[which(outliers_0.15$Qstat > 10),] 
  outlier_Q_10_0.20 <- outliers_0.20[which(outliers_0.20$Qstat > 10),] 
  outlier_Q_10_0.25 <- outliers_0.25[which(outliers_0.25$Qstat > 10),] 
  
  #########################################################
  #STEP 2: get the covariance matrices and format the data#
  #########################################################
  
  #First we need to remove the outlier from the exposure and the outcome:
  
  dat_q_10_0 <- dat[-which(dat$SNP%in%outlier_Q_10_0$SNP),]
  dat_q_10_0.05 <- dat[-which(dat$SNP%in%outlier_Q_10_0.05$SNP),]
  dat_q_10_0.10 <- dat[-which(dat$SNP%in%outlier_Q_10_0.10$SNP),]
  dat_q_10_0.15 <- dat[-which(dat$SNP%in%outlier_Q_10_0.15$SNP),]
  dat_q_10_0.20 <- dat[-which(dat$SNP%in%outlier_Q_10_0.20$SNP),]
  dat_q_10_0.25 <- dat[-which(dat$SNP%in%outlier_Q_10_0.25$SNP),]
  
  outcome_q_10_0 <- outcome[-which(outcome$MarkerName%in%outlier_Q_10_0$SNP),]
  outcome_q_10_0.05 <- outcome[-which(outcome$MarkerName%in%outlier_Q_10_0.05$SNP),]
  outcome_q_10_0.10 <- outcome[-which(outcome$MarkerName%in%outlier_Q_10_0.10$SNP),]
  outcome_q_10_0.15 <- outcome[-which(outcome$MarkerName%in%outlier_Q_10_0.15$SNP),]
  outcome_q_10_0.20 <- outcome[-which(outcome$MarkerName%in%outlier_Q_10_0.20$SNP),]
  outcome_q_10_0.25 <- outcome[-which(outcome$MarkerName%in%outlier_Q_10_0.25$SNP),]
  
  #Now calculate the covariance matrices again:
  
  seBXGs_0 = cbind(dat_q_10_0$se.smk, dat_q_10_0$se.edu)
  seBXGs_0.05 = cbind(dat_q_10_0.05$se.smk, dat_q_10_0.05$se.edu)
  seBXGs_0.10 = cbind(dat_q_10_0.10$se.smk, dat_q_10_0.10$se.edu)
  seBXGs_0.15 = cbind(dat_q_10_0.15$se.smk, dat_q_10_0.15$se.edu)
  seBXGs_0.20 = cbind(dat_q_10_0.20$se.smk, dat_q_10_0.20$se.edu)
  seBXGs_0.25 = cbind(dat_q_10_0.25$se.smk, dat_q_10_0.25$se.edu)
  
  colnames(seBXGs_0) <- c("1", "2")
  colnames(seBXGs_0.05) <- c("1", "2")
  colnames(seBXGs_0.10) <- c("1", "2")
  colnames(seBXGs_0.15) <- c("1", "2")
  colnames(seBXGs_0.20) <- c("1", "2")
  colnames(seBXGs_0.25) <- c("1", "2")
  
  #And now we run the covariance matrix:
  
  covariance_matrix_0_q_10 <- phenocov_mvmr(Pcov = Pcov_0, seBXGs = seBXGs_0)
  covariance_matrix_0.05_q_10 <- phenocov_mvmr(Pcov = Pcov_0.05, seBXGs = seBXGs_0.05)
  covariance_matrix_0.10_q_10 <- phenocov_mvmr(Pcov = Pcov_0.10, seBXGs = seBXGs_0.10)
  covariance_matrix_0.15_q_10 <- phenocov_mvmr(Pcov = Pcov_0.15, seBXGs = seBXGs_0.15)
  covariance_matrix_0.20_q_10 <- phenocov_mvmr(Pcov = Pcov_0.20, seBXGs = seBXGs_0.20)
  covariance_matrix_0.25_q_10 <- phenocov_mvmr(Pcov = Pcov_0.25, seBXGs = seBXGs_0.25)
  
  #And finally we format the data:
  
  #For 0:
  
  BXGs_0 <- cbind(dat_q_10_0$beta.smk, dat_q_10_0$beta.edu)
  colnames(BXGs_0) <- c("1", "2")
  BYG_0 <-outcome_q_10_0$final_beta
  seBYG_0 <- outcome_q_10_0$se
  final_df_0 <- format_mvmr(BXGs = BXGs_0, seBXGs = seBXGs_0, BYG = BYG_0, seBYG = seBYG_0, RSID = dat_q_10_0$SNP)
  
  #For 0.05:
  
  BXGs_0.05 <- cbind(dat_q_10_0.05$beta.smk, dat_q_10_0.05$beta.edu)
  colnames(BXGs_0.05) <- c("1", "2")
  BYG_0.05 <-outcome_q_10_0.05$final_beta
  seBYG_0.05 <- outcome_q_10_0.05$se
  final_df_0.05 <- format_mvmr(BXGs = BXGs_0.05, seBXGs = seBXGs_0.05, BYG = BYG_0.05, seBYG = seBYG_0.05, RSID = dat_q_10_0.05$SNP)
  
  #For 0.10:
  
  BXGs_0.10 <- cbind(dat_q_10_0.10$beta.smk, dat_q_10_0.10$beta.edu)
  colnames(BXGs_0.10) <- c("1", "2")
  BYG_0.10 <-outcome_q_10_0.10$final_beta
  seBYG_0.10 <- outcome_q_10_0.10$se
  final_df_0.10 <- format_mvmr(BXGs = BXGs_0.10, seBXGs = seBXGs_0.10, BYG = BYG_0.10, seBYG = seBYG_0.10, RSID = dat_q_10_0.10$SNP)
  
  #For 0.15:
  
  BXGs_0.15 <- cbind(dat_q_10_0.15$beta.smk, dat_q_10_0.15$beta.edu)
  colnames(BXGs_0.15) <- c("1", "2")
  BYG_0.15 <-outcome_q_10_0.15$final_beta
  seBYG_0.15 <- outcome_q_10_0.15$se
  final_df_0.15 <- format_mvmr(BXGs = BXGs_0.15, seBXGs = seBXGs_0.15, BYG = BYG_0.15, seBYG = seBYG_0.15, RSID = dat_q_10_0.15$SNP)
  
  #For 0.20:
  
  BXGs_0.20 <- cbind(dat_q_10_0.20$beta.smk, dat_q_10_0.20$beta.edu)
  colnames(BXGs_0.20) <- c("1", "2")
  BYG_0.20 <-outcome_q_10_0.20$final_beta
  seBYG_0.20 <- outcome_q_10_0.20$se
  final_df_0.20 <- format_mvmr(BXGs = BXGs_0.20, seBXGs = seBXGs_0.20, BYG = BYG_0.20, seBYG = seBYG_0.20, RSID = dat_q_10_0.20$SNP)
  
  #For 0.25:
  
  BXGs_0.25 <- cbind(dat_q_10_0.25$beta.smk, dat_q_10_0.25$beta.edu)
  colnames(BXGs_0.25) <- c("1", "2")
  BYG_0.25 <-outcome_q_10_0.25$final_beta
  seBYG_0.25 <- outcome_q_10_0.25$se
  final_df_0.25 <- format_mvmr(BXGs = BXGs_0.25, seBXGs = seBXGs_0.25, BYG = BYG_0.25, seBYG = seBYG_0.25, RSID = dat_q_10_0.25$SNP)
  
  ##################################################################
  #Calculating conditional mF and Cochran's Q after outlier removal#
  ##################################################################
  
  IV_0_q_10 <- strength_mvmr(final_df_0, covariance_matrix_0_q_10) 
  IV_0.05_q_10 <- strength_mvmr(final_df_0.05, covariance_matrix_0.05_q_10)
  IV_0.10_q_10 <- strength_mvmr(final_df_0.10, covariance_matrix_0.10_q_10)
  IV_0.15_q_10 <- strength_mvmr(final_df_0.15, covariance_matrix_0.15_q_10)
  IV_0.20_q_10 <- strength_mvmr(final_df_0.20, covariance_matrix_0.20_q_10)
  IV_0.25_q_10 <- strength_mvmr(final_df_0.25, covariance_matrix_0.25_q_10)
  
  #Now let's get the df with the info:
  
  IV_df_q_10 <- rbind(IV_0_q_10, IV_0.05_q_10, IV_0.10_q_10, IV_0.15_q_10, IV_0.20_q_10, IV_0.25_q_10)
  
  IV_df_q_10$SNP <- c(length(covariance_matrix_0_q_10), length(covariance_matrix_0.05_q_10), length(covariance_matrix_0.10_q_10), length(covariance_matrix_0.15_q_10), length(covariance_matrix_0.20_q_10), length(covariance_matrix_0.25_q_10))
  
  rownames(IV_df_q_10) <- c("Correlation_0", "Correlation_0.05", "Correlation_0.10", "Correlation_0.15", "Correlation_0.20", "Correlation_0.25")
  colnames(IV_df_q_10) <- c("Smoking_mF_lax_outlier_out", "Edu_mF_lax_outliers_out", "Number_of_SNPs_lax_outliers_out")
  
  #Now we do the same for pleiotropy:
  
  pleiotropy_0_q_10 <- pleiotropy_mvmr(final_df_0, covariance_matrix_0_q_10) 
  pleiotropy_0.05_q_10 <- pleiotropy_mvmr(final_df_0.05, covariance_matrix_0.05_q_10) 
  pleiotropy_0.10_q_10 <- pleiotropy_mvmr(final_df_0.10, covariance_matrix_0.10_q_10) 
  pleiotropy_0.15_q_10 <- pleiotropy_mvmr(final_df_0.15, covariance_matrix_0.15_q_10) 
  pleiotropy_0.20_q_10 <- pleiotropy_mvmr(final_df_0.20, covariance_matrix_0.20_q_10)
  pleiotropy_0.25_q_10 <- pleiotropy_mvmr(final_df_0.25, covariance_matrix_0.25_q_10)
  
  pleiotropy_df_q_10 <- rbind(pleiotropy_0_q_10, pleiotropy_0.05_q_10, pleiotropy_0.10_q_10, pleiotropy_0.15_q_10, pleiotropy_0.20_q_10, pleiotropy_0.25_q_10)
  
  rownames(pleiotropy_df_q_10) <- c("Correlation_0", "Correlation_0.05", "Correlation_0.10", "Correlation_0.15", "Correlation_0.20", "Correlation_0.25")
  colnames(pleiotropy_df_q_10) <- c("Qstat_lax_outlier_out", "Qpval_lax_outlier_out")
  
  ##########################
  #Final step: run IVW MVMR#
  ##########################
  
  results_ivw_q_10_0 <- ivw_mvmr(final_df_0, gencov = covariance_matrix_0_q_10)
  results_ivw_q_10_0.05 <- ivw_mvmr(final_df_0.05, gencov = covariance_matrix_0.05_q_10)
  results_ivw_q_10_0.10 <- ivw_mvmr(final_df_0.10, gencov = covariance_matrix_0.10_q_10)
  results_ivw_q_10_0.15 <- ivw_mvmr(final_df_0.15, gencov = covariance_matrix_0.15_q_10)
  results_ivw_q_10_0.20 <- ivw_mvmr(final_df_0.20, gencov = covariance_matrix_0.20_q_10)
  results_ivw_q_10_0.25 <- ivw_mvmr(final_df_0.25, gencov = covariance_matrix_0.25_q_10)
  
  results_ivw_q_10_all <- rbind(results_ivw_q_10_0, results_ivw_q_10_0.05, results_ivw_q_10_0.10, results_ivw_q_10_0.15, results_ivw_q_10_0.20, results_ivw_q_10_0.25)
  
  colnames(results_ivw_q_10_all) <-  c("Estimate_lax_outlier_out", "Std. Error_lax_outlier_out", "t value_lax_outlier_out", "Pr(>|t|)_lax_outlier_edu")
  rownames(results_ivw_q_10_all) <-  c("SMK_0", "EDU_0", "SMK_0.05", "EDU_0.05", "SMK_0.10", "EDU_0.10", "SMK_0.15", "EDU_0.15", "SMK_0.20", "EDU_0.20", "SMK_0.25", "EDU_0.25")
  
  ##########################################################################
  # REDOING ANALYSIS WITHOUT EXTREME OUTLIERS THAT HAVE COCHRAN'S Q > 3.84 #
  ##########################################################################
  
  ##################################################################
  #STEP 1: get the outliers for each type of phenotypic correlation#
  ##################################################################
  
  outlier_q_3.84_0 <- outliers_0[which(outliers_0$Qstat > 3.84),] 
  outlier_q_3.84_0.05 <- outliers_0.05[which(outliers_0.05$Qstat > 3.84),] 
  outlier_q_3.84_0.10 <- outliers_0.10[which(outliers_0.10$Qstat > 3.84),] 
  outlier_q_3.84_0.15 <- outliers_0.15[which(outliers_0.15$Qstat > 3.84),] 
  outlier_q_3.84_0.20 <- outliers_0.20[which(outliers_0.20$Qstat > 3.84),] 
  outlier_q_3.84_0.25 <- outliers_0.25[which(outliers_0.25$Qstat > 3.84),] 
  
  #########################################################
  #STEP 2: get the covariance matrices and format the data#
  #########################################################
  
  #First we need to remove the outlier from the exposure and the outcome:
  
  dat_q_3.84_0 <- dat[-which(dat$SNP%in%outlier_q_3.84_0$SNP),]
  dat_q_3.84_0.05 <- dat[-which(dat$SNP%in%outlier_q_3.84_0.05$SNP),]
  dat_q_3.84_0.10 <- dat[-which(dat$SNP%in%outlier_q_3.84_0.10$SNP),]
  dat_q_3.84_0.15 <- dat[-which(dat$SNP%in%outlier_q_3.84_0.15$SNP),]
  dat_q_3.84_0.20 <- dat[-which(dat$SNP%in%outlier_q_3.84_0.20$SNP),]
  dat_q_3.84_0.25 <- dat[-which(dat$SNP%in%outlier_q_3.84_0.25$SNP),]
  
  outcome_q_3.84_0 <- outcome[-which(outcome$MarkerName%in%outlier_q_3.84_0$SNP),]
  outcome_q_3.84_0.05 <- outcome[-which(outcome$MarkerName%in%outlier_q_3.84_0.05$SNP),]
  outcome_q_3.84_0.10 <- outcome[-which(outcome$MarkerName%in%outlier_q_3.84_0.10$SNP),]
  outcome_q_3.84_0.15 <- outcome[-which(outcome$MarkerName%in%outlier_q_3.84_0.15$SNP),]
  outcome_q_3.84_0.20 <- outcome[-which(outcome$MarkerName%in%outlier_q_3.84_0.20$SNP),]
  outcome_q_3.84_0.25 <- outcome[-which(outcome$MarkerName%in%outlier_q_3.84_0.25$SNP),]
  
  #Now calculate the covariance matrices again:
  
  seBXGs_0 = cbind(dat_q_3.84_0$se.smk, dat_q_3.84_0$se.edu)
  seBXGs_0.05 = cbind(dat_q_3.84_0.05$se.smk, dat_q_3.84_0.05$se.edu)
  seBXGs_0.10 = cbind(dat_q_3.84_0.10$se.smk, dat_q_3.84_0.10$se.edu)
  seBXGs_0.15 = cbind(dat_q_3.84_0.15$se.smk, dat_q_3.84_0.15$se.edu)
  seBXGs_0.20 = cbind(dat_q_3.84_0.20$se.smk, dat_q_3.84_0.20$se.edu)
  seBXGs_0.25 = cbind(dat_q_3.84_0.25$se.smk, dat_q_3.84_0.25$se.edu)
  
  colnames(seBXGs_0) <- c("1", "2")
  colnames(seBXGs_0.05) <- c("1", "2")
  colnames(seBXGs_0.10) <- c("1", "2")
  colnames(seBXGs_0.15) <- c("1", "2")
  colnames(seBXGs_0.20) <- c("1", "2")
  colnames(seBXGs_0.25) <- c("1", "2")
  
  #And now we run the covariance matrix:
  
  covariance_matrix_0_q_3.84 <- phenocov_mvmr(Pcov = Pcov_0, seBXGs = seBXGs_0)
  covariance_matrix_0.05_q_3.84 <- phenocov_mvmr(Pcov = Pcov_0.05, seBXGs = seBXGs_0.05)
  covariance_matrix_0.10_q_3.84 <- phenocov_mvmr(Pcov = Pcov_0.10, seBXGs = seBXGs_0.10)
  covariance_matrix_0.15_q_3.84 <- phenocov_mvmr(Pcov = Pcov_0.15, seBXGs = seBXGs_0.15)
  covariance_matrix_0.20_q_3.84 <- phenocov_mvmr(Pcov = Pcov_0.20, seBXGs = seBXGs_0.20)
  covariance_matrix_0.25_q_3.84 <- phenocov_mvmr(Pcov = Pcov_0.25, seBXGs = seBXGs_0.25)
  
  #And finally we format the data:
  
  #For 0:
  
  BXGs_0 <- cbind(dat_q_3.84_0$beta.smk, dat_q_3.84_0$beta.edu)
  colnames(BXGs_0) <- c("1", "2")
  BYG_0 <-outcome_q_3.84_0$final_beta
  seBYG_0 <- outcome_q_3.84_0$se
  final_df_0 <- format_mvmr(BXGs = BXGs_0, seBXGs = seBXGs_0, BYG = BYG_0, seBYG = seBYG_0, RSID = dat_q_3.84_0$SNP)
  
  #For 0.05:
  
  BXGs_0.05 <- cbind(dat_q_3.84_0.05$beta.smk, dat_q_3.84_0.05$beta.edu)
  colnames(BXGs_0.05) <- c("1", "2")
  BYG_0.05 <-outcome_q_3.84_0.05$final_beta
  seBYG_0.05 <- outcome_q_3.84_0.05$se
  final_df_0.05 <- format_mvmr(BXGs = BXGs_0.05, seBXGs = seBXGs_0.05, BYG = BYG_0.05, seBYG = seBYG_0.05, RSID = dat_q_3.84_0.05$SNP)
  
  #For 0.10:
  
  BXGs_0.10 <- cbind(dat_q_3.84_0.10$beta.smk, dat_q_3.84_0.10$beta.edu)
  colnames(BXGs_0.10) <- c("1", "2")
  BYG_0.10 <-outcome_q_3.84_0.10$final_beta
  seBYG_0.10 <- outcome_q_3.84_0.10$se
  final_df_0.10 <- format_mvmr(BXGs = BXGs_0.10, seBXGs = seBXGs_0.10, BYG = BYG_0.10, seBYG = seBYG_0.10, RSID = dat_q_3.84_0.10$SNP)
  
  #For 0.15:
  
  BXGs_0.15 <- cbind(dat_q_3.84_0.15$beta.smk, dat_q_3.84_0.15$beta.edu)
  colnames(BXGs_0.15) <- c("1", "2")
  BYG_0.15 <-outcome_q_3.84_0.15$final_beta
  seBYG_0.15 <- outcome_q_3.84_0.15$se
  final_df_0.15 <- format_mvmr(BXGs = BXGs_0.15, seBXGs = seBXGs_0.15, BYG = BYG_0.15, seBYG = seBYG_0.15, RSID = dat_q_3.84_0.15$SNP)
  
  #For 0.20:
  
  BXGs_0.20 <- cbind(dat_q_3.84_0.20$beta.smk, dat_q_3.84_0.20$beta.edu)
  colnames(BXGs_0.20) <- c("1", "2")
  BYG_0.20 <-outcome_q_3.84_0.20$final_beta
  seBYG_0.20 <- outcome_q_3.84_0.20$se
  final_df_0.20 <- format_mvmr(BXGs = BXGs_0.20, seBXGs = seBXGs_0.20, BYG = BYG_0.20, seBYG = seBYG_0.20, RSID = dat_q_3.84_0.20$SNP)
  
  #For 0.25:
  
  BXGs_0.25 <- cbind(dat_q_3.84_0.25$beta.smk, dat_q_3.84_0.25$beta.edu)
  colnames(BXGs_0.25) <- c("1", "2")
  BYG_0.25 <-outcome_q_3.84_0.25$final_beta
  seBYG_0.25 <- outcome_q_3.84_0.25$se
  final_df_0.25 <- format_mvmr(BXGs = BXGs_0.25, seBXGs = seBXGs_0.25, BYG = BYG_0.25, seBYG = seBYG_0.25, RSID = dat_q_3.84_0.25$SNP)
  
  
  ##################################################################
  #Calculating conditional mF and Cochran's Q after outlier removal#
  ##################################################################
  
  IV_0_q_3.84 <- strength_mvmr(final_df_0, covariance_matrix_0_q_3.84) 
  IV_0.05_q_3.84 <- strength_mvmr(final_df_0.05, covariance_matrix_0.05_q_3.84)
  IV_0.10_q_3.84 <- strength_mvmr(final_df_0.10, covariance_matrix_0.10_q_3.84)
  IV_0.15_q_3.84 <- strength_mvmr(final_df_0.15, covariance_matrix_0.15_q_3.84)
  IV_0.20_q_3.84 <- strength_mvmr(final_df_0.20, covariance_matrix_0.20_q_3.84)
  IV_0.25_q_3.84 <- strength_mvmr(final_df_0.25, covariance_matrix_0.25_q_3.84)
  
  #Now let's get the df with the info:
  
  IV_df_q_3.84 <- rbind(IV_0_q_3.84, IV_0.05_q_3.84, IV_0.10_q_3.84, IV_0.15_q_3.84, IV_0.20_q_3.84, IV_0.25_q_3.84)
  
  IV_df_q_3.84$SNP <- c(length(covariance_matrix_0_q_3.84), length(covariance_matrix_0.05_q_3.84), length(covariance_matrix_0.10_q_3.84), length(covariance_matrix_0.15_q_3.84), length(covariance_matrix_0.20_q_3.84), length(covariance_matrix_0.25_q_3.84))
  
  rownames(IV_df_q_3.84) <- c("Correlation_0", "Correlation_0.05", "Correlation_0.10", "Correlation_0.15", "Correlation_0.20", "Correlation_0.25")
  colnames(IV_df_q_3.84) <- c("Smoking_mF_strict_outlier_out", "Edu_mF_strict_outliers_out", "Number_of_SNPs_strict_outliers_out")
  
  #Now we do the same for pleiotropy:
  
  pleiotropy_0_q_3.84 <- pleiotropy_mvmr(final_df_0, covariance_matrix_0_q_3.84) 
  pleiotropy_0.05_q_3.84 <- pleiotropy_mvmr(final_df_0.05, covariance_matrix_0.05_q_3.84) 
  pleiotropy_0.10_q_3.84 <- pleiotropy_mvmr(final_df_0.10, covariance_matrix_0.10_q_3.84) 
  pleiotropy_0.15_q_3.84 <- pleiotropy_mvmr(final_df_0.15, covariance_matrix_0.15_q_3.84) 
  pleiotropy_0.20_q_3.84 <- pleiotropy_mvmr(final_df_0.20, covariance_matrix_0.20_q_3.84)
  pleiotropy_0.25_q_3.84 <- pleiotropy_mvmr(final_df_0.25, covariance_matrix_0.25_q_3.84)
  
  pleiotropy_df_q_3.84 <- rbind(pleiotropy_0_q_3.84, pleiotropy_0.05_q_3.84, pleiotropy_0.10_q_3.84, pleiotropy_0.15_q_3.84, pleiotropy_0.20_q_3.84, pleiotropy_0.25_q_3.84)
  
  rownames(pleiotropy_df_q_3.84) <- c("Correlation_0", "Correlation_0.05", "Correlation_0.10", "Correlation_0.15", "Correlation_0.20", "Correlation_0.25")
  colnames(pleiotropy_df_q_3.84) <- c("Qstat_strict_outlier_out", "Qpval_strict_outlier_out")
  
  ##########################
  #Final step: run IVW MVMR#
  ##########################
  
  results_ivw_q_3.84_0 <- ivw_mvmr(final_df_0, gencov = covariance_matrix_0_q_3.84)
  results_ivw_q_3.84_0.05 <- ivw_mvmr(final_df_0.05, gencov = covariance_matrix_0.05_q_3.84)
  results_ivw_q_3.84_0.10 <- ivw_mvmr(final_df_0.10, gencov = covariance_matrix_0.10_q_3.84)
  results_ivw_q_3.84_0.15 <- ivw_mvmr(final_df_0.15, gencov = covariance_matrix_0.15_q_3.84)
  results_ivw_q_3.84_0.20 <- ivw_mvmr(final_df_0.20, gencov = covariance_matrix_0.20_q_3.84)
  results_ivw_q_3.84_0.25 <- ivw_mvmr(final_df_0.25, gencov = covariance_matrix_0.25_q_3.84)
  
  results_ivw_q_3.84_all <- rbind(results_ivw_q_3.84_0, results_ivw_q_3.84_0.05, results_ivw_q_3.84_0.10, results_ivw_q_3.84_0.15, results_ivw_q_3.84_0.20, results_ivw_q_3.84_0.25)
  
  colnames(results_ivw_q_3.84_all) <-  c("Estimate_strict_outlier_out", "Std. Error_strict_outlier_out", "t value_strict_outlier_out", "Pr(>|t|)_strict_outlier_edu")
  rownames(results_ivw_q_3.84_all) <-  c("SMK_0", "EDU_0", "SMK_0.05", "EDU_0.05", "SMK_0.10", "EDU_0.10", "SMK_0.15", "EDU_0.15", "SMK_0.20", "EDU_0.20", "SMK_0.25", "EDU_0.25")
  
  #######################################################################################################################
  #Finally let's merge ALL the data into one CSV file that we are gonna save in the path that we specify in the function#
  #######################################################################################################################
  
  #First let's merge the mF data:
  
  IV_pleiotropy_df_all <- cbind(IV_df, pleiotropy_df, IV_df_q_10, pleiotropy_df_q_10, IV_df_q_3.84, pleiotropy_df_q_3.84)
  
  IV_pleiotropy_df_all <- as.data.frame(IV_pleiotropy_df_all)
  
  #And now merge the results data:
  
  ivw_results_all <- cbind(results_ivw_original_all, results_ivw_q_10_all, results_ivw_q_3.84_all)
  
  #Let's use the path to save them. The first one is a bit tricky since it has a weird format.
  
  path_iv_pleiotropy <- paste(path_, "conditional_mF_cochrans_Q", sep = "")
  write_rds(IV_pleiotropy_df_all, path_iv_pleiotropy)
  
  #Now the same with the results, which are easier to do:
  
  path_ivw_results <- paste(path_, "mvmr_ivw_results.csv", sep = "")
  write.csv(ivw_results_all, path_ivw_results, row.names = FALSE, quote = FALSE)
  
}

mvmr_cochran_Q <-function(r_input,gencov){
  
  # convert MRMVInput object to mvmr_format
  #if ("MRMVInput" %in% class(r_input)) {
  #  r_input <- mrmvinput_to_mvmr_format(r_input)
  #}
  
  # Perform check that r_input has been formatted using format_mvmr function
  #if(!("mvmr_format" %in%
  #     class(r_input))) {
  #  stop('The class of the data object must be "mvmr_format", please resave the object with the output of format_mvmr().')
  #}
  
  #gencov is the covariance between the effect of the genetic variants on each exposure.
  #By default it is set to 0.
  
  if(missing(gencov)) {
    gencov<-0
    warning("Covariance between effect of genetic variants on each exposure not specified. Fixing covariance at 0.")
  }
  
  # Inverse variance weighting is used.
  
  Wj<-1/r_input[,3]^2
  
  #Determine the number of exposures included in the model
  
  exp.number<-length(names(r_input)[-c(1,2,3)])/2
  
  #Fit the IVW MVMR model
  
  A_sum<-summary(lm(as.formula(paste("betaYG~ -1 +", paste(names(r_input)[
    seq(4,3+exp.number,by=1)], collapse="+")))
    ,weights=Wj,data=r_input))
  
  A<-summary(lm(as.formula(paste("betaYG~ -1 +", paste(names(r_input)[
    seq(4,3+exp.number,by=1)], collapse="+")))
    ,weights=Wj,data=r_input))$coef
  
  #Rename the regressors for ease of interpretation
  for(i in 1:exp.number){
    dimnames(A)[[1]][i]<- paste0("exposure",i,collapse="")
  }
  
  
  #Create a subset containing only standard errors for exposure effect estimates
  sebetas<-r_input[,(exp.number + 4):length(r_input)]
  
  
  ########################
  ## Instrument Validity #
  ########################
  
  if(length(gencov) <2){
    
    # Generate Sigma^2_A values
    sigma2A<-r_input[,3]^2
    for(i in 1:exp.number){
      sigma2A<-sigma2A + (A[i]^2 * sebetas[,i]^2)
    }
    
    #Create a subset of exposure effect estimates
    betas<-r_input[,c(4:(3+exp.number))]
    
    #Generates the component of the Q statistic to be subtracted from the outcome estimates
    temp.sub2<-0
    for(i in 1:exp.number){
      temp.sub2<-temp.sub2 + (betas[,i] * A[i])
    }
    
    #Calculates Q statistic for instrument validity
    Q_valid<- sum ((1/sigma2A)*(r_input[,2]-temp.sub2)^2)
    
    #Calculates p_value for instrument validity
    Q_chiValid<-pchisq(Q_valid,length(r_input[,2])-exp.number-1,lower.tail = FALSE)
    
    
  }
  
  if(length(gencov) >2){
    
    # Generate Sigma^2_A values
    sigma2A<-r_input[,3]^2
    for(i in 1:length(r_input[,3])){
      sigma2A[i]<-sigma2A[i] + (t(as.matrix(A[,1])) %*% gencov[[i]]%*% as.matrix(A[,1]))
    }
    
    #Create a subset of exposure effect estimates
    betas<-r_input[,c(4:(3+exp.number))]
    
    #Generates the component of the Q statistic to be subtracted from the outcome estimates
    temp.sub2<-0
    for(i in 1:exp.number){
      temp.sub2<-temp.sub2 + (betas[,i] * A[i])
    }
    
    #Calculates Q statistic for instrument validity
    Q_valid<- ((1/sigma2A)*(r_input[,2]-temp.sub2)^2)
    
  }
  
  
  ##########
  # Output #
  ##########
  
  cat("Q-Statistic for instrument validity:")
  
  multi_return <- function() {
    Out_list <- list("Qstat" = Q_valid)
    
    #Defines class of output object
    
    return(Out_list)
  }
  OUT<-multi_return()
}

##############
#Loading data#
##############

###############################################
#First we are gonna load the data that we need#
###############################################

lifetime_smk <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/Lifetime_smk/Lifetime_Smoking_Curated.txt")

lifetime_edu_df <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/Curated_data/LifeTime_edu_Risk_Factors.txt") #241

##############################
#Let's do some checks here...#
##############################

length(which(lifetime_edu_df$pval.smk < 0.00000005)) #72
length(which(lifetime_edu_df$pval.edu < 0.00000005)) #185

#this seems suspicious...

lifetime_edu_df[which.min(lifetime_edu_df$pval.smk),]

lifetime_smk[which(lifetime_smk$SNP == "rs8042849"),]

##############################
#It seems that we are fine???#
##############################

#Now let's add wc.

wc <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/2SMR/Curating_data/Curated_data/WC/WC_combined_Curated.txt")

######################
#Let's match the data#
######################

wc_match <- wc[which(wc$chr_pos_37%in%lifetime_edu_df$chr_pos),] #118
check <- wc[which(wc$MarkerName%in%lifetime_edu_df$SNP),] #118 too.

#That means that we need to run proxies for about 130 SNPs.
#Let's go:

missings_snps <- lifetime_edu_df[which(!(lifetime_edu_df$chr_pos%in%wc$chr_pos_37)),]

#Let's check if these SNPs are the same as in WHRAdjBMI..., 
#since they come from the same study, the possibility is high.

check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/Lifetime/proxies/storage/Lifetime_EDU/WHRAdjBMI/missing_snps.txt")

length(which(check$SNP%in%missings_snps$SNP)) #133/133.

#All of them there!! 
#Thus we can retrieve them.

proxies_df <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/Lifetime/proxies/storage/Lifetime_edu/WHRAdjBMI/combined_query_snp_list.txt", fill = TRUE)

#Perfect, now we have the data for the SNP aligned and ready to go.
#We have our own SNP as a proxy but with r2 == 1, distance = 0. We need to remove those.
#Or we will end up with the same results all the time.

proxies_df_clean <- proxies_df[which(!(proxies_df$query_snp == proxies_df$Coord & proxies_df$Distance == 0 & proxies_df$R2 == 1)),]

head(proxies_df_clean)

#Now we are only interested in those with a r2 > 0.8

proxies_df_clean_curated <- proxies_df_clean[which(proxies_df_clean$R2 > 0.8),] #6678

#Just in case, let's check the MAF:

summary(proxies_df_clean_curated$MAF) #perfect. No need to remove any.

#Cool, then we can take the proxies in two forms:

proxies_rsid <- proxies_df_clean_curated$RS_Number

proxies_chr_pos <- proxies_df_clean_curated$Coord

###############################################
#Let's get those that are present in wc#
###############################################

wc_proxies_match <- wc[which(wc$chr_pos_37%in%proxies_chr_pos),] #1838

#This is going to be always the final threshold since wc is a bottleneck. Hence...
#we are gonna use these SNPs to match the rest:

lifetime_proxies_match <- lifetime_smk[which(lifetime_smk$chr_pos%in%wc_proxies_match$chr_pos_37),] #1800/1838

summary(lifetime_proxies_match$EAF) #perfect

#And now we retrieve the edu proxies.
#REALLY IMPORTANT: proxies are turned off here. 
#What if we get, as a proxy, the original SNP we were using?
#It would be repeating all again.

edu_proxies_match <- extract_outcome_data(
  snps = lifetime_proxies_match$SNP,
  outcomes = 'ukb-b-16489',
  proxies = FALSE
)

#The number of proxies and non-proxies match..., so I am just gonna stick with these.

edu_proxies_match <- edu_proxies_match[which(edu_proxies_match$eaf.outcome > 0.01),] #1757
edu_proxies_match <- edu_proxies_match[which(edu_proxies_match$eaf.outcome < 0.99),] #1757

edu_proxies_match$chr_pos <- paste(edu_proxies_match$chr, edu_proxies_match$pos, sep = ":")
edu_proxies_match$chr_pos <- paste("chr", edu_proxies_match$chr_pos, sep = "")

#############################################
#Check alleles on the df that have less data#
#############################################

proxies_check <- proxies_df_clean_curated[which(proxies_df_clean_curated$Coord%in%edu_proxies_match$chr_pos),] #it seems that we have a duplicate, but that is expected.
proxies_check <- proxies_check[order(proxies_check$Alleles),]

head(proxies_check) #PERFECT
tail(proxies_check) #PERFECT

##########################
#We got all the we wanted#
##########################

#Let's summarize: 

#1. We got the proxies for the missing SNPs.
#2. Then we checked how many of those are in wc because it is the bottleneck.
#3. Then, we checked how many of those are present in lifetime smoking
#4. And how many of those are in edu.
#5. Now we need to merge the data between lifetime and edu.

#The next step is to obtain a dataframe with all the data for the proxies.
#And combine it with the original lifetime_edu data.

colnames(lifetime_proxies_match) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
                               "other_allele.exposure", "eaf.exposure", "info.exposure", "beta.exposure",
                               "se.exposure", "pval.exposure", "chr_pos")

lifetime_proxies_match$id.exposure <- "lifetime_smk"
lifetime_proxies_match$exposure <- "lifetime_smk"

dat_1_smk2edu_proxies <- harmonise_data(lifetime_proxies_match, edu_proxies_match, action = 3) #Note it removes the wrong triallelic combinations.

#We end up with 1792!

#Many of them are palindromic.

dat_1_smk2edu_proxies <- dat_1_smk2edu_proxies[which(dat_1_smk2edu_proxies$palindromic == FALSE),]
dat_1_smk2edu_proxies <- dat_1_smk2edu_proxies[which(dat_1_smk2edu_proxies$remove == FALSE),] #it removes the incompatible triallelic combination. But keeps the good ones. That is our boy!!

#In the last version we ended up with 1757, so we have a bit more... 
#The thing is... how is it going to be clumped?

#Now change the dataframe so it matches the original:

dat_1_smk2edu_proxies$samplesize.exposure <- 462690

dat_1_smk2edu_proxies_clean <- dat_1_smk2edu_proxies %>%
  select(chr_pos.x, SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure, samplesize.exposure, eaf.outcome, beta.outcome, se.outcome, pval.outcome, samplesize.outcome)

colnames(dat_1_smk2edu_proxies_clean) <- c("chr_pos", "SNP", "effect_allele", "other_allele", "eaf.smk", "beta.smk", "se.smk", "pval.smk", "samplesize.smk", "eaf.edu", "beta.edu", "se.edu", "pval.edu", "samplesize.edu")

lifetime_edu_df_not_missing <- lifetime_edu_df[which(lifetime_edu_df$chr_pos%in%wc$chr_pos_37),]

lifetime_edu_df_post_proxies <- rbind(lifetime_edu_df_not_missing, dat_1_smk2edu_proxies_clean)

#Only one thing lacking: the minimum p-value to get the proper lead SNP according to each of the phenotypes:

pval_min <- c()

for(i in seq(1, length(lifetime_edu_df_post_proxies$SNP))){
  
  pval_vect <- c(lifetime_edu_df_post_proxies$pval.smk[i], lifetime_edu_df_post_proxies$pval.edu[i])
  
  index <- which.min(pval_vect)
  
  final_pval <- pval_vect[index]
  
  pval_min <- c(pval_min, final_pval)
  
}

lifetime_edu_df_post_proxies$pval <- pval_min #worked perfectly.

head(lifetime_edu_df_post_proxies)

#And finally we get the RSIDs to prune them:

lifetime_edu_df_post_proxies$rsid <- lifetime_edu_df_post_proxies$SNP

#Let's do a quick check on the RSIDs, just to check that they are OK.

head(lifetime_edu_df_post_proxies$rsid[order(lifetime_edu_df_post_proxies$rsid)])

#Since we have duplicates, I will keep out one of them and check what happens:

which(duplicated(lifetime_edu_df_post_proxies$SNP) == TRUE) #there are always two associations with this SNP. Why?

#In the same study we have the same association... 
#Thus we are gonna remove the one that is not genome-wide for education
#Since we already have on that is:

lifetime_edu_df_post_proxies_single <- lifetime_edu_df_post_proxies[-which(lifetime_edu_df_post_proxies$SNP == "rs2668692" & lifetime_edu_df_post_proxies$pval.edu > 0.00000005),]

#########################
#Let's get the lead SNPs#
#########################

lifetime_edu_independent <- ieugwasr::ld_clump_local(lifetime_edu_df_post_proxies, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 1) #no need to do anything. LD r2 < 0.001 and kb = 10000.

#We end up with 212/251. Not bad!!

test <- ieugwasr::ld_clump_local(lifetime_edu_df_post_proxies_single, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 1) #no need to do anything. LD r2 < 0.001 and kb = 10000.

#If we take into account the independent one... we get the same amount of SNPs.
#Are they the same??

length(which(test$SNP%in%lifetime_edu_independent$SNP)) #They all match whether we remove the weird SNP or not. That is fine by me, then.

#########################################################
#Let's do a quick exploration of the snps in the proxies#
#########################################################

View(proxies_df_clean_curated[which(proxies_df_clean_curated$Coord%in%lifetime_edu_df_post_proxies$chr_pos),])

#############################################################
#Now we need to get the data from wc aligned and curated too#
#############################################################

#############################
#Let the harmonization begin#
#############################

#Since we remove the palindromic SNPs, this should be a problem.
#We are gonna match the alleles, the betas and the allele frequences in the old-fashioned way.

#1. We match the snps:

wc_outcome <- wc[which(wc$MarkerName%in%lifetime_edu_independent$SNP),] #206/206

#1. we order the SNPs:

wc_outcome <- wc_outcome[order(match(wc_outcome$MarkerName, lifetime_edu_independent$SNP)),]

#Let's check this:

length(which(wc_outcome$MarkerName != lifetime_edu_independent$SNP)) #0
length(which(wc_outcome$MarkerName == lifetime_edu_independent$SNP)) #206

#Good.

#2. let's switch the alleles, the beta and the eaf:

wc_outcome$final_a1 <- ifelse(wc_outcome$Allele1 != lifetime_edu_independent$effect_allele, wc_outcome$Allele2, wc_outcome$Allele1)
wc_outcome$final_a2 <- ifelse(wc_outcome$Allele1 != lifetime_edu_independent$effect_allele, wc_outcome$Allele1, wc_outcome$Allele2)
wc_outcome$final_beta <- ifelse(wc_outcome$Allele1 != lifetime_edu_independent$effect_allele, wc_outcome$b*(-1), wc_outcome$b)
wc_outcome$final_eaf <- ifelse(wc_outcome$Allele1 != lifetime_edu_independent$effect_allele, 1-wc_outcome$FreqAllele1HapMapCEU, wc_outcome$FreqAllele1HapMapCEU)

#Let's check:

length(which(wc_outcome$final_a1 != lifetime_edu_independent$effect_allele)) #0
length(which(wc_outcome$final_a1 == lifetime_edu_independent$effect_allele)) #206

length(which(wc_outcome$final_a2 != lifetime_edu_independent$other_allele)) #0
length(which(wc_outcome$final_a2 == lifetime_edu_independent$other_allele)) #206

head(wc_outcome) #worked like a charm. 

#WE CAN FINALLY RUN THIS.

################################################
#Now starts the fun. Let's follow MVMR pipeline#
################################################

library(MVMR)

path_ <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/Lifetime/proxies/storage/Lifetime_edu/WC/"

mvmr_Mario(lifetime_edu_independent, wc_outcome, path_)

#DONE! My functions does everything already :)

################################################################################
################################################################################

      # DEPRECTATED OLD CODE USED TO BUILD THE REST OF THE DATA #

################################################################################
################################################################################

#########################################
#STEP 1: construct the covariance matrix#
#########################################

#We first need a covariance matrix for each type of phenotypic correlation.

Pcov_0 <- matrix(data = c(1, 0, 0 , 1), nrow = 2, ncol = 2) #we make the fake phenotypic correlations
Pcov_0.25 <- matrix(data = c(1, 0.25, 0.25 , 1), nrow = 2, ncol = 2) #we make the fake phenotypic correlations
Pcov_0.5 <- matrix(data = c(1, 0.5, 0.5 , 1), nrow = 2, ncol = 2) #we make the fake phenotypic correlations
Pcov_0.75 <- matrix(data = c(1, 0.75, 0.75 , 1), nrow = 2, ncol = 2) #we make the fake phenotypic correlations
Pcov_1 <- matrix(data = c(1, 1, 1 , 1), nrow = 2, ncol = 2) #we make the fake phenotypic correlations

#Let's change the names of the matrices:

colnames(Pcov_0) <- c("1", "2")
rownames(Pcov_0) <- c("1", "2")

colnames(Pcov_0.25) <- c("1", "2")
rownames(Pcov_0.25) <- c("1", "2")

colnames(Pcov_0.5) <- c("1", "2")
rownames(Pcov_0.5) <- c("1", "2")

colnames(Pcov_0.75) <- c("1", "2")
rownames(Pcov_0.75) <- c("1", "2")

colnames(Pcov_1) <- c("1", "2")
rownames(Pcov_1) <- c("1", "2")

#Now we need a dataframe with the SEs.
#It needs to have the double of rows because MVMR is only made for <1 SNP.

seBXGs = cbind(lifetime_edu_independent$se.smk, lifetime_edu_independent$se.edu)

colnames(seBXGs) <- c("1", "2")

#And now we run the covariance matrix:

covariance_matrix_0 <- phenocov_mvmr(Pcov = Pcov_0, seBXGs = seBXGs)
covariance_matrix_0.25 <- phenocov_mvmr(Pcov = Pcov_0.25, seBXGs = seBXGs)
covariance_matrix_0.50 <- phenocov_mvmr(Pcov = Pcov_0.5, seBXGs = seBXGs)
covariance_matrix_0.75 <- phenocov_mvmr(Pcov = Pcov_0.75, seBXGs = seBXGs)
covariance_matrix_1 <- phenocov_mvmr(Pcov = Pcov_1, seBXGs = seBXGs)

#Maybe there is something wrong with my matrix. 
#I can't be 100% sure SINCE THERE ARE NO EXAMPLES ANYWHERE

##########################################
#Let's format the summary statistics data#
##########################################

BXGs <- cbind(lifetime_edu_independent$beta.smk, lifetime_edu_independent$beta.edu)

colnames(BXGs) <- c("1", "2")

BYG <- wc_outcome$final_beta

seBYG <- wc_outcome$se

final_df <- format_mvmr(BXGs = BXGs, seBXGs = seBXGs, BYG = BYG, seBYG = seBYG, RSID = lifetime_edu_independent$SNP)

##########################################
#Let's check for the strength of the SNPs#
##########################################

IV_0_test <- strength_mvmr(final_df) 

#exposure1 exposure2
#F-statistic  18.14334  21.88765

IV_0 <- strength_mvmr(final_df, covariance_matrix_0) 

#Conditional F-statistics for instrument strength

#exposure1 exposure2
#F-statistic  18.14334  21.88765

IV_0.25 <- strength_mvmr(final_df, covariance_matrix_0.25)

Conditional F-statistics for instrument strength

#exposure1 exposure2
#F-statistic  18.03631  21.76074

#IV_0.50 <- strength_mvmr(final_df, covariance_matrix_0.50)

Conditional F-statistics for instrument strength

#exposure1 exposure2
#F-statistic  17.93069  21.63533

IV_0.75 <- strength_mvmr(final_df, covariance_matrix_0.75)

#Conditional F-statistics for instrument strength

#exposure1 exposure2
#F-statistic  17.82645  21.51141

IV_1 <- strength_mvmr(final_df, covariance_matrix_1)

#Conditional F-statistics for instrument strength

#exposure1 exposure2
#F-statistic  17.72357  21.38895

#The strength of the F-statistic increased so much for smoking when we increased the correlation.
#We do not know how much the correlation is, but it is only over 25% that it started working.

################################
#This part is not working so...#
################################

pleoitropy_test <- pleiotropy_mvmr(final_df) 

#Q-Statistic for instrument validity:
#  297.9312 on 204 DF , p-value: 1.940013e-05

pleoitropy_0 <- pleiotropy_mvmr(final_df, covariance_matrix_0) #0.02. Pvalue = 0

#Q-Statistic for instrument validity:
#  297.9312 on 204 DF , p-value: 1.940013e-05

pleoitropy_0.25 <- pleiotropy_mvmr(final_df, covariance_matrix_0.25) #same

#Q-Statistic for instrument validity:
#  298.7749 on 204 DF , p-value: 1.683962e-05

pleoitropy_0.50 <- pleiotropy_mvmr(final_df, covariance_matrix_0.50) #same

#Q-Statistic for instrument validity:
#  299.6242 on 204 DF , p-value: 1.459209e-05

pleoitropy_0.75 <- pleiotropy_mvmr(final_df, covariance_matrix_0.75) #same

#Q-Statistic for instrument validity:
#  300.4792 on 204 DF , p-value: 1.262265e-05

pleoitropy_1 <- pleiotropy_mvmr(final_df, covariance_matrix_1) #same.

#Q-Statistic for instrument validity:
#  301.34 on 204 DF , p-value: 1.089987e-05

#Cochran's Q is strong.

############################
#Finally, let's run the ivw#
############################

ivw_mvmr(final_df)

#Multivariable MR

#Estimate Std. Error   t value     Pr(>|t|)
#exposure1  0.2838924 0.05872361  4.834382 2.616691e-06
#exposure2 -0.1550645 0.06727150 -2.305054 2.216439e-02

#Residual standard error: 1.214 on 205 degrees of freedom

#Both of them are significant despite all.

##############################
#Let's try the next one, then#
##############################

qhet_mvmr(final_df, pcor = Pcov_0)

#Effect Estimates
#Exposure 1        0.2981290
#Exposure 2       -0.1539708

qhet_mvmr(final_df, pcor = Pcov_0.25)

qhet_mvmr(final_df, pcor = Pcov_0.50)

qhet_mvmr(final_df, pcor = Pcov_0.75)

qhet_mvmr(final_df, pcor = Pcov_1)

###################
#Checking outliers#
###################

#First we generate the outliers:

snps <- final_df$SNP

outliers_0 <- mvmr_cochran_Q(final_df, covariance_matrix_0) 
outliers_0.25 <- mvmr_cochran_Q(final_df, covariance_matrix_0.25)
outliers_0.5 <- mvmr_cochran_Q(final_df, covariance_matrix_0.50) 
outliers_0.75 <- mvmr_cochran_Q(final_df, covariance_matrix_0.75) 
outliers_1 <- mvmr_cochran_Q(final_df, covariance_matrix_1) 

#Now let's filter those SNPs:

snps_cochran_0 <- as.data.frame(cbind(snps, unlist(outliers_0)))
snps_cochran_0 <- snps_cochran_0[order(as.numeric(snps_cochran_0$V2), decreasing = TRUE),]
View(snps_cochran)
snps_cochran_0 <- snps_cochran_0[which(as.numeric(snps_cochran_0$V2) > 5),]

final_df_filt_0 <- final_df[-which(final_df$SNP%in%snps_cochran_0$snps),]

####################
#Now let's run MVMR#
####################

seBXGs_filt_0

covariance_matrix_0 <- phenocov_mvmr(Pcov = Pcov_0, seBXGs = seBXGs)
covariance_matrix_0.25 <- phenocov_mvmr(Pcov = Pcov_0.25, seBXGs = seBXGs)
covariance_matrix_0.50 <- phenocov_mvmr(Pcov = Pcov_0.5, seBXGs = seBXGs)
covariance_matrix_0.75 <- phenocov_mvmr(Pcov = Pcov_0.75, seBXGs = seBXGs)
covariance_matrix_1 <- phenocov_mvmr(Pcov = Pcov_1, seBXGs = seBXGs)

ivw_mvmr(final_df_filt_0, covariance_matrix_0)

strength_mvmr(final_df_filt_0, covariance_matrix_0)

#Conditional F-statistics for instrument strength

#exposure1 exposure2
#F-statistic   19.1634  23.28395

pleoitropy_0 <- pleiotropy_mvmr(final_df_filt_0, covariance_matrix_0) #same


pleoitropy_0.25 <- pleiotropy_mvmr(final_df, covariance_matrix_0.25) #same

#Q-Statistic for instrument validity:
#  298.7749 on 204 DF , p-value: 1.683962e-05

pleoitropy_0.50 <- pleiotropy_mvmr(final_df, covariance_matrix_0.50) #same

#Q-Statistic for instrument validity:
#  299.6242 on 204 DF , p-value: 1.459209e-05

pleoitropy_0.75 <- pleiotropy_mvmr(final_df, covariance_matrix_0.75) #same

#Q-Statistic for instrument validity:
#  300.4792 on 204 DF , p-value: 1.262265e-05

pleoitropy_1 <- pleiotropy_mvmr(final_df, covariance_matrix_1) #same.

#Q-Statistic for instrument validity:
#  301.34 on 204 DF , p-value: 1.089987e-05
