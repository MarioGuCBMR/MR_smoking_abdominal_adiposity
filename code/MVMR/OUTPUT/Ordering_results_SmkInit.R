##############
#INTRODUCTION#
##############

#This is a code to open data results for MVMR for each combination of trait and see wassup with dat.

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)

###################
#Loading functions#
###################

saving_original <- function(data, path_){
  
  mF_ordered <- c()
  
  for(i in seq(1, length(data$Smoking_mF))){
    
    mF_ordered <- c(mF_ordered, data$Smoking_mF[i]) #first Smoking...
    mF_ordered <- c(mF_ordered, data$Edu_mF[i]) #Then edu...
    
  }
  
  #Now we do the same with Cochran's Q data:
  
  Q_ordered <- c()
  
  for(i in seq(1, length(data$Smoking_mF))){
    
    Q_ordered <- c(Q_ordered, unlist(data$Qstat_original[i])) #first qval...
    Q_ordered <- c(Q_ordered, unlist(data$Qpval_original[i])) #Then qstat...
    
  }
  
  #####################################################
  #Saving the data in the path and naming it correctly#
  #####################################################
  
  final_path_mf <- paste(path_, "mF_original.txt", sep = "")
  final_path_Q <- paste(path_, "Cochrans_Q_test_original.txt", sep = "")
  
  ####################
  #Let's save this...#
  ####################
  
  write.table(mF_ordered, final_path_mf, row.names = FALSE, quote = FALSE, col.names = FALSE)
  write.table(Q_ordered, final_path_Q, row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  
}

saving_lax <- function(data, path_){
  
  mF_ordered <- c()
  
  for(i in seq(1, length(data$Smoking_mF_lax_outlier_out))){
    
    mF_ordered <- c(mF_ordered, data$Smoking_mF_lax_outlier_out[i]) #first Smoking...
    mF_ordered <- c(mF_ordered, data$Edu_mF_lax_outliers_out[i]) #Then edu...
    
  }
  
  #Now we do the same with Cochran's Q data:
  
  Q_ordered <- c()
  
  for(i in seq(1, length(data$Smoking_mF_lax_outlier_out))){
    
    Q_ordered <- c(Q_ordered, unlist(data$Qstat_lax_outlier_out[i])) #first qval...
    Q_ordered <- c(Q_ordered, unlist(data$Qpval_lax_outlier_out[i])) #Then qstat...
    
  }
  
  #####################################################
  #Saving the data in the path and naming it correctly#
  #####################################################
  
  final_path_mf <- paste(path_, "mF_lax.txt", sep = "")
  final_path_Q <- paste(path_, "Cochrans_Q_test_lax.txt", sep = "")
  
  ####################
  #Let's save this...#
  ####################
  
  write.table(mF_ordered, final_path_mf, row.names = FALSE, quote = FALSE, col.names = FALSE)
  write.table(Q_ordered, final_path_Q, row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  
}

saving_strict <- function(data, path_){
  
  mF_ordered <- c()
  
  for(i in seq(1, length(data$Smoking_mF_strict_outlier_out))){
    
    mF_ordered <- c(mF_ordered, data$Smoking_mF_strict_outlier_out[i]) #first Smoking...
    mF_ordered <- c(mF_ordered, data$Edu_mF_strict_outliers_out[i]) #Then edu...
    
  }
  
  #Now we do the same with Cochran's Q data:
  
  Q_ordered <- c()
  
  for(i in seq(1, length(data$Smoking_mF_strict_outlier_out))){
    
    Q_ordered <- c(Q_ordered, unlist(data$Qstat_strict_outlier_out[i])) #first qval...
    Q_ordered <- c(Q_ordered, unlist(data$Qpval_strict_outlier_out[i])) #Then qstat...
    
  }
  
  #####################################################
  #Saving the data in the path and naming it correctly#
  #####################################################
  
  final_path_mf <- paste(path_, "mF_strict.txt", sep = "")
  final_path_Q <- paste(path_, "Cochrans_Q_test_strict.txt", sep = "")
  
  ####################
  #Let's save this...#
  ####################
  
  write.table(mF_ordered, final_path_mf, row.names = FALSE, quote = FALSE, col.names = FALSE)
  write.table(Q_ordered, final_path_Q, row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  
}

csv_save_original <- function(csv_dat, path_){
  
  #################################
  #Selecting the data that we want#
  #################################
  
  final_data <- csv_dat %>%
    select(Estimate_original, `Std. Error_original`, `Pr(>|t|)_original`)
  
  final_path_data <- paste(path_, "results_original.txt", sep = "")

  ####################
  #Let's save this...#
  ####################
  
  write.table(final_data, final_path_data, row.names = FALSE, quote = FALSE, col.names = FALSE)
  
}

csv_save_lax <- function(csv_dat, path_){
  
  #################################
  #Selecting the data that we want#
  #################################
  
  final_data <- csv_dat %>%
    select(Estimate_lax_outlier_out, `Std. Error_lax_outlier_out`, `Pr(>|t|)_lax_outlier_edu`)
  
  final_path_data <- paste(path_, "results_lax.txt", sep = "")
  
  ####################
  #Let's save this...#
  ####################
  
  write.table(final_data, final_path_data, row.names = FALSE, quote = FALSE, col.names = FALSE)
  
}

csv_save_strict <- function(csv_dat, path_){
  
  #################################
  #Selecting the data that we want#
  #################################
  
  final_data <- csv_dat %>%
    select(Estimate_strict_outlier_out, `Std. Error_strict_outlier_out`, `Pr(>|t|)_strict_outlier_edu`)
  
  final_path_data <- paste(path_, "results_strict.txt", sep = "")
  
  ####################
  #Let's save this...#
  ####################
  
  write.table(final_data, final_path_data, row.names = FALSE, quote = FALSE, col.names = FALSE)
  
}

##############
#Loading data#
##############

pleiotropy_data_whr <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/MVMR/Smoking_Initiation/OUTPUT/WHR/conditional_mF_cochrans_Q")
pleiotropy_data_whradjbmi <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/MVMR/Smoking_Initiation/OUTPUT/WHRAdjBMI/conditional_mF_cochrans_Q")
pleiotropy_data_wc <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/MVMR/Smoking_Initiation/OUTPUT/WC/conditional_mF_cochrans_Q")
pleiotropy_data_wcadjbmi <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/MVMR/Smoking_Initiation/OUTPUT/WCAdjBMI/conditional_mF_cochrans_Q")
pleiotropy_data_hip <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/MVMR/Smoking_Initiation/OUTPUT/HIP/conditional_mF_cochrans_Q")
pleiotropy_data_hipadjbmi <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/MVMR/Smoking_Initiation/OUTPUT/HIPAdjBMI/conditional_mF_cochrans_Q")


#########################
#Getting the paths ready#
#########################

path_whr <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/OUTPUT/SmkInit_results/WHR/WHR_"
path_wc <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/OUTPUT/SmkInit_results/WC/WC_"
path_whradjbmi <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/OUTPUT/SmkInit_results/WHRAdjBMI/WHRAdjBMI_"
path_wcadjbmi <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/OUTPUT/SmkInit_results/WCAdjBMI/WCAdjBMI_"
path_hip <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/OUTPUT/SmkInit_results/HIP/HIP_"
path_hipadjbmi <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/OUTPUT/SmkInit_results/HipAdjBMI/HipAdjBMI_"

#####################
#Let's run this baby#
#####################

#First the originals:

saving_original(pleiotropy_data_whr, path_ = path_whr)
saving_original(pleiotropy_data_whradjbmi, path_ = path_whradjbmi)
saving_original(pleiotropy_data_wc, path_ = path_wc)
saving_original(pleiotropy_data_wcadjbmi, path_ = path_wcadjbmi)
saving_original(pleiotropy_data_hip, path_ = path_hip)
saving_original(pleiotropy_data_hipadjbmi, path_ = path_hipadjbmi)

#Now the lax version:

saving_lax(pleiotropy_data_whr, path_ = path_whr)
saving_lax(pleiotropy_data_whradjbmi, path_ = path_whradjbmi)
saving_lax(pleiotropy_data_wc, path_ = path_wc)
saving_lax(pleiotropy_data_wcadjbmi, path_ = path_wcadjbmi)
saving_lax(pleiotropy_data_hip, path_ = path_hip)
saving_lax(pleiotropy_data_hipadjbmi, path_ = path_hipadjbmi)

#Finally the strict:

saving_strict(pleiotropy_data_whr, path_ = path_whr)
saving_strict(pleiotropy_data_whradjbmi, path_ = path_whradjbmi)
saving_strict(pleiotropy_data_wc, path_ = path_wc)
saving_strict(pleiotropy_data_wcadjbmi, path_ = path_wcadjbmi)
saving_strict(pleiotropy_data_hip, path_ = path_hip)
saving_strict(pleiotropy_data_hipadjbmi, path_ = path_hipadjbmi)

###################################################
#Let's see if I can save the data from the csvs...#
###################################################

whr_csv <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/MVMR/Smoking_Initiation/OUTPUT/WHR/mvmr_ivw_results.csv")

#YES!! The data might be weird when seen in excel, but R reads it properly...

#we will save it again to avoid dramas.

whradjbmi_csv <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/MVMR/Smoking_Initiation/OUTPUT/WHRAdjBMI/mvmr_ivw_results.csv")
wc_csv <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/MVMR/Smoking_Initiation/OUTPUT/WC/mvmr_ivw_results.csv")
wcadjbmi_csv <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/MVMR/Smoking_Initiation/OUTPUT/WCAdjBMI/mvmr_ivw_results.csv")

hip_csv <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/MVMR/Smoking_Initiation/OUTPUT/HIP/mvmr_ivw_results.csv")
hipadjbmi_csv <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/MVMR/CODE/MVMR/Smoking_Initiation/OUTPUT/HIPAdjBMI/mvmr_ivw_results.csv")


#######################
#Saving the data again#
#######################

fwrite(whr_csv, paste(path_whr, "mvmr_ivw_results.csv", sep = ""))
fwrite(whradjbmi_csv, paste(path_whradjbmi, "mvmr_ivw_results.csv", sep = ""))
fwrite(wc_csv, paste(path_wc, "mvmr_ivw_results.csv", sep = ""))
fwrite(wcadjbmi_csv, paste(path_wcadjbmi, "mvmr_ivw_results.csv", sep = ""))
fwrite(hip_csv, paste(path_hip, "mvmr_ivw_results.csv", sep = ""))
fwrite(hipadjbmi_csv, paste(path_hipadjbmi, "mvmr_ivw_results.csv", sep = ""))

###########################################################################################################################
# Let's save the data in the prefered column format so we can copy paste it in the same way as we did with the pleiotropy #
###########################################################################################################################

#First the originals:

csv_save_original(whr_csv, path_ = path_whr)
csv_save_original(whradjbmi_csv, path_ = path_whradjbmi)
csv_save_original(wc_csv, path_ = path_wc)
csv_save_original(wcadjbmi_csv, path_ = path_wcadjbmi)
csv_save_original(hip_csv, path_ = path_hip)
csv_save_original(hipadjbmi_csv, path_ = path_hipadjbmi)

#Now the lax:

csv_save_lax(whr_csv, path_ = path_whr)
csv_save_lax(whradjbmi_csv, path_ = path_whradjbmi)
csv_save_lax(wc_csv, path_ = path_wc)
csv_save_lax(wcadjbmi_csv, path_ = path_wcadjbmi)
csv_save_lax(hip_csv, path_ = path_hip)
csv_save_lax(hipadjbmi_csv, path_ = path_hipadjbmi)

#Now the strict:

csv_save_strict(whr_csv, path_ = path_whr)
csv_save_strict(whradjbmi_csv, path_ = path_whradjbmi)
csv_save_strict(wc_csv, path_ = path_wc)
csv_save_strict(wcadjbmi_csv, path_ = path_wcadjbmi)
csv_save_strict(hip_csv, path_ = path_hip)
csv_save_strict(hipadjbmi_csv, path_ = path_hipadjbmi)
