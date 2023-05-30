##############
#INTRODUCTION#
##############

#In this very small code we are gonna check the number of cases and controls for Smoking Initiation from Liu et al 2019.
#This info has been extracted from supplementary table number 7: description of binary traits from Liu et al 2019.

#Obviously, the column taken is that for smoking initiation.

###########
#Libraries#
###########

library(data.table)
library(tidyverse)
library(readxl)

##############
#Loading data#
##############

liu_data <- readxl::read_xlsx("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/RAW_Data/Cases_Controls_4_Smoking_Initiation_Liu_et_al_2019.xlsx")

liu_data_df <- as.data.frame(liu_data)

##########################
#QC the data that we have#
##########################

#We have two columns: the total samplesize for the data and the proportion of Ever Smoker. 
#The objective is calculate the amount of cases (ever smoker) and controls (never smoker).

#Let's go.

liu_data_df_with_sample <- liu_data_df[which(!(liu_data_df$N == "-")),] #we remove the studies that do not have data for our trait.

#Now we are gonna make sure that the total sum of the N column matches with the total sample size that is described in the paper.

sum(as.numeric(liu_data_df_with_sample$N)) #1232091 -> matches perfectly the paper: (SmkInit, N=1,232,091; 378 associated variants).

###################
#Getting the cases#
###################

#This is the easiest thing in the world. In Ever Smoker we have the proportions of cases. They are percentages/100. 
#Thus, we just need to do the product between them.

liu_data_df_with_sample$cases <- as.numeric(liu_data_df_with_sample$N)*as.numeric(liu_data_df_with_sample$`Ever Smoker`)

#Let's check:

head(liu_data_df_with_sample)

#N Ever Smoker     cases
#1 599289        0.41 245708.49
#2  11345        0.42   4764.90
#4   4293        0.65   2790.45
#5   1238        0.46    569.48
#6   1192        0.85   1013.20
#9  57097         0.7  39967.90

#Makes all the sense in the world.

#########################
#Generating the controls#
#########################

liu_data_df_with_sample$controls <- as.numeric(liu_data_df_with_sample$N)-as.numeric(liu_data_df_with_sample$cases)

#Let's check:

head(liu_data_df_with_sample)

#N Ever Smoker     cases  controls
#1 599289        0.41 245708.49 353580.51
#2  11345        0.42   4764.90   6580.10
#4   4293        0.65   2790.45   1502.55
#5   1238        0.46    569.48    668.52
#6   1192        0.85   1013.20    178.80
#9  57097         0.7  39967.90  17129.10

#############
#FINAL CHECK#
#############

#The final check is producing the percentages ourselves:

liu_data_df_with_sample$perc_cases <- as.numeric(liu_data_df_with_sample$cases)/as.numeric(liu_data_df_with_sample$N)
liu_data_df_with_sample$perc_control <- as.numeric(liu_data_df_with_sample$controls)/as.numeric(liu_data_df_with_sample$N)

head(liu_data_df_with_sample)

#N Ever Smoker     cases  controls perc_cases perc_control
#1 599289        0.41 245708.49 353580.51       0.41         0.59
#2  11345        0.42   4764.90   6580.10       0.42         0.58
#4   4293        0.65   2790.45   1502.55       0.65         0.35
#5   1238        0.46    569.48    668.52       0.46         0.54
#6   1192        0.85   1013.20    178.80       0.85         0.15
#9  57097         0.7  39967.90  17129.10       0.70         0.30

##################
#FREAKING PERFECT#
##################

################################################
#CALCULATING TOTAL NUMBER OF CASES AND CONTROLS#
################################################

sum(as.numeric(liu_data_df_with_sample$cases)) #557337.1
sum(as.numeric(liu_data_df_with_sample$controls)) #674753.9

#We can approximate this to 557337 cases and 674754 controls.

#############
#WE ARE DONE#
#############
