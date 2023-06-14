##############
#INTRODUCTION#
##############

#This code will generate files for the loci of interest for the visual co-localization of the FIadjBMI-HDL-TG triangulation loci.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

start_parser <- function(loci){
  #This function takes the starting location of a loci defined as chr_start_end
  
  start_pos <- as.numeric(as.character(strsplit(loci, "_")[[1]][2]))
  
  return(start_pos)
  
}

end_parser <- function(loci){
  #This function takes the ending location of a loci defined as chr_start_end
  
  end_pos <- as.numeric(as.character(strsplit(loci, "_")[[1]][3]))
  
  return(end_pos)
  
}


loci_df_parser <- function(loci_df){
  #With this function we are going to parse the data of each loci.
  
  loci_df$start_pos <- as.numeric(as.character(unlist(sapply(loci_df$loci, start_parser))))
  loci_df$end_pos <- as.numeric(as.character(unlist(sapply(loci_df$loci, end_parser))))
  
  loci_tmp <- loci_df %>%
    dplyr::select(chromosome, start_pos, end_pos, chr_pos)
  
  return(loci_tmp)
  
}

##############
#Loading data#
##############

loci_information <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Colocalization/Colocalization/info_4_loci/1_curated_data_4_locus_zoom/cigday_4_locus_zoom.txt")

loci_information <- loci_information[which(loci_information$variant == "rs1317286"),]

loci_information$loci <- paste(loci_information$chromosome, "_", loci_information$base_pair_location-500000, "_", loci_information$base_pair_location+500000, sep="")

######################################################################
#STEP 1: now let's make new datasets to make them go through the data#
######################################################################

#This might work in all our SNPs, because in our data the chr_pos of the lead SNPs to direct to any other variant,
#but we will make sure of that.
#for the reference allele, we will take into consideration the alleles.

loci_information_4_lz <- loci_df_parser(loci_information)

#############################
#STEP 2: let's save the data#
#############################

write.table(loci_information_4_lz, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Colocalization/Colocalization/info_4_loci/2_input_commands/loci_of_interest.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


