
## DATA CURATION (FAT DEPOTS)

library(vroom)
library(data.table)

fat_data_curation <- function(input_path, N, output_path){
  data <- as.data.frame(vroom(input_path))
  print(dim(data))
  print(summary(data$A1FREQ))
  
  # 1. Remove rare variants
  data_maf <- data[which(data$A1FREQ > 0.01),]
  data_maf <- data_maf[which(data_maf$A1FREQ < 0.99),]
  
  print(summary(data_maf$A1FREQ))
  
  # 2. Remove the MHC region:
  data_maf$chr_pos <- paste("chr", data_maf$CHR, ":", data_maf$BP, sep = "") # add chr:pos column
  data_maf_mhc <- data_maf[which(as.numeric(data_maf$CHR) == 6 & as.numeric(data_maf$BP) >= 26000000 & as.numeric(data_maf$BP) <= 34000000),]
  
  print(summary(as.numeric(data_maf_mhc$CHR))) #perfect.
  print(summary(as.numeric(data_maf_mhc$BP))) #perfect.
  
  #Now let's check if we had any interesting variants there:
  #data_maf_no_mhc <- data_maf[which(!(data_maf$chr_pos%in%data_maf_mhc$chr_pos)),]
  
  data_maf_no_mhc$sample_size <- N
  
  ################################################################
  #Let's make chromosome and position and then rearrange the data#
  ################################################################
  
  data_end <- data_maf_no_mhc %>%
    select(SNP, CHR, BP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF, sample_size, chr_pos)
  
  colnames(data_end) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", 
                          "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")
  
  #In the whole process we will be very careful if one of the SNPs RSID is not there.
  #Because if that is the case we will have a problem when clumping.
  print(dim(data_end))
  fwrite(data_end, output_path)
}

for (file in list.files("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/raw_summary_statistics",
                        pattern = "0321", full.names = T)){
  name = gsub("_bgen_stats.gz", "", gsub("0321_", "", strsplit(file, "/")[[1]][8]))
  fat_data_curation(file, 39076, paste0("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/", name, ".txt"))
}
