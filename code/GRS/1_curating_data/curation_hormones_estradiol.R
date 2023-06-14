
library(vroom)
library(data.table)
library(tidyverse)

hormones_data_curation <- function(input_path, N, output_path){
  data <- as.data.frame(vroom(input_path))
  print(dim(data))
  print(summary(data$effect_allele_frequency))
  
  # 1. Remove rare variants
  data_maf <- data[which(data$effect_allele_frequency > 0.01),]
  data_maf <- data_maf[which(data_maf$effect_allele_frequency < 0.99),]
  
  print(summary(data_maf$effect_allele_frequency))
  
  # 2. Remove the MHC region:
  data_maf$chr_pos <- paste("chr", data_maf$chromosome, ":", data_maf$base_pair_location, sep = "") # add chr:pos column
  data_maf_mhc <- data_maf[which(as.numeric(data_maf$chromosome) == 6 & as.numeric(data_maf$base_pair_location) >= 26000000 & as.numeric(data_maf$base_pair_location) <= 34000000),]
  
  print(summary(as.numeric(data_maf_mhc$chromosome))) #perfect.
  print(summary(as.numeric(data_maf_mhc$base_pair_location))) #perfect.
  
  #Now let's check if we had any interesting variants there:
  data_maf_no_mhc <- data_maf[which(!(data_maf$chr_pos%in%data_maf_mhc$chr_pos)),]
  
  data_maf_no_mhc$sample_size <- N
  
  ################################################################
  #Let's make chromosome and position and then rearrange the data#
  ################################################################
  
  data_end <- data_maf_no_mhc %>%
    select(variant_id, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, chr_pos)
  
  colnames(data_end) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", 
                          "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")
  
  #In the whole process we will be very careful if one of the SNPs RSID is not there.
  #Because if that is the case we will have a problem when clumping.
  print(dim(data_end))
  print(head(data_end))
  fwrite(data_end, output_path)
}

# GCST90020091 = male ; 147690
input_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/raw_summary_statistics/GCST90020091_buildGRCh38.tsv"
# variant_id  p_value chromosome base_pair_location effect_allele other_allele effect_allele_frequency odds_ration ci_lower ci_upper        beta
# 1 rs556128788 0.750698          1              99388             T            G               0.0125806    0.975707 0.838362  1.13555 -0.02459294
# standard_error
# 1      0.0774054
output_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/estradiol_male.txt"
hormones_data_curation(input_path, 147690, output_path)

# [1] 7871694      12
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001153 0.037808 0.127634 0.172665 0.288100 0.516357 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01000 0.03858 0.12883 0.17350 0.28910 0.51636 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6       6       6       6       6       6 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 26000393 29770528 31243184 30872844 32587716 33999992 
# [1] 7783170      11

# GCST90020092 = female ; 163985
input_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/raw_summary_statistics/GCST90020092_buildGRCh38.tsv"
output_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/estradiol_female.txt"
hormones_data_curation(input_path, 163985, output_path)
# [1] 7870546      12
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001252 0.037813 0.127618 0.172669 0.288090 0.516368 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0100  0.0386  0.1289  0.1735  0.2891  0.5164 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6       6       6       6       6       6 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 26000393 29770079 31242695 30872468 32587702 33999992 
# [1] 7782036      11
