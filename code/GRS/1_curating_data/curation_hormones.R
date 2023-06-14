
## DATA CURATION (Hormones)
### Files to curate: 
# Cortisol - NEEDS ANOTHER FUNCTION
# Testosterone
# SHBG

library(vroom)
library(data.table)

# TESTOSTERONE AND SHBG CURATION
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

# TESTOSTERONE
# > head(data, 1)
# variant_id chromosome base_pair_location effect_allele other_allele effect_allele_frequency imputation_quality      beta standard_error p_value
# 1 rs533090414          1              18849             C            G               0.0243876           0.536412 -0.014725     0.00582251   0.009
input_path <- "/emc/cbmr/users/wfs758/smoking_project/data/raw_data/testosterone_hormones.tsv.gz"
output_path <- "/emc/cbmr/users/wfs758/smoking_project/data/curated_data/testosterone.txt"
hormones_data_curation(input_path, 425097, output_path)
# >   print(summary(data$effect_allele_frequency))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0010  0.7664  0.9707  0.8366  0.9966  0.9990 
# >   print(summary(data_maf$effect_allele_frequency))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0100  0.5841  0.8423  0.7390  0.9567  0.9900 
# >   print(summary(as.numeric(data_maf_mhc$chromosome))) #perfect.
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6       6       6       6       6       6 
# >   print(summary(as.numeric(data_maf_mhc$base_pair_location))) #perfect.
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 26000393 29854418 31304538 30991232 32575070 33999992 
# >   print(head(data_end))
# variant chromosome base_pair_location effect_allele other_allele effect_allele_frequency         beta standard_error p_value sample_size     chr_pos
# 1  rs533090414          1              18849             C            G               0.0243876 -0.014725000     0.00582251   0.009      425097  chr1:18849
# 4  rs568927205          1              54712             T        TTTTC               0.4144420  0.000484583     0.00178834   0.870      425097  chr1:54712
# 9         <NA>          1              82133           CAA            C               0.0182434 -0.003394390     0.00706713   0.710      425097  chr1:82133
# 13  rs28619159          1              91421             T            C               0.9829950  0.005176830     0.00688747   0.540      425097  chr1:91421
# 16 rs377100675          1             125271             C            T               0.0379787 -0.001196640     0.00476360   0.880      425097 chr1:125271
# 18  rs62642131          1             135982             A            G               0.0126978  0.005685250     0.00849939   0.460      425097 chr1:135982

# SHBG
input_path <- "/emc/cbmr/users/wfs758/smoking_project/data/raw_data/shbg_hormones.tsv.gz"
output_path <- "/emc/cbmr/users/wfs758/smoking_project/data/curated_data/shbg.txt"
hormones_data_curation(input_path, 370125, output_path)
# [1] 16582964       10
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0010  0.7665  0.9707  0.8366  0.9966  0.9990 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0100  0.5841  0.8423  0.7390  0.9567  0.9900 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6       6       6       6       6       6 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 26000393 29854577 31304726 30991856 32575167 33999992 
# [1] 10014857       11
# variant chromosome base_pair_location effect_allele other_allele effect_allele_frequency         beta standard_error p_value sample_size     chr_pos
# 1  rs533090414          1              18849             C            G               0.0243762  0.000100134     0.00420384   0.890      370125  chr1:18849
# 4  rs568927205          1              54712             T        TTTTC               0.4142650  0.002382490     0.00129031   0.099      370125  chr1:54712
# 9         <NA>          1              82133           CAA            C               0.0182309  0.001597830     0.00510779   0.750      370125  chr1:82133
# 13  rs28619159          1              91421             T            C               0.9829760 -0.003305900     0.00497020   0.370      370125  chr1:91421
# 16 rs377100675          1             125271             C            T               0.0379328  0.002309270     0.00344015   0.380      370125 chr1:125271
# 18  rs62642131          1             135982             A            G               0.0126889 -0.006014340     0.00612645   0.300      370125 chr1:135982

## CORTISOL CURATION
data <- as.data.frame(vroom("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/raw_summary_statistics/cortisol.txt"))
output_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_data/cortisol.txt"
print(head(data))
# > print(head(data))
# SNP Chromosome Position Allele1 Allele2  Freq1 FreqSE MinFreq MaxFreq  Effect StdErr      Pval         Direction HetISq HetChiSq HetDf   HetPVal     N
# rs9989237         14 94795202       t       c 0.2104 0.0259  0.0629  0.2300  0.0857 0.0095 2.157e-19 +++++++++++++++++   60.4   40.410    16 0.0006781 25314
print(summary(data$Freq1))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0101  0.0958  0.4091  0.4589  0.8317  0.9899 
# 1. Remove rare variants
data_maf <- data[which(data$Freq1 > 0.01),]
data_maf <- data_maf[which(data_maf$Freq1 < 0.99),]
print(summary(data_maf$Freq1))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0101  0.0958  0.4091  0.4589  0.8317  0.9899
# 2. Remove the MHC region:
data_maf$chr_pos <- paste("chr", data_maf$Chromosome, ":", data_maf$Position, sep = "") # add chr:pos column
data_maf_mhc <- data_maf[which(as.numeric(data_maf$Chromosome) == 6 & as.numeric(data_maf$Position) >= 26000000 & as.numeric(data_maf$Position) <= 34000000),]
print(summary(as.numeric(data_maf_mhc$Chromosome))) #perfect.
print(summary(as.numeric(data_maf_mhc$Position))) #perfect.
# > print(summary(as.numeric(data_maf_mhc$Chromosome))) #perfect.
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6       6       6       6       6       6 
# > print(summary(as.numeric(data_maf_mhc$Position))) #perfect.
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 26000393 29870326 31311242 31016976 32559683 33999992
#Now let's check if we had any interesting variants there:
data_maf_no_mhc <- data_maf[which(!(data_maf$chr_pos%in%data_maf_mhc$chr_pos)),]
data_maf_no_mhc$sample_size <- 25314
data_end <- data_maf_no_mhc %>%
  select(SNP, Chromosome, Position, Allele1, Allele2, Freq1, Effect, StdErr, Pval, sample_size, chr_pos)
colnames(data_end) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", 
                        "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")
print(dim(data_end))
print(head(data_end))
# variant chromosome base_pair_location effect_allele other_allele effect_allele_frequency    beta standard_error   p_value sample_size        chr_pos
# 1  rs9989237         14           94795202             t            c                  0.2104  0.0857         0.0095 2.157e-19       25314 chr14:94795202
# 2  rs6575415         14           94791601             a            t                  0.2104  0.0855         0.0095 2.970e-19       25314 chr14:94791601
# 3  rs7161521         14           94787288             t            c                  0.2104  0.0855         0.0095 3.073e-19       25314 chr14:94787288
# 4 rs12589136         14           94793686             t            g                  0.2101  0.0853         0.0095 3.226e-19       25314 chr14:94793686
# 5   rs941599         14           94788341             t            c                  0.2105  0.0852         0.0095 4.406e-19       25314 chr14:94788341
# 6   rs718187         14           94801860             t            c                  0.7897 -0.0849         0.0095 4.520e-19       25314 chr14:94801860
fwrite(data_end, output_path)
any(is.na(data_end$variant)) # FALSE
grepl("NA", data_end$variant)
df_sorted <- data_end[order(data_end$variant),]
head(df_sorted)
tail(df_sorted)
