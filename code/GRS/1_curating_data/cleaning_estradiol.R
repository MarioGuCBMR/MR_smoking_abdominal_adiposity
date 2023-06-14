library(vroom)
library(tidyverse)
library(data.table)

male <- as.data.frame(vroom("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/estradiol_male.txt")) # 7783170
female <- as.data.frame(vroom("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/estradiol_female.txt")) # 7782036 
estradiol <- as.data.frame(vroom("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/metaanalysis_data/METAANALYSIS_estradiol2_1.tbl")) # 7464078
colnames(estradiol) <- c("variant", "effect_allele", "other_allele", "effect_allele_frequency", "freqse", "beta", "standard_error", "p_value", "direction") 
# Missing - chr, pos, sample_size
estradiol$effect_allele <- toupper(estradiol$effect_allele)
estradiol$other_allele <- toupper(estradiol$other_allele)
estradiol$chr_pos <- NA
estradiol$base_pair_location <- NA
estradiol$chromosome <- NA
estradiol$sample_size <- 311675

matching_rows <- estradiol$variant %in% male$variant # 7464078

estradiol$chromosome[matching_rows] <- male$chromosome[match(estradiol$variant[matching_rows], male$variant)]
estradiol$base_pair_location[matching_rows] <- male$base_pair_location[match(estradiol$variant[matching_rows], male$variant)]
estradiol$chr_pos[matching_rows] <- male$chr_pos[match(estradiol$variant[matching_rows], male$variant)]

subsetted_df <- select(estradiol, variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, chr_pos)
colnames(subsetted_df) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")

# I have realized there are X chromosomes... 
str(subsetted_df)
# 'data.frame':	7464078 obs. of  11 variables:
#   $ variant                : chr  "rs2326918" "rs56252526" "rs977590" "rs66941928" ...
# $ chromosome             : chr  "6" "6" "7" "3" ...
# $ base_pair_location     : num  1.31e+08 1.32e+08 3.44e+07 1.77e+08 2.57e+07 ...
# $ effect_allele          : chr  "A" "A" "A" "T" ...
# $ other_allele           : chr  "G" "G" "G" "C" ...
# $ effect_allele_frequency: num  0.846 0.988 0.936 0.802 0.764 ...
# $ beta                   : num  -0.0166 0.0466 0.0067 -0.0018 -0.0128 0.0116 -0.0096 0.0109 0.0106 0.0031 ...
# $ standard_error         : num  0.0122 0.0425 0.0181 0.011 0.0103 0.0154 0.0092 0.0132 0.0125 0.0092 ...
# $ p_value                : num  0.173 0.272 0.712 0.873 0.215 ...
# $ sample_size            : num  311675 311675 311675 311675 311675 ...
# $ chr_pos                : chr  "chr6:130840091" "chr6:132270085" "chr7:34448765" "chr3:176666749" ...

fwrite(subsetted_df, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/curated_estradiol_includesX.txt")

noX_subsetted_df = subsetted_df
df <- noX_subsetted_df[grepl("^[0-9]+$", noX_subsetted_df$chromosome), ]
fwrite(df, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/curated_estradiol.txt")
