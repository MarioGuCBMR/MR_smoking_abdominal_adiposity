library(vroom)
library(tidyverse)
library(dplyr)
library(data.table)

estradiol <- as.data.frame(vroom("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/metaanalysis_data/METAANALYSIS_estradiol2_1.tbl")) # 7464078
colnames(estradiol) <- c("variant", "effect_allele", "other_allele", "effect_allele_frequency", "freqse", "beta", "standard_error", "p_value", "direction") 
# Missing - chr, pos, sample_size
estradiol$effect_allele <- toupper(estradiol$effect_allele)
estradiol$other_allele <- toupper(estradiol$other_allele)
estradiol$sample_size <- 311675

# > head(estradiol)
# variant effect_allele other_allele effect_allele_frequency freqse    beta standard_error p_value direction
# 1  rs2326918             A            G                  0.8459  6e-04 -0.0166         0.0122  0.1728        -+
#   2 rs56252526             A            G                  0.9882  1e-04  0.0466         0.0425  0.2724        ++
#   3   rs977590             A            G                  0.9361  4e-04  0.0067         0.0181  0.7118        +-
#   4 rs66941928             T            C                  0.8025  5e-04 -0.0018         0.0110  0.8729        +-
#   5 rs17278013             T            G                  0.7644  0e+00 -0.0128         0.0103  0.2153        --
#   6 rs62234673             T            C                  0.0885  5e-04  0.0116         0.0154  0.4508        ++

estradiol <- estradiol[, c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele",           
                           "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size",            
                           "chr_pos")]

male <- as.data.frame(vroom("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/estradiol_male.txt")) # 7783170
# female <- as.data.frame(vroom("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/estradiol_female.txt")) # 7782036 

male_sel <- male[, c("variant", "chromosome", "base_pair_location", "chr_pos")]

new_estradiol <- merge(estradiol, male_sel, by = "variant")
head(new_estradiol)
tail(new_estradiol)
any(is.na(new_estradiol$chromosome))
any(is.na(new_estradiol$base_pair_location))
any(is.na(new_estradiol$chr_pos))

# estradiol dim - 7464078
# new estradiol dim - 7783170

noNA_estradiol <- new_estradiol[complete.cases(new_estradiol$variant), ] # 7463479
# this means that there is another variant that has the same rsID

### SHOULD I REMOVE THE NA - IT IS ADDING A LOT OF NAs FROM MALE THAT WERE NOT IN THE METAANALYSIS
# > tail(new_estradiol)
# variant effect_allele other_allele effect_allele_frequency freqse    beta standard_error p_value direction sample_size
# 7783165    <NA>           CAG            C                  0.9883  2e-04 -0.0446         0.0496  0.3685        -+      311675
# 7783166    <NA>           CAG            C                  0.9883  2e-04 -0.0446         0.0496  0.3685        -+      311675
# 7783167    <NA>           CAG            C                  0.9883  2e-04 -0.0446         0.0496  0.3685        -+      311675
# 7783168    <NA>           CAG            C                  0.9883  2e-04 -0.0446         0.0496  0.3685        -+      311675
# 7783169    <NA>           CAG            C                  0.9883  2e-04 -0.0446         0.0496  0.3685        -+      311675
# 7783170    <NA>           CAG            C                  0.9883  2e-04 -0.0446         0.0496  0.3685        -+      311675
# chromosome base_pair_location        chr_pos
# 7783165          2          194803889 chr2:194803889
# 7783166         22           35986067 chr22:35986067
# 7783167          7          155392635 chr7:155392635
# 7783168          8          137212751 chr8:137212751
# 7783169          8          139173628 chr8:139173628
# 7783170         18           69048914 chr18:69048914
# 
# > estradiol[which(is.na(estradiol$variant) == TRUE),]
# variant effect_allele other_allele effect_allele_frequency freqse    beta standard_error p_value direction sample_size
# 2152957    <NA>           CAG            C                  0.9883  2e-04 -0.0446         0.0496  0.3685        -+      311675

# I have realized there are X chromosomes... 
str(new_estradiol)
# 'data.frame':	7783170 obs. of  13 variables:
#   $ variant                : chr  "rs10" "rs1000000" "rs10000003" "rs10000005" ...
# $ effect_allele          : chr  "A" "A" "A" "A" ...
# $ other_allele           : chr  "C" "G" "G" "G" ...
# $ effect_allele_frequency: num  0.0594 0.2271 0.2988 0.5316 0.5168 ...
# $ freqse                 : num  1e-04 1e-04 9e-04 1e-04 9e-04 5e-04 9e-04 5e-04 2e-04 7e-04 ...
# $ beta                   : num  -0.0075 0.0008 0.021 0.0125 -0.0126 -0.0004 0.0069 -0.0018 0.0064 -0.0061 ...
# $ standard_error         : num  0.0191 0.0105 0.0096 0.0088 0.0088 0.0127 0.0107 0.0092 0.0194 0.0106 ...
# $ p_value                : num  0.693 0.942 0.0283 0.1543 0.1519 ...
# $ direction              : chr  "--" "+-" "++" "++" ...
# $ sample_size            : num  311675 311675 311675 311675 311675 ...
# $ chromosome             : chr  "7" "12" "4" "4" ...
# $ base_pair_location     : num  9.24e+07 1.27e+08 5.76e+07 8.52e+07 2.16e+07 ...
# $ chr_pos                : chr  "chr7:92383888" "chr12:126890980" "chr4:57561647" "chr4:85161558" ...

noX_subsetted_df = new_estradiol # 7783170
df <- noX_subsetted_df[grepl("^[0-9]+$", noX_subsetted_df$chromosome), ] # 7564117
fwrite(df, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/curated_estradiol_chrpos_noX.txt")

### CHECKING THE BUILD...
df_ordered <- df[order(df$p_value), ]
df_ordered
# variant effect_allele other_allele effect_allele_frequency freqse    beta standard_error   p_value direction sample_size
# 5787271  rs727479             A            C                  0.6491  7e-04  0.1250         0.0094 9.654e-41        ++      311675
# 3346220 rs2414095             A            G                  0.3500  7e-04 -0.1246         0.0093 1.304e-40        --      311675
# 5674295 rs7173595             T            C                  0.6496  7e-04  0.1247         0.0094 1.348e-40        ++      311675
# 5675172 rs7175531             T            C                  0.3503  7e-04 -0.1247         0.0094 1.364e-40        --      311675
# 3346222 rs2414097             A            G                  0.6495  7e-04  0.1245         0.0093 1.593e-40        ++      311675
# 3346223 rs2414098             T            C                  0.3860  8e-04 -0.1222         0.0092 1.861e-40        --      311675
# chromosome base_pair_location        chr_pos
# 5787271         15           51534547 chr15:51534547 --> 37
# 3346220         15           51524292 chr15:51524292 --> 37
# 5674295         15           51533736 chr15:51533736 --> 37
# 5675172         15           51534055 chr15:51534055
# 3346222         15           51529835 chr15:51529835
# 3346223         15           51537806 chr15:51537806
