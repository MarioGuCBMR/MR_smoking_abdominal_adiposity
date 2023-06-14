library(vroom)
library(data.table)

oestradiolB37 <- as.data.frame(vroom("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/estradiol_B37.txt", delim = "/t"))
colnames(oestradiolB37)
# > head(oestradiolB37)
# group group_name seqnames     start       end width strand    variant chromosome base_pair_location effect_allele other_allele
# 1     1         NA     chr6 131161231 131161231     1      *  rs2326918          6          130840091             A            G
# 2     2         NA     chr6 132591224 132591224     1      * rs56252526          6          132270085             A            G
# 3     3         NA     chr7  34488377  34488377     1      *   rs977590          7           34448765             A            G
# 4     4         NA     chr3 176384537 176384537     1      * rs66941928          3          176666749             T            C
# 5     5         NA    chr14  26182670  26182670     1      * rs17278013         14           25713464             T            G
# 6     6         NA     chr3  31528328  31528328     1      * rs62234673          3           31486836             T            C
# effect_allele_frequency    beta standard_error p_value sample_size        chr_pos old_hg19_POS new_build37_bp
# 1                  0.8459 -0.0166         0.0122  0.1728      311675 chr6:130840091    130840091      131161231
# 2                  0.9882  0.0466         0.0425  0.2724      311675 chr6:132270085    132270085      132591224
# 3                  0.9361  0.0067         0.0181  0.7118      311675  chr7:34448765     34448765       34488377
# 4                  0.8025 -0.0018         0.0110  0.8729      311675 chr3:176666749    176666749      176384537
# 5                  0.7644 -0.0128         0.0103  0.2153      311675 chr14:25713464     25713464       26182670
# 6                  0.0885  0.0116         0.0154  0.4508      311675  chr3:31486836     31486836       31528328

oestradiol_sel <- oestradiolB37[, c("variant", "chromosome", "new_build37_bp", "effect_allele", "other_allele", 
                                 "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", 
                                 "chr_pos")]

colnames(oestradiol_sel) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", 
                              "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", 
                              "chr_pos")

oestradiol_sel$chr_pos <- paste("chr", oestradiol_sel$chromosome, ":", oestradiol_sel$base_pair_location, sep = "") # add chr:pos column

fwrite(oestradiol_sel, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/final_estradiol_B37.txt")
