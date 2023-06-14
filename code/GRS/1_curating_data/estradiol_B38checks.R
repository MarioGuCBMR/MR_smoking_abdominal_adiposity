
oestradiol <- as.data.frame(vroom("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/metaanalysis_data/METAANALYSIS_estradiol2_1.tbl"))
head(oestradiol)

df_ordered <- oestradiol[order(oestradiol$`P-value`), ]
head(df_ordered)
# MarkerName Allele1 Allele2  Freq1 FreqSE  Effect StdErr   P-value Direction
# 4245379   rs727479       a       c 0.6491  7e-04  0.1250 0.0094 9.654e-41        ++     
#   3911704  rs2414095       a       g 0.3500  7e-04 -0.1246 0.0093 1.304e-40        --
#   1737044  rs7173595       t       c 0.6496  7e-04  0.1247 0.0094 1.348e-40        ++
#   1812268  rs7175531       t       c 0.3503  7e-04 -0.1247 0.0094 1.364e-40        --
#   3193709  rs2414097       a       g 0.6495  7e-04  0.1245 0.0093 1.593e-40        ++
#   468448   rs2414098       t       c 0.3860  8e-04 -0.1222 0.0092 1.861e-40        --


oestradiol_male <- as.data.frame(vroom("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/estradiol_male.txt"))
# > head(oestradiol_male)
# variant chromosome base_pair_location effect_allele other_allele effect_allele_frequency        beta standard_error  p_value sample_size     chr_pos
# 1 rs556128788          1              99388             T            G               0.0125875 -0.05114699      0.0844242 0.544623      147690  chr1:99388
# 2 rs146209971          1              99671             T            A               0.0104829  0.07647018      0.1052350 0.467436      147690  chr1:99671
# 3   rs9726668          1             108230             T            C               0.0118564  0.03660190      0.0877148 0.676508      147690 chr1:108230
# 4 rs576317820          1             324822             T            A               0.0133099  0.09787054      0.0873551 0.262565      147690 chr1:324822 --> GRCh37
# 5  rs41485244          1             565343             G            A               0.0125681  0.02808198      0.0824964 0.733525      147690 chr1:565343 --> GRCh37
# 6   rs9285835          1             569004             C            T               0.0157125  0.10050575      0.0785332 0.200622      147690 chr1:569004 --> GRCh37
df_ordered_male <- oestradiol_male[order(oestradiol_male$p_value), ]
head(df_ordered_male)
# variant chromosome base_pair_location effect_allele other_allele effect_allele_frequency       beta standard_error     p_value sample_size        chr_pos
# 2325824 rs28892005         15           51519945             A         AAAG                0.350362 -0.2289491      0.0139272 1.00598e-60      147690 chr15:51519945 --> NONE
# 2325859  rs7173595         15           51533736             C            T                0.351149 -0.2287404      0.0139251 1.23745e-60      147690 chr15:51533736 --> GRCh37
# 2325836 rs12592697         15           51525173             T            C                0.351046 -0.2284614      0.0139092 1.26097e-60      147690 chr15:51525173 --> GRCh37
# 2325860  rs7175531         15           51534055             T            C                0.351077 -0.2287065      0.0139263 1.31677e-60      147690 chr15:51534055 --> GRCh37
# 2325861   rs727479         15           51534547             C            A                0.351669 -0.2286839      0.0139285 1.41295e-60      147690 chr15:51534547 --> GRCh37
# 2325846  rs2414097         15           51529835             G            A                0.351335 -0.2282842      0.0139060 1.46381e-60      147690 chr15:51529835 --> GRCh37
oestradiol_female <- as.data.frame(vroom("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/estradiol_female.txt"))
df_ordered_female <- oestradiol_female[order(oestradiol_female$p_value), ]
head(df_ordered_female)
# variant chromosome base_pair_location effect_allele other_allele effect_allele_frequency       beta standard_error     p_value sample_size       chr_pos
# 6660628  rs45446698          7           99332948             G            T               0.0427387 -0.2064413      0.0301572 7.61947e-12      163985 chr7:99332948 --> GRCh37
# 6660393 rs148982377          7           99075038             C            T               0.0413413 -0.1931213      0.0308472 3.83601e-10      163985 chr7:99075038 --> GRCh37
# 6660448  rs34670419          7           99130834             T            G               0.0398637 -0.1889088      0.0312199 1.44019e-09      163985 chr7:99130834 --> GRCh37
# 6660348 rs574637649          7           99025328             T            C               0.0390673 -0.1891951      0.0316999 2.39705e-09      163985 chr7:99025328 --> GRCh37
# 6660630  rs45575938          7           99332998             G            A               0.0287107 -0.2254061      0.0389497 7.16037e-09      163985 chr7:99332998 --> GRCh37
# 6660629  rs45467892          7           99332997             A            T               0.0286697 -0.2252696      0.0389482 7.30148e-09      163985 chr7:99332997 --> GRCh37

