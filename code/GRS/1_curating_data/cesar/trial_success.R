
#INSTALL "rtracklayer" in BiocManager
# install.packages("BiocManager")
#BiocManager::install("rtracklayer")

#LOAD "rtracklayer"
library(BiocManager)
library("rtracklayer")
library(vroom)

#ESTABLISH FILE PATH TO THE GWAS DATA TO BE CONVERTED
#file_path <- "C:/Users/msd830/Desktop/Triangulation/EAGWAS2022/EA4_additive_excl_23andMe.txt"
#file_path <- "C:/Users/msd830/Desktop/Triangulation/BMIGWAS2018/Meta-analysis_Locke_et_al_UKBiobank_2018_UPDATED.txt"
file_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/curated_estradiol.txt"

#READ THE FILE
snpsFile <- as.data.frame(vroom(file_path))
head(snpsFile)

#ADD 'CHR' TO THE NUMBER FOR IT TO WORK CORRECTLY WITH 'rtracklayer'
snpsFile$newChromosomeColumn = paste0("chr", snpsFile$chromosome)
head(snpsFile)

#CREATE A COLUMNS WITH THE OLD BP COLUMN NAME 
snpsFile$old_hg19_POS = snpsFile$base_pair_location
head(snpsFile)

snpsFile$startField = snpsFile$base_pair_location
snpsFile$endField = snpsFile$base_pair_location

head(snpsFile)

#IMPORT THE CHAIN FILE (TO DO THE CONVERSION)

chainobject <- import.chain("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/hg38ToHg19.over.chain") 

#CREATE A GENOMIC RANGES OBJECT (FROM DATA.FRAME TO GRANGES)

grGWAS_SNPs <- makeGRangesFromDataFrame(
  snpsFile,
  seqnames.field = "newChromosomeColumn",
  start.field = "startField",
  end.field = "endField",
  keep.extra.columns = TRUE
  
)

#RUN LIFTOVER
results <- as.data.frame(liftOver(grGWAS_SNPs, chainobject))
head(results)

results$new_build37_bp = results$start
head(results)

write.table(results, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/curated_summary_statistics/estradiol_B37.txt", sep = "/t", row.names = FALSE)
