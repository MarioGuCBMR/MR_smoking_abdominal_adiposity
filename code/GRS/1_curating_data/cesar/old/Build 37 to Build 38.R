
#INSTALL "rtracklayer" in BiocManager
#BiocManager::install("rtracklayer")

#LOAD "rtracklayer"
library("rtracklayer")

#ESTABLISH FILE PATH TO THE GWAS DATA TO BE CONVERTED
#file_path <- "C:/Users/msd830/Desktop/Triangulation/EAGWAS2022/EA4_additive_excl_23andMe.txt"
#file_path <- "C:/Users/msd830/Desktop/Triangulation/BMIGWAS2018/Meta-analysis_Locke_et_al_UKBiobank_2018_UPDATED.txt"
file_path <- "C:/Users/msd830/Desktop/Triangulation/T2DGWAS2022/DIAMANTE-EUR.sumstat.txt"

#READ THE FILE
snpsFile <- read.delim(file_path, header = TRUE, sep = " ")
head(snpsFile)

#ADD 'CHR' TO THE NUMBER FOR IT TO WORK CORRECTLY WITH 'rtracklayer'
snpsFile$newChromosomeColumn = paste0("chr", snpsFile$chromosome.b37.)
head(snpsFile)

#CREATE A COLUMNS WITH THE OLD BP COLUMN NAME 
snpsFile$old_hg19_POS = snpsFile$position.b37.
head(snpsFile)

snpsFile$startField = snpsFile$position.b37.
snpsFile$endField = snpsFile$position.b37.

head(snpsFile)

#IMPORT THE CHAIN FILE (TO DO THE CONVERSION)

chainobject <- import.chain("C://Users//msd830//Desktop//Original_Build_Conversion_Files//hg19ToHg38.over.chain") 

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
results$new_hg_38_bp = results$start
head(results)

write.table(results, "T2D_B38.txt", sep = "\t", row.names = FALSE)
