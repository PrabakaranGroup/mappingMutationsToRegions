
altorfs <- read.csv("/Users/naru/Downloads/altorf_database.txt",sep = "\t",header = T)

head(altorfs)


codingCosmic <- read.csv("/Users/naru/Downloads/CosmicCodingMuts.vcf-14-12-17-SP", sep = "\t", header = T,skip = 13)

codingCosmic$X.CHROM <- paste("chr", codingCosmic$X.CHROM, sep="")

head(codingCosmic)

library(dplyr)
x = altorfs %>%
  left_join(codingCosmic, by = c("chr_no." = "X.CHROM")) %>%
  filter( (start_coordinate <= POS & POS <= end_coordinate)
          