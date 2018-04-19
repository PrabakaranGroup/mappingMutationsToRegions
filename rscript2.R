
altorfs <- read.csv("/Users/naru/Downloads/Pseudo_Coord_NonTranscribed.csv",sep = ",",header = T)

headaltorfs <-altorfs

head(headaltorfs )
hgmd <- read.csv("/Users/naru/Downloads/HGMD_PRO_2016.4_hg38.vcf", sep = "\t", header = T,skip = 13)

head(hgmd)
#don't click on it twice 
#hgmd$X.CHROM <- paste("chr", hgmd$X.CHROM, sep="")

headcodingCosmic <- hgmd
combined <- sort(union(levels(headcodingCosmic$X.CHROM), levels(headaltorfs$ChrN)))
head(headcodingCosmic)
library(dplyr)
x = mutate(headaltorfs,ChrN=factor(ChrN, levels=combined)) %>%
  inner_join(mutate(headcodingCosmic,X.CHROM=factor(X.CHROM, levels=combined)), by = c("ChrN" = "X.CHROM")) %>%
  filter( Start <= POS & POS <= Stop)


write.csv(x,file = "/Users/naru/Downloads/dPseudo_Coord_TranslatedL_mapping_of_HGMD.csv")
