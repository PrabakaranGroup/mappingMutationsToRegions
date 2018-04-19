
sorfs <- read.csv("/Users/naru/Downloads/all_identified_6248_sorfs_info.hg38_OverlapedsORfsoftblastnAndliftOver.csv",sep = ",",header = F)

headaltorfs <-sorfs

head(headaltorfs )
codingCosmic <- read.csv("/Users/naru/Downloads/HGMD_PRO_2016.4_hg38.vcf", sep = "\t", header = T,skip = 13)

##don't click on it twice 
codingCosmic$X.CHROM <- paste("chr", codingCosmic$X.CHROM, sep="")

headcodingCosmic <- codingCosmic
#combined <- sort(union(levels(headcodingCosmic$X.CHROM), levels(headaltorfs$V10)))
head(headcodingCosmic)
library(dplyr)
#x = mutate(headaltorfs,V10=factor(V10, levels=combined)) %>%
#  inner_join(mutate(headcodingCosmic,X.CHROM=factor(X.CHROM, levels=combined)), by = c("V10" = "X.CHROM")) %>%
#  filter(  V16 <= POS & POS <= V17)

colnames(headaltorfs)[10] <- "X.CHROM"

y<- inner_join(headaltorfs,headcodingCosmic, by = "X.CHROM")  %>% filter(  V16 <= POS & POS <= V17)

write.csv(y,file="/Users/naru/Downloads/HGMD_mutations_mapped_to_all_identified_6248_sorfs_info.hg38_OverlapedsORfsoftblastnAndliftOver_inner_join.csv")
