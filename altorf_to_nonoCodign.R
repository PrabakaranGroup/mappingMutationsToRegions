





altorfs <- read.csv("altorfs/HS_genbank_GRCh38_v1.2_altorf_filtered.tsv",sep = "\t",header = T)


headaltorfs <-altorfs


head(headaltorfs )
length(unique(headaltorfs$Chromosome))


#cleaning chromosome columne 
library(plyr)
library(dplyr)

library(tidyr)
#update.packages()
headaltorfs$Chromosome <- gsub("\\.", ',', headaltorfs$Chromosome)

newData <- headaltorfs %>% 
  mutate(Chromosome = strsplit(as.character(Chromosome), ",")) %>% 
  unnest(Chromosome)

newData$Chromosome <- gsub("23", 'X',newData$Chromosome)

unique(newData$Chromosome)

#write.csv(newData, file="/Users/naru/Downloads/human_altorf_database.csv")
#spliting 
colnames(newData)[8] <- "X.CHROM"
#newData
#newData[, 'X.CHROM'] <- as.factor(newData[, 'X.CHROM'])
#altOrfsPart <- split( newData  , f = newData$X.CHROM )
#head(newData)

head(newData)





#files <- list.files(path="/home/ngs/narumeena/Documents/arkarplication/CosmicNoCodingDataSetPart", pattern="*", full.names=T, recursive=FALSE)
#codingCosmic <- read.csv("CosmicNoCodingDataSetPart/CosmicNonCodingVariants_aa", sep = "\t", header = T,skip = 13)
#colnames(codingCosmic)[1] <- "X.CHROM"
#newData[, 'X.CHROM'] <- as.character(newData[, 'X.CHROM'])
#codingCosmic[, 'X.CHROM'] <- as.character(codingCosmic[, 'X.CHROM'])
#result <- inner_join(newData,codingCosmic, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)
#write.csv(result,file="altorfs/noncoding/altorf_mapped_by_cosmic_noncoding_mutations_part_aa.csv"))
files <- list.files(path="/home/ngs/narumeena/Documents/arkarplication/CosmicNoCodingDataSetPart", pattern="*", full.names=T, recursive=FALSE)
for(i in 1:length(files)){

codingCosmic <- read.csv(files[i], sep = "\t", header = T,skip = 13)
colnames(codingCosmic)[1] <- "X.CHROM"
newData[, 'X.CHROM'] <- as.character(newData[, 'X.CHROM'])
codingCosmic[, 'X.CHROM'] <- as.character(codingCosmic[, 'X.CHROM'])
result <- inner_join(newData,codingCosmic, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)
write.csv(result,file=paste0("altorfs/noncoding/altorf_mapped_by_cosmic_noncoding_mutations_part_",i,".csv"))
  #  t <- read.table(x, header=T) # load file
    # apply function
 #   out <- function(t)
    # write to file
#    write.table(out, "path/to/output", sep="\t", quote=F, row.names=F, col.names=T)
}




#codingCosmic <- read.csv("CosmicNonCodingVariantsab_abaa", sep = "\t", header = T,skip = 13)
#colnames(codingCosmic)[1] <- "X.CHROM"
#head(codingCosmic)
#mutate(altorfs, X.CHROM=factor(X.CHROM, levels=combined))
#mutate(codingCosmic, X.CHROM=factor(X.CHROM, levels=combined))
#newData[, 'X.CHROM'] <- as.character(newData[, 'X.CHROM'])
#codingCosmic[, 'X.CHROM'] <- as.character(codingCosmic[, 'X.CHROM'])
#head(codingCosmic)
#combined <- sort(union(levels(altorfs$X.CHROM), levels(codingCosmic$X.CHROM)))
#result <- inner_join(newData,codingCosmic, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)
#write.csv(result,file="altorfs/noncoding/altorf_mapped_by_cosmic_noncoding_mutations_part_ab_ab_aa.csv")
