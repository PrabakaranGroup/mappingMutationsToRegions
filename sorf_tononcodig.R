# sorfs coding 
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
#codingCosmic <- read.csv("CosmicNonCodingVariantsad_new", sep = "\t", header = T,skip = 13)

altorfs <- read.csv("human_PLsorf_database.txt",sep = "\t",header = T)
#codingCosmic$X.CHROM <- paste("chr", codingCosmic$X.CHROM, sep="")
colnames(altorfs)[4] <- "X.CHROM"
#combined <- sort(union(levels(altorfs$X.CHROM), levels(codingCosmic$X.CHROM)))

#mutate(altorfs, X.CHROM=factor(X.CHROM, levels=combined)

#result <- inner_join(mutate(altorfs, X.CHROM=factor(X.CHROM, levels=combined)),mutate(codingCosmic, X.CHROM=factor(X.CHROM, levels=combined)), by = "X.CHROM") %>% filter(   Start_position <= POS & POS <= End_position)

#partitions 

#cosmicNonCodingPart <- split( codingCosmic  , f = codingCosmic$X.CHROM )

#altOrfsPart <- split( altorfs , f = altorfs$X.CHROM )
newData <- altorfs
#newData[, 'X.CHROM'] <- as.character(newData[, 'X.CHROM'])
#codingCosmic[, 'X.CHROM'] <- as.character(codingCosmic[, 'X.CHROM'])
#combined <- sort(union(levels(altorfs$X.CHROM), levels(codingCosmic$X.CHROM)))
#result <- inner_join(newData,codingCosmic, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

files <- list.files(path="/home/ngs/narumeena/Documents/arkarplication/cosmicPart", pattern="*", full.names=T, recursive=FALSE)
for(i in 1:length(files)){

codingCosmic <- read.csv(files[i], sep = "\t", header = T,skip = 13)
codingCosmic$X.CHROM <- paste("chr", codingCosmic$X.CHROM, sep="")
colnames(codingCosmic)[1] <- "X.CHROM"
newData[, 'X.CHROM'] <- as.character(newData[, 'X.CHROM'])
codingCosmic[, 'X.CHROM'] <- as.character(codingCosmic[, 'X.CHROM'])
result <- inner_join(newData,codingCosmic, by = "X.CHROM") %>% filter(   Start_position <= POS & POS <= End_position )
write.csv(result,file=paste0("sorfs/noncoding/sorfs_mapped_by_cosmic_noncoding_mutations_part_",i,".csv"))
  #  t <- read.table(x, header=T) # load file
    # apply function
 #   out <- function(t)
    # write to file
#    write.table(out, "path/to/output", sep="\t", quote=F, row.names=F, col.names=T)
}

#write.csv(result,file="sorfs/noncoding/sorf_non_coding_part_ad_inner_join.csv")


