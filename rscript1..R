
altorfs <- read.csv("/Users/naru/Downloads/human_genbank_GRCh38/HS_genbank_GRCh38_v1.2_altorf_filtered.tsv",sep = "\t",header = T)
#altorfs <- read.csv("altorfs/HS_genbank_GRCh38_v1.2_altorf_filtered.tsv",sep = "\t",header = T)


headaltorfs <-altorfs

head(headaltorfs )
length(unique(headaltorfs$Chromosome))


#cleaning chromosome columne 
library(plyr)
library(dplyr)

library(tidyr)
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
newData[, 'X.CHROM'] <- as.factor(newData[, 'X.CHROM'])
altOrfsPart <- split( newData  , f = newData$X.CHROM )
head(newData)

head(newData)
#cosmic mutations 
codingCosmic <- read.csv("/Users/naru/Downloads/CosmicNonCodingVariants-14-12-17-SP.vcf", sep = "\t", header = T,skip = 13)
#codingCosmic <- read.csv("CosmicNonCodingVariants-14-12-17-SP.vcf", sep = "\t", header = T,skip = 13)



#don't click on it twice 
#codingCosmic$X.CHROM <- paste("chr", codingCosmic$X.CHROM, sep="")

headcodingCosmic <- codingCosmic
head(headcodingCosmic)

#split data frame based on numbers of rows 
chunk <- 10000
n <- nrow(headcodingCosmic)
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
d <- split(headcodingCosmic,r)

##loop 

result <- inner_join(newData,d[[1]], by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

# for (i in 2:18359){
#   result <-  rbind(result,inner_join(newData,d[[i]], by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate))
# 
# }

# loop in parrallel
library(foreach)
library(doParallel)
#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


finalMatrix <- foreach(i=2:1836, .combine=cbind) %dopar% {
 result <-  rbind(result,inner_join(newData,d[[i]], by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate))
 #calling a function
  #do other things if you want
  
 result #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
#stop cluster
stopCluster(cl)




#combine and inner group 
results <- mapply(function(model, df) {
   cbind(df[rep(1, 72), c("lat", "long")], resid(model))
   }, models, pieces)

###spliting cosmic data frame 

cosmicNonCodingPart <- split( headcodingCosmic  , f = headcodingCosmic$X.CHROM )

#

head(cosmicNonCodingPart$`1`)
head(altOrfsPart$`1`)
#table join 

dt1 <- data.table(altOrfsPart$`1`, key='X.CHROM')
dt2 <- data.table(cosmicNonCodingPart$`1`, key='X.CHROM')

## Do an INNER JOIN-like operation, where non-matching rows are removed
chr1 <- dt1[dt2, nomatch=0, by=.EACHI ]


#inner join 
chr1 <- inner_join(altOrfsPart$`1`,cosmicNonCodingPart$`1`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr2 <- inner_join(altOrfsPart$`2`,cosmicNonCodingPart$`2`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr3 <- inner_join(altOrfsPart$`3`,cosmicNonCodingPart$`3`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr4 <- inner_join(altOrfsPart$`4`,cosmicNonCodingPart$`4`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr5 <- inner_join(altOrfsPart$`5`,cosmicNonCodingPart$`5`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr6 <- inner_join(altOrfsPart$`6`,cosmicNonCodingPart$`6`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr7 <- inner_join(altOrfsPart$`7`,cosmicNonCodingPart$`7`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr8 <- inner_join(altOrfsPart$`8`,cosmicNonCodingPart$`8`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr9 <- inner_join(altOrfsPart$`9`,cosmicNonCodingPart$`9`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr10 <- inner_join(altOrfsPart$`10`,cosmicNonCodingPart$`10`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr11 <- inner_join(altOrfsPart$`11`,cosmicNonCodingPart$`11`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr12 <- inner_join(altOrfsPart$`12`,cosmicNonCodingPart$`12`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr13 <- inner_join(altOrfsPart$`13`,cosmicNonCodingPart$`13`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr14 <- inner_join(altOrfsPart$`14`,cosmicNonCodingPart$`14`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr15 <- inner_join(altOrfsPart$`15`,cosmicNonCodingPart$`15`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr16 <- inner_join(altOrfsPart$`16`,cosmicNonCodingPart$`16`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr17 <- inner_join(altOrfsPart$`17`,cosmicNonCodingPart$`17`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr18 <- inner_join(altOrfsPart$`18`,cosmicNonCodingPart$`18`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr19 <- inner_join(altOrfsPart$`19`,cosmicNonCodingPart$`19`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr20 <- inner_join(altOrfsPart$`20`,cosmicNonCodingPart$`20`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr21 <- inner_join(altOrfsPart$`21`,cosmicNonCodingPart$`21`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chr22 <- inner_join(altOrfsPart$`22`,cosmicNonCodingPart$`22`, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chrX <- inner_join(altOrfsPart$X,cosmicNonCodingPart$X, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chrY <- inner_join(altOrfsPart$Y,cosmicNonCodingPart$Y, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)

chrMT <- inner_join(altOrfsPart$MT,cosmicNonCodingPart$MT, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)


result <- bind_rows(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrMT)



#result <- ldply(cosmicNonCodingPart, function(x){ inner_join(newData,x, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)}) 

#write.csv(result, file = "altorf_mapped_by_noncoding_cosmic_mutation.csv")
# 
# unique(newData$Chromosome)
# combined <- sort(union(levels(headcodingCosmic$X.CHROM), levels(headaltorfs$Chromosome)))
# 
# 
# x = mutate(headaltorfs,Chromosome=factor(Chromosome, levels=combined)) %>%
#   inner_join(mutate(headcodingCosmic,X.CHROM=factor(X.CHROM, levels=combined)), by = c("Chromosome" = "X.CHROM")) %>%
#   filter( Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)
#           
# #tmp <- merge(data_A, data_B, by = c("USER_A","USER_B"), all.x = TRUE)
# write.csv(x,file="altorf_mapped_by_cosmic_noncoding_mutations.csv")
# 
# 
# colnames(headaltorfs)[5] <- "X.CHROM"
# 
# y <- inner_join(headaltorfs,headcodingCosmic, by = "X.CHROM") %>% filter(   Start.genomic.coordinate <= POS & POS <= End.genomic.coordinate)
