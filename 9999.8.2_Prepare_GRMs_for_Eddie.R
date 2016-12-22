
library(GenABEL)
library(plyr)

AnalysisSuffix <- "a"
AnalysisSuffix2 <- "a_cel"
AnalysisSuffix3 <- "a_skel"

out.prefix <- "data/Deer31_QC"

mapdata <- read.table("results/6_Linkage_Map_Positions_CEL_run5_a.txt", header = T)

#~~ Read & format genetic data

load("data/Deer31_QC.RData", verbose = T)

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")


#~~ Create the PLINK file input

write.table(mapdata$SNP.Name, "data/Deer31.recoded.v2.snplist", row.names = F, col.names = F, quote = F)

system("cmd", input = "plink --file data/Deer31.recoded --extract data/Deer31.recoded.v2.snplist --recode  --cow --out data/Deer31.recoded.v2")

mapfile <- read.table("data/Deer31.recoded.v2.map")
head(mapfile)

mapdata2 <- subset(mapdata, select = c(SNP.Name, cMPosition.run5, CEL.LG))
names(mapdata2)[1] <- c("V2")

mapfile <- join(mapfile, mapdata2)

mapfile <- mapfile[,c("CEL.LG", "V2", "cMPosition.run5")]
mapfile$Sequence <- c(0, rep(NA, times = nrow(mapfile) - 1))
for(i in 2:nrow(mapfile)) mapfile$Sequence[i] <- ifelse(mapfile$cMPosition.run5[i] == mapfile$cMPosition.run5[i-1], mapfile$Sequence[i-1] + 1, 0)

head(mapfile)

mapfile$Mb <- (mapfile$cMPosition.run5 * 1e6) + (mapfile$Sequence * 100)
mapfile <- subset(mapfile, select = -Sequence)

write.table(mapfile, "data/Deer31.recoded.v2.map", row.names = F, col.names = F, quote = F)

#~~ Create BED file

system("cmd", input = "plink --file data/Deer31.recoded.v2 --recode --make-bed --out gcta/Deer31.recoded.v2")

#system("cmd", input = "plink --bfile data/Deer31.recoded.v2 --recode --out gcta/TEST")

#~~ Create GRM snplists

for(i in 1:34){
  snplist <- subset(mapdata, CEL.LG == i)
  write.table(snplist$SNP.Name, paste0("gcta/chr", i, "snplist.txt"), row.names = F, col.names = F, quote = F)
}

write.table(subset(mapdata, CEL.LG != 34)$SNP.Name, "gcta/chrALLsnplist.txt", row.names = F, col.names = F, quote = F)





