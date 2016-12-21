
library(GenABEL)
library(plyr)

AnalysisSuffix <- "a"
AnalysisSuffix2 <- "a_cel"
AnalysisSuffix3 <- "a_skel"

out.prefix <- "data/Deer31_QC"

mapdata <- read.table("results/5_Linkage_Map_Positions_CEL_run4_a.txt", header = T)

load("data/check_marker.Rdata")

abeldata <- load.gwaa.data(phenofile = paste0(out.prefix, ".pheno"),
                           genofile  = paste0(out.prefix, ".gen"))

abeldata <- abeldata[qc.abel$idok, qc.abel$snpok]

rm(qc.abel)


#~~ Create the PLINK file input

nsnps(abeldata)
nrow(mapdata)


write.table(snp.names(abeldata), "data/Deer31_QC.recoded.v2.snplist", row.names = F, col.names = F, quote = F)

system("cmd", input = "plink --file data/Deer31_QC.recoded --extract data/Deer31_QC.recoded.v2.snplist --recode  --cow --out data/Deer31_QC.recoded.v2")

mapfile <- read.table("data/Deer31_QC.recoded.v2.map")
head(mapfile)

mapdata2 <- subset(mapdata, select = c(SNP.Name, cMPosition.run4, CEL.LG))
names(mapdata2)[1] <- c("V2")

mapfile <- join(mapfile, mapdata2)

mapfile <- mapfile[,c("CEL.LG", "V2", "cMPosition.run4")]
mapfile$Sequence <- c(0, rep(NA, times = nrow(mapfile) - 1))
for(i in 2:nrow(mapfile)) mapfile$Sequence[i] <- ifelse(mapfile$cMPosition.run4[i] == mapfile$cMPosition.run4[i-1], mapfile$Sequence[i-1] + 1, 0)

head(mapfile)

mapfile$Mb <- (mapfile$cMPosition.run4 * 1e6) + (mapfile$Sequence * 100)
mapfile <- subset(mapfile, select = -Sequence)

write.table(mapfile, "data/Deer31_QC.recoded.v2.map", row.names = F, col.names = F, quote = F)

#~~ Create BED file

system("cmd", input = "plink --file data/Deer31_QC.recoded.v2 --recode --make-bed --out gcta/Deer31_QC.recoded.v2")

system("cmd", input = "plink --bfile data/Deer31_QC.recoded.v2 --recode --out gcta/TEST")

#~~ Create GRM snplists

for(i in 1:33){
  snplist <- subset(mapdata, CEL.LG == i)
  write.table(snplist$SNP.Name, paste0("gcta/chr", i, "snplist.txt"), row.names = F, col.names = F, quote = F)
}

write.table(mapdata$SNP.Name, "gcta/chrALLsnplist.txt", row.names = F, col.names = F, quote = F)





