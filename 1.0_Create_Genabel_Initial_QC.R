#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   File to create working genetic files                        #
#   Susan Johnston Nov 2016                                     #
#                                                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# For CRIMAP, ID's must be coded as numbers.
# In conversion from bed format to ped/map format, ID's have been recoded
# in the pedigree and new PLINK files.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set up working environment and load in data         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Load libraries

library(ggplot2)
library(plyr)
library(GenABEL)
library(crimaptools)
library(beepr)
source("r/recoderFunc.R")

#~~ read in data and format

pedigree   <- read.table("data/Pedigree_16-05-02.txt", header = T, stringsAsFactors = F)
plink.file <- "data/Deer31"
out.prefix <- "data/Deer31"

firstRun <- TRUE

pedigree <- format_pedigree(pedigree)

#~~ Read in phenotypic data

basedata <- read.table("data/20160725_Deer_BasicData.txt", header = T, sep = "\t")
head(basedata)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Recode IDs for CRIMAP                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(firstRun == TRUE){
  
  #~~ convert PLINK to normal ped/map format and remove duplicated IDs
  
  system("cmd", input = paste0("plink --bfile ", plink.file, " --remove data/Deer31duplicateIDs.txt --recode --cow --out ", plink.file))
  
  #~~ extract information from the PLINK file
  
  system("cmd", input = paste0("gawk < ", plink.file, ".ped \"{print $1,$2,$3,$4,$5,$6}\" > plink.temp"))
  famfile <- read.table("plink.temp", stringsAsFactors = F)
  
  #~~ Create a new vector of IDs
  
  new.ids <- data.frame(oldID = c(pedigree[,1], pedigree[,2], pedigree[,3], famfile$V2), stringsAsFactors = F)
  new.ids <- unique(new.ids)
  new.ids <- subset(new.ids, oldID != "0")
  new.ids$newID <- 1:nrow(new.ids)
  head(new.ids)
  
  #~~ Create file for PLINK to change
  
  famfile.new <- famfile[,1:2]
  famfile.new$new.fam <- famfile.new[,1]
  famfile.new$new.id  <- recoderFunc(famfile.new$V2, new.ids$oldID, new.ids$newID)
  
  head(famfile.new)
  
  write.table(famfile.new, "data/plinkrecode.txt", col.names = F, row.names = F, quote = F)
  
  rm(famfile, famfile.new)
  
  system("cmd", input = paste0("plink --file ", plink.file, " --update-ids data/plinkrecode.txt --recode  --cow --out ", plink.file, ".recoded"))
  
  #~~ recode pedigree
  
  pedigree.new <- pedigree
  for(i in 1:3) pedigree.new[,i] <- recoderFunc(pedigree.new[,i], new.ids$oldID, new.ids$newID)
  
  head(pedigree.new)
  
  write.table(pedigree.new, "data/Pedigree_16-05-02.recoded.txt", row.names = F, quote = F)
  write.table(new.ids, "data/DeerRecodedIDs.txt", row.names = F, sep = "\t", quote = F)
  
  rm(pedigree)
}
 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Recode IDs for CRIMAP                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
if(firstRun == TRUE){    

  pedigree.new <- read.table("data/Pedigree_16-05-02.recoded.txt", header = T)
  new.ids <- read.table("data/DeerRecodedIDs.txt", header = T)
  
  new.ids$Code <- new.ids$oldID
  basedata <- join(basedata, new.ids)
  
  system("cmd", input = paste0("gawk < ", plink.file, ".recoded.ped \"{print $2,$5}\" > plink.temp"))
  famfile <- read.table("plink.temp", stringsAsFactors = F)
  
  #~~ Check for duplicate IDs
  
  dupids <- data.frame(table(famfile$V1))
  head(dupids)
  dupids <- subset(dupids, Freq > 1)

  #~~ Define sex from PLINK where males = 1, females = 2 to genabel males = 1,
  #    females = 0. Keep unknowns as female at present.
  
  names(famfile)[1:2] <- c("id", "plinksex")   
  famfile$sex <- famfile$plinksex
  famfile$sex[which(famfile$sex == 2)] <- 0
  
  pedfile <- subset(basedata, select = c(newID, Sex))
  pedfile <- subset(pedfile, !is.na(newID))
  pedfile$Sex[which(pedfile$Sex == 1)] <- "F"
  pedfile$Sex[which(pedfile$Sex == 2)] <- "M"
  names(pedfile) <- c("id", "datasex")
  
  famfile <- join(famfile, pedfile)
  
  head(famfile)
  
  table(famfile$plinksex)
  table(famfile$sex)
  table(famfile$datasex)
  
  subset(famfile, id == 2510)
  
  #~~ create pheno file for GenABEL   
  
  write.table(famfile, paste0(out.prefix, ".pheno"), sep = "\t", row.names = F, quote = F)
  
  #~~ Create map file by extracting first 4 columns of tped file
  
  mapfile <- read.table(paste0(plink.file,".recoded.map"))
  names(mapfile)<- c("Chromosome", "SNP.Name", "MapDist", "Position")
  mapfile <- mapfile[,-3]
  mapfile$Chromosome[which(mapfile$Chromosome == 30)] <- "X"
  
  #~~ Create GenABEL File
  
  write.table(mapfile, paste0(out.prefix, ".genabelmap"), col.names  = F, row.names = F, quote = F)
  
  convert.snp.ped(pedfile = paste0(plink.file, ".recoded.ped"), 
                  mapfile = paste0(out.prefix, ".genabelmap"),
                  outfile = paste0(out.prefix, ".gen"),
                  strand = "u", bcast = 1000, traits = 1, mapHasHeaderLine = F)   #added traits = 1
  
  rm(famfile, pedigree.new)
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Open GenABEL files and check sex                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

abeldata <- load.gwaa.data(phe = paste0(out.prefix, ".pheno"), 
                           gen = paste0(out.prefix, ".gen"))

#~~ Check sex

abelx <- abeldata[,chromosome(abeldata) == "X"]

table(chromosome(abeldata))

raw.test <- summary.snp.data(gtdata(abeldata))

ggplot(raw.test, aes(CallRate)) + geom_histogram(binwidth = 0.001)
ggplot(raw.test, aes(Q.2))      + geom_histogram(binwidth = 0.005)

raw.id   <- perid.summary(abeldata)
raw.id$id <- row.names(raw.id)
raw.id <- join(raw.id, phdata(abeldata))

ggplot(raw.id, aes(F)) + geom_histogram(binwidth = 0.01)
ggplot(raw.id, aes(CallPP)) + geom_histogram(binwidth = 0.01)
ggplot(raw.id, aes(Het)) + geom_histogram(binwidth = 0.01)

raw.id$Ffull <- raw.id$F

qc.initial <- check.marker(abeldata, callrate = 0.99, perid.call = 0.99, ibs.mrk = 0, maf = 0.01)
qc.initial$Xmrkfail

raw.id.x   <- perid.summary(abelx[, snpnames(abelx)[which(!snpnames(abelx) %in% qc.initial$Xmrkfail)]])
raw.id.x$id <- row.names(raw.id.x)
raw.id.x <- join(raw.id.x, phdata(abelx))
raw.id.x <- join(raw.id.x, raw.id[,c("id", "Ffull")])

raw.id.x <- arrange(raw.id.x, F)
raw.id.x$F.Order <- 1:nrow(raw.id.x)
ggplot(raw.id.x, aes(F.Order, F, col = factor(datasex))) + geom_point()
ggplot(raw.id.x, aes(F.Order, F, col = factor(sex))) + geom_point()
ggplot(raw.id.x, aes(F.Order, F, col = factor(plinksex))) + geom_point()

ggplot(raw.id.x, aes(Ffull, F, col = factor(datasex))) + geom_point() + geom_hline(yintercept = 0.90)
ggplot(raw.id.x, aes(Ffull, F, col = factor(sex))) + geom_point() + geom_hline(yintercept = 0.90)
ggplot(raw.id.x, aes(Ffull, F, col = factor(plinksex))) + geom_point() + geom_hline(yintercept = 0.90)


newphenofile <- phdata(abeldata)
newphenofile$sex <- newphenofile$datasex
newphenofile$sex[which(newphenofile$sex == "M")] <- 1
newphenofile$sex[which(newphenofile$sex == "F")] <- 0
newphenofile$sex[which(newphenofile$sex == 3)] <- 0
newphenofile$sex[which(is.na(newphenofile$sex))] <- 0


newphenofile$sex[newphenofile$id %in% raw.id.x$id[which(raw.id.x$F < 0.9)]] <- 0

subset(newphenofile, id == 2510)

not.males   <- subset(newphenofile, sex == 1 & id %in% pedigree.new$MOTHER)$id
not.females <- subset(newphenofile, sex == 0 & id %in% pedigree.new$FATHER)$id


if(length(not.males) > 0) 
  newphenofile$sex[which(newphenofile$id %in% not.males)] <- 0
if(length(not.females) > 0) 
  newphenofile$sex[which(newphenofile$id %in% not.females)] <- 1


write.table(newphenofile, paste0(out.prefix, ".pheno"), row.names = F, sep = "\t", quote = F)

abeldata <- load.gwaa.data(phe = paste0(out.prefix, ".pheno"), 
                           gen = paste0(out.prefix, ".gen"))

qc.initial <- check.marker(abeldata, callrate = 0.95, perid.call = 0.99, ibs.mrk = 0, maf = 0.01)

writeLines(qc.initial$Xmrkfail, "data/Pseudoautosomal_SNPs.txt")

#~~ redo the map file

mapfile <- read.table(paste0(plink.file,".recoded.map"))
names(mapfile)<- c("Chromosome", "SNP.Name", "MapDist", "Position")
mapfile <- mapfile[,-3]
mapfile$Chromosome[which(mapfile$Chromosome == 30)] <- "X"
mapfile$Chromosome[which(mapfile$SNP.Name %in% qc.initial$Xmrkfail)] <- 30

#~~ Create GenABEL File

write.table(mapfile, paste0(out.prefix, ".genabelmap"), col.names  = F, row.names = F, quote = F)

convert.snp.ped(pedfile = paste0(plink.file, ".recoded.ped"), 
                mapfile = paste0(out.prefix, ".genabelmap"),
                outfile = paste0(out.prefix, ".gen"),
                strand = "u", bcast = 1000, traits = 1, mapHasHeaderLine = F)   #added traits = 1


abeldata <- load.gwaa.data(phe = paste0(out.prefix, ".pheno"), 
                           gen = paste0(out.prefix, ".gen"))

qc.initial <- check.marker(abeldata, callrate = 0.95, perid.call = 0.99, ibs.mrk = 0, maf = 0.01, fdrate = 0.1)
qc.initial$isfemale %in% qc.initial$idok

abeldata <- Xfix(abeldata)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Conduct quality control on GenABEL files            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# We will use a perid.call of 0.95 for the time-being, and a minor allele frequency of 0.01, and an ibs.threshold 
# of 0.9. THe following is the output from R, showing the rounds of QC followed until all arguments are fulfilled:

qc.abel <- check.marker(abeldata, callrate = 0.95, perid.call = 0.98, maf = 0.01, ibs.threshold = 0.9, fdrate = 0.1)


#~~ get and check QC'ed data

abeldata <- abeldata[qc.abel$idok, qc.abel$snpok]


nids(abeldata)
nsnps(abeldata)

save(qc.abel, file = "data/Deer31_QC.check_marker.Rdata")

save(abeldata, file = "data/Deer31_QC.RData")

save.gwaa.data(abeldata, phenofile = "data/Deer31_QC.phe", genofile = "data/Deer31_QC.gen")

