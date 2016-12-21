
library(plyr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Read in and format data                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

basedata <- read.table("data/20160725_Deer_BasicData.txt", header = T, sep = "\t", stringsAsFactors = F)

popsize <- read.table("data/20160725_Deer_PopSize.txt", header = T, stringsAsFactors = F)

#fitness <- read.table("data/AnnualData2b.txt", header = T, stringsAsFactors = F)

pedigree <- read.table("data/Pedigree_16-05-02.recoded.txt", header = T, stringsAsFactors = F)

recodes <- read.table("data/DeerRecodedIDs.txt", header = T, stringsAsFactors = F)

recsumm <- read.table("results/8_Per_ID_Recomb.txt", header = T, sep = "\t", stringsAsFactors = F)

#~~ format the basic data

head(basedata)

basedata$Sex[which(basedata$Sex == 1)] <- "F"
basedata$Sex[which(basedata$Sex == 2)] <- "M"
basedata$Sex[which(basedata$Sex == 3)] <- NA

names(basedata)[which(names(basedata) == "CaptureWt")]   <- "BirthWt"
names(basedata)[which(names(basedata) == "HindLeg")]     <- "BirthHindLeg"
names(basedata)[which(names(basedata) == "CaptureDate")] <- "BirthCapDate"

head(basedata)


#~~ Add pop size information

head(popsize)
names(popsize) <- c("BirthYear", "BirthDensity")

basedata <- join(basedata, popsize)


#~~ add code and pedigree information

head(basedata)
head(recodes)

names(recodes) <- c("Code", "ANIMAL")
basedata <- join(basedata, recodes)

basedata$ANIMAL[which(is.na(basedata$ANIMAL))] <- basedata$Code[which(is.na(basedata$ANIMAL))]
head(basedata)


#~~ add missing values to pedigree

pedigree <- rbind(pedigree,
                  cbind(ANIMAL = basedata$ANIMAL[which(!basedata$ANIMAL %in% pedigree$ANIMAL)],
                        MOTHER = 0, FATHER = 0))

head(pedigree)
tail(pedigree)


basedata <- join(basedata, pedigree)


#~~ Create the recsumm pheno object

head(recsumm)
names(recsumm)[which(names(recsumm) == "ANIMAL")] <- "RRID.Offspring"

recsumm <- recsumm[,c("Family", "RRID", "RRID.Offspring", "RRID.Sex",
                      "TotalInfLoci", "TotalChrCount", "TotalRecombCount")]

head(recsumm)

head(basedata)

temp.base <- subset(basedata, select = c(Code, ANIMAL, MOTHER, FATHER, BirthYear,
                                         BirthMonth, BirthWt, BirthDensity))
names(temp.base) <- c("RRID.Code", "RRID", "RRID.MOTHER", "RRID.FATHER", "RRID.BirthYear",
                      "RRID.BirthMonth", "RRID.BirthWt", "RRID.BirthDensity")

recsumm <- join(recsumm, temp.base)
head(recsumm)

rm(temp.base)

temp.base <- subset(basedata, select = c(Code, ANIMAL, BirthYear, BirthMonth, BirthWt, BirthDensity))
names(temp.base) <- c("Offspring.Code", "RRID.Offspring", "Offspring.BirthYear",
                      "Offspring.BirthMonth", "Offspring.BirthWt", "Offspring.BirthDensity")

recsumm <- join(recsumm, temp.base)
head(recsumm)

recsumm$Age <- recsumm$Offspring.BirthYear - recsumm$RRID.BirthYear

rm(temp.base)

save(list=ls(), file = "results/8_Pheno_Recsumm_Data_for_Analysis.RData")






