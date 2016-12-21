# Much of this manipulation has been specified manually in code after reading 
# Slate et al 2002 and Ensembl Cattle BTA and any reference to Chromosome are to
# the cattle genome. It should be run line by line and checked. CEL is the
# cervus elaphus information.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set Working Environment and Load in Data  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(ggplot2)
library(reshape)
library(beepr)
library(plyr)
library(dplyr)
library(GenABEL)
library(crimaptools)

source("r/countIF.R")

AnalysisSuffix <- "a"
out.prefix <- "data/Deer31_QC"
chr.vec <- 1:29    # 1:29
runLDplots <- FALSE
cMDiffCutoff <- 3

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")


load("data/Deer31_QC.RData", verbose = T)


famped <- read.table(paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"),
                     header = T, stringsAsFactors = F)

mapdata <- read.table(paste0("results/2_Linkage_Map_Positions_run2_", AnalysisSuffix, ".txt"),
                      header = T, stringsAsFactors = F)

# According to Slate et al 2002, the following rearrangements have occurred:

cattle.v.deer <- read.table("data/CattleDeerComparison.txt", header = T, stringsAsFactors = F)


#~~ format map data ~~~~~~~~~~~~~~~~~~~~~

mapdata <- arrange(mapdata, Chr, Position)

mapdata$BTA.Order <- NA
mapdata$BTA.Order[1] <- 1
for(i in 2:nrow(mapdata)) mapdata$BTA.Order[i] <- ifelse(mapdata$Chr[i] > mapdata$Chr[(i-1)], 1, mapdata$BTA.Order[(i-1)] + 1)

mapdata$PosDiff <- c(diff(mapdata$Position), NA)
mapdata$PosDiff[which(mapdata$PosDiff < 0)] <- NA

mapdata$PosDiffb4 <- c(NA, diff(mapdata$Position))
mapdata$cMdiffb4  <- c(NA, diff(mapdata$cMPosition))

mapdata$PosDiffb4[which(mapdata$PosDiffb4 < 0)] <- NA
mapdata$cMdiffb4[which(mapdata$cMdiffb4 < 0)] <- NA

mapdata$cMStretch  <- mapdata$cMdiff + mapdata$cMdiffb4
mapdata$PosStretch <- mapdata$PosDiff + mapdata$PosDiffb4


mapdata$chunk <- NA
mapdata$chunk[1] <- 1

for(i in 2:nrow(mapdata)){
  if(i %in% seq(1, nrow(mapdata), 1000)) print(paste("Analysing row", i, "of", nrow(mapdata)))
  mapdata$chunk[i] <- ifelse(is.na(mapdata$cMdiffb4[i]) || mapdata$cMdiffb4[i] > cMDiffCutoff, mapdata$chunk[i-1] + 1,  mapdata$chunk[i-1])
}

mapdata <- subset(mapdata, select = -c(cMdiffb4, PosDiffb4, PosStretch, cMStretch))

head(mapdata)

table(mapdata$chunk)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. LD clustering by BTA Chromosome           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

table(chromosome(abeldata))
table(mapdata$Chr)

#~~ LD Plots per Chromosome 

if(runLDplots == TRUE){
  
  ld.mats <- list()
  
  for(i in c(chr.vec, 30)){
    
    print(paste("Running chromosome", i, "of", length(chr.vec)))
    
    chrNo <- abeldata[,mapdata$SNP.Name[which(mapdata$Chr == i)]]
    test <- t(r2fast(chrNo))
    test[upper.tri(test)] <- NA
    
    ld.mats[[i]] <- test
    
    flat.test <- melt(test)
    flat.test$X1.Order <- rep(1:nrow(test), times = nrow(test))
    flat.test$X2.Order <- rep(1:nrow(test), each = nrow(test))
    flat.test <- subset(flat.test, !is.na(value))
    
    head(flat.test)
    
    ggplot(flat.test, aes(X1.Order, X2.Order, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red")
    
    ggsave(filename = paste0("figs/LDmap_chr", i, ".png"), device = "png", width = 10, height = 10, units = "in")
    
    rm(test, flat.test)
  }
  
  save(ld.mats, file = paste0("results/3_LD_matrices_before_rearrange_", AnalysisSuffix, ".RData"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Deal with re-arrangements                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

mapdata$CEL.order <- NA
mapdata$CEL.LG    <- NA

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2a. CHROMOSOME 1                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Middle section of chromosome 1 is at the opposite end of the third segment (fission and fusion).

ggplot(subset(mapdata, Chr == 1), aes(Position, cMPosition)) + geom_point()
data.frame(table(subset(mapdata, Chr == 1)$Chr, subset(mapdata, Chr == 1)$chunk))

chr1.newmap <- NULL
for(i in c(1, 3, 2)) chr1.newmap <- rbind(chr1.newmap, subset(mapdata, chunk == i))
chr1.newmap$test.order <- 1:nrow(chr1.newmap)
chr1.newmap <- droplevels(chr1.newmap)

chrNo <- abeldata[,which(snpnames(abeldata) %in% chr1.newmap$SNP.Name)]
analysisID <- paste0(1, AnalysisSuffix, "_re")

# create_crimap_input(gwaa.data = chrNo,
#                     familyPedigree = famped,
#                     analysisID = analysisID,
#                     snplist = chr1.newmap$SNP.Name,
#                     outdir = paste0("crimap/crimap_", AnalysisSuffix),
#                     use.specific.mnd = paste0("crimap/crimap_", AnalysisSuffix, "/chr1", AnalysisSuffix, ".mnd"),
#                     clear.existing.analysisID = TRUE)
# 
# run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"))
# 
# parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".pre"),
#                genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"),
#                familyPedigree = famped,
#                save.mendfile = TRUE)
# 
# 
# 
# run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"))

#~~ read in recombination map

recombmap <- parse_map_chrompic(paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".cmp"))
head(recombmap)
recombmap <- join(recombmap, subset(chr1.newmap, select = c(SNP.Name, test.order)))

ggplot(recombmap, aes(test.order, cMPosition)) + geom_point()

#~~ Edit the master map

mapdata$CEL.LG[which(mapdata$chunk == 1)] <- 31
mapdata$CEL.order[which(mapdata$chunk == 1)] <- rev(1:length(which(mapdata$chunk == 1)))

mapdata$CEL.LG[which(mapdata$chunk %in% c(2, 3))] <- 19
mapdata$CEL.order[which(mapdata$chunk == 3)] <- 1:length(which(mapdata$chunk == 3))
mapdata$CEL.order[which(mapdata$chunk == 2)] <- (length(which(mapdata$chunk == 3))+1):length(which(mapdata$chunk %in% c(2, 3)))

table(table(mapdata$CEL.order[which(mapdata$chunk %in% c(2, 3))]))

rm(recombmap, chr1.newmap, analysisID)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2b. CHROMOSOMES 2, 5, 6, 8, 9                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ details from investgation of Slate et al

chr.splits <- subset(mapdata, Chr %in% c(2, 5, 6, 8, 9))
chr.splits.tab <- data.frame(table(chr.splits$Chr, chr.splits$chunk))
chr.splits.tab <- subset(chr.splits.tab, Freq != 0)
names(chr.splits.tab) <- c("Chr", "chunk", "Freq")

subset(cattle.v.deer, Cattle %in% chr.splits.tab$Chr)

ggplot(subset(mapdata, Chr == 9), aes(BTA.Order, cMPosition, col = chunk)) + geom_point()

mapdata$chunk[which(mapdata$chunk == 16)] <- 15

chr.splits <- subset(mapdata, Chr %in% c(2, 5, 6, 8, 9))
chr.splits.tab <- data.frame(table(chr.splits$Chr, chr.splits$chunk))
chr.splits.tab <- subset(chr.splits.tab, Freq != 0)
names(chr.splits.tab) <- c("Chr", "chunk", "Freq")

subset(cattle.v.deer, Cattle %in% chr.splits.tab$Chr)


chr.splits.tab$CEL.LG <- c(33, 8, 3, 22, 17, 6, 29, 16, 28, 26)
chr.splits.tab$direction <- c("fwd", "rev", "fwd", "rev", "fwd", "rev", "fwd", "rev", "fwd", "fwd")


for(i in 1:nrow(chr.splits.tab)){
  
  temp.which <- which(mapdata$chunk == chr.splits.tab$chunk[i])
  mapdata$CEL.LG[temp.which] <- chr.splits.tab$CEL.LG[i]
  
  if(chr.splits.tab$direction[i] == "fwd") mapdata$CEL.order[temp.which] <- 1:length(temp.which)
  if(chr.splits.tab$direction[i] == "rev") mapdata$CEL.order[temp.which] <- rev(1:length(temp.which))
  
  rm(temp.which)
  
}

table(mapdata$CEL.LG)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2c. Single Chromosomes                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# AT THE MOMENT ASSUME SAME DIRECTION AS CATTLE MAPS

cattle.v.deer$CattleCount <- countIF(cattle.v.deer$Cattle)
cattle.v.deer.singletons <- subset(cattle.v.deer, CattleCount ==1)

#~~ Deal with simple cases where BTA corresponds to a single deer chromosome

for(i in 1:nrow(cattle.v.deer.singletons)){
  mapdata$CEL.LG[which(mapdata$Chr == cattle.v.deer.singletons$Cattle[i])] <- cattle.v.deer.singletons$Deer[i]
  mapdata$CEL.order[which(mapdata$Chr == cattle.v.deer.singletons$Cattle[i])] <- 1:length(which(mapdata$Chr == cattle.v.deer.singletons$Cattle[i]))
  
}

table(mapdata$CEL.LG, useNA = "always")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2d. BTA17 and BTA19 are CEL5                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

chrNo <- abeldata[,which(chromosome(abeldata) %in%  c(17, 19))]

test <- t(r2fast(chrNo))
test[upper.tri(test)] <- NA

flat.test <- melt(test)
flat.test$X1.Order <- rep(1:nrow(test), times = nrow(test))
flat.test$X2.Order <- rep(1:nrow(test), each = nrow(test))
flat.test <- subset(flat.test, !is.na(value))

head(flat.test)

ggplot(flat.test, aes(X1.Order, X2.Order, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")


#~~ BTA 17 needs to be reversed

mapdata$CEL.order[which(mapdata$Chr == 17)] <- rev(1:length(which(mapdata$Chr == 17)))
mapdata$CEL.order[which(mapdata$Chr == 19)] <- (length(which(mapdata$Chr == 17)) + 1):length(which(mapdata$Chr %in% c(17, 19)))

#~~ make subset

chr17.19.subset <- subset(mapdata, Chr %in% c(17, 19))
table(chr17.19.subset$CEL.LG)

chr17.19.subset <- arrange(chr17.19.subset, CEL.order)

table(table(chr17.19.subset$CEL.order))

analysisID <- paste0("17.19", AnalysisSuffix, "_re")

# create_crimap_input(gwaa.data = chrNo,
#                     familyPedigree = famped,
#                     analysisID = analysisID,
#                     snplist = chr17.19.subset$SNP.Name,
#                     outdir = paste0("crimap/crimap_", AnalysisSuffix),
#                     use.specific.mnd = "crimap/crimap_a/chrafull.mnd",
#                     clear.existing.analysisID = TRUE)
# 
# run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"))
# 
# parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".pre"),
#                genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"),
#                familyPedigree = famped,
#                save.mendfile = TRUE)
# 
# run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"))




#~~ read in recombination map


recombmap <- parse_map_chrompic(paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".cmp"))
head(recombmap)
recombmap <- join(recombmap, subset(chr17.19.subset, select = c(SNP.Name, CEL.order)))

ggplot(recombmap, aes(CEL.order, cMPosition)) + geom_point()

rm(chrNo, test, flat.test)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2e. BTA26 and BTA28 are CEL15                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

chrNo <- abeldata[,which(chromosome(abeldata) %in%  c(26, 28))]

# test <- t(r2fast(chrNo))
# test[upper.tri(test)] <- NA
# 
# flat.test <- melt(test)
# flat.test$X1.Order <- rep(1:nrow(test), times = nrow(test))
# flat.test$X2.Order <- rep(1:nrow(test), each = nrow(test))
# flat.test <- subset(flat.test, !is.na(value))
# 
# head(flat.test)
# 
# ggplot(flat.test, aes(X1.Order, X2.Order, fill = value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "red")


#~~ BTA 28, then 26

mapdata <- arrange(mapdata, Chr, Position)

mapdata$CEL.order[which(mapdata$Chr == 28)] <- 1:length(which(mapdata$Chr == 28))
mapdata$CEL.order[which(mapdata$Chr == 26)] <- (length(which(mapdata$Chr == 28)) + 1):length(which(mapdata$Chr %in% c(26, 28)))

#~~ deal with small rearrangement on Chr 28
temp <- table(mapdata$chunk[which(mapdata$Chr == 28)])
mapdata$CEL.order[which(mapdata$chunk == names(temp)[1])] <- rev(mapdata$CEL.order[which(mapdata$chunk == names(temp)[1])])

#~~ make subset

chr28.26.subset <- subset(mapdata, Chr %in% c(28, 26))

table(chr28.26.subset$chunk)

chr28.26.subset <- arrange(chr28.26.subset, CEL.order)

analysisID <- paste0("28.26", AnalysisSuffix, "_re")

create_crimap_input(gwaa.data = chrNo,
                    familyPedigree = famped,
                    analysisID = analysisID,
                    snplist = chr28.26.subset$SNP.Name,
                    outdir = paste0("crimap/crimap_", AnalysisSuffix),
                    use.specific.mnd = "crimap/crimap_a/chrafull.mnd",
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"))

parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".pre"),
               genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"),
               familyPedigree = famped,
               save.mendfile = TRUE)


run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"))

beep()
#~~ read in recombination map

recombmap <- parse_map_chrompic(paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".cmp"))
head(recombmap)
recombmap <- join(recombmap, subset(chr28.26.subset, select = c(SNP.Name, CEL.order)))

ggplot(recombmap, aes(CEL.order, cMPosition)) + geom_point()

rm(chr28.26.subset, chrNo, test, flat.test)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2f. Inversion on Chromosome 13               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ggplot(subset(mapdata, Chr == 13), aes(Position, cMPosition)) + geom_point()
inversionchunk <- data.frame(table(subset(mapdata, Chr == 13)$Chr, subset(mapdata, Chr == 13)$chunk))

mapdata$CEL.order[which(mapdata$chunk == inversionchunk$Var2[2])] <- rev(mapdata$CEL.order[which(mapdata$chunk == inversionchunk$Var2[2])])

chr13.subset <- subset(mapdata, Chr == 13)

chrNo <- abeldata[,which(chromosome(abeldata) %in%  c(13))]

#~~ make subset

chr13.subset <- arrange(chr13.subset, CEL.order)

analysisID <- paste0("13", AnalysisSuffix, "_re")

create_crimap_input(gwaa.data = chrNo,
                    familyPedigree = famped,
                    analysisID = analysisID,
                    snplist = chr13.subset$SNP.Name,
                    use.specific.mnd = "crimap/crimap_a/chrafull.mnd",
                    outdir = paste0("crimap/crimap_", AnalysisSuffix),
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"))

parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".pre"),
               genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"),
               familyPedigree = famped,
               save.mendfile = TRUE)

run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"))

beep()
#~~ read in recombination map

recombmap <- parse_map_chrompic(paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".cmp"))
head(recombmap)
recombmap <- join(recombmap, subset(chr13.subset, select = c(SNP.Name, CEL.order)))

ggplot(recombmap, aes(Order, cMPosition)) + geom_point()

rm(chr13.subset, chrNo)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. WRITE TO FILE                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

names(mapdata)[which(names(mapdata) == "cMPosition")] <- "cMPosition.run2"
names(mapdata)[which(names(mapdata) == "Position")] <- "BTA.Position"
names(mapdata)[which(names(mapdata) == "Chr")] <- "BTA.Chr"

head(mapdata)
str(mapdata)

ggplot(mapdata, aes(CEL.order, cMPosition.run2)) + geom_point() + facet_wrap(~CEL.LG)


write.table(mapdata, paste0("results/3_Linkage_Map_Positions_CEL.order_", AnalysisSuffix, "_woX.txt"),
            row.names = F, quote = F)

save.image("3.0_Large_Scale_Temp.Rdata")
