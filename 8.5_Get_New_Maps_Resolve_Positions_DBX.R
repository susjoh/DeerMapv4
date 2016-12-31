
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

AnalysisSuffix <- "a"
AnalysisSuffix2 <- "a_cel"
out.prefix <- "data/Deer31_QC"
firstRun <- FALSE

#~~ Read & format genetic data

load("data/Deer31_QC.RData", verbose = T)

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

#~~ Read family data

famped <- read.table(paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"),
                     header = T, stringsAsFactors = F)



lg.vec <- sort(unique(mapdata$CEL.LG))
lg.sex <- 34

analysis.vec <- unique(mapdata$analysisID)

#~~ Load new map

load("results/8_Remap_missing_loci_new_orders.RData")

AnalysisSuffix2 <- "a_cel_rebuild"
newmaporder$SNP.Name <- as.character(newmaporder$SNP.Name)

mapdata <- read.table("results/6_Linkage_Map_Positions_CEL_run5_dummyBTApos_a.txt", header = T, stringsAsFactors = F)
mapdata_run1 <- read.table("results/2_Linkage_Map_Positions_run1_a.txt", header = T, stringsAsFactors = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Deal with anything informative on male X         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



lg <- 34
  
create_crimap_input(gwaa.data = abeldata,
                      familyPedigree = famped,
                      analysisID = lg,
                      snplist = subset(newmaporder, CEL.LG == lg)$SNP.Name,
                      is.X = ifelse(lg == lg.sex, TRUE, FALSE),
                      pseudoautoSNPs = pseudoautoSNPs,
                      use.specific.mnd = "crimap/crimap_a/chrafull.mnd",
                      outdir = paste0("crimap/crimap_", AnalysisSuffix2),
                      clear.existing.analysisID = TRUE)
  
run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".gen"))
  
run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".gen"))


x.map <- parse_map_chrompic(paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".cmp"))
x.cmp <- parse_crossovers(paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".cmp"), familyPedigree = famped)

phase.x <- lapply(x.cmp$data, function(x) data.frame(InfMarker = strsplit(x, split = "")[[1]],
                                                     Order = 1:nchar(x)))
for(i in 1:nrow(x.cmp)) phase.x[[i]] <- cbind(phase.x[[i]], UniqueID   = x.cmp$UniqueID[i], RRID = x.cmp$RRID[i])

phase.x <- do.call(rbind, phase.x)

head(phase.x)
head(x.map)

phase.x <- subset(phase.x, InfMarker != "-")
phase.x <- join(phase.x, subset(x.map, select = c(Order, SNP.Name)))
phase.x$PAR.SNP <- ifelse(phase.x$SNP.Name %in% pseudoautoSNPs, "yes", "no")

phase.x$SEX <- NA
phase.x$SEX[grep("FATHER", phase.x$UniqueID)] <- "Male"
phase.x$SEX[grep("MOTHER", phase.x$UniqueID)] <- "Female"

head(phase.x)

table(phase.x$PAR.SNP, phase.x$SEX)

bad.snps <- subset(phase.x, PAR.SNP == "no" & SEX == "Male")

#~~ Get rid of these dodgy SNPs and rerun the map   

bad.snps$Offspring <- NA
for(i in 1:nrow(bad.snps)) bad.snps$Offspring[i] <- as.numeric(strsplit(as.character(bad.snps$UniqueID[i]), split = "_")[[1]][4])

bad.mnd <- subset(bad.snps, select = c(SNP.Name, Offspring, RRID))
bad.mnd <- melt(bad.mnd, id.vars = "SNP.Name")
bad.mnd$variable <- NULL
names(bad.mnd)[2] <- "ANIMAL"

master.mnd <- read.table("crimap/crimap_a/chrafull.mnd", header = T, stringsAsFactors = F)
head(master.mnd)
master.mnd <- rbind(master.mnd, bad.mnd)
master.mnd <- unique(master.mnd)
write.table(master.mnd, "crimap/crimap_a_cel_rebuild/chrafull_doubxover.mnd", row.names = F, sep = "\t", quote = F)


lg <- 34

create_crimap_input(gwaa.data = abeldata,
                    familyPedigree = famped,
                    analysisID = lg,
                    snplist = subset(newmaporder, CEL.LG == lg)$SNP.Name,
                    is.X = ifelse(lg == lg.sex, TRUE, FALSE),
                    pseudoautoSNPs = pseudoautoSNPs,
                    use.specific.mnd = "crimap/crimap_a_cel_rebuild/chrafull_doubxover.mnd",
                    outdir = paste0("crimap/crimap_", AnalysisSuffix2),
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".gen"))

run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".gen"))

run_crimap_map(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".gen"))

mapx <- parse_map(paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".map"))

sex.mapx <- melt(subset(mapx, select = c(Order, SNP.Name, cMPosition.Female, cMPosition.Male)),
                 id.vars = c("SNP.Name", "Order"))

ggplot(sex.mapx, aes(Order, value, col = variable)) + geom_point() + scale_color_brewer(palette = "Set1")

rm(sex.mapx, x.cmp, x.map, phase.x)

beep()

#~~ Rerun

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Create a new crimap build                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


for(lg in 1:34){

      print(paste("Running linkage group", lg, "of", 34))
      
      create_crimap_input(gwaa.data = abeldata,
                          familyPedigree = famped,
                          analysisID = lg,
                          snplist = subset(newmaporder, CEL.LG == lg)$SNP.Name,
                          is.X = ifelse(lg == lg.sex, TRUE, FALSE),
                          pseudoautoSNPs = pseudoautoSNPs,
                          use.specific.mnd = "crimap/crimap_a_cel_rebuild/chrafull_doubxover.mnd",
                          outdir = paste0("crimap/crimap_", AnalysisSuffix2),
                          clear.existing.analysisID = TRUE)
      
      run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".gen"))
      
      run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".gen"))
      
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Parse double crossovers                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rectab <- NULL
doub.xovers <- NULL
cm.map <- NULL

for(i in 1:34){
  
  test <- parse_crossovers(paste0("crimap/crimap_a_cel_rebuild/chr", i, ".cmp"), familyPedigree = famped)
  rectab <- rbind(rectab, test)
  
  tempmap <- parse_map_chrompic(paste0("crimap/crimap_a_cel_rebuild/chr", i, ".cmp"))
  tempmap <- subset(tempmap, select = c(SNP.Name, cMPosition, Order, analysisID))
  names(tempmap)[2] <- "Position"
  
  
  test2 <- check_double_crossovers(test, physical.map = tempmap)
  doub.xovers <- rbind(doub.xovers, test2)
  
  tempmap$analysisID <- i
  cm.map <- rbind(cm.map, tempmap)
  
}

ggplot(rectab, aes(RecombCount)) + geom_histogram()
rectemp <- data.frame(TotalRecombCount = tapply(rectab$RecombCount, rectab$Family, sum))
rectemp$Family<- row.names(rectemp)

ggplot(rectemp, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)


head(rectemp)


#~ Deal with the crossovers

doub.xovers$is.X <- ifelse(doub.xovers$analysisID == "34a_cel", "yes", "no")

doub.xovers$SpanCountDist <- doub.xovers$StopSpan - doub.xovers$StartSpan
doub.xovers <- subset(doub.xovers, Type == "Mid")

ggplot(doub.xovers, aes(InfCount, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(SpanCountDist, fill = Singleton)) + geom_histogram(binwidth = 5) + scale_fill_brewer(palette = "Set1") 
ggplot(doub.xovers, aes(log10(SpanLength), fill = Singleton)) + geom_histogram(binwidth = 0.01) + scale_fill_brewer(palette = "Set1") + facet_wrap(~is.X)

ggplot(doub.xovers, aes(InfCount, SpanCountDist)) + geom_point(alpha = 0.2)
ggplot(doub.xovers, aes(SpanCountDist, fill = Singleton)) + geom_histogram(binwidth = 5) + scale_fill_brewer(palette = "Set1") + facet_wrap(~analysisID)

# Some short double crossover regions are caused because they have a singleton
# within a long region. At this point remove singletons.

singletons <- subset(doub.xovers, Singleton == "yes")

singleton.analysisIDs <- unique(singletons$analysisID)
singleton.info <- cbind(singletons, counter = 0)


singletons <- subset(singletons, select = c(StartPos, analysisID, RRID, UniqueID))
names(singletons)[1] <- "Order"
singletons <- join(singletons, cm.map)
table(table(singletons$SNP.Name, singletons$RRID))
singletons$Offspring <- sapply(singletons$UniqueID, function(x) strsplit(x, split = "_")[[1]][4])
table(singletons$analysisID)
singletons <- subset(singletons, select = c(RRID, Offspring, SNP.Name))
singletons <- melt(singletons, id.vars = "SNP.Name")

head(singletons)
singletons <- subset(singletons, select = c(value, SNP.Name))
names(singletons)[1] <- "ANIMAL"

nrow(singletons)
nrow(doub.xovers)

#~~ Create a master .mnd file and add the singleton information

master.mnd <- read.table("crimap/crimap_a_cel_rebuild/chrafull_doubxover.mnd", header = T, stringsAsFactors = F)  
master.mnd <- rbind(master.mnd, singletons)
master.mnd <- unique(master.mnd) # 2543

write.table(master.mnd, "crimap/crimap_a_cel_rebuild/chrafull_doubxover.mnd", row.names = F, quote = F)

write.table(rectab     , "results/6_Per_Chromosome_Recomb_raw_a_rebuild.txt", row.names = F, sep = "\t", quote = F)
write.table(doub.xovers, "results/6_Double_Xovers_raw_a_rebuild.txt", row.names = F, sep = "\t", quote = F)
write.table(cm.map, "results/6_cM_Map_raw_a_rebuild.txt", row.names = F, sep = "\t", quote = F)








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Rerun with masked double crossovers              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


for(lg in 1:34){
  
  print(paste("Running linkage group", lg, "of", 34))
  
  create_crimap_input(gwaa.data = abeldata,
                      familyPedigree = famped,
                      analysisID = paste0(lg, "_dbx"),
                      snplist = subset(newmaporder, CEL.LG == lg)$SNP.Name,
                      is.X = ifelse(lg == lg.sex, TRUE, FALSE),
                      pseudoautoSNPs = pseudoautoSNPs,
                      use.specific.mnd = "crimap/crimap_a_cel_rebuild/chrafull_doubxover.mnd",
                      outdir = paste0("crimap/crimap_", AnalysisSuffix2),
                      clear.existing.analysisID = TRUE)

  run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, "_dbx.gen"))

  run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, "_dbx.gen"))
  
  run_crimap_map(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, "_dbx.gen"))
  
}
















rectab <- NULL
doub.xovers <- NULL
cm.map <- NULL

for(i in 1:34){
  
  test <- parse_crossovers(paste0("crimap/crimap_a_cel_rebuild/chr", i, "_dbx.cmp"), familyPedigree = famped)
  rectab <- rbind(rectab, test)
  
  tempmap <- parse_map_chrompic(paste0("crimap/crimap_a_cel_rebuild/chr", i, "_dbx.cmp"))
  tempmap <- subset(tempmap, select = c(SNP.Name, cMPosition, Order, analysisID))
  names(tempmap)[2] <- "Position"
  
  
  test2 <- check_double_crossovers(test, physical.map = tempmap)
  doub.xovers <- rbind(doub.xovers, test2)
  
  tempmap$analysisID <- i
  cm.map <- rbind(cm.map, tempmap)
  
}

ggplot(rectab, aes(RecombCount)) + geom_histogram()
rectemp <- data.frame(TotalRecombCount = tapply(rectab$RecombCount, rectab$Family, sum))
rectemp$Family<- row.names(rectemp)

ggplot(rectemp, aes(TotalRecombCount)) + geom_histogram(binwidth = 1)


head(rectemp)


#~ Deal with the crossovers

doub.xovers$is.X <- ifelse(doub.xovers$analysisID == "34a_cel", "yes", "no")

doub.xovers$SpanCountDist <- doub.xovers$StopSpan - doub.xovers$StartSpan
doub.xovers <- subset(doub.xovers, Type == "Mid")

ggplot(doub.xovers, aes(InfCount, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(SpanCountDist, fill = Singleton)) + geom_histogram(binwidth = 5) + scale_fill_brewer(palette = "Set1") 
ggplot(doub.xovers, aes(log10(SpanLength), fill = Singleton)) + geom_histogram(binwidth = 0.01) + scale_fill_brewer(palette = "Set1") + facet_wrap(~is.X)

ggplot(doub.xovers, aes(InfCount, SpanCountDist)) + geom_point(alpha = 0.2)
ggplot(doub.xovers, aes(SpanCountDist, fill = Singleton)) + geom_histogram(binwidth = 5) + scale_fill_brewer(palette = "Set1") + facet_wrap(~analysisID)


dodgy.xovers <- subset(doub.xovers, log10(doub.xovers$SpanLength) < 0.5 & Type == "Mid")

table(dodgy.xovers$InfCount)
dodgy.xovers <- arrange(dodgy.xovers, Family)

#~~ What are the dodgy SNPs?

dodgy.xovers <- join(dodgy.xovers, subset(rectab, select = c(UniqueID, data)))
head(dodgy.xovers)

dodgy.markers <- NULL

for(i in 1:nrow(dodgy.xovers)){
  
  dodgy.markers <- rbind(dodgy.markers, data.frame(InfMarker  = strsplit(substring(dodgy.xovers$data[i], dodgy.xovers$StartPos[i], dodgy.xovers$StopPos[i]), split = "")[[1]],
                                                   Order      = dodgy.xovers$StartPos[i]:dodgy.xovers$StopPos[i],
                                                   analysisID = dodgy.xovers$analysisID[i],
                                                   UniqueID   = dodgy.xovers$UniqueID[i],
                                                   RRID       = dodgy.xovers$RRID[i]))
}

dodgy.markers <- subset(dodgy.markers, InfMarker != "-")
dodgy.markers$analysisID <- gsub("_dbx", "", dodgy.markers$analysisID)

dodgy.markers <- join(dodgy.markers, cm.map)

table(table(dodgy.markers$SNP.Name))

head(dodgy.markers)

dodgy.markers$Offspring <- sapply(dodgy.markers$UniqueID, function(x){
  x <- strsplit(as.character(x), split = "_")[[1]]
  x <- x[length(x)-2]
  x
})

dodgy.markers <- melt(subset(dodgy.markers, select = c(RRID, Offspring, SNP.Name)), id.vars = "SNP.Name")
head(dodgy.markers)
dodgy.markers <- subset(dodgy.markers, select = -variable)
names(dodgy.markers)[2] <- "ANIMAL"


#~~ Create a master .mnd file and add the dodgy marker information

master.mnd <- read.table("crimap/crimap_a_cel_rebuild/chrafull_doubxover.mnd", header = T, stringsAsFactors = F)  
master.mnd <- rbind(master.mnd, dodgy.markers)
master.mnd <- unique(master.mnd) # 2543

write.table(master.mnd, "crimap/crimap_a_cel_rebuild/chrafull_doubxover.mnd", row.names = F, quote = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# RUN AGAIN                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Deal with X
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


x.map <- parse_map_chrompic(paste0("crimap/crimap_a_cel_rebuild/chr34_dbx.cmp"))
x.cmp <- parse_crossovers(paste0("crimap/crimap_a_cel_rebuild/chr34_dbx.cmp"), familyPedigree = famped)

phase.x <- lapply(x.cmp$data, function(x) data.frame(InfMarker = strsplit(x, split = "")[[1]],
                                                     Order = 1:nchar(x)))
for(i in 1:nrow(x.cmp)) phase.x[[i]] <- cbind(phase.x[[i]], UniqueID   = x.cmp$UniqueID[i], RRID = x.cmp$RRID[i])

phase.x <- do.call(rbind, phase.x)

head(phase.x)
head(x.map)

phase.x <- subset(phase.x, InfMarker != "-")
phase.x <- join(phase.x, subset(x.map, select = c(Order, SNP.Name)))
phase.x$PAR.SNP <- ifelse(phase.x$SNP.Name %in% pseudoautoSNPs, "yes", "no")

phase.x$SEX <- NA
phase.x$SEX[grep("FATHER", phase.x$UniqueID)] <- "Male"
phase.x$SEX[grep("MOTHER", phase.x$UniqueID)] <- "Female"

head(phase.x)

table(phase.x$PAR.SNP, phase.x$SEX)

bad.snps <- subset(phase.x, PAR.SNP == "no" & SEX == "Male")

#~~ Get rid of these dodgy SNPs and rerun the map   

bad.snps$Offspring <- NA
for(i in 1:nrow(bad.snps)) bad.snps$Offspring[i] <- as.numeric(strsplit(as.character(bad.snps$UniqueID[i]), split = "_")[[1]][5])

bad.mnd <- subset(bad.snps, select = c(SNP.Name, Offspring, RRID))
bad.mnd <- melt(bad.mnd, id.vars = "SNP.Name")
bad.mnd$variable <- NULL
names(bad.mnd)[2] <- "ANIMAL"

master.mnd <- read.table("crimap/crimap_a_cel_rebuild/chrafull_doubxover.mnd", header = T, stringsAsFactors = F)
head(master.mnd)
master.mnd <- rbind(master.mnd, bad.mnd)
master.mnd <- unique(master.mnd)
write.table(master.mnd, "crimap/crimap_a_cel_rebuild/chrafull_doubxover.mnd", row.names = F, sep = "\t", quote = F)


# 
# create_crimap_input(gwaa.data = abeldata,
#                     familyPedigree = famped,
#                     analysisID = paste0(lg),
#                     snplist = x.map$SNP.Name[which(x.map$analysisID == lg)],
#                     is.X = TRUE,
#                     pseudoautoSNPs = pseudoautoSNPs,
#                     use.specific.mnd = "crimap/crimap_a_cel/chrafull_doubxover.mnd",
#                     outdir = paste0("crimap/crimap_", AnalysisSuffix2),
#                     clear.existing.analysisID = TRUE)
# 
# run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".gen"))
# 
# run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".gen"))
# 
# run_crimap_map(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".gen"))
# 
# mapx <- parse_map(paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".map"))
# 
# sex.mapx <- melt(subset(mapx, select = c(Order, SNP.Name, cMPosition.Female, cMPosition.Male)),
#                  id.vars = c("SNP.Name", "Order"))
# 
# ggplot(sex.mapx, aes(Order, value, col = variable)) + geom_point() + scale_color_brewer(palette = "Set1")
# 






for(lg in 1:34){
  
  print(paste("Running linkage group", lg, "of", 34))
  
  create_crimap_input(gwaa.data = abeldata,
                      familyPedigree = famped,
                      analysisID = paste0(lg, "_dbx"),
                      snplist = subset(newmaporder, CEL.LG == lg)$SNP.Name,
                      is.X = ifelse(lg == lg.sex, TRUE, FALSE),
                      pseudoautoSNPs = pseudoautoSNPs,
                      use.specific.mnd = "crimap/crimap_a_cel_rebuild/chrafull_doubxover.mnd",
                      outdir = paste0("crimap/crimap_", AnalysisSuffix2),
                      clear.existing.analysisID = TRUE)
  
  run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, "_dbx.gen"))
  
  run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, "_dbx.gen"))
  
  run_crimap_map(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, "_dbx.gen"))
  
}



#~~ Parse double crossovers

rectab <- NULL
doub.xovers <- NULL
cm.map <- NULL
sexsp.map <- NULL

for(i in 1:34){
  
  test <- parse_crossovers(paste0("crimap/crimap_a_cel_rebuild/chr", i, "_dbx.cmp"), familyPedigree = famped)
  rectab <- rbind(rectab, test)
  
  tempmap <- parse_map_chrompic(paste0("crimap/crimap_a_cel_rebuild/chr", i, "_dbx.cmp"))
  tempmap <- subset(tempmap, select = c(SNP.Name, cMPosition, Order, analysisID))
  names(tempmap)[2] <- "Position"
  
  
  test2 <- check_double_crossovers(test, physical.map = tempmap)
  doub.xovers <- rbind(doub.xovers, test2)
  

  
}


for(i in 1:34){
  
  
  tempmap <- parse_map_chrompic(paste0("crimap/crimap_a_cel_rebuild/chr", i, "_dbx.cmp"))
  tempmap$analysisID <- i
  cm.map <- rbind(cm.map, tempmap)
  
  sexsp.map <- rbind(sexsp.map, parse_map(paste0("crimap/crimap_a_cel_rebuild/chr", i, "_dbx.map")))
  
  
}


doub.xovers$is.X <- ifelse(doub.xovers$analysisID == "34", "yes", "no")

doub.xovers$SpanCountDist <- doub.xovers$StopSpan - doub.xovers$StartSpan
doub.xovers <- subset(doub.xovers, Type == "Mid")

ggplot(doub.xovers, aes(InfCount, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(SpanCountDist, fill = Singleton)) + geom_histogram(binwidth = 5) + scale_fill_brewer(palette = "Set1") 
ggplot(doub.xovers, aes(log10(SpanLength), fill = Singleton)) + geom_histogram(binwidth = 0.01) + scale_fill_brewer(palette = "Set1") + facet_wrap(~is.X)

ggplot(doub.xovers, aes(InfCount, SpanCountDist)) + geom_point(alpha = 0.2)
ggplot(doub.xovers, aes(SpanCountDist, fill = Singleton)) + geom_histogram(binwidth = 5) + scale_fill_brewer(palette = "Set1") + facet_wrap(~analysisID)

ggplot(doub.xovers, aes(analysisID, fill = Singleton)) + geom_bar() + scale_fill_brewer(palette = "Set1")

# Some short double crossover regions are caused because they have a singleton
# within a long region. At this point remove singletons and rerun the maps.

singletons <- subset(doub.xovers, Singleton == "yes")
singletons$analysisID <- gsub("_dbx", "", singletons$analysisID)
singletons$UniqueID <- gsub("_dbx", "", singletons$UniqueID)

print(singletons)


#~~ Sort out the map

head(cm.map)
head(sexsp.map)
names(cm.map)[c(1, 6)] <- c("CEL.order", "CEL.LG")
sexsp.map <- subset(sexsp.map, select = -c(Order, analysisID))

cm.map <- join(cm.map, sexsp.map)
names(cm.map)[3] <- "cMPosition.run5"

cm.map$Order <- NULL

mapdata_run1 <- subset(mapdata_run1, select = c(SNP.Name, Chr, Position))
names(mapdata_run1) <- c("SNP.Name", "BTA.Chr", "BTA.Position")

cm.map <- join(cm.map, mapdata_run1)

ggplot(cm.map, aes(BTA.Position, cMPosition.run5, col = factor(BTA.Chr))) +
  geom_point(alpha = 0.2) +
  facet_wrap(~CEL.LG, ncol = 5, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "BTA Position (cM)",
       y = "Linkage Map Length (cM)")


write.table(rectab, "results/8_Per_Chromosome_Recomb_final_a.txt", row.names = F, sep = "\t", quote = F)
write.table(cm.map, "results/8_Linkage_Map_Positions_CEL_run5_a.txt", row.names = F, sep = "\t", quote = F)
write.table(doub.xovers, "results/8_Double_Crossovers_Final_a.txt", row.names = F, sep = "\t", quote = F)






