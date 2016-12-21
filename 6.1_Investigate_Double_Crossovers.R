
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

#~~ Read & format linkage map data

mapdata <- read.table(paste0("results/5_Linkage_Map_Positions_CEL_run4_", AnalysisSuffix, ".txt"),
                      header = T, stringsAsFactors = F)

table(mapdata$cMdiff == c(diff(mapdata$cMPosition.run4), NA))

mapdata$cMdiff[which(mapdata$cMdiff < 0)] <- NA

mapdata <- arrange(mapdata, CEL.LG, CEL.order)


lg.vec <- sort(unique(mapdata$CEL.LG))
lg.sex <- 34

analysis.vec <- unique(mapdata$analysisID)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Deal with informative SNPs that occur on non-PAR X between fathers and sons  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

lg <- "34a_cel"

create_crimap_input(gwaa.data = abeldata,
                    familyPedigree = famped,
                    analysisID = paste0(lg),
                    snplist = mapdata$SNP.Name[which(mapdata$CEL.LG == 34)],
                    is.X = TRUE,
                    pseudoautoSNPs = pseudoautoSNPs,
                    use.specific.mnd = "crimap/crimap_a/chrafull.mnd",
                    outdir = paste0("crimap/crimap_", AnalysisSuffix2),
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".gen"))

run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".gen"))


x.map <- parse_map_chrompic(paste0("crimap/crimap_a_cel/chr34a_cel.cmp"))
x.cmp <- parse_crossovers(paste0("crimap/crimap_a_cel/chr34a_cel.cmp"), familyPedigree = famped)

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

subset(x.map, SNP.Name == "cela1_red_x_75365938")

bad.snps$Offspring <- NA
for(i in 1:nrow(bad.snps)) bad.snps$Offspring[i] <- as.numeric(strsplit(as.character(bad.snps$UniqueID[i]), split = "_")[[1]][5])

bad.mnd <- subset(bad.snps, select = c(SNP.Name, Offspring, RRID))
bad.mnd <- melt(bad.mnd, id.vars = "SNP.Name")
bad.mnd$variable <- NULL
names(bad.mnd)[2] <- "ANIMAL"

master.mnd <- read.table("crimap/crimap_a/chrafull.mnd", header = T, stringsAsFactors = F)
head(master.mnd)
master.mnd <- rbind(master.mnd, bad.mnd)
master.mnd <- unique(master.mnd)
write.table(master.mnd, "crimap/crimap_a_cel/chrafull_doubxover.mnd", row.names = F, sep = "\t", quote = F)


lg <- "34a_cel"

create_crimap_input(gwaa.data = abeldata,
                    familyPedigree = famped,
                    analysisID = paste0(lg),
                    snplist = x.map$SNP.Name[which(x.map$analysisID == lg)],
                    is.X = TRUE,
                    pseudoautoSNPs = pseudoautoSNPs,
                    use.specific.mnd = "crimap/crimap_a_cel/chrafull_doubxover.mnd",
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Parse crossovers                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


rectab <- NULL
doub.xovers <- NULL
cm.map <- NULL

for(i in analysis.vec){
  
  test <- parse_crossovers(paste0("crimap/crimap_a_cel/chr", i, ".cmp"), familyPedigree = famped)
  rectab <- rbind(rectab, test)
  
  tempmap <- parse_map_chrompic(paste0("crimap/crimap_a_cel/chr", i, ".cmp"))
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

write.table(rectab     , "results/6_Per_Chromosome_Recomb_raw_a.txt", row.names = F, sep = "\t", quote = F)
write.table(doub.xovers, "results/6_Double_Xovers_raw_a.txt", row.names = F, sep = "\t", quote = F)
write.table(cm.map, "results/6_cM_Map_raw_a.txt", row.names = F, sep = "\t", quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Deal with double crossovers           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

doub.xovers <- read.table("results/6_Double_Xovers_raw_a.txt"        , header = T, stringsAsFactors = F)
rectab      <- read.table("results/6_Per_Chromosome_Recomb_raw_a.txt", header = T, stringsAsFactors = F)
cm.map      <- read.table("results/6_cM_Map_raw_a.txt"                   , header = T, stringsAsFactors = F)

doub.xovers$is.X <- ifelse(doub.xovers$analysisID == "34a_cel", "yes", "no")

doub.xovers$SpanCountDist <- doub.xovers$StopSpan - doub.xovers$StartSpan
doub.xovers <- subset(doub.xovers, Type == "Mid")

ggplot(doub.xovers, aes(InfCount, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(SpanCountDist, fill = Singleton)) + geom_histogram(binwidth = 5) + scale_fill_brewer(palette = "Set1") 
ggplot(doub.xovers, aes(log10(SpanLength), fill = Singleton)) + geom_histogram(binwidth = 0.01) + scale_fill_brewer(palette = "Set1") + facet_wrap(~is.X)

ggplot(doub.xovers, aes(InfCount, SpanCountDist)) + geom_point(alpha = 0.2)
ggplot(doub.xovers, aes(SpanCountDist, fill = Singleton)) + geom_histogram(binwidth = 5) + scale_fill_brewer(palette = "Set1") + facet_wrap(~analysisID)

# Some short double crossover regions are caused because they have a singleton
# within a long region. At this point remove singletons and runs of <10cM and rerun the maps.

singletons <- subset(doub.xovers, Singleton == "yes")

singleton.analysisIDs <- unique(singletons$analysisID)
singleton.info <- cbind(singletons, counter = 0)

table(singletons$analysisID)

head(singletons)

singletons <- subset(singletons, select = c(StartPos, analysisID, RRID, UniqueID))
head(singletons)
names(singletons)[1] <- "Order"
singletons <- join(singletons, cm.map)
table(table(singletons$SNP.Name, singletons$RRID))
singletons$Offspring <- sapply(singletons$UniqueID, function(x) strsplit(x, split = "_")[[1]][5])
head(singletons)
table(singletons$analysisID)
singletons <- subset(singletons, select = c(RRID, Offspring, SNP.Name))
singletons <- melt(singletons, id.vars = "SNP.Name")

head(singletons)
singletons <- subset(singletons, select = c(value, SNP.Name))
names(singletons)[1] <- "ANIMAL"


#~~ Create a master .mnd file and add the singleton information

master.mnd <- read.table("crimap/crimap_a_cel/chrafull_doubxover.mnd", header = T, stringsAsFactors = F)  
master.mnd <- rbind(master.mnd, singletons)
master.mnd <- unique(master.mnd) # 2543

write.table(master.mnd, "crimap/crimap_a_cel/chrafull_doubxover.mnd", row.names = F, quote = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Re-run crimap and remove singletons   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

lg.sex <- "34a_cel"

counter <- 0

singleton.analysisIDs <- unique(gsub("_dbx", "", singleton.analysisIDs))


system.time({
  for(lg in singleton.analysisIDs){
    
    print(paste("Running linkage group", lg, "of", length(analysis.vec)))
    
    create_crimap_input(gwaa.data = abeldata,
                        familyPedigree = famped,
                        analysisID = paste0(lg, "_dbx"),
                        snplist = mapdata$SNP.Name[which(mapdata$analysisID == lg)],
                        is.X = ifelse(lg == lg.sex, TRUE, FALSE),
                        pseudoautoSNPs = pseudoautoSNPs,
                        use.specific.mnd = "crimap/crimap_a_cel/chrafull_doubxover.mnd",
                        outdir = paste0("crimap/crimap_", AnalysisSuffix2),
                        clear.existing.analysisID = TRUE)
    
    run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, "_dbx.gen"))
    
    run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, "_dbx.gen"))
    
  }
})

#~~ Parse double crossovers

rectab <- NULL
doub.xovers <- NULL
cm.map <- NULL

for(i in analysis.vec){
  
  test <- parse_crossovers(paste0("crimap/crimap_a_cel/chr", i, "_dbx.cmp"), familyPedigree = famped)
  rectab <- rbind(rectab, test)
  
  tempmap <- parse_map_chrompic(paste0("crimap/crimap_a_cel/chr", i, "_dbx.cmp"))
  tempmap <- subset(tempmap, select = c(SNP.Name, cMPosition, Order, analysisID))
  names(tempmap)[2] <- "Position"
  
  
  test2 <- check_double_crossovers(test, physical.map = tempmap)
  doub.xovers <- rbind(doub.xovers, test2)
  
  tempmap$analysisID <- i
  cm.map <- rbind(cm.map, tempmap)
  
}

rectemp <- data.frame(TotalRecombCount = tapply(rectab$RecombCount, rectab$Family, sum))
rectemp$Family<- row.names(rectemp)

ggplot(rectemp, aes(TotalRecombCount)) + geom_histogram()

head(rectemp)
subset(rectemp,TotalRecombCount > 50)


write.table(rectab     , "results/6_Per_Chromosome_Recomb_dbx_a.txt", row.names = F, sep = "\t", quote = F)
write.table(doub.xovers, "results/6_Double_Xovers_dbx_a.txt", row.names = F, sep = "\t", quote = F)
write.table(cm.map,      "results/6_cM_Map_dbx_a.txt", row.names = F, sep = "\t", quote = F)

doub.xovers$is.X <- ifelse(doub.xovers$analysisID == "34a_cel", "yes", "no")

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

# singleton.analysisIDs <- unique(singletons$analysisID)
# head(singletons)
# 
# if(length(singleton.analysisIDs) > 0){
#   singleton.info <- rbind(singleton.info, cbind(singletons, counter = counter))
#   
#   head(singleton.info)
#   
#   table(singleton.info$counter, singleton.info$analysisID)
#   
#   write.table(singleton.info, paste0("results/6_singletoninfo_iteration", counter, ".txt"), row.names = F, sep = "\t", quote = F) 
#   
#   singletons <- subset(singletons, select = c(StartPos, analysisID, RRID, UniqueID, Family))
#   head(singletons)
#   names(singletons)[1] <- "Order"
#   
#   singletons <- join(singletons, cm.map)
#   singletons$Offspring <- sapply(singletons$Family, function(x) strsplit(x, split = "_")[[1]][3])
#   singletons <- subset(singletons, select = c(RRID, Offspring, SNP.Name))
#   singletons <- melt(singletons, id.vars = "SNP.Name")
#   
#   head(singletons)
#   singletons <- subset(singletons, select = c(value, SNP.Name))
#   names(singletons)[1] <- "ANIMAL"
#   
#   #~~ Create a master .mnd file and add the singleton information
#   
#   master.mnd <- read.table("crimap/crimap_a_cel/chrafull_doubxover.mnd", header = T, stringsAsFactors = F)
#   master.mnd <- rbind(master.mnd, singletons)
#   master.mnd <- unique(master.mnd)
#   
#   write.table(master.mnd, "crimap/crimap_a_cel/chrafull_doubxover.mnd", row.names = F, quote = F)
# }
# 
# }
# 
# }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Check short stretches                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rectab <- read.table("results/6_Per_Chromosome_Recomb_dbx_a.txt", header = T, stringsAsFactors = F)
doub.xovers <- read.table("results/6_Double_Xovers_dbx_a.txt", header = T, stringsAsFactors = F)
cm.map <- read.table("results/6_cM_Map_dbx_a.txt", header = T, stringsAsFactors = F)

dodgy.xovers <- subset(doub.xovers, log10(doub.xovers$SpanLength) < 0.5 & Type == "Mid")
table(dodgy.xovers$InfCount)
dodgy.xovers <- arrange(dodgy.xovers, InfCount)

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

master.mnd <- read.table("crimap/crimap_a_cel/chrafull_doubxover.mnd", header = T)
head(master.mnd)

master.mnd <- rbind(master.mnd, dodgy.markers)
master.mnd <- unique(master.mnd)
write.table(master.mnd, "crimap/crimap_a_cel/chrafull_doubxover.mnd", row.names = F, sep = "\t", quote = F)

singleton.analysisIDs <- unique(dodgy.xovers$analysisID)


print(paste0("Running iteration ", counter))

singleton.analysisIDs <- unique(gsub("_dbx", "", singleton.analysisIDs))


system.time({
  for(lg in singleton.analysisIDs){
    
    print(paste("Running linkage group", lg, "of", length(singleton.analysisIDs)))
    
    create_crimap_input(gwaa.data = abeldata,
                        familyPedigree = famped,
                        analysisID = paste0(lg, "_dbx"),
                        snplist = mapdata$SNP.Name[which(mapdata$analysisID == lg)],
                        is.X = ifelse(lg == lg.sex, TRUE, FALSE),
                        pseudoautoSNPs = pseudoautoSNPs,
                        use.specific.mnd = "crimap/crimap_a_cel/chrafull_doubxover.mnd",
                        outdir = paste0("crimap/crimap_", AnalysisSuffix2),
                        clear.existing.analysisID = TRUE)
    
    run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, "_dbx.gen"))
    

    run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, "_dbx.gen"))
    
  }
})


#~~ Parse double crossovers

rectab <- NULL
doub.xovers <- NULL
cm.map <- NULL

for(i in analysis.vec){
  
  test <- parse_crossovers(paste0("crimap/crimap_a_cel/chr", i, "_dbx.cmp"), familyPedigree = famped)
  rectab <- rbind(rectab, test)
  
  tempmap <- parse_map_chrompic(paste0("crimap/crimap_a_cel/chr", i, "_dbx.cmp"))
  tempmap <- subset(tempmap, select = c(SNP.Name, cMPosition, Order, analysisID))
  names(tempmap)[2] <- "Position"
  
  
  test2 <- check_double_crossovers(test, physical.map = tempmap)
  doub.xovers <- rbind(doub.xovers, test2)
  
  tempmap$analysisID <- i
  cm.map <- rbind(cm.map, tempmap)
  
}


doub.xovers$is.X <- ifelse(doub.xovers$analysisID == "34a_cel", "yes", "no")

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

write.table(rectab     , "results/6_Per_Chromosome_Recomb_dbx_a.txt", row.names = F, sep = "\t", quote = F)
write.table(doub.xovers, "results/6_Double_Xovers_dbx_a.txt", row.names = F, sep = "\t", quote = F)
write.table(cm.map,      "results/6_cM_Map_dbx_a.txt", row.names = F, sep = "\t", quote = F)


# 
# singleton.analysisIDs <- unique(singletons$analysisID)
# head(singletons)
# 
# if(length(singleton.analysisIDs) > 0){
#   singleton.info <- rbind(singleton.info, cbind(singletons, counter = counter))
#   
#   head(singleton.info)
#   
#   table(singleton.info$counter, singleton.info$analysisID)
#   
#   write.table(singleton.info, paste0("results/6_singletoninfo_iteration", counter, ".txt"), row.names = F, sep = "\t", quote = F) 
#   
#   singletons <- subset(singletons, select = c(StartPos, analysisID, RRID, UniqueID, Family))
#   head(singletons)
#   names(singletons)[1] <- "Order"
#   
#   singletons <- join(singletons, cm.map)
#   singletons$Offspring <- sapply(singletons$Family, function(x) strsplit(x, split = "_")[[1]][3])
#   singletons <- subset(singletons, select = c(RRID, Offspring, SNP.Name))
#   singletons <- melt(singletons, id.vars = "SNP.Name")
#   
#   head(singletons)
#   singletons <- subset(singletons, select = c(value, SNP.Name))
#   names(singletons)[1] <- "ANIMAL"
#   
#   #~~ Create a master .mnd file and add the singleton information
#   
#   master.mnd <- read.table("crimap/crimap_a_cel/chrafull_doubxover.mnd", header = T, stringsAsFactors = F)
#   master.mnd <- rbind(master.mnd, singletons)
#   master.mnd <- unique(master.mnd)
#   
#   write.table(master.mnd, "crimap/crimap_a_cel/chrafull_doubxover.mnd", row.names = F, quote = F)
# }
# 
# }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Get the final map information            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rectab <- read.table("results/6_Per_Chromosome_Recomb_dbx_a.txt", header = T, stringsAsFactors = F)
doub.xovers <- read.table("results/6_Double_Xovers_dbx_a.txt", header = T, stringsAsFactors = F)
cm.map <- read.table("results/6_cM_Map_dbx_a.txt", header = T, stringsAsFactors = F)


#~~ Get sex specific maps and save map information

head(cm.map)
head(mapdata)

table(mapdata$analysisID)

cm.map <- subset(cm.map, select = c(SNP.Name, Position))
names(cm.map)[2] <- "cMPosition.run5"
mapdata <- join(mapdata, cm.map)
mapdata <- subset(mapdata, select = -cMPosition.run4)

which(is.na(mapdata$cMPosition.run5))

analysis.vec <- unique(mapdata$analysisID)

sexsp.map <- NULL

for(i in analysis.vec){
  print(paste("Running ", i))
  run_crimap_map(paste0("crimap/crimap_a_cel/chr", i, "_dbx.gen"))
  sexsp.map <- rbind(sexsp.map, parse_map(paste0("crimap/crimap_a_cel/chr", i, "_dbx.map")))
  
}

sexsp.map <- subset(sexsp.map, select = -c(analysisID, Order))
mapdata <- join(mapdata, sexsp.map)

head(mapdata)
sex.map <- melt(subset(mapdata, select = c(SNP.Name, CEL.order, CEL.LG, cMPosition.Female, cMPosition.Male)),
                id.vars = c("SNP.Name", "CEL.order", "CEL.LG"))

head(sex.map)
sex.map$variable <- gsub("cMPosition.", "", sex.map$variable)

ggplot(sex.map, aes(CEL.order, value, col = variable)) +
  geom_point() +
  facet_wrap(~CEL.LG, scales = "free") +
  scale_color_brewer(palette = "Set1")


write.table(mapdata, "results/6_Linkage_Map_Positions_CEL_run5_a_HOLD.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Deal with X again
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ggplot(subset(sex.map, CEL.LG == 34), aes(CEL.order, value, col = variable)) +
  geom_point() +
  facet_wrap(~CEL.LG, scales = "free") +
  scale_color_brewer(palette = "Set1")


x.map <- parse_map_chrompic(paste0("crimap/crimap_a_cel/chr34a_cel_dbx.cmp"))
x.cmp <- parse_crossovers(paste0("crimap/crimap_a_cel/chr34a_cel_dbx.cmp"), familyPedigree = famped)

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
for(i in 1:nrow(bad.snps)) bad.snps$Offspring[i] <- as.numeric(strsplit(as.character(bad.snps$UniqueID[i]), split = "_")[[1]][6])

bad.mnd <- subset(bad.snps, select = c(SNP.Name, Offspring, RRID))
bad.mnd <- melt(bad.mnd, id.vars = "SNP.Name")
bad.mnd$variable <- NULL
names(bad.mnd)[2] <- "ANIMAL"

master.mnd <- read.table("crimap/crimap_a_cel/chrafull_doubxover.mnd", header = T, stringsAsFactors = F)
head(master.mnd)
master.mnd <- rbind(master.mnd, bad.mnd)
master.mnd <- unique(master.mnd)
write.table(master.mnd, "crimap/crimap_a_cel/chrafull_doubxover.mnd", row.names = F, sep = "\t", quote = F)


lg <- "34a_cel_dbx"

create_crimap_input(gwaa.data = abeldata,
                    familyPedigree = famped,
                    analysisID = paste0(lg),
                    snplist = x.map$SNP.Name[which(x.map$analysisID == lg)],
                    is.X = TRUE,
                    pseudoautoSNPs = pseudoautoSNPs,
                    use.specific.mnd = "crimap/crimap_a_cel/chrafull_doubxover.mnd",
                    outdir = paste0("crimap/crimap_", AnalysisSuffix2),
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".gen"))

run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".gen"))

run_crimap_map(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".gen"))

mapx <- parse_map(paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, ".map"))

sex.mapx <- melt(subset(mapx, select = c(Order, SNP.Name, cMPosition.Female, cMPosition.Male)),
                 id.vars = c("SNP.Name", "Order"))

ggplot(sex.mapx, aes(Order, value, col = variable)) + geom_point() + scale_color_brewer(palette = "Set1")


beep()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Make final files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


rectab <- NULL
cm.map <- NULL

sexsp.map <- NULL

for(i in analysis.vec){
  
  test <- parse_crossovers(paste0("crimap/crimap_a_cel/chr", i, "_dbx.cmp"), familyPedigree = famped)
  rectab <- rbind(rectab, test)
  
  tempmap <- parse_map_chrompic(paste0("crimap/crimap_a_cel/chr", i, "_dbx.cmp"))
  tempmap <- subset(tempmap, select = c(SNP.Name, cMPosition, Order, analysisID))
  names(tempmap)[2] <- "Position"
  
  tempmap$analysisID <- i
  cm.map <- rbind(cm.map, tempmap)
  
  sexsp.map <- rbind(sexsp.map, parse_map(paste0("crimap/crimap_a_cel/chr", i, "_dbx.map")))
  
  
}

head(mapdata)
mapdata <- subset(mapdata, select = -c(cMPosition.Female, cMPosition.Male, Female.r, Male.r, cMdiff.Female, cMdiff.Male, cMPosition.run5))
table(mapdata$analysisID)

cm.map <- subset(cm.map, select = c(SNP.Name, Position))
names(cm.map)[2] <- "cMPosition.run5"
mapdata <- join(mapdata, cm.map)



sexsp.map <- subset(sexsp.map, select = -c(analysisID, Order))
mapdata <- join(mapdata, sexsp.map)

head(mapdata)
sex.map <- melt(subset(mapdata, select = c(SNP.Name, CEL.order, CEL.LG, cMPosition.Female, cMPosition.Male)),
                id.vars = c("SNP.Name", "CEL.order", "CEL.LG"))

head(sex.map)
sex.map$variable <- gsub("cMPosition.", "", sex.map$variable)

ggplot(sex.map, aes(CEL.order, value, col = variable)) +
  geom_point() +
  facet_wrap(~CEL.LG, scales = "free") +
  scale_color_brewer(palette = "Set1")

ggsave(filename = paste0("figs/Map_run5_", AnalysisSuffix, ".png"), device = "png", width = 15, height = 10, units = "in")


write.table(rectab, "results/6_Per_Chromosome_Recomb_final_a.txt")
write.table(mapdata, "results/6_Linkage_Map_Positions_CEL_run5_a.txt")


