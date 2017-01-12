
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Calculate recombination fractions between unmapped and mapped SNPs  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(plyr)
library(GenABEL)
library(crimaptools)
library(reshape)
library(ggplot2)


AnalysisSuffix <- "a"
AnalysisSuffix2 <- "a_cel"
AnalysisSuffix3 <- "a_skel"

out.prefix <- "data/Deer31_QC"
runMapPhasing <- FALSE

runLDplots <- FALSE
runLDCalcs <- FALSE
runRecComparison <- TRUE

#~~ Read & format linkage map data

mapdata_run1 <- read.table("results/2_Linkage_Map_Positions_run1_a.txt", header = T, stringsAsFactors = F)
mapdata      <- read.table("results/8_Linkage_Map_Positions_CEL_run5_dumpos_a.txt", header = T, stringsAsFactors = F)
rectab       <- read.table("results/8_Per_Chromosome_Recomb_final_a.txt", header = T, stringsAsFactors = F)

max.vals <- read.table("results/8_Predicted_Physical_Size_run5_a.txt", header = T)

#~~ Read & format genetic data

load("data/Deer31_QC.RData", verbose = T)

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

#~~ Read family data

famped <- read.table(paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"),
                     header = T, stringsAsFactors = F)

lg.vec <- sort(unique(mapdata$CEL.LG))

#~~ determine unmapped SNPs and split into autosomal/X

unsnp <- mapdata_run1$SNP.Name[which(!mapdata_run1$SNP.Name %in% mapdata$SNP.Name)]
unsnp.sex <- unsnp[grep("_x_", unsnp)]
unsnp.sex <- unsnp.sex[-which(unsnp.sex %in% pseudoautoSNPs)]
unsnp.nonsex <- unsnp[which(!unsnp %in% unsnp.sex)]


#~~ Determine phase with crimap. Do this with females only to bypass sex-chromosome links

head(famped)
famped.f <- famped[grep("Mum", famped$Family),]

lg <- 99


create_crimap_input(gwaa.data = abeldata,
                    familyPedigree = famped.f,
                    analysisID = paste0(lg, AnalysisSuffix2),
                    snplist = unsnp,
                    is.X = TRUE,
                    pseudoautoSNPs = unsnp.nonsex,
                    outdir = paste0("crimap/crimap_", AnalysisSuffix2),
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".gen"))

parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".pre"),
               genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".gen"),
               familyPedigree = famped.f,
               save.mendfile = TRUE)

mend.err <- read.table(paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".mndverbose"), header = T, sep = "\t")

create_crimap_input(gwaa.data = abeldata,
                    familyPedigree = famped.f,
                    analysisID = paste0(lg, AnalysisSuffix2),
                    snplist = unsnp,
                    is.X = TRUE,
                    pseudoautoSNPs = unsnp.nonsex,
                    outdir = paste0("crimap/crimap_", AnalysisSuffix2),
                    use.mnd = TRUE,
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".gen"))


run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".gen"))

unmap.rec <- parse_crossovers( paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".cmp"), familyPedigree = famped.f)


#~~ extract the phase information for X

head(rectab)
head(unmap.rec)

unmap.list <- lapply(unmap.rec$data, function(x){
  y <- data.frame(Phase = strsplit(x, split = "")[[1]],
                  SNP.Name = unsnp)
  y <- subset(y, Phase != "-")
  y
})

for(i in 1:length(unmap.list)) unmap.list[[i]]$Family <- unmap.rec$Family[i]

unmap.phase <- do.call(rbind, unmap.list)

#~~ Now for all mapped markers

x <- data.frame(analysisID = unique(rectab$analysisID),
                CEL.LG = 1:34)
rectab <- join(rectab, x)

if(runMapPhasing == TRUE){
  
  map.phase <- list()
  
  for(i in 1:34){
    
    print(paste("Running chromosome", i, "of", 34))
    
    temptab <- subset(rectab, CEL.LG == i)
    
    temp.list <- lapply(temptab$data, function(x){
      y <- data.frame(Phase = strsplit(x, split = "")[[1]],
                      SNP.Name = subset(mapdata, CEL.LG == i)$SNP.Name)
      y <- subset(y, Phase != "-")
      y
    })
    
    for(j in 1:length(temp.list)) temp.list[[j]]$Family <- temptab$Family[j]
    
    list.phase <- data.table::rbindlist(temp.list)
    
    map.phase[[i]] <- list.phase
    
    rm(temp.list, temptab, list.phase)
    
  }
  
  map.phase <- data.table::rbindlist(map.phase)
  
  save(map.phase, file = "results/8_Map_Phase_Information_v2.RData")
  
} else {
  load("results/8_Map_Phase_Information_v2.RData")
}

rm(x, AnalysisSuffix, AnalysisSuffix2, AnalysisSuffix3, i, lg, j, out.prefix, unmap.list)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Look at LD with females only (deals with X)      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

abeldata.fem <- abeldata[phdata(abeldata)$sex == 0,]
nids(abeldata.fem)
head(phdata(abeldata.fem))
abeldata.fem <- recodeChromosome(abeldata.fem, list("X" = 34))
table(chromosome(abeldata.fem))

#~~ Look 

if(runLDCalcs == TRUE){
  
  best.postab <- NULL
  
  for(i in 1:length(unsnp)){
    
    if(i %in% seq(1, length(unsnp), 10))
      print(paste("Calculating LD position for SNP", i, "of", length(unsnp)))
    
    
    x <- GenABEL::r2fast(abeldata.fem,
                         snpsubset = unsnp[i],
                         cross.snpsubset = snpnames(abeldata.fem[,snpnames(abeldata.fem) %in% mapdata$SNP.Name]))
    
    
    attributes(x$r2) <- NULL
    y <- data.frame(Count = x$num,
                    R2 = x$r2)
    
    y$SNP.Name <- row.names(y)
    
    suppressMessages(y <- join(y, mapdata))
    
    y <- arrange(y, -R2)[1:10, c("R2", "SNP.Name", "CEL.LG", "cMPosition.run5")]
    
    y$Unmapped.SNP <- unsnp[i]
    
    best.postab <- rbind(best.postab, y)
    
  }
  
  write.table(best.postab, "results/8_Best_LD_pos_for_unmapped_SNPs_v2.txt", row.names = F, sep = "\t", quote = F)
  
} else {
  
  best.postab <- read.table("results/8_Best_LD_pos_for_unmapped_SNPs_v2.txt", header = T)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Calculate best genomic position                         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

unmap.tab <- data.frame(LG.Count     = tapply(best.postab$CEL.LG,
                                              best.postab$Unmapped.SNP,
                                              function (x) length(table(x))),
                        Top.LG.Count = tapply(best.postab$CEL.LG,
                                              best.postab$Unmapped.SNP,
                                              function (x) ifelse(!all(is.na(x)),
                                                                  sort(table(x), decreasing = T)[[1]],
                                                                  NA)),
                        LG.NA    = tapply(best.postab$CEL.LG,
                                          best.postab$Unmapped.SNP,
                                          function (x) length(which(is.na(x)))),
                        Pos.Min  = tapply(best.postab$cMPosition.run5,
                                          best.postab$Unmapped.SNP,
                                          min, na.rm = T),
                        Pos.Max  = tapply(best.postab$cMPosition.run5,
                                          best.postab$Unmapped.SNP,
                                          max, na.rm = T),
                        Mean.R2  = tapply(best.postab$R2,
                                          best.postab$Unmapped.SNP,
                                          mean, na.rm = T),
                        Max.R2   = tapply(best.postab$R2,
                                          best.postab$Unmapped.SNP,
                                          max, na.rm = T))

best.postab$SNP.Rank <- rep(1:10, length.out = nrow(best.postab))

unmap.tab$Pos.Min[which(unmap.tab$LG.Count != 1)] <- NA
unmap.tab$Pos.Max[which(unmap.tab$LG.Count != 1)] <- NA
unmap.tab$Pos.Range <- unmap.tab$Pos.Max - unmap.tab$Pos.Min
unmap.tab$SNP.Name <- row.names(unmap.tab)

unmap.tab <- join(unmap.tab, mapdata_run1[,c("SNP.Name", "chunk")])

temp <- subset(best.postab, SNP.Rank == 1)
temp <- subset(temp, select = c(SNP.Name, Unmapped.SNP))
names(temp) <- c("Top.SNP", "SNP.Name")

unmap.tab <- join(unmap.tab, temp)
rm(temp)

temp <- subset(mapdata, select = c(SNP.Name, CEL.LG, cMPosition.run5))
names(temp) <- c("Top.SNP", "Top.CEL.LG", "Top.cM")

unmap.tab <- join(unmap.tab, temp)
rm(temp)

ggplot(unmap.tab, aes(Max.R2)) + geom_histogram() + facet_wrap(~LG.Count)
ggplot(unmap.tab, aes(Top.LG.Count)) + geom_histogram()
ggplot(unmap.tab, aes(Pos.Range)) + geom_histogram()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Get out chunks                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

names(map.phase)[1] <- "Phase.Map"

map.phase$Phase.Map <- as.numeric(as.character(map.phase$Phase.Map))
unmap.phase$Phase <- as.numeric(as.character(unmap.phase$Phase))

mapdata_run1.removed <- subset(mapdata_run1, SNP.Name %in% unsnp)
head(mapdata_run1.removed)

mapdata_run1.removed$Likely.LG <- NA
mapdata_run1.removed$MinLDPos <- NA
mapdata_run1.removed$MaxLDPos <- NA


#~~ Get the most likely CEL.LG

for(i in 1:nrow(mapdata_run1.removed)){
  
  chunk.id <- mapdata_run1.removed$chunk[i]
  
  chunk.tab <- subset(mapdata_run1, chunk == chunk.id)
  chunk.ld <- subset(unmap.tab, SNP.Name %in% chunk.tab$SNP.Name)
  chunk.LG <- names(sort(table(chunk.ld$Top.CEL.LG), decreasing = T))[1]
  
  mapdata_run1.removed$Likely.LG[i] <- chunk.LG
  
  map.vec <- which(mapdata$SNP.Name %in% chunk.ld$Top.SNP & mapdata$CEL.LG == chunk.LG)
  mapdata_run1.removed$MinLDPos[i] = min(map.vec) - 60
  mapdata_run1.removed$MaxLDPos[i] = max(map.vec) + 60
  
  
}


#~~ Make a table of all comparisons to be run

if(runRecComparison == TRUE){
  
  comparison.tab <- NULL
  
  for(i in 1:nrow(mapdata_run1.removed)){
    x <- data.frame(Map.Marker = mapdata$SNP.Name[mapdata_run1.removed$MinLDPos[i]:mapdata_run1.removed$MaxLDPos[i]],
                    Unmap.Marker = mapdata_run1.removed$SNP.Name[i])
    head(x)
    x <- subset(x, Map.Marker %in% subset(mapdata, CEL.LG == mapdata_run1.removed$Likely.LG[i])$SNP.Name)
    
    comparison.tab <- rbind(comparison.tab, x)
    rm(x)
  }
  
  comparison.tab$Map.Marker <- as.character(comparison.tab$Map.Marker)
  comparison.tab$Unmap.Marker <- as.character(comparison.tab$Unmap.Marker)
  
  comparison.tab$R2 <- NA
  comparison.tab$N  <- NA
  
  map.phase.1 <- subset(map.phase, SNP.Name %in% comparison.tab$Map.Marker)
  
  
  for(i in 1:nrow(comparison.tab)){
    if(i %in% seq(1, nrow(comparison.tab), 100)) print(paste("Running line", i, "of", nrow(comparison.tab)))
    suppressMessages(x <- join(map.phase.1[which(map.phase.1$SNP.Name == comparison.tab$Map.Marker[i]),c(1, 3)],
                               unmap.phase[which(unmap.phase$SNP.Name == comparison.tab$Unmap.Marker[i]),c(1, 3)],
                               type = "inner"))
    if(nrow(x) > 20) comparison.tab$R2[i] <- cor.test(x$Phase, x$Phase.Map)$estimate
    comparison.tab$N[i] <- nrow(x)
  }
  
  
  head(comparison.tab)
  head(mapdata_run1.removed)
  
  temp <- subset(mapdata_run1.removed, select = c(SNP.Name, chunk, Likely.LG))
  names(temp)[1] <- "Unmap.Marker"
  
  comparison.tab <- join(comparison.tab, temp)
  
  temp <- subset(mapdata, select = c(SNP.Name, BTA.Position, Dummy.Position))
  names(temp)[1] <- "Map.Marker"
  
  comparison.tab <- join(comparison.tab, temp)
  
  write.table(comparison.tab, "results/8_Rec_Fractions_Unplaced_Loci_v2.txt", row.names = F, sep = "\t", quote = F)  
  
} else {
  
  comparison.tab <- read.table("results/8_Rec_Fractions_Unplaced_Loci_v2.txt", header = T, stringsAsFactors = F)
  
}



head(comparison.tab)
table(comparison.tab$chunk)

# for(i in unique(comparison.tab$chunk)){
#   
#   x <- unique(subset(comparison.tab, chunk == i)$Likely.LG) 
#   y <- unique(subset(comparison.tab, chunk == i)$Unmap.Marker) 
#   
#   print(ggplot(subset(comparison.tab, chunk == i), aes(Dummy.Position, R2)) +
#           geom_point() +
#           stat_smooth() +
#           facet_wrap(~chunk, scales = "free") +
#           ggtitle(paste("Chunk", i, "with", length(y), "SNPs on linkage group", x)))
#   
# }


rm(runLDCalcs, runLDplots, map.vec, i, chunk.id, chunk.LG, runMapPhasing, runRecComparison, x, y)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Getting the approximate positions of all chunks       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(mapdata_run1.removed)
table(mapdata_run1.removed$chunk)
source("r/countIF.R")

mapdata_run1.removed$Freq <- countIF(mapdata_run1.removed$chunk)

#all.chunk.tab <- subset(mapdata_run1.removed, Freq < 5)
all.comp <- subset(comparison.tab, Unmap.Marker %in% mapdata_run1.removed$SNP.Name)
head(all.comp)

#~~ Choose a threshold of 0.9?

all.comp <- subset(all.comp, R2 >= 0.9)
head(all.comp)
names(all.comp)[1] <- "SNP.Name"
all.comp <- join(all.comp, mapdata[,c("SNP.Name", "cMPosition.run5", "CEL.LG")])




all.chunk.summ <- data.frame(Window.Start = tapply(all.comp$cMPosition.run5, all.comp$chunk, min),
                             Window.Stop  = tapply(all.comp$cMPosition.run5, all.comp$chunk, max),
                             CEL.LG       = tapply(as.numeric(all.comp$Likely.LG), all.comp$chunk, mean),
                             Unmap.SNP.vec = tapply(all.comp$Unmap.Marker, all.comp$chunk, function(x) paste(unique(x), collapse = ",")))

all.chunk.summ$Unmap.SNP.vec <- as.character(all.chunk.summ$Unmap.SNP.vec)
all.chunk.summ$chunk <- row.names(all.chunk.summ$chunk)
all.chunk.summ$SNP.Start <- NA
all.chunk.summ$SNP.Stop  <- NA

head(all.chunk.summ)

for(i in 1:nrow(all.chunk.summ)){
  x <- subset(all.comp, cMPosition.run5 == all.chunk.summ$Window.Start[i] & CEL.LG == all.chunk.summ$CEL.LG[i])
  all.chunk.summ$SNP.Start[i] <- x$SNP.Name[1]
  rm(x)
  x <- subset(all.comp, cMPosition.run5 == all.chunk.summ$Window.Stop[i] & CEL.LG == all.chunk.summ$CEL.LG[i])
  all.chunk.summ$SNP.Stop[i] <- x$SNP.Name[nrow(x)]
  rm(x)
  
}

head(all.chunk.summ)

all.chunk.summ$chunk <- row.names(all.chunk.summ)

all.chunk.summ$SNP.Start.Of.Chromosome <- NA
all.chunk.summ$SNP.End.Of.Chromosome <- NA

for(i in 1:nrow(all.chunk.summ)){

  all.chunk.summ$SNP.Start.Of.Chromosome[i] <-  ifelse(all.chunk.summ$SNP.Start[i] == subset(mapdata, CEL.LG == all.chunk.summ$CEL.LG[i])[1,"SNP.Name"], 
                                                     "yes", "no")
  all.chunk.summ$SNP.End.Of.Chromosome[i] <-  ifelse(all.chunk.summ$SNP.Stop[i] == subset(mapdata, CEL.LG == all.chunk.summ$CEL.LG[i])[nrow(subset(mapdata, CEL.LG == all.chunk.summ$CEL.LG[i])),"SNP.Name"], 
                                                     "yes", "no")

}

str(all.chunk.summ)

all.chunk.summ <- cbind(subset(all.chunk.summ, select = -Unmap.SNP.vec), subset(all.chunk.summ, select = Unmap.SNP.vec))

write.table(all.chunk.summ, "results/8_Predicted_Positions_of_Unmapped_Loci.txt", row.names = F, sep = "\t", quote = F)
write.table(all.chunk.summ, "doc/Table_Predicted_Positions_of_Unmapped_Loci.txt", row.names = F, sep = "\t", quote = F)
