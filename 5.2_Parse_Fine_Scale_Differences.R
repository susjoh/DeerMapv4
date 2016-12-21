
library(ggplot2)
library(reshape)
library(beepr)
library(plyr)
library(dplyr)
library(GenABEL)
library(crimaptools)

load("finescale/5_Data_Fine_Scale_Chunk_Diffs.RData")


chunk.res <- NULL

for(i in 1:nrow(map.chunk.50)){
  
  print(paste("Running problem", i, "of", nrow(map.chunk.50)))
  
  test <- subset(mapdata, CEL.LG == map.chunk.50$CEL.LG[i])
  
  invmap <- read.table(paste0("finescale/chr", map.chunk.50$CEL.LG[i], AnalysisSuffix, "_cel_ckP", map.chunk.50$chunk[i], "_inv.parsemap"), header = T, stringsAsFactors = F)
  delmap <- read.table(paste0("finescale/chr", map.chunk.50$CEL.LG[i], AnalysisSuffix, "_cel_ckP", map.chunk.50$chunk[i], "_del.parsemap"), header = T, stringsAsFactors = F)
  
  map.chunk.50$Inversion.Len[i] <- invmap$cMPosition[nrow(invmap)]
  map.chunk.50$Initial.Len[i] <- test$cMPosition.run3[nrow(test)]
  map.chunk.50$Deletion.Len[i] <- delmap$cMPosition[nrow(delmap)]
  

  chunk.res <- rbind(chunk.res, invmap, delmap)
  
}

write.table(map.chunk.50, paste("results/5_MapChunk_DelInvCheck_Summary_", AnalysisSuffix2, ".txt"),
            row.names = F, quote = F)

write.table(chunk.res, paste("results/5_MapChunk_DelInvCheck_Full_", AnalysisSuffix2, ".txt"),
            row.names = F, quote = F)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. For small chunks, try exclusion/inversion #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

map.chunk.50 <- read.table(paste("results/5_MapChunk_DelInvCheck_Summary_", AnalysisSuffix2, ".txt"),
                           header = T, stringsAsFactors = F)


map.chunk.50$Inversion.Diff <-   map.chunk.50$Initial.Len - map.chunk.50$Inversion.Len
map.chunk.50$Deletion.Diff <-   map.chunk.50$Initial.Len - map.chunk.50$Deletion.Len
map.chunk.50$Order <- 1:nrow(map.chunk.50)

chunk.res <- read.table(paste("results/5_MapChunk_DelInvCheck_Full_", AnalysisSuffix2, ".txt"),
                        header = T, stringsAsFactors = F)

chunk.res$chunk <- unlist(lapply(chunk.res$analysisID,
                                 function(x) gsub("ck", "",  strsplit(x, split = "_")[[1]][3])))

names(chunk.res)[which(names(chunk.res) == "chunk")] <- "chunk.edit"
chunk.res <- subset(chunk.res, select = -Order)

chunk.res <- join(chunk.res, mapdata[,c("SNP.Name", "chunk", "cMPosition.run3", "Order")])

chunk.res$chunk.focal <- ifelse(chunk.res$chunk.edit == chunk.res$chunk, "focal", "nf")



ggplot(map.chunk.50, aes(Initial.Len, Inversion.Len)) +
  geom_text(aes(label = Freq)) +
  geom_abline(slope = 1, intercept = 0)

ggplot(map.chunk.50, aes(Initial.Len, Deletion.Len)) +
  geom_text(aes(label = Freq)) +
  geom_abline(slope = 1, intercept = 0)

ggplot(map.chunk.50, aes(Deletion.Diff, Freq)) +
  geom_point() 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Remove small chunks with large deletion difference # 
#    and inversions that result in > 0.5cM shorter map  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

edit.tab <- rbind(cbind(map.chunk.50[which(map.chunk.50$Freq < 10 & map.chunk.50$Deletion.Diff > 0.5),], edit = "del"),
                  cbind(map.chunk.50[which(map.chunk.50$Inversion.Diff > 0.5),], edit = "inv"))

arrange(edit.tab, edit, Freq)

table(edit.tab$chunk)

arrange(edit.tab, chunk)

#~~ If deletion and inversion, remove deletion

doub.chunk <- data.frame(table(edit.tab$chunk))
doub.chunk <- subset(doub.chunk, Freq > 1)

edit.tab <- edit.tab[-which(edit.tab$chunk %in% doub.chunk$Var1 & edit.tab$edit == "del"),]

head(mapdata)

#~~ create mapdata for LGs that were fine

newmap <- mapdata

for(i in 1:nrow(edit.tab)){

  if(edit.tab$edit[i] == "del"){
    
    newmap$CEL.order[which(newmap$chunk == edit.tab$chunk[i])] <- NA
    
  }
  
  if(edit.tab$edit[i] == "inv"){
    
    newmap$CEL.order[which(newmap$chunk == edit.tab$chunk[i])] <- rev(newmap$CEL.order[which(newmap$chunk == edit.tab$chunk[i])])
    
  }
  
}

newmap <- arrange(newmap, CEL.LG, CEL.order)
newmap <- subset(newmap, !is.na(CEL.order))

#~~ rerun crimap

counter <- 1

for(i in sort(unique(edit.tab$CEL.LG))){
  
  print(paste0("Running chromosome ", counter, " of ", length(unique(edit.tab$CEL.LG))))
  
  create_crimap_input(gwaa.data = abeldata,
                      familyPedigree = famped,
                      analysisID = paste0(i, AnalysisSuffix2),
                      snplist = subset(newmap, CEL.LG == i)$SNP.Name,
                      pseudoautoSNPs = pseudoautoSNPs,
                      is.X = ifelse(i == lg.sex, T, F),
                      outdir = paste0("crimap/crimap_", AnalysisSuffix2),
                      use.specific.mnd = paste0("crimap/crimap_", AnalysisSuffix, "/chr", AnalysisSuffix, "full.mnd"),
                      clear.existing.analysisID = TRUE)
  
  run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", i, AnalysisSuffix2, ".gen"))
  
  run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", i, AnalysisSuffix2, ".gen"))
  
  counter <- counter + 1
  
}


#~~ Parse Chrompic Files

fullmap <- NULL

for(lg in c(lg.vec)){
  
  recombmap <- parse_map_chrompic(paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".cmp"))
  recombmap$Chr <- lg
  
  fullmap <- rbind(fullmap, recombmap)
  rm(recombmap)
}

names(fullmap)[which(names(fullmap) == "cMPosition")] <- "cMPosition.run4"

head(fullmap)

fullmap <- join(fullmap, subset(newmap, select = c(SNP.Name, BTA.Chr, BTA.Position, BTA.Order, CEL.order, CEL.LG, chunk)))
table(is.na(fullmap$CEL.LG))


tapply(mapdata$cMPosition.run3, mapdata$CEL.LG, max) - tapply(fullmap$cMPosition.run4, fullmap$CEL.LG, max)
tapply(mapdata$cMPosition.run3, mapdata$CEL.LG, length) - tapply(fullmap$cMPosition.run4, fullmap$CEL.LG, length)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. write everything to file and update unmapped snps  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


ggplot(fullmap, aes(CEL.order, cMPosition.run4)) +
  geom_point() +
  facet_wrap(~Chr, scales= "free")

table(is.na(fullmap$cMPosition.run4))
table(is.na(fullmap$CEL.order))

ggsave(filename = paste0("figs/Linkage_Map_run4_", AnalysisSuffix, ".png"), device = "png", width = 15, height = 10, units = "in")


ggplot(subset(fullmap, CEL.LG == 34), aes(CEL.order, cMPosition.run4)) +
  geom_point() +
  facet_wrap(~Chr, scales= "free")

write.table(fullmap, paste0("results/5_Linkage_Map_Positions_CEL_run4_", AnalysisSuffix, ".txt"), row.names = F, quote = F, sep = "\t")

