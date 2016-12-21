# First linkage map for the deer, before dealing with fine scale rearrangements.
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


load("finescale/5_Data_Fine_Scale_Chunk_Diffs.RData")

#~~ Read & format linkage map data

mapdata <- read.table(paste0("results/5_Linkage_Map_Positions_CEL_run4_", AnalysisSuffix, ".txt"),
                      header = T, stringsAsFactors = F)
head(mapdata)

table(mapdata$cMdiff == c(diff(mapdata$cMPosition.run4), NA))

mapdata$cMdiff[which(mapdata$cMdiff < 0)] <- NA

mapdata$BTA.Pos.Diff <- c(diff(mapdata$BTA.Position), NA)
mapdata$BTA.Pos.Diff[which(mapdata$BTA.Pos.Diff < 0)] <- NA
mapdata$BTA.Pos.Diff[which(mapdata$BTA.Pos.Diff > 2e6)] <- NA

mapdata <- arrange(mapdata, CEL.LG, CEL.order)

tapply(mapdata$cMPosition.run4, mapdata$CEL.LG, max)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Locate gaps in the assembly           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ggplot(mapdata, aes(cMdiff)) + geom_histogram(binwidth = 0.05)
ggplot(mapdata, aes(BTA.Pos.Diff)) + geom_histogram()
ggplot(mapdata, aes(BTA.Pos.Diff, cMdiff)) + geom_point()
ggplot(snpmap, aes(Pos.Diff)) + geom_histogram()
ggplot(subset(mapdata, cMdiff > 0), aes(cMdiff)) + geom_histogram(binwidth = 0.05)


mapdata <- arrange(mapdata, Chr, CEL.order)
mapdata$cMdiffb4  <- c(NA, diff(mapdata$cMPosition))
mapdata$cMdiffb4[which(mapdata$cMdiffb4 < 0)] <- NA
mapdata$cMStretch  <- mapdata$cMdiff + mapdata$cMdiffb4

cMDiffCutoff <- 1

mapdata$chunk <- NA
mapdata$chunk[1] <- 1


for(i in 2:nrow(mapdata)){
  if(i %in% seq(1, nrow(mapdata), 1000)) print(paste("Analysing row", i, "of", nrow(mapdata)))
  mapdata$chunk[i] <- ifelse(is.na(mapdata$cMdiffb4[i]) || mapdata$cMdiffb4[i] > cMDiffCutoff, mapdata$chunk[i-1] + 1,  mapdata$chunk[i-1])
}

mapdata <- subset(mapdata, select = -c(cMdiffb4, cMStretch))

head(mapdata)
tail(mapdata)


map.chunk <- data.frame(table(mapdata$chunk, mapdata$CEL.LG))
map.chunk <- subset(map.chunk, Freq > 0)
names(map.chunk) <- c("chunk", "CEL.LG", "Freq")

map.chunk$test <- c(-999, map.chunk$CEL.LG[-nrow(map.chunk)])
map.chunk$test2 <- c(map.chunk$CEL.LG[-1], 999)

map.chunk$Type <- ifelse(map.chunk$test < as.numeric(map.chunk$CEL.LG), "First",
                         ifelse(map.chunk$test2 > as.numeric(map.chunk$CEL.LG), "Last", "Mid"))
map.chunk <- subset(map.chunk, select = -c(test, test2))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Check possible fine-scale reaarangements  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(mapdata)
mapdata$chunk.col <- ifelse(mapdata$chunk %in% seq(1, nrow(map.chunk), 2), 1, 2)

for(lg in c(lg.vec, lg.sex)){
  
  sub.map <- subset(mapdata, CEL.LG == lg)
  ggplot(sub.map, aes(CEL.order, cMPosition.run4, col = as.character(chunk.col))) +
    geom_point(alpha = 0.5) +
    scale_color_brewer(palette = "Set1")
  ggsave(filename = paste0("figs/Linkage_Map_run4_", lg, "_",  AnalysisSuffix2, ".png"),
         device = "png", width = 10, height = 10, units = "in")
}

rm(lg)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. For small chunks, try exclusion/inversion #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ggplot(map.chunk, aes(Freq)) + geom_histogram(binwidth = 20)





map.chunk.50 <- subset(map.chunk, Freq < 50)
map.chunk.50$Inversion.Len <- NA
map.chunk.50$Deletion.Len <- NA
map.chunk.50$Initial.Len <- NA

save.image(file = "finescale/5_Recheck_Data_Fine_Scale_Chunk_Diffs.RData")

