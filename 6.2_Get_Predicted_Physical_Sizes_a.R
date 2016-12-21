#~~ deer cattle comparison
source("r/countIF.R")


cat.deer <- read.table("data/CattleDeerComparison.txt", header = T)
cat.deer$CattleCount <- countIF(cat.deer$Cattle)
cat.deer$DeerCount <- countIF(cat.deer$Deer)


#~~ new CEL map

mapdata <- read.table("results/6_Linkage_Map_Positions_CEL_run5_a.txt", header = T)

max.map <- data.frame(max.cM = tapply(mapdata$cMPosition.run5, mapdata$CEL.LG, max),
                      LocusCount = tapply(mapdata$SNP.Name, mapdata$CEL.LG, length))
max.map$CEL.LG <- row.names(max.map)
max.map$Est.Length <- NA

#~~ old bovine map

snpmap <- read.table("data/Deer31.map")
names(snpmap) <- c("BTA.Chr", "SNP.Name", "X", "BTA.Position")
snpmap <- subset(snpmap, select = -X)

max.bta <- data.frame(Max.Position = tapply(snpmap$BTA.Position, snpmap$BTA.Chr, max))
max.bta$BTA.Chr <- row.names(max.bta)

#~~ single chromosomes

sing.cat <- cat.deer[which(cat.deer$CattleCount == 1 & cat.deer$DeerCount == 1),]
cat.deer <- cat.deer[-which(cat.deer$CattleCount == 1 & cat.deer$DeerCount == 1),]

for(i in 1:nrow(sing.cat)){
  max.map$Est.Length[which(max.map$CEL.LG == sing.cat$Deer[i])] <- max.bta$Max.Position[which(max.bta$BTA.Chr == sing.cat$Cattle[i])]
}

#~~ fused chromosomes

max.map$Est.Length[which(max.map$CEL.LG == 5)] <- sum(max.bta$Max.Position[which(max.bta$BTA.Chr %in% c(17, 19))])
max.map$Est.Length[which(max.map$CEL.LG == 15)] <- sum(max.bta$Max.Position[which(max.bta$BTA.Chr %in% c(26, 28))])

cat.deer <- subset(cat.deer, DeerCount == 1)

#~~ split chromosomes

split.map <- data.frame(Max.BTA = tapply(mapdata$BTA.Position, mapdata$CEL.LG, max),
                        Min.BTA = tapply(mapdata$BTA.Position, mapdata$CEL.LG, min))
split.map$CEL.LG <- row.names(split.map)
split.map$Est.Length <- split.map$Max.BTA - split.map$Min.BTA
split.map

split.map <- subset(split.map, CEL.LG %in% max.map$CEL.LG[is.na(max.map$Est.Length)])

for(i in 1:nrow(split.map)){
  max.map$Est.Length[which(max.map$CEL.LG == split.map$CEL.LG[i])] <- split.map$Est.Length[i]
}


library(ggplot2)
max.map$CEL.LG2 <- max.map$CEL.LG
max.map$CEL.LG2[34] <- "X"

ggplot(max.map, aes(Est.Length/1e6, max.cM)) +
  stat_smooth(method = "lm") +
  geom_text(aes(label = CEL.LG2)) +
  labs(x = "Estimated Chromosome Length (MB)", y = "Linkage Map Length (cM)")

write.table(max.map, "results/6_Predicted_Physical_Size_run5_a.txt", row.names = F, sep = "\t", quote = F)
