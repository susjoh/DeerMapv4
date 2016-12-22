library(ggplot2)
pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

mapdata <- read.table("results/6_Linkage_Map_Positions_CEL_run5_reorient_a.txt", header = T, stringsAsFactors = F)
head(mapdata)
mapdata$PAR <- ifelse(mapdata$SNP.Name %in% pseudoautoSNPs, "yes", "no")


mapdata$BTA.Diff <- c(diff(mapdata$BTA.Position), NA)
hist(mapdata$BTA.Diff)

mapdata$BTA.Diff <- ifelse(mapdata$BTA.Diff < 0, mapdata$BTA.Diff * -1, mapdata$BTA.Diff)
mapdata$BTA.Diff <- ifelse(is.na(mapdata$r), NA, mapdata$BTA.Diff)
hist(mapdata$BTA.Diff)

badgaps <- subset(mapdata, BTA.Diff > 1e6)
mapdata$BTA.Diff2 <- ifelse(mapdata$BTA.Diff > 1e6,
                            ifelse(mapdata$cMdiff < 0.01, 0.5e6,
                                   ifelse(mapdata$cMdiff*1e6 < mapdata$BTA.Diff, mapdata$cMdiff*1e6 , mapdata$BTA.Diff)),
                                   mapdata$BTA.Diff)

mean(tapply(mapdata$BTA.Position, mapdata$BTA.Chr, min))

hist(mapdata$BTA.Diff2)

mapdata$Dummy.Position <- NA
mapdata$Dummy.Position[1] <- mean(tapply(mapdata$BTA.Position, mapdata$BTA.Chr, min))      

for(i in 2:nrow(mapdata)){
  if(i %in% seq(1, nrow(mapdata), 10000)) print(paste("Running row", i, "of", nrow(mapdata)))
  mapdata$Dummy.Position[i] <- mapdata$Dummy.Position[i-1] + mapdata$BTA.Diff2[i-1]
  if(mapdata$CEL.LG[i] != mapdata$CEL.LG[i-1]) mapdata$Dummy.Position[i] <- mean(tapply(mapdata$BTA.Position, mapdata$BTA.Chr, min))
}

ggplot(mapdata, aes(Dummy.Position, CEL.order)) +
  geom_point() +
  facet_wrap(~CEL.LG, scales = "free")


max.map <- data.frame(max.cM = tapply(mapdata$cMPosition.run5, mapdata$CEL.LG, max),
                      max.cM.male = tapply(mapdata$cMPosition.Male, mapdata$CEL.LG, max),
                      max.cM.female = tapply(mapdata$cMPosition.Female, mapdata$CEL.LG, max),
                      Est.Length = tapply(mapdata$Dummy.Position, mapdata$CEL.LG, max),
                      LocusCount = tapply(mapdata$SNP.Name, mapdata$CEL.LG, length))

max.map$CEL.LG <- row.names(max.map)

max.map

max.map$CEL.LG2 <- max.map$CEL.LG
max.map$CEL.LG2[34] <- "X"

ggplot(max.map, aes(Est.Length/1e6, max.cM)) +
  stat_smooth(method = "lm") +
  geom_text(aes(label = CEL.LG2)) +
  labs(x = "Estimated Chromosome Length (MB)", y = "Linkage Map Length (cM)")

ggplot(max.map, aes(Est.Length/1e6, max.cM.male)) +
  stat_smooth(method = "lm") +
  geom_text(aes(label = CEL.LG2)) +
  labs(x = "Estimated Chromosome Length (MB)", y = "Linkage Map Length (cM)")

ggplot(max.map, aes(Est.Length/1e6, max.cM.female)) +
  stat_smooth(method = "lm") +
  geom_text(aes(label = CEL.LG2)) +
  labs(x = "Estimated Chromosome Length (MB)", y = "Linkage Map Length (cM)")

ggplot(max.map, aes(LocusCount, max.cM)) +
  stat_smooth(method = "lm") +
  geom_text(aes(label = CEL.LG2)) +
  labs(x = "Number of Markers", y = "Linkage Map Length (cM)")

ggplot(max.map, aes(LocusCount, max.cM.male)) +
  stat_smooth(method = "lm") +
  geom_text(aes(label = CEL.LG2)) +
  labs(x = "Number of Markers", y = "Linkage Map Length (cM)")

ggplot(max.map, aes(LocusCount, max.cM.female)) +
  stat_smooth(method = "lm") +
  geom_text(aes(label = CEL.LG2)) +
  labs(x = "Number of Markers", y = "Linkage Map Length (cM)")


write.table(max.map, "results/6_Predicted_Physical_Size_run5_reorient_a.txt", row.names = F, sep = "\t", quote = F)

write.table(mapdata, "results/6_Linkage_Map_Positions_CEL_run5_reorient_dumpos_a.txt", row.names = F, sep = "\t", quote = F)
