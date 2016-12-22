library(ggplot2)
library(plyr)
library(reshape)
source("r/multiplot.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Sex-averaged maps                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Run 1

run1 <- read.table("results/2_Linkage_Map_Positions_run1_a.txt", header = T)

head(run1)
run1$BTA.Chr <- paste0("BTA", run1$Chr)
run1$BTA.Chr <- factor(run1$BTA.Chr, levels = paste0("BTA", sort(unique(run1$Chr))))

ggplot(run1, aes(Position/1e6, cMPosition)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~BTA.Chr, ncol = 5, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Predicted Physical Position (Mb)",
       y = "Linkage Map Length (cM)")

ggsave("figs/LinkageMapRun1.png", width = 10, height = 14, device = "png")


#~~ Run 2

run2 <- read.table("results/2_Linkage_Map_Positions_run2_a.txt", header = T)

head(run2)
run2$BTA.Chr <- paste0("BTA", run2$Chr)
run2$BTA.Chr <- factor(run2$BTA.Chr, levels = paste0("BTA", sort(unique(run2$Chr))))

ggplot(run2, aes(Position/1e6, cMPosition)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~BTA.Chr, ncol = 5, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Predicted Physical Position (Mb)",
       y = "Linkage Map Length (cM)")

ggsave("figs/LinkageMapRun2.png", width = 10, height = 14, device = "png")


#~~ Run 3

run3 <- read.table("results/4_Linkage_Map_Positions_CEL_run3_a.txt", header = T)

head(run3)
run3$CEL.LG.lab <- paste0("CEL", run3$CEL.LG)
run3$CEL.LG.lab <- factor(run3$CEL.LG.lab, levels = paste0("CEL", sort(unique(run3$CEL.LG))))

ggplot(run3, aes(CEL.order, cMPosition.run3)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~CEL.LG.lab, ncol = 5, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "CEL Marker Order",
       y = "Linkage Map Length (cM)")

ggsave("figs/LinkageMapRun3.png", width = 10, height = 14, device = "png")

#~~ Run 4

run4 <- read.table("results/5_Linkage_Map_Positions_CEL_run4_a.txt", header = T)

head(run4)
run4$CEL.LG.lab <- paste0("CEL", run4$CEL.LG)
run4$CEL.LG.lab <- factor(run4$CEL.LG.lab, levels = paste0("CEL", sort(unique(run4$CEL.LG))))

ggplot(run4, aes(CEL.order, cMPosition.run4)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~CEL.LG.lab, ncol = 5, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "CEL Marker Order",
       y = "Linkage Map Length (cM)")

ggsave("figs/LinkageMapRun4.png", width = 10, height = 14, device = "png")



#~~ Run 5

run5 <- read.table("results/6_Linkage_Map_Positions_CEL_run5_reorient_a.txt", header = T)

head(run5)
run5$CEL.LG.lab <- paste0("CEL", run5$CEL.LG)
run5$CEL.LG.lab <- factor(run5$CEL.LG.lab, levels = paste0("CEL", sort(unique(run5$CEL.LG))))



ggplot(run5, aes(CEL.order, cMPosition.run5)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~CEL.LG.lab, ncol = 5, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "CEL Marker Order",
       y = "Linkage Map Length (cM)")

ggsave("figs/LinkageMapRun5.png", width = 10, height = 14, device = "png")

ggplot(run5, aes(BTA.Position, cMPosition.run5, col = factor(BTA.Chr))) +
  geom_point(alpha = 0.2) +
  facet_wrap(~CEL.LG.lab, ncol = 5, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "BTA Position (cM)",
       y = "Linkage Map Length (cM)",
       colour = "Cattle\nchromosome")
ggsave("figs/Run5_Compare_with_Cattle.png", width = 10, height = 14, device = "png")

ggplot(run5, aes(BTA.Position, cMPosition.run5, col = factor(CEL.LG))) +
  geom_point(alpha = 0.2) +
  facet_wrap(~BTA.Chr, ncol = 5, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "BTA Position (cM)",
       y = "Linkage Map Length (cM)")


ggplot(run5, aes(cMdiff)) +
  geom_histogram(binwidth = 0.05) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Intra-marker distance (cM)",
       y = "Frequency")

ggplot(subset(run5, cMdiff > 0), aes(cMdiff)) +
  geom_histogram(binwidth = 0.05) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Intra-marker distance (cM)",
       y = "Frequency")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Sex-specific maps                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(reshape)
sex.map <- melt(subset(run5, select = c(SNP.Name, CEL.order, CEL.LG, cMPosition.Female, cMPosition.Male)),
                id.vars = c("SNP.Name", "CEL.order", "CEL.LG"))

head(sex.map)
sex.map$variable <- gsub("cMPosition.", "", sex.map$variable)

ggplot(sex.map, aes(CEL.order, value, col = variable)) +
  geom_point() +
  facet_wrap(~CEL.LG, ncol = 5, scales = "free") +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "CEL Marker Order",
       y = "Linkage Map Length (cM)")

ggsave(filename = paste0("figs/LinkageMapRun5_sexspecific.png"), device = "png", width = 10, height = 14, units = "in")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# X chromosome rearrangements from run3       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")


x.map <- subset(run2, Chr == 30) 
head(x.map)

x.map$Diff <- c(NA, diff(x.map$cMPosition))
x.map$chunk <- NA

x.map$chunk[1] <- 1

for(i in 2:nrow(x.map)){
  x.map$chunk[i] <- ifelse(x.map$Diff[i] > 3,
                           x.map$chunk[i-1] + 1,
                           x.map$chunk[i-1])
}
x.map$SNP.Type <- ifelse(x.map$SNP.Name %in% pseudoautoSNPs, "PAR", "Sex-linked")

ggplot(x.map, aes(Order, cMPosition, col = factor(chunk))) +
  geom_point() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "BTA Marker Order",
       y = "Linkage Map Length (cM)",
       colour = "Chunk")

ggsave("figs/X_chromosome_chunks.png", width = 6, height = 5, device = "png")

x.map.hold <- x.map

#~~ How does it relate to the bovine genome?

x.map <- subset(run5, CEL.LG == 34)
head(x.map)

x.map$newbovmap <- x.map$BTA.Position/1e6

x.map2 <- melt(subset(x.map, select = c(SNP.Name, cMPosition.run5, newbovmap)), id.vars = "SNP.Name")
tail(x.map2)
x.map2$variable <- as.character(x.map2$variable)

x.map2$variable[which(x.map2$variable == "cMPosition.run5")] <- "Deer Position (cM)"
x.map2$variable[which(x.map2$variable == "newbovmap")] <- "Cattle Position (Mb)"

x.map2 <- join(x.map2, subset(x.map.hold, select = c(SNP.Name, chunk)))
x.map2$SNP.Type <- ifelse(x.map2$SNP.Name %in% pseudoautoSNPs, "PAR", "Sex-linked")

# temp <- subset(x.map2, SNP.Name %in% pseudoautoSNPs)

ggplot(x.map2, aes(variable, value, group = SNP.Name, col = factor(chunk))) +
  geom_line(alpha = 0.1) +
  geom_point() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Map",
       y = "Map Length (cM or Mb)",
       colour = "Chunk") +
  coord_cartesian(xlim = c(1.3, 1.7))
  

ggsave("figs/X_CEL_vs_BTA.png", width = 6, height = 8, device = "png")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Overall length                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

max.vals <- read.table("results/6_Predicted_Physical_Size_run5_a.txt", header = T)
max.vals$CEL.LG2 <- max.vals$CEL.LG
max.vals$CEL.LG2[34]<- "X"

ggplot(max.vals, aes(Est.Length/1e6, max.cM)) +
  stat_smooth(method = "lm") +
  geom_text(aes(label = CEL.LG2)) +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Estimated chromosome size (MB)",
       y = "Linkage Map Length (cM)")

summary(lm(Est.Length ~ max.cM, data = max.vals))

ggsave("figs/ChromosomeSizeVsMapLength", width = 6, height = 5, device = "png")

max.vals.sex <- melt(subset(max.vals, select = c(CEL.LG2, Est.Length, max.cM.male, max.cM.female)),
                     id.vars = c("CEL.LG2", "Est.Length"))
max.vals.sex$variable <- as.character(max.vals.sex$variable)
max.vals.sex$variable[which(max.vals.sex$variable == "max.cM.male")] <- "Male"
max.vals.sex$variable[which(max.vals.sex$variable == "max.cM.female")] <- "Female"

max.vals.sex.x <- subset(max.vals.sex, CEL.LG2 == "X" & variable == "Male")
max.vals.sex <- max.vals.sex[-which(max.vals.sex$CEL.LG2 == "X" & max.vals.sex$variable == "Male"),]

ggplot(max.vals.sex, aes(Est.Length/1e6, value, colour = variable)) +
  stat_smooth(method = "lm") +
  geom_text(aes(label = CEL.LG2)) +
  geom_text(data = max.vals.sex.x, aes(x = Est.Length/1e6, y = value, colour = variable, label = CEL.LG2)) +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Estimated chromosome size (MB)",
       y = "Linkage Map Length (cM)",
       colour = "Sex")

ggsave("figs/ChromosomeSizeVsMapLengthBySex", width = 6, height = 5, device = "png")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~ Linkage disequilibrium                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load("results/3_LD_matrices_before_rearrange_a.RData")

# BTA 13
i = 13
flat.test <- melt(ld.mats[[i]])
flat.test$X1.Order <- rep(1:nrow(ld.mats[[i]]), times = nrow(ld.mats[[i]]))
flat.test$X2.Order <- rep(1:nrow(ld.mats[[i]]), each = nrow(ld.mats[[i]]))
flat.test <- subset(flat.test, !is.na(value))

flat.test <- subset(flat.test, X1.Order < 500 & X2.Order < 500  & X1.Order > 0 & X2.Order > 0)

head(flat.test)
flat.test$Strip <- "(a) BTA13"

ld1 <- ggplot(flat.test, aes(X1.Order, X2.Order, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "red", high = "red", midpoint = 0.5) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank(),
        legend.position = "none") +
  facet_wrap(~Strip) +
  labs(x = "BTA13 Marker Order",
       y = "BTA13 Marker Order")

# BTA 28
i = 28
flat.test <- melt(ld.mats[[i]])
flat.test$X1.Order <- rep(1:nrow(ld.mats[[i]]), times = nrow(ld.mats[[i]]))
flat.test$X2.Order <- rep(1:nrow(ld.mats[[i]]), each = nrow(ld.mats[[i]]))
flat.test <- subset(flat.test, !is.na(value))

flat.test <- subset(flat.test, X1.Order < 400 & X2.Order < 400  & X1.Order > 0 & X2.Order > 0)

head(flat.test)
flat.test$Strip <- "(b) BTA28"

ld2 <- ggplot(flat.test, aes(X1.Order, X2.Order, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "red", high = "red", midpoint = 0.5) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank(),
        legend.position = "none") +
  facet_wrap(~Strip) +
  labs(x = "BTA28 Marker Order",
       y = "BTA28 Marker Order")



# BTA 28
i = 1
flat.test <- melt(ld.mats[[i]])
flat.test$X1.Order <- rep(1:nrow(ld.mats[[i]]), times = nrow(ld.mats[[i]]))
flat.test$X2.Order <- rep(1:nrow(ld.mats[[i]]), each = nrow(ld.mats[[i]]))
flat.test <- subset(flat.test, !is.na(value))

#flat.test <- subset(flat.test, X1.Order < 250 & X2.Order < 250  & X1.Order > 0 & X2.Order > 0)

head(flat.test)
flat.test$Strip <- "(c) BTA1"

#flat.test <- subset(flat.test,X1.Order > 926 & X2.Order > 926)

ld3 <- ggplot(flat.test, aes(X1.Order, X2.Order, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "red", high = "red", midpoint = 0.5) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank(),
        legend.position = "top") +
  facet_wrap(~Strip) +
  labs(x = "BTA1 Marker Order",
       y = "BTA1 Marker Order",
       fill = "R2")

m <- matrix(c(1, 2, 3, 3, 3, 3), ncol  = 2, byrow = TRUE)



png("figs/LD_Patterns_Inversions.png", height = 12, width = 9, units = "in", res = 300)
multiplot(ld1, ld2, ld3, cols = 1, layout = m)
dev.off()
beepr::beep()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~ Double crossovers                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dbxtab <- read.table("results/6_Double_Xovers_raw_a.txt", header = T)
head(dbxtab)

dbxtab <- subset(dbxtab, Type == "Mid")

ggplot(dbxtab, aes(SpanLength, fill = Singleton)) +
  geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = 10, linetype = "dashed") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Span distance (cM)",
       y = "Count")
  
ggsave("figs/DoubleCrossovers.png", width = 6, height = 5, device = "png")

singtab <- subset(dbxtab, Singleton == "yes")

dbx.rec <- read.table("results/6_Per_Chromosome_Recomb_dbx_a.txt", header = T, stringsAsFactors = F)
library(crimaptools)

distab <- subset(dbx.rec, select = c(Family, RRID, 
                                           parent, UniqueID, data, analysisID))
message("Splitting chromosome into segments of shared grandparental origin")
test <- mapply(switch_position_cmpstring, distab$data, distab$UniqueID, 
               SIMPLIFY = F)
switchtab <- data.frame(data.table::rbindlist(test))
suppressMessages(switchtab <- join(switchtab, distab))

library(plyr)
physical.map <- subset(run4, select = c(SNP.Name, cMPosition.run4, CEL.order, analysisID))
names(physical.map) <- c("SNP.Name", "Position", "Order", "analysisID")
physical.map$analysisID <- paste0(physical.map$analysisID, "_dbx")

switchtab <- data.frame(data.table::rbindlist(test))
suppressMessages(switchtab <- plyr::join(switchtab, distab))
switchtab <- subset(switchtab, select = -data)
if (!is.null(physical.map)) {
  suppressMessages({
    names(physical.map)[2:3] <- c("StartPos.GenomePos", 
                                  "StartPos")
    switchtab <- join(switchtab, physical.map[, 2:4])
    names(physical.map)[2:3] <- c("StopPos.GenomePos", 
                                  "StopPos")
    switchtab <- join(switchtab, physical.map[, 2:4])
    names(physical.map)[2:3] <- c("StartSpan.GenomePos", 
                                  "StartSpan")
    switchtab <- join(switchtab, physical.map[, 2:4])
    names(physical.map)[2:3] <- c("StopSpan.GenomePos", 
                                  "StopSpan")
    switchtab <- join(switchtab, physical.map[, 2:4])
  })
  switchtab$PosLength <- switchtab$StopPos.GenomePos - 
    switchtab$StartPos.GenomePos
  switchtab$SpanLength <- switchtab$StopSpan.GenomePos - 
    switchtab$StartSpan.GenomePos
}

switchtab$Singleton <- ifelse(switchtab$InfCount == 1, "yes", "no")
switchtab <- subset(switchtab, Type == "Mid")
table(switchtab$Singleton, switchtab$Type)

table(switchtab$SpanLength < 10)

ggplot(switchtab, aes(SpanLength)) + geom_histogram(binwidth = 0.5) + geom_vline(xintercept =10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# How informative are loci?               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load("data/Deer31_QC.RData", verbose = T)

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

snpinfo <- summary.snp.data(gtdata(abeldata))
snpinfo$SNP.Name <- row.names(snpinfo)
head(snpinfo)

infloci <- NULL
for(i in 1:34){
  x <- crimaptools::parse_loc(paste0("crimap/crimap_a_cel/chr", i, "a_cel_dbx.loc"))
  x$CEL.LG <- i
  infloci <- rbind(infloci, x)
  rm(x)
}

infloci <- join(infloci, snpinfo)
head(infloci)

ggplot(infloci, aes(factor(CEL.LG), inf.mei)) +
  geom_boxplot()

infmelt <- melt(subset(infloci, select = c(SNP.Name, CEL.LG, tot_m, tot_f)), id.vars = c("SNP.Name", "CEL.LG"))
infmelt$variable <- as.character(infmelt$variable)

infmelt$variable[which(infmelt$variable == "tot_m")] <- "Male"
infmelt$variable[which(infmelt$variable == "tot_f")] <- "Female"

infmelt <- infmelt[-which(infmelt$CEL.LG == 34 & infmelt$variable == "Male" & !infmelt$SNP.Name %in% pseudoautoSNPs),]



ggplot(infmelt, aes(factor(CEL.LG), value, colour = variable)) +
  geom_boxplot() +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~variable)


head(infloci)

ggplot(infloci, aes(Q.2, inf.mei)) + geom_point(alpha = 0.2) + stat_smooth()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Skeleton map                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


load("flipstest/6_run5_Deer_Data_for_Eddie_a.RData", verbose = T)

head(mapdata)
sum(table(table(mapdata$chunk)))

1 - table(mapdata$cMdiff == 0)[["FALSE"]]/sum(table(mapdata$cMdiff))
mapdata.skel <- subset(mapdata, Skeleton == 1)
table(table(mapdata.skel$chunk))

source("r/countIF.R")
head(skel.map)
skel.map$chunk.count <- countIF(skel.map$chunk)

head(subset(skel.map, chunk.count > 1))
skel.map$chunk.order <- NA
skel.map$chunk.order[1] <- 1
for(i in 2:nrow(skel.map)){
  skel.map$chunk.order[i] <- ifelse(skel.map$chunk[i] == skel.map$chunk[i-1],
                                    skel.map$chunk.order[i-1] + 1,
                                    1)
}

head(skel.map)

skel.map <- subset(skel.map, chunk.order == 1)
skel.map$chunk.order <- NULL
skel.map$chunk.count <- NULL 
skel.map$Skeleton <- NULL

head(run5)

#~~~~~~~~~~~~~~~~~ Make Final Table

finalmap <- run5
finalmap <- subset(finalmap, select = -c(Order, analysisID, Chr, chunk, cMdiff, r, cMdiff.Female, Female.r, Male.r, cMdiff.Male, CEL.LG.lab))
head(finalmap)
finalmap$Skeleton.SNP <- ifelse(finalmap$SNP.Name %in% skel.map$SNP.Name, "yes", "no")

head(infloci)
temp <- infloci
temp <- subset(temp, select = -c(CrimapOrder, CEL.LG, Chromosome, Position, Strand, P.11, P.12, P.22, NoMeasured, Pexact, Fmax, Plrt))
head(temp)

finalmap <- join(finalmap, temp)
head(finalmap)
names(finalmap)[which(names(finalmap) == "cMPosition.run5")] <- "cMPosition.SexAveraged"

head(finalmap)

finalmap$PseudoAutosomalSNP <- ifelse(finalmap$CEL.LG != 34, NA, ifelse(finalmap$SNP.Name %in% pseudoautoSNPs, "yes", "no"))

write.table(finalmap, "doc/Table_CervusElaphus_Final_Linkage_Map.txt", row.names = F, sep = "\t", quote = F)
names(finalmap)
