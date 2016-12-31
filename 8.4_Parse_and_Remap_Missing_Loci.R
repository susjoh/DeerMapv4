library(plyr)
library(GenABEL)
library(crimaptools)
library(reshape)
library(ggplot2)

load("results/8.3_Data_Output_for_Build6.RData")

chr.insert.vec <- sort(unique(as.numeric(long.chunk.tab$Likely.LG)))

source("runCrimapFunction.R")
results.tab <- NULL
new.marker.order <- list()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. CHROMOSOME 1                                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

lg <- 1
x <- subset(long.chunk.tab, Likely.LG == lg)
j = unique(x$chunk)
z <- subset(long.comp, chunk == j)
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1)
marker.gap <- c(z1$Map.Marker[which(z1$cMdiff == max(z1$cMdiff, na.rm = T))],
                z1$Map.Marker[which(z1$cMdiff == max(z1$cMdiff, na.rm = T)) + 1])

#~~ Try putting the chunk in the gap between markers

snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name

snpfwd <- c(snpvec[1:which(snpvec == marker.gap[1])], x$SNP.Name, snpvec[which(snpvec == marker.gap[2]):length(snpvec)])
snprev <- c(snpvec[1:which(snpvec == marker.gap[1])], rev(x$SNP.Name), snpvec[which(snpvec == marker.gap[2]):length(snpvec)])


results.tab <- rbind(results.tab,
                     runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                       abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                     runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                       abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))


new.marker.order[[length(new.marker.order) + 1]] <- snpfwd

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. CHROMOSOME 2                                                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


lg <- 2
x <- subset(long.chunk.tab, Likely.LG == lg)
j = unique(x$chunk)
z <- subset(long.comp, chunk == j)
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1)
marker.gap <- c(z1$Map.Marker[which(z1$cMdiff == max(z1$cMdiff, na.rm = T))],
                z1$Map.Marker[which(z1$cMdiff == max(z1$cMdiff, na.rm = T)) + 1])


#~~ Try putting the chunk in the gap between markers and also at the end of the chromosome.

snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name

snpfwd <- c(snpvec[1:which(snpvec == marker.gap[1])], x$SNP.Name, snpvec[which(snpvec == marker.gap[2]):length(snpvec)])
snprev <- c(snpvec[1:which(snpvec == marker.gap[1])], rev(x$SNP.Name), snpvec[which(snpvec == marker.gap[2]):length(snpvec)])
rm(snpvec)

results.tab <- rbind(results.tab,
                     runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd1",
                                       abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                     runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev1",
                                       abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))

#~~ Now try the end of the chromosome

snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name
marker.gap <- c(snpvec[length(snpvec)], NA)

snpfwd <- c(snpvec, x$SNP.Name)
snprev <- c(snpvec, rev(x$SNP.Name))
rm(snpvec)

results.tab <- rbind(results.tab,
                     runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                       abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                     runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                       abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))


plotLDcMMapRegions(j = j, direction = "fwd", markerrange = c(600, length(snpfwd)), snpvec = snpfwd, abeldata = abeldata)
plotLDcMMapRegions(j = j, direction = "rev", markerrange = c(600, length(snprev)), snpvec = snprev, abeldata = abeldata)


new.marker.order[[length(new.marker.order) + 1]] <- snpfwd


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. CHROMOSOME 6                                                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


lg <- 6
x <- subset(long.chunk.tab, Likely.LG == lg)
unique(x$chunk)

results.tab.6 <- NULL


for(j in unique(x$chunk)){
  
  z <- subset(long.comp, chunk == j)
  z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
  z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)
  
  print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
          geom_point() +
          stat_smooth(se = F) +
          ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))
  
  print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
          geom_point() +
          stat_smooth(se = F) +
          ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))
  
  print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))
  
  #~~ Put the chunks at the end of the chromosome
  
  snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name
  marker.gap <- c(snpvec[length(snpvec)], NA)
  
  snpfwd <- c(snpvec, subset(x, chunk == j)$SNP.Name)
  snprev <- c(snpvec, rev( subset(x, chunk == j)$SNP.Name))
  rm(snpvec)
  
  results.tab.6 <- rbind(results.tab.6,
                         runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                           abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                         runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                           abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))
  
  
}

#~~ Try chunks in order at the end of the chromosome.


x <- subset(long.chunk.tab, Likely.LG == lg)
unique(x$chunk)

j = "45_47_49"



z <- subset(long.comp, chunk %in% unique(x$chunk))
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

# print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
#         geom_point() +
#         stat_smooth(se = F) +
#         ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))
# 
# print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
#         geom_point() +
#         stat_smooth(se = F) +
#         ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))
# 
# print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))

#~~ Put the chunks at the end of the chromosome

snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name
marker.gap <- c(snpvec[length(snpvec)], NA)

snpfwd <- c(snpvec, subset(x, chunk %in% unique(x$chunk))$SNP.Name)
snprev <- c(snpvec, rev( subset(x, chunk %in% unique(x$chunk))$SNP.Name))
rm(snpvec)

results.tab.6 <- rbind(results.tab.6,
                       runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                         abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                       runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                         abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))



# ld.measure <- r2fast(abeldata[,snpfwd])
# ld.measure <- melt(ld.measure)
# ld.measure <- subset(ld.measure, value < 1.01)
# 
# 
# snpld <- data.frame(Order = 1:length(snpfwd),
#                     SNP.Name = snpfwd)
# names(snpld) <- c("X1.Order", "X1")
# ld.measure <- join(ld.measure, snpld)
# names(snpld) <- c("X2.Order", "X2")
# ld.measure <- join(ld.measure, snpld)
# 
# y <- parse_map_chrompic("crimap/crimap_build6/chr45_47_49chunk_fwd.cmp")
# head(y)
# y <- join(y, x[,c("SNP.Name", "chunk")])
# ggplot(subset(y, Order > 600), aes(Order, cMPosition, colour = factor(chunk))) + geom_point()
# 
# 
# ld.measure <- subset(ld.measure, X1.Order > 600 & X2.Order > 600)
# 
# print(ggplot(ld.measure, aes(X2.Order, X1.Order, fill = value)) +
#         geom_tile() +
#         scale_fill_gradient(low = "white", high = "red") +
#         geom_vline(xintercept = unique(ld.measure$X2.Order[which(ld.measure$X2 %in% c("cela1_red_6_106556351", "cela1_red_6_107749392", "cela1_red_6_108844094"))])))

results.tab <- rbind(results.tab, results.tab.6)

new.marker.order[[length(new.marker.order) + 1]] <- snpfwd


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. CHROMOSOME 9                                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

lg <- 9
x <- subset(long.chunk.tab, Likely.LG == lg)
unique(x$chunk)
x

j = unique(x$chunk)

z <- subset(long.comp, chunk == j)
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1)


#~~ try the end of the chromosome

snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name
marker.gap <- c(snpvec[length(snpvec)], NA)

snpfwd <- c(snpvec, x$SNP.Name)
snprev <- c(snpvec, rev(x$SNP.Name))
rm(snpvec)

results.tab.9<- rbind(runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                        abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                      runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                        abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))


plotLDcMMapRegions(j = j, direction = "fwd", markerrange = c(1650, length(snpfwd)), snpvec = snpfwd, abeldata = abeldata)
plotLDcMMapRegions(j = j, direction = "rev", markerrange = c(1650, length(snprev)), snpvec = snprev, abeldata = abeldata)


results.tab <- rbind(results.tab, results.tab.9)

new.marker.order[[length(new.marker.order) + 1]] <- snpfwd


beepr::beep()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. CHROMOSOME 14                                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


lg <- 14
x <- subset(long.chunk.tab, Likely.LG == lg)
unique(x$chunk)
x


results.tab.14 <- NULL


for(j in unique(x$chunk)){
  
  z <- subset(long.comp, chunk == j)
  z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
  z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)
  
  print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
          geom_point() +
          stat_smooth(se = F) +
          ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))
  
  print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
          geom_point() +
          stat_smooth(se = F) +
          ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))
  
  print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))
  
}

#~~ Put the chunks at the end of the chromosome: 111 at the start, 112 at the end.

j = 111
snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name
marker.gap <- c(NA, snpvec[length(snpvec)])

snpfwd <- c(subset(x, chunk == j)$SNP.Name, snpvec)
snprev <- c(rev( subset(x, chunk == j)$SNP.Name), snpvec)
rm(snpvec)

results.tab.14 <- rbind(results.tab.14,
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))


plotLDcMMapRegions(j = j, direction = "fwd", markerrange = c(1, 50), snpvec = snpfwd, abeldata = abeldata)
plotLDcMMapRegions(j = j, direction = "rev", markerrange = c(1, 50), snpvec = snprev, abeldata = abeldata)

j = 112
snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name
marker.gap <- c(snpvec[length(snpvec)], NA)

snpfwd <- c(snpvec, subset(x, chunk == j)$SNP.Name)
snprev <- c(snpvec, rev( subset(x, chunk == j)$SNP.Name))
rm(snpvec)

results.tab.14 <- rbind(results.tab.14,
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))


plotLDcMMapRegions(j = j, direction = "fwd", markerrange = c(1080, length(snpfwd)), snpvec = snpfwd, abeldata = abeldata)
plotLDcMMapRegions(j = j, direction = "rev", markerrange = c(1080, length(snprev)), snpvec = snprev, abeldata = abeldata)

# ACCEPT FORWARD

results.tab <- rbind(results.tab, results.tab.14)

snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name
snpvec <- c(subset(x, chunk == 111)$SNP.Name, snpvec, subset(x, chunk == 112)$SNP.Name)

new.marker.order[[length(new.marker.order) + 1]] <- snpvec


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. CHROMOSOME 15                                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

chr.insert.vec

lg <- 15
x <- subset(long.chunk.tab, Likely.LG == lg)
j = unique(x$chunk)
x

#~~ separate into finer chunks

x$chunk2 <- c(x$chunk[1:10] + 0.1, x$chunk[11:15] + 0.2, x$chunk[16:17] + 0.3)


results.tab.15 <- NULL

for(j in unique(x$chunk2)){
  z <- subset(long.comp, Unmap.Marker %in% subset(x, chunk2 == j)$SNP.Name)
  z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
  z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)
  
  print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
          geom_point() +
          stat_smooth(se = F) +
          ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))
  
  print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
          geom_point() +
          stat_smooth(se = F) +
          ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))
  
  print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))
}


#~~ remove last two SNPs and Try putting the chunk in the gap between markers
x <- x[1:14,]
j = unique(x$chunk)

marker.gap <- c(z1$Map.Marker[which(z1$cMdiff == max(z1$cMdiff, na.rm = T))],
                z1$Map.Marker[which(z1$cMdiff == max(z1$cMdiff, na.rm = T)) + 1])

snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name

snpfwd <- c(snpvec[1:which(snpvec == marker.gap[1])], x$SNP.Name, snpvec[which(snpvec == marker.gap[2]):length(snpvec)])
snprev <- c(snpvec[1:which(snpvec == marker.gap[1])], rev(x$SNP.Name), snpvec[which(snpvec == marker.gap[2]):length(snpvec)])


results.tab.15 <- rbind(results.tab.15,
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))

plotLDcMMapRegions(j = j, direction = "fwd", markerrange = c(150, 230), snpvec = snpfwd, abeldata = abeldata)
plotLDcMMapRegions(j = j, direction = "rev", markerrange = c(150, 230), snpvec = snprev, abeldata = abeldata)

results.tab <- rbind(results.tab, results.tab.15)

new.marker.order[[length(new.marker.order) + 1]] <- snprev



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 7. CHROMOSOME 19                                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

chr.insert.vec

lg <- 19
x <- subset(long.chunk.tab, Likely.LG == lg)
unique(x$chunk)
subset(x, select = -c(PosDiff, PosDiffb4, cMdiffb4, PosStretch))

#~~ separate into finer chunks

results.tab.19 <- NULL

#~~~~ First chunk #########################################

j = 2
z <- subset(long.comp, Unmap.Marker %in% subset(x, chunk == j)$SNP.Name)
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))

#~~ remove last two SNPs and Try putting the chunk in the gap between markers

marker.gap <- c(z1$Map.Marker[which(z1$cMdiff == max(z1$cMdiff, na.rm = T))],
                z1$Map.Marker[which(z1$cMdiff == max(z1$cMdiff, na.rm = T)) + 1])

snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name

snpfwd <- c(snpvec[1:which(snpvec == marker.gap[1])], subset(x, chunk == j)$SNP.Name, snpvec[which(snpvec == marker.gap[2]):length(snpvec)])
snprev <- c(snpvec[1:which(snpvec == marker.gap[1])], rev(subset(x, chunk == j)$SNP.Name), snpvec[which(snpvec == marker.gap[2]):length(snpvec)])


results.tab.19 <- rbind(results.tab.19,
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))

plotLDcMMapRegions(j = j, direction = "fwd", markerrange = c(500, 600), snpvec = snpfwd, abeldata = abeldata)
plotLDcMMapRegions(j = j, direction = "rev", markerrange = c(500, 600), snpvec = snprev, abeldata = abeldata)

subset(all.chunk.summ, CEL.LG == 19)

#~~~~ Second chunk #########################################

unique(x$chunk)
j = 10
z <- subset(long.comp, Unmap.Marker %in% subset(x, chunk == j)$SNP.Name)
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))

#~~~~ Third chunk #########################################

unique(x$chunk)
j = 18
z <- subset(long.comp, Unmap.Marker %in% subset(x, chunk == j)$SNP.Name)
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))

new.marker.order[[length(new.marker.order) + 1]] <- subset(mapdata, CEL.LG == lg)$SNP.Name

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 8. CHROMOSOME 20                                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

chr.insert.vec

lg <- 20
x <- subset(long.chunk.tab, Likely.LG == lg)
unique(x$chunk)
subset(x, select = -c(PosDiff, PosDiffb4, cMdiffb4, PosStretch))


results.tab.20 <- NULL

j = unique(x$chunk)
z <- subset(long.comp, Unmap.Marker %in% subset(x, chunk == j)$SNP.Name)
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))


snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name

marker.gap <- c(NA, snpvec[1])

snpfwd <- c(subset(x, chunk == j)$SNP.Name, snpvec)
snprev <- c(rev( subset(x, chunk == j)$SNP.Name), snpvec)


results.tab.20 <- rbind(results.tab.20,
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))

plotLDcMMapRegions(j = j, direction = "fwd", markerrange = c(1, 50), snpvec = snpfwd, abeldata = abeldata)
plotLDcMMapRegions(j = j, direction = "rev", markerrange = c(1, 50), snpvec = snprev, abeldata = abeldata)

results.tab <- rbind(results.tab, results.tab.20)
subset(all.chunk.summ, CEL.LG == 20)

new.marker.order[[length(new.marker.order) + 1]] <- snpfwd


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 9. CHROMOSOME 22                                                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

chr.insert.vec

lg <- 22
x <- subset(long.chunk.tab, Likely.LG == lg)
unique(x$chunk)
subset(x, select = -c(PosDiff, PosDiffb4, cMdiffb4, PosStretch))

#~~ separate into finer chunks

results.tab.22 <- NULL

j = unique(x$chunk)
z <- subset(long.comp, Unmap.Marker %in% subset(x, chunk == j)$SNP.Name)
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))

#~~ Put at start

snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name

marker.gap <- c(NA, snpvec[1])

snpfwd <- c(subset(x, chunk == j)$SNP.Name, snpvec)
snprev <- c(rev( subset(x, chunk == j)$SNP.Name), snpvec)


results.tab.22 <- rbind(results.tab.22,
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))

plotLDcMMapRegions(j = j, direction = "fwd", markerrange = c(1, 50), snpvec = snpfwd, abeldata = abeldata)
plotLDcMMapRegions(j = j, direction = "rev", markerrange = c(1, 50), snpvec = snprev, abeldata = abeldata)

results.tab <- rbind(results.tab, results.tab.22)
subset(all.chunk.summ, CEL.LG == 22)

new.marker.order[[length(new.marker.order) + 1]] <- snpfwd


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 10. CHROMOSOME 32                                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

chr.insert.vec

lg <- 32
x <- subset(long.chunk.tab, Likely.LG == lg)
unique(x$chunk)
subset(x, select = -c(PosDiff, PosDiffb4, cMdiffb4, PosStretch))

#~~ separate into finer chunks

results.tab.32 <- NULL

j = unique(x$chunk)
z <- subset(long.comp, Unmap.Marker %in% subset(x, chunk == j)$SNP.Name)
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))

#~~ Put at start

snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name

marker.gap <- c(NA, snpvec[1])

snpfwd <- c(subset(x, chunk == j)$SNP.Name, snpvec)
snprev <- c(rev( subset(x, chunk == j)$SNP.Name), snpvec)


results.tab.32 <- rbind(results.tab.32,
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))

plotLDcMMapRegions(j = j, direction = "fwd", markerrange = c(1, 50), snpvec = snpfwd, abeldata = abeldata)
plotLDcMMapRegions(j = j, direction = "rev", markerrange = c(1, 50), snpvec = snprev, abeldata = abeldata)

results.tab <- rbind(results.tab, results.tab.32)

new.marker.order[[length(new.marker.order) + 1]] <- snprev

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 11. CHROMOSOME 34                                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

chr.insert.vec

lg <- 34
x <- subset(long.chunk.tab, Likely.LG == lg)
unique(x$chunk)
subset(x, select = -c(PosDiff, PosDiffb4, cMdiffb4, PosStretch))
x$PAR <- ifelse(x$SNP.Name %in% pseudoautoSNPs, "yes", "no")
results.tab.34 <- NULL


j = 176
subset(x, chunk == j)
x <- subset(x, SNP.Name != "cela1_red_x_57663610") # remove a dodgy SNP

z <- subset(long.comp, Unmap.Marker %in% subset(x, chunk == j)$SNP.Name)
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))


#~~ Put in gap

snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name

marker.gap <- c(z1$Map.Marker[which(z1$cMdiff == max(z1$cMdiff, na.rm = T))],
                z1$Map.Marker[which(z1$cMdiff == max(z1$cMdiff, na.rm = T)) + 1])

snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name

snpfwd <- c(snpvec[1:which(snpvec == marker.gap[1])], subset(x, chunk == j)$SNP.Name, snpvec[which(snpvec == marker.gap[2]):length(snpvec)])
snprev <- c(snpvec[1:which(snpvec == marker.gap[1])], rev(subset(x, chunk == j)$SNP.Name), snpvec[which(snpvec == marker.gap[2]):length(snpvec)])



results.tab.34 <- rbind(results.tab.34,
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))

plotLDcMMapRegions(j = j, direction = "fwd", markerrange = c(1, 100), snpvec = snpfwd, abeldata = abeldata.fem)
plotLDcMMapRegions(j = j, direction = "rev", markerrange = c(1, 100), snpvec = snprev, abeldata = abeldata.fem) # REVERSE


#~~ NExt chunk
unique(x$chunk)
j = 182
subset(x, chunk == j)


z <- subset(long.comp, Unmap.Marker %in% subset(x, chunk == j)$SNP.Name)
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))


temp <- read.table("results/4_Linkage_Map_Positions_CEL_run3_a.txt", header = T)
which(temp$SNP.Name %in% subset(x, chunk == j)$SNP.Name)
temp[36275:36290,]

#~~ remove SNP cela1_sika_x_68995376 from mapdata
mapdata.34 <- subset(mapdata, !SNP.Name %in% "cela1_sika_x_68995376" & CEL.LG == 34)

#~~ Put chunk between cela1_red_x_68207800 & cela1_red_x_1474252   
subset(x, chunk == j)

snpvec <- mapdata.34$SNP.Name

marker.gap <- c("cela1_red_x_68207800", "cela1_red_x_1474252")

snpfwd <- c(snpvec[1:which(snpvec == marker.gap[1])], subset(x, chunk == j)$SNP.Name, snpvec[which(snpvec == marker.gap[2]):length(snpvec)])
snprev <- c(snpvec[1:which(snpvec == marker.gap[1])], rev(subset(x, chunk == j)$SNP.Name), snpvec[which(snpvec == marker.gap[2]):length(snpvec)])



results.tab.34 <- rbind(results.tab.34,
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))

plotLDcMMapRegions(j = j, direction = "fwd", markerrange = c(150, 200), snpvec = snpfwd, abeldata = abeldata.fem)
plotLDcMMapRegions(j = j, direction = "rev", markerrange = c(150, 200), snpvec = snprev, abeldata = abeldata.fem) # FORWARD



#~~ NExt chunk
unique(x$chunk)

j = 183
subset(x, chunk == j)
subset(mapdata_run1, chunk %in% 183:185)

z <- subset(long.comp, Unmap.Marker %in% subset(x, chunk == j)$SNP.Name)
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F, alpha = 0.1) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))


#~~ TRy between cela1_red_x_41783448 & cela1_red_x_72708173

snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name

marker.gap <- c("cela1_red_x_41783448", "cela1_red_x_72708173")


snpfwd <- c(snpvec[1:which(snpvec == marker.gap[1])], subset(x, chunk == j)$SNP.Name, snpvec[which(snpvec == marker.gap[2]):length(snpvec)])
snprev <- c(snpvec[1:which(snpvec == marker.gap[1])], rev(subset(x, chunk == j)$SNP.Name), snpvec[which(snpvec == marker.gap[2]):length(snpvec)])



results.tab.34 <- rbind(results.tab.34,
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))

plotLDcMMapRegions(j = j, direction = "fwd", markerrange = c(800, 900), snpvec = snpfwd, abeldata = abeldata.fem) # FORWARD
plotLDcMMapRegions(j = j, direction = "rev", markerrange = c(800, 900), snpvec = snprev, abeldata = abeldata.fem) 

unique(x$chunk)

#~~ Next chunk

j = 184
subset(x, chunk == j)
subset(mapdata_run1, chunk %in% 183:185)
subset(mapdata.34, cMPosition.run5 > 17 & cMPosition.run5 < 19.5)

z <- subset(long.comp, Unmap.Marker %in% subset(x, chunk == j)$SNP.Name)
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F, alpha = 0.1) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))


#~~ TRy between cela1_red_x_68207800 & cela1_red_x_1474252

snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name

marker.gap <- c("cela1_red_x_68207800", "cela1_red_x_1474252")


snpfwd <- c(snpvec[1:which(snpvec == marker.gap[1])], subset(x, chunk == j)$SNP.Name, snpvec[which(snpvec == marker.gap[2]):length(snpvec)])
snprev <- c(snpvec[1:which(snpvec == marker.gap[1])], rev(subset(x, chunk == j)$SNP.Name), snpvec[which(snpvec == marker.gap[2]):length(snpvec)])



results.tab.34 <- rbind(results.tab.34,
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))

plotLDcMMapRegions(j = j, direction = "fwd", markerrange = c(125, 225), snpvec = snpfwd, abeldata = abeldata.fem) # FORWARD
plotLDcMMapRegions(j = j, direction = "rev", markerrange = c(125, 225), snpvec = snprev, abeldata = abeldata.fem) 


#~~ Next chunk
unique(x$chunk)


j = 189
subset(x, chunk == j)
subset(mapdata_run1, chunk %in% 189:191)
mapdata.34$PAR <- ifelse(mapdata.34$SNP.Name %in% pseudoautoSNPs, "yes", "no")

#~~ remove the non-PAR SNPs
mapdata.34 <- subset(mapdata.34, !SNP.Name %in% x[which(x$PAR == "no" & x$chunk == j),"SNP.Name"])
x <- x[-which(x$PAR == "no" & x$chunk == j),]

subset(x, chunk == j)


z <- subset(long.comp, Unmap.Marker %in% subset(x, chunk == j)$SNP.Name)
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F, alpha = 0.1) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))


#~~ TRy at end

snpvec <- subset(mapdata.34, CEL.LG == lg)$SNP.Name
marker.gap <- c(snpvec[length(snpvec)], NA)

snpfwd <- c(snpvec, subset(x, chunk == j)$SNP.Name)
snprev <- c(snpvec, rev( subset(x, chunk == j)$SNP.Name))
rm(snpvec)


results.tab.34 <- rbind(results.tab.34,
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))

plotLDcMMapRegions(j = j, direction = "fwd", markerrange = c(1750, length(snpfwd)), snpvec = snpfwd, abeldata = abeldata.fem) # FORWARD
plotLDcMMapRegions(j = j, direction = "rev", markerrange = c(1750, length(snprev)), snpvec = snprev, abeldata = abeldata.fem) 



#~~ Next chunk
unique(x$chunk)


j = 190
subset(x, chunk == j)

subset(x, chunk == j)$SNP.Name

#~~ remove the PAR SNPs

x <- x[-which(x$PAR == "yes" & x$chunk == j),]

z <- subset(long.comp, Unmap.Marker %in% subset(x, chunk == j)$SNP.Name)
z1 <- unique(subset(z, select = c(Map.Marker, cMPosition.run5)))
z1$cMdiff <- c(diff(z1$cMPosition.run5), NA)

print(ggplot(z, aes(Dummy.Position, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z, aes(cMPosition.run5, R2, colour = Unmap.Marker)) +
        geom_point() +
        stat_smooth(se = F, alpha = 0.1) +
        ggtitle(paste("Chunk", j, "with", length(unique(z$Unmap.Marker)), "SNPs on linkage group", z$Likely.LG[1])))

print(ggplot(z1, aes(cMdiff)) + geom_histogram(binwidth = 0.1))


#~~ TRy between  cela1_red_x_143680253 and cela1_red_x_144801845     

snpvec <- subset(mapdata, CEL.LG == lg)$SNP.Name

marker.gap <- c("cela1_red_x_143680253", "cela1_red_x_144801845")


snpfwd <- c(snpvec[1:which(snpvec == marker.gap[1])], subset(x, chunk == j)$SNP.Name, snpvec[which(snpvec == marker.gap[2]):length(snpvec)])
snprev <- c(snpvec[1:which(snpvec == marker.gap[1])], rev(subset(x, chunk == j)$SNP.Name), snpvec[which(snpvec == marker.gap[2]):length(snpvec)])


results.tab.34 <- rbind(results.tab.34,
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snpfwd, direction = "fwd",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F),
                        runCrimapFunction(j = j, lg = lg, marker.gap = marker.gap, snpvec = snprev, direction = "rev",
                                          abeldata = abeldata, famped = famped, pseudoautoSNPs = pseudoautoSNPs, runCrimap = F))

plotLDcMMapRegions(j = j, direction = "fwd", markerrange = c(1700, 1800), snpvec = snpfwd, abeldata = abeldata.fem) # FORWARD
plotLDcMMapRegions(j = j, direction = "rev", markerrange = c(1700, 1800), snpvec = snprev, abeldata = abeldata.fem) 

head(z)

results.tab <- rbind(results.tab, results.tab.34)

temptab <- results.tab.34[c(2, 3, 5, 7, 9, 11),]
tempmarker <- subset(mapdata.34, select = c(SNP.Name, cMPosition.run5))
names(tempmarker)[1] <- "MarkerLeft"
temptab <- join(temptab, tempmarker)
temptab
rm(tempmarker)

temptab <- arrange(temptab, cMPosition.run5)



snpvec <- mapdata.34$SNP.Name
snpvec2 <- 
  c(snpvec[1:which(snpvec == temptab$MarkerLeft[1])],
    rev(subset(x, chunk == temptab$chunk[1])$SNP.Name),
    snpvec[which(snpvec == temptab$MarkerRight[1]):which(snpvec == temptab$MarkerLeft[2])],
    subset(x, chunk == temptab$chunk[2])$SNP.Name,
    subset(x, chunk == temptab$chunk[3])$SNP.Name,
    snpvec[which(snpvec == temptab$MarkerRight[3]):which(snpvec == temptab$MarkerLeft[4])],
    subset(x, chunk == temptab$chunk[4])$SNP.Name,
    snpvec[which(snpvec == temptab$MarkerRight[4]):which(snpvec == temptab$MarkerLeft[5])],
    #subset(x, chunk == temptab$chunk[5])$SNP.Name,
    snpvec[which(snpvec == temptab$MarkerRight[5]):which(snpvec == temptab$MarkerLeft[6])],
    subset(x, chunk == temptab$chunk[6])$SNP.Name)
    
    length(unique(snpvec2))
    


#new.marker.order[[length(new.marker.order) + 1]] <- snpvec2

#~~ check X

j = 34
lg = 34
direction = "RERUN"
# create_crimap_input(gwaa.data = abeldata,
#                     familyPedigree = famped,
#                     analysisID = paste0(j, "chunk_", direction),
#                     snplist = snpvec2,
#                     is.X = ifelse(lg == 34, TRUE, FALSE),
#                     pseudoautoSNPs = pseudoautoSNPs,
#                     outdir = paste0("crimap/crimap_build6"),
#                     clear.existing.analysisID = TRUE)
# 
# run_crimap_prepare(genfile = paste0("crimap/crimap_build6/chr", j, "chunk_", direction, ".gen"))
# 
# parse_mend_err(prefile = paste0("crimap/crimap_build6/chr", j, "chunk_", direction, ".pre"),
#                genfile = paste0("crimap/crimap_build6/chr", j, "chunk_", direction, ".gen"),
#                familyPedigree = famped,
#                is.X = TRUE,
#                pseudoautoSNPs = pseudoautoSNPs,
#                genabel.phdata = phdata(abeldata),
#                save.mendfile = TRUE)

create_crimap_input(gwaa.data = abeldata,
                    familyPedigree = famped,
                    analysisID = paste0(j, "chunk_", direction),
                    snplist = snpvec2,
                    is.X = ifelse(lg == 34, TRUE, FALSE),
                    pseudoautoSNPs = pseudoautoSNPs,
                    outdir = paste0("crimap/crimap_build6"),
                    use.mnd = TRUE,
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = paste0("crimap/crimap_build6/chr", j, "chunk_", direction, ".gen"))

run_crimap_chrompic(genfile = paste0("crimap/crimap_build6/chr", j, "chunk_", direction, ".gen"))

y <- parse_map_chrompic(paste0("crimap/crimap_build6/chr", j, "chunk_", direction, ".cmp"))
ggplot(y, aes(Order, cMPosition)) + geom_point()

new.marker.order[[length(new.marker.order) + 1]] <- snpvec2




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CREATE NEW MAP ORDERS                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

newmaporder <- NULL

for(i in 1:34){
  
  if(i %in% chr.insert.vec){
    
    x <- data.frame(SNP.Name = new.marker.order[[which(chr.insert.vec == i)]],
                    CEL.LG = i)
    x$CEL.order <- 1:nrow(x)
    
  } else {
    x <- subset(mapdata, CEL.LG == i)
    x <- subset(x, select = c(SNP.Name, CEL.LG, CEL.order))
    
  }
  
  newmaporder <- rbind(newmaporder, x)
  
}

save(results.tab, temptab, new.marker.order, newmaporder, file = "results/8_Remap_missing_loci_new_orders.RData")


