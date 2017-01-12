library(plyr)
library(RCircos)


map.tab <- read.table("results/8_Linkage_Map_Positions_CEL_run5_dumpos_a.txt", header = T)

head(map.tab)
map.tab$CEL.Position <- map.tab$cMPosition.run5*1e6

max.tab <- data.frame(ChromStart = 1,
                      ChromEnd = tapply(map.tab$CEL.Position, map.tab$CEL.LG, max)+2)

max.tab$Chromosome <- paste0("CEL", row.names(max.tab))

max.tab2 <- data.frame(ChromStart = 1,
                       ChromEnd = tapply(map.tab$BTA.Position, map.tab$BTA.Chr, max)+2)


max.tab2$Chromosome <- paste0("BTA", row.names(max.tab2))

unique(map.tab$BTA.Chr)

max.tab2$Chromosome <- factor(max.tab2$Chromosome, levels = paste0("BTA",unique(rev(unique(map.tab$BTA.Chr)))))
max.tab2 <- arrange(max.tab2, Chromosome)
max.tab2$Chromosome <- as.character(max.tab2$Chromosome)

max.tab <- rbind(max.tab, max.tab2)
max.tab$Band = "x"
max.tab$Stain = "gvar"

max.tab <- max.tab[,c("Chromosome", "ChromStart", "ChromEnd", "Band", "Stain")]

max.tab$Chromosome <- paste0("chr", max.tab$Chromosome)



link.data2 <- subset(map.tab, select = c(CEL.LG, CEL.Position, BTA.Chr, BTA.Position))
link.data2$CEL.Position.Stop <- link.data2$CEL.Position + 1
link.data2$BTA.Position.Stop <- link.data2$BTA.Position + 1

head(link.data2)
link.data2 <- link.data2[,c("CEL.LG", "CEL.Position", "CEL.Position.Stop", "BTA.Chr", "BTA.Position", "BTA.Position.Stop")]
link.data2$CEL.LG <- paste0("CEL", link.data2$CEL.LG)
link.data2$BTA.Chr <- paste0("BTA", link.data2$BTA.Chr)

names(link.data2) <- c("Chromosome", "chromStart", "chromEnd", "Chromosome.1", "chromStart.1", "chromEnd.1")
str(link.data2)

#link.data2 <- link.data2[sample(1:nrow(link.data2), 20),]
link.data2$Chromosome <- as.factor(link.data2$Chromosome)
link.data2$Chromosome.1 <- as.factor(link.data2$Chromosome.1)


link.data2$chromStart[which(link.data2$chromStart == 0)] <- 1
link.data2$chromEnd[which(link.data2$chromStart == 1)] <- 2


#link.data <- RCircos.Link.Data;
link.colors <- rep("blue", nrow(link.data2));
rows <- seq(1, nrow(link.data2), by=5);
link.colors[rows] <- "red";
link.data2$Number  <- as.factor(as.numeric(link.data2$Chromosome));

colourtab <- data.frame(Number = 1:34,
                        PlotColor = rep(c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                                          "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                                          "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"),  length.out = 34))
link.data2 <- join(link.data2, colourtab)

link.data2 <- subset(link.data2, select = -Number)

head(link.data2)

for(i in 35:nrow(max.tab)){
  x.chromstart.which <- which(link.data2$Chromosome.1 == gsub("chr", "", max.tab$Chromosome[i]))
  link.data2$chromStart.1[x.chromstart.which] <- max.tab$ChromEnd[i] - link.data2$chromStart.1[x.chromstart.which] -1
  link.data2$chromEnd.1[x.chromstart.which] <- max.tab$ChromEnd[i] - link.data2$chromEnd.1[x.chromstart.which] +1
  
}


track.num <- 2;

link.data3 <- link.data2[sort(sample(1:nrow(link.data2), 2000)),]


rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$chrom.paddings <- 1
RCircos.Reset.Plot.Parameters(new.params = rcircos.params)
RCircos.List.Parameters();

# pdf(file="figs/Circos_Plot_2000_sampled.pdf", height=8, width=8);

png(file="figs/Circos_Plot_2000_sampled.png", units = "in", res = 300, height=8, width=8);


RCircos.Set.Core.Components(max.tab, NULL, 0, 0)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Link.Plot(link.data3, track.num, FALSE)
dev.off()


stains <- c("gneg", "acen", "stalk", "gvar", "gpos", "gpos100", 
            "gpos75", "gpos66", "gpos50", "gpos33", "gpos25");
color.index <- c(1, 552, 615, 418, 24, 24, 193, 203, 213, 223, 233);

colors()[color.index]



# install.packages("circlize")
# library(circlize)
# bed1 = generateRandomBed(nr = 100)
# bed1 = bed1[sample(nrow(bed1), 20), ]
# bed2 = generateRandomBed(nr = 100)
# bed2 = bed2[sample(nrow(bed2), 20), ]
# 
# circos.initializeWithIdeogram(plotType = c("axis", "labels"))
# circos.genomicLink(bed1, bed2)
# circos.clear()
# circos.initializeWithIdeogram(plotType = c("axis", "labels"))
# circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5),
#                    border = NA)


