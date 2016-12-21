#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set up working environment and load in data    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(plyr)
library(GenABEL)
library(crimaptools)
library(reshape)
library(data.table)
library(ggplot2)
library(beepr)

source("r/multiplot.R")
source("r/countIF.R")

parse_map_li <- function(mapfile){
  foo <- system(paste0("grep log10_like ", mapfile), intern = T)
  as.numeric(strsplit(foo[length(foo)], split = " ")[[1]][3])
}

load("flipstest/6_Input_Flips_Analysis_window100_5.RData")

rm(flips.vals, flips.val, i, iteration.max, iteration.no.improve, lg, lg.vec, Li.threshold, outdir, snplist, test.run.no, window_size)

flips.vec <- c(2:5)


skel.map <- arrange(skel.map, CEL.LG, cMPosition.run5)
skel.map$Skeleton.Order <- 1
for(i in 2:nrow(skel.map)) skel.map$Skeleton.Order[i] <- ifelse(skel.map$CEL.LG[i] == skel.map$CEL.LG[i-1],
                                                                skel.map$Skeleton.Order[i-1] + 1,
                                                                1)

flips.temp <- flips.tab
flips.tab <- NULL
for(i in 1:length(flips.vec)){
  flips.tab <- rbind(flips.tab, cbind(flips.temp, flips = flips.vec[i]))
}

flips.tab$window.name.stem <- sapply(flips.tab$window.name, function(x) substr(x, 1, nchar(x)-2))
flips.tab$window.name <- paste0(flips.tab$window.name.stem, "_", flips.tab$flips)

flips.final     <- read.table("results/7_Full_Flips_Output_NegLi.txt", header = T)
map.stats.final <- read.table("results/7_Map_Statistics_per_Flips_Run.txt", header = T)
snp.order.final <- read.table("results/7_SNP_Order_per_Flips_Run.txt", header = T)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Extract problem regions               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

map.threshold <- -1

map.changes.info <- NULL
map.changes.flips <- NULL
map.changes.snps <- NULL

#~~ Find the difference between the largest and smallest map length and likelihood changes (map.changes.info)
#   and pull out the flips regions (map.changes.flips)

for(h in 1:nrow(flips.tab)){
  x <- subset(map.stats.final, window.name == flips.tab$window.name[h])
  x.passli <- which(x$Min.Li <= map.threshold)
  
  if(length(x.passli) > 0){
    
    x <- x[1:(max(x.passli) + 1),]
    
    map.changes.info <- rbind(map.changes.info,
                              data.frame(window.name = flips.tab$window.name[h],
                                         window.name.stem = flips.tab$window.name.stem[h],
                                         Li.Change = max(x$Min.Li) - min(x$Min.Li),
                                         Map.Change = max(x$Map.Len) - min(x$Map.Len)))
    
    map.changes.flips <- rbind(map.changes.flips,
                               subset(flips.final, window.name == flips.tab$window.name[h] & Min.Li <= map.threshold))
    
  }
  
  rm(x, x.passli)
  
}

#~~ Pull out the specific SNPs that were in the changed order part of the flips, and which position they were in.
# (map.changes.snps)

for(i in as.character(unique(map.changes.info$window.name.stem))){
  
  x <- droplevels(subset(map.changes.flips, window.name.stem == i))
  
  x <- melt(x, id.vars = c("Li.Change", "Map.Li", "Min.Li", "Iteration", "window.name", "window.name.stem", "flips"))
  head(x)
  
  x <- subset(x, value != "-")
  x$SNP.Order <- as.numeric(as.character(x$value)) + 1
  
  snp.list <- droplevels(subset(snp.order.final, window.name.stem == i))
  
  snp.list <- melt(snp.list, id.vars = c("window.name", "window.name.stem", "flips", "Iteration"))
  snp.list$variable <- as.numeric(gsub("X", "", snp.list$variable))
  names(snp.list)[which(names(snp.list) == "variable")] <- "SNP.Order"
  names(snp.list)[which(names(snp.list) == "value")] <- "SNP.Name"
  
  
  suppressMessages(x <- droplevels(join(x, snp.list)))
  
  table(x$SNP.Name)
  
  map.changes.snps <- rbind(map.changes.snps, x)
  
  rm(x, snp.list)
  
  
}

#~~ What are these SNPs? Group them into regions

snps.to.change <- data.frame(table(map.changes.snps$SNP.Name))
names(snps.to.change)[1] <- "SNP.Name"

snps.to.change$SNP.Name <- as.character(snps.to.change$SNP.Name)

snps.to.change <- join(snps.to.change, skel.map)
snps.to.change <- arrange(snps.to.change, CEL.LG, cMPosition.run5)

head(snps.to.change)

snps.to.change$SNP.Group <- 1

for(i in 2:nrow(snps.to.change)) snps.to.change$SNP.Group[i] <- ifelse(snps.to.change$CEL.LG[i-1] != snps.to.change$CEL.LG[i],
                                                                       snps.to.change$SNP.Group[i-1] + 1,
                                                                       ifelse((snps.to.change$Skeleton.Order[i] - snps.to.change$Skeleton.Order[i-1]) < 5,
                                                                              snps.to.change$SNP.Group[i-1],
                                                                              snps.to.change$SNP.Group[i-1] + 1))

table(snps.to.change$SNP.Group)

map.changes.info$CEL.LG <- sapply(map.changes.info$window.name.stem, function(x) as.numeric(gsub("a", "", strsplit(as.character(x), split = "_")[[1]][1])))

str(map.changes.info)
map.changes.info$window.name      <- as.character(map.changes.info$window.name)
map.changes.info$window.name.stem <- as.character(map.changes.info$window.name.stem)

map.changes.info <- arrange(map.changes.info, CEL.LG, window.name)

ggplot(map.changes.info, aes(Map.Change, Li.Change)) + geom_point() + stat_smooth()

skel.map$Confidence <- ifelse(as.character(skel.map$SNP.Name) %in% as.character(snps.to.change$SNP.Name), "low", "high")


save(snps.to.change, skel.map, abeldata, famped, flips.tab, file = "flipstest/7.4_data_for_flips_test.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Find proper order                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



system("cmd", input = "mkdir flipstest")

#~~ extract map segment for that section, and remove any other problem regions that are present.

for(i in 1:max(snps.to.change$SNP.Group)){
  
  print(paste("Running group", i, "of", max(snps.to.change$SNP.Group)))
  
  x <- subset(snps.to.change, SNP.Group == i)
  x.map <- subset(skel.map, CEL.LG == x$CEL.LG[1] & Skeleton.Order > (min(x$Skeleton.Order) - 20) & Skeleton.Order < (max(x$Skeleton.Order) + 20))
  x.omitted <- x.map$SNP.Name[which(x.map$Confidence == "low" & !x.map$SNP.Name %in% x$SNP.Name)]
  x.map <- subset(x.map, !SNP.Name %in% x.omitted)
  
  #~~ create input, run chrompic and flips5
  
  create_crimap_input(abeldata,
                      famped,
                      analysisID = paste0(i, "gptest"),
                      snplist = x.map$SNP.Name,
                      outdir = "flipstest",
                      clear.existing.analysisID = TRUE) 
  
  run_crimap_prepare(paste0("flipstest/chr", i, "gptest.gen"))
  run_crimap_chrompic(paste0("flipstest/chr", i, "gptest.gen"))
  x.cmp <- parse_map_chrompic(paste0("flipstest/chr", i, "gptest.cmp"))
  parse_map_li(paste0("flipstest/chr", i, "gptest.cmp"))
  
  system.time(run_crimap_flips(paste0("flipstest/chr", i, "gptest.gen"), 5))
  beep()
  
}



for(i in 1:max(snps.to.change$SNP.Group)){
  
  x.flips <- parse_flips(paste0("flipstest/chr", i, "gptest.fl5"))[-1,]
  
  head(x.flips)
  
  x.flips <- subset(x.flips, LogLi < 0)                     
  x.flips <- arrange(x.flips, LogLi)
  
  write.table(x.flips, paste0("flipstest/chr", i, "gptest.fl5ordered"), row.names = F, quote = F)
  
  ggplot(x.flips, aes(LogLi)) + geom_histogram()
  
}


# for(i in 1:max(snps.to.change$SNP.Group)){
#   
#   x <- subset(snps.to.change, SNP.Group == i)
#   x.map <- subset(skel.map, CEL.LG == x$CEL.LG[1] & Skeleton.Order > (min(x$Skeleton.Order) - 20) & Skeleton.Order < (max(x$Skeleton.Order) + 20))
#   x.omitted <- x.map$SNP.Name[which(x.map$Confidence == "low" & !x.map$SNP.Name %in% x$SNP.Name)]
#   x.map <- subset(x.map, !SNP.Name %in% x.omitted)
#   
#   snp.1 <- x.map$SNP.Name[which(!x.map$SNP.Name %in% x$SNP.Name)]
#   snp.2 <- x$SNP.Name
#   
#   create_crimap_input(abeldata,
#                       famped,
#                       analysisID = paste0(1, "gptest"),
#                       snplist = x.map$SNP.Name,
#                       outdir = "flipstest",
#                       clear.existing.analysisID = TRUE) 
#   
#   run_crimap_prepare(paste0("flipstest/chr", 1, "gptest.gen"), snplist = snp.1, snpinsert = snp.2)
#   
#   system.time(run_crimap_build(paste0("flipstest/chr", 1, "gptest.gen")))
#   
#   parse_map_build(paste0("flipstest/chr", 1, "gptest.bld1"))
#   
# }



