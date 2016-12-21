library(plyr)
library(GenABEL)
library(crimaptools)
library(reshape)
library(data.table)
library(ggplot2)

source("r/multiplot.R")
source("r/countIF.R")


load("flipstest/6_Input_Flips_Analysis_window100_5.RData")

rm(flips.vals, flips.val, i, iteration.max, iteration.no.improve, lg, lg.vec, Li.threshold, outdir, snplist, test.run.no, window_size)

flips.vec <- c(2:5)

parseFlipsRun <- TRUE

Li.threshold <- 0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Prepare for the analysis              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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


flips.tab$FinishedRunning <- NA
flips.tab$ProblemSNPCount <- NA

flips.final <- NULL
map.stats.final <- NULL
snp.order.final <- NULL   # add SNPs to this

rm(flips.temp, i)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Extract Data from the Eddie Run       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(parseFlipsRun == TRUE){
  
  for(h in 1:nrow(flips.tab)){
    
    
    
    if(h %in% seq(1, nrow(flips.tab), 20)) print(paste("Run", h, "of", nrow(flips.tab), ": Analysis ID", flips.tab$window.name[h]))
    
    
    
    if(paste0(flips.tab$window.name[h], ".Rdata") %in% dir("flipstest/")){
      
      load(paste0("flipstest/", flips.tab$window.name[h], ".Rdata"))
      
      #~~ did the analysis finish running?
      
      x <- readLines(paste0("flipstest/flipsrun", flips.tab$window.name[h], ".Rout"))
      
      if(length(grep("proc.time()", x)) > 0){
        flips.tab$FinishedRunning[h] <- "yes"
      } else {
        print(paste(h, "Did not finish running."))
        flips.tab$FinishedRunning[h] <- "no"
      }
      
      #~~ Concatenate the map stats for each window
      
      map.stats$window.name <- flips.tab$window.name[h]
      map.stats$window.name.stem <- flips.tab$window.name.stem[h]
      map.stats$flips <- flips.tab$flips[h]
      map.stats$last.run <- as.character(map.stats$last.run)
      map.stats$last.run[nrow(map.stats)] <- "yes"
      map.stats.final <- rbind(map.stats.final, map.stats)
      
      
      #~~ Make list of the SNP orders for each iteration
      
      snp.order.final <- rbind(snp.order.final, 
                               cbind(window.name = flips.tab$window.name[h],
                                     window.name.stem = flips.tab$window.name.stem[h],
                                     flips = flips.tab$flips[h],
                                     Iteration = 1:length(snp.order.list),
                                     data.frame(do.call(rbind, snp.order.list))))
      
      #~~ Concatenate the flips information if a new order was needed
      
      try({
        for(i in 1:length(flips.list)){
          flips.list[[i]] <- subset(flips.list[[i]], Li.Change < Li.threshold)
          flips.list[[i]] <- subset(flips.list[[i]], Li.Change == Min.Li)
          if(nrow(flips.list[[i]]) > 0){
            flips.list[[i]]$Iteration <- i
          } else {
            flips.list[[i]] <- cbind(flips.list[[i]], Iteration = numeric(0))
          }
        }
        
        flips.rbind <- rbindlist(flips.list)
        flips.rbind$window.name <- flips.tab$window.name[h]
        flips.rbind$window.name.stem <- flips.tab$window.name.stem[h]
        flips.rbind$flips <- flips.tab$flips[h]
        
        flips.rbind <- data.frame(flips.rbind)
        
        flips.final <- rbind(flips.final, flips.rbind)
        
        y <- NULL
        for(i in 1:ncol(flips.rbind)) if(all(flips.rbind[,i] == "-")) y <- c(y, i)
        
        flips.rbind <- flips.rbind[,-y]
        flips.melt <- melt(flips.rbind, id.vars = c("Li.Change", "Map.Li", "Min.Li", "Iteration", "window.name", "window.name.stem", "flips"))
        flips.melt$variable <- as.numeric(as.character(gsub("X", "", flips.melt$variable)))
        flips.melt$value[which(flips.melt$value == "-")] <- NA
        flips.melt$value <- as.numeric(as.character(flips.melt$value)) + 1
        flips.melt$value[which(is.na(flips.melt$value))] <- flips.melt$variable[which(is.na(flips.melt$value))]
        flips.tab$ProblemSNPCount[h] <- length(table(flips.melt$value))
        
      })
      
      rm(flips.list, flips.melt, flips.rbind, i, map.list, 
         map.stats, snp.order.list, y)
      
    } else {
      print(paste(h, "No data file found"))
      flips.tab$FinishedRunning[h] <- "no"
    }
    
    rm()
  }
  
  #~~ check the tables
  
  flips.final <- unique(flips.final)
  
  head(flips.final)
  head(flips.tab)
  head(map.stats.final)
  head(snp.order.final)
  
  #~~ did everything finish running?
  
  table(flips.tab$FinishedRunning)
  
  #~~ output a shell script to rerun anything that hung up on eddie:
  
  if(length(which(flips.tab$FinishedRunning == "no")) > 0){

    flips.tab[which(flips.tab$FinishedRunning == "no"),]
    flips.tab$window.name[which(flips.tab$window.name %in% flips.final$window.name)]

    writeLines(paste0("qsub flipsrun", flips.tab[which(flips.tab$FinishedRunning == "no"),"window.name"], ".sh"),
               "flipstest/6_rerun_delvec_error.sh")
  }
  
  #flips.tab <- subset(flips.tab, select = -FinishedRunning)
  
  flips.hold <- flips.tab
  
  write.table(flips.final    , "results/7_Full_Flips_Output_NegLi.txt", row.names = F, sep = "\t", quote = F)
  write.table(map.stats.final, "results/7_Map_Statistics_per_Flips_Run.txt", row.names = F, sep = "\t", quote = F)
  write.table(snp.order.final, "results/7_SNP_Order_per_Flips_Run.txt", row.names = F, sep = "\t", quote = F)
  
} 


