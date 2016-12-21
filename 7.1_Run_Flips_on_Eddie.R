rm(list=ls())

library(plyr)
library(GenABEL)
library(crimaptools)
library(reshape)


#load("eddie/6_run5_Deer_Data_for_Eddie_a.RData")
load("6_run5_Deer_Data_for_Eddie_a.RData")
#reruns <- read.table("6_Reruns_2.txt", header = T)


window_size <- 100

iteration.no.improve <- 15
iteration.max <- 30

test.run.no <- NULL

flips.vals <- 2:5

Li.threshold <- 0

mapdata <- arrange(mapdata, CEL.LG, cMPosition.run5)

#outdir <- paste0("crimap/crimap_", AnalysisSuffix3)

outdir <- ""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Test on LG 10                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Creat a master data frame for a particular window size.


for(flips.val in flips.vals){

flips.tab <- NULL

for(lg in lg.vec){
  
  snplist <- skel.map$SNP.Name[which(skel.map$CEL.LG == lg)]
  
  #~~ Make a data.frame with window information.
  
  window.frame <- data.frame(Start = sort(c(seq(1                , length(snplist), window_size),
                                            seq((window_size/2)+1, length(snplist), window_size))))
  window.frame$Start <- as.integer(window.frame$Start)
  window.frame$Stop <- window.frame$Start + window_size - 1
  window.frame <- subset(window.frame, Stop <= length(snplist))
  window.frame <- rbind(window.frame, c(length(snplist) - window_size + 1, length(snplist)))
  window.frame <- unique(window.frame)
  window.frame$CEL.LG <- lg
  
  flips.tab <- rbind(flips.tab, window.frame)
  
  rm(window.frame)
}

flips.tab$window.name <- paste0(flips.tab$CEL.LG, AnalysisSuffix3,"_", flips.tab$Start, "_", window_size, "_", flips.val)
if(!is.null(test.run.no)) flips.tab <- flips.tab[1:test.run.no,]

save.image(file = paste0("6_Input_Flips_Analysis_window", window_size, "_", flips.val, ".RData"))


for(i in 1:nrow(flips.tab)){
  
  writeLines(paste0("load(\"6_Input_Flips_Analysis_window", window_size, "_", flips.val, ".RData\")\ni = ", i),
             con = paste0("flipsrun", flips.tab$window.name[i], ".R"))
  
  if(Sys.info()["sysname"] == "Windows") {
    
    system("cmd", input = paste0("copy /b flipsrun", flips.tab$window.name[i], ".R + 6.2_Run_Flips_on_Eddie_Stem.R flipsrun", flips.tab$window.name[i], ".R"))
    
  } else {
    
    system(paste0("cat 6.2_Run_Flips_on_Eddie_Stem.R >> flipsrun", flips.tab$window.name[i], ".R"))

  }
  
  writeLines(paste0("#!/bin/sh
#$ -N fl_run_", i, "_", flips.val, "              
#$ -cwd                  
#$ -l h_rt=24:00:00 
#$ -l h_vmem=2G
                    
# Initialise the environment modules
. /etc/profile.d/modules.sh
                    
# Load R
                    
module load R
                    
R CMD BATCH flipsrun", flips.tab$window.name[i], ".R"), paste0("flipsrun", flips.tab$window.name[i], ".sh"))
  
  if("reruns" %in% ls()){
    if(flips.tab$window.name[i] %in% droplevels(reruns$window.name)) system(paste0("qsub flipsrun", flips.tab$window.name[i], ".sh"))
  } else {
    system(paste0("qsub flipsrun", flips.tab$window.name[i], ".sh"))
  }
  
  
}

}