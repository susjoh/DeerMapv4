rm(list=ls())

library(plyr)
library(GenABEL)
library(crimaptools)
library(reshape)


load("5_Recheck_Data_Fine_Scale_Chunk_Diffs.RData")

#for(i in 1:nrow(map.chunk.50)){

for(i in 93:109){  
  writeLines(paste0("load(\"5_Recheck_Data_Fine_Scale_Chunk_Diffs.RData\")\ni = ", i),
             con = paste0("mapchunk50_", map.chunk.50$chunk[i], ".R"))
  
      
    system(paste0("cat 5.4_Fine_Scale_Differences_Stem.R >> mapchunk50_Re_", map.chunk.50$chunk[i], ".R"))
      
    
    writeLines(paste0("#!/bin/sh
#$ -N fl_run_", map.chunk.50$chunk[i], "              
#$ -cwd                  
#$ -l h_rt=24:00:00 
#$ -l h_vmem=2G
                      
# Initialise the environment modules
. /etc/profile.d/modules.sh
                      
# Load R
                      
module load R
                      
R CMD BATCH mapchunk50_Re_", map.chunk.50$chunk[i], ".R"), paste0("mapchunk50_Re_", map.chunk.50$chunk[i], ".sh"))
    
system(paste0("qsub mapchunk50_Re_", map.chunk.50$chunk[i], ".sh"))
   
  
}