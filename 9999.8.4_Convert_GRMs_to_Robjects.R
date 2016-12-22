source("r/makeGRM.R")
library(beepr)

memory.limit(size = 800000000) #ESSENTIAL OTHERWISE YOU RUN INTO MEMORY ISSUES


#~~ read the GRM file from GCTA

system.time({
  
  grm.auto <- read.table("gcta/Deer_autoGRM_adj.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
  ids.auto <- read.table("gcta/Deer_autoGRM_adj.grm.id")  # CONTAINS ID LIST
  
  ids.auto$Order <- 1:nrow(ids.auto)
  
  #~~ remove IDs not used in the analysis
  
#   ids.auto <- ids.auto[which(ids.auto$V2 %in% famped[,1]),]
#   grm.auto <- grm.auto[which(grm.auto$V1 %in% ids.auto$Order & grm.auto$V2 %in% ids.auto$Order),-3]
  
  save(grm.auto, ids.auto, file = "gcta/Deer_autoGRM_adj.Rdata")
  
  rm(ids.auto, grm.auto)
  
})


beep()

system.time({
  for(i in 1:33){
    
    print(paste("Running chromosome", i, "of 33"))
    
    grm.chr <- read.table(paste0("gcta/Deer_chr", i, "_GRM_adj.grm.gz"))
    ids.chr <- read.table(paste0("gcta/Deer_chr", i, "_GRM_adj.grm.id"))
    
    ids.chr$Order <- 1:nrow(ids.chr)
    
#     ids.chr <- ids.chr[which(ids.chr$V2 %in% famped[,1]),]
#     grm.chr <- grm.chr[which(grm.chr$V1 %in% ids.chr$Order & grm.chr$V2 %in% ids.chr$Order),-3]
    
    grm.wochr <- read.table(paste0("gcta/Deer_WOchr", i, "_GRM_adj.grm.gz"))
    ids.wochr <- read.table(paste0("gcta/Deer_WOchr", i, "_GRM_adj.grm.id"))
    
    ids.wochr$Order <- 1:nrow(ids.wochr)
    
#     ids.wochr <- ids.wochr[which(ids.wochr$V2 %in% famped[,1]),]
#     grm.wochr <- grm.wochr[which(grm.wochr$V1 %in% ids.wochr$Order & grm.wochr$V2 %in% ids.wochr$Order),-3]
    
    save(grm.chr, ids.chr, grm.wochr, ids.wochr, file = paste0("gcta/Deer_Chr_", i, "_GRMs.Rdata"))
    
    rm(grm.chr, ids.chr, grm.wochr, ids.wochr)
  }
})

beep()
