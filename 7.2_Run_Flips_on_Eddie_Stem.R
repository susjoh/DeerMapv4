
library(crimaptools)
library(reshape)

print(paste("Running window", flips.tab$window.name[i]))

#~~ take the original SNP order

snp.order <- subset(skel.map, CEL.LG == flips.tab$CEL.LG[i])[flips.tab$Start[i]:flips.tab$Stop[i],"SNP.Name"]

#~~ create lists adata frame to store information

map.list   <- list()
flips.list <- list()
snp.order.list <- list()
reordered.snps <- list()

map.stats <- NULL

system(paste0("mkdir ", flips.tab$window.name[i]))

setwd(flips.tab$window.name[i])

iteration = 1

#~~ Run flips analysis

for(loop.it in 1:iteration.max){
  
  print(paste("Running iteration", iteration))
  
  #~~ Create crimap files and run flips
  snp.order.list[[iteration]] <- snp.order
  
  create_crimap_input(gwaa.data = abeldata,
                      familyPedigree = famped,
                      analysisID = flips.tab$window.name[i],
                      snplist = snp.order,
                      use.specific.mnd = "6_Mendelian_Errors_DoubXovers_Master.txt")
  
  run_crimap_prepare(genfile = paste0("chr", flips.tab$window.name[i]))
  run_crimap_chrompic(genfile = paste0("chr", flips.tab$window.name[i]))
  run_crimap_flips(genfile = paste0("chr", flips.tab$window.name[i]), flips = flips.val)
  
  #~~ parse map, flips and stats, and save to file
  
  map.list[[iteration]] <- parse_map_chrompic(chrompicfile = paste0("chr", flips.tab$window.name[i], ".cmp"))
  
  flip.res <- parse_flips(paste0("chr", flips.tab$window.name[i], ".fl", flips.val))
  names(flip.res)[window_size+1] <- "Li.Change"
  flip.res$Li.Change <- as.numeric(as.character(flip.res$Li.Change))
  flip.res$Map.Li <- flip.res[1, "Li.Change"]
  flip.res <- flip.res[-1,]
  flip.res$Min.Li <- min(flip.res$Li.Change)
  flips.list[[iteration]] <- flip.res
  
  map.stats <- rbind(map.stats,
                     data.frame(Iteration = iteration,
                                Map.Li    = flip.res$Map.Li[1],
                                Min.Li    = flip.res$Min.Li[1],
                                Map.Len   = map.list[[iteration]]$cMPosition[window_size],
                                last.run  = ifelse(nrow(flip.res) == 0, "yes", "no")))
  
  flip.res <- subset(flip.res, Li.Change < Li.threshold)
  
  save(snp.order.list, map.list, flips.list, map.stats, file = paste0("../", flips.tab$window.name[i], ".Rdata"))
  
  
  #~~ delete the input files
  
  del.vec <- dir()
  del.vec <- del.vec[grep(paste0("chr", flips.tab$window.name[i], "."), del.vec, fixed = T)]
  if(length(grep(".R", del.vec, fixed = T)) > 0) del.vec <- del.vec[-grep(".R", del.vec, fixed = T)]
  if(length(grep("flipsrun", del.vec, fixed = T)) > 0) del.vec <- del.vec[-grep("flipsrun", del.vec, fixed = T)]



  if(Sys.info()["sysname"] == "Windows") {
    for(j in del.vec){
      system("cmd", input = paste0("del ", j))
    }
  } else {
    for(j in del.vec){
      system(paste0("rm ", j))
    }
  }
  
  #~~ determine if the order is the best.
  
  print(map.stats)
  
  if(nrow(flip.res) == 0){
    message(paste0("Best SNP order found on iteration ", iteration, ": no improvement obtained from flips analysis. Analysis terminated."))
    break
  } 
  
  if(iteration > iteration.no.improve & map.stats$Min.Li[iteration] != max(map.stats$Min.Li)){
    
    message(paste0("No improvement after ", iteration, " iterations, best SNP order found on iteration(s) ",
                   paste(which(map.stats$Min.Li == max(map.stats$Min.Li)), collapse = ", "), ". Analysis terminated."))
    break
  
    
  } else {
    
    snp.order.hold <- NULL
    
    for(k in 1:nrow(flip.res)){
      
      flip.res2 <- subset(flip.res, Li.Change == Min.Li)[k,]
      
      x.melt <- melt(flip.res2, id.vars = c("Li.Change", "Map.Li", "Min.Li"))
      x.melt$value <- as.numeric(as.character(x.melt$value))
      x.melt$new.order <- 1:nrow(x.melt)
      x.melt$new.order[which(!is.na(x.melt$value))] <- x.melt$value[which(!is.na(x.melt$value))]+1
      
      x.melt$SNP.Name <- snp.order[x.melt$new.order]
      
      reordered.snps[[iteration]] <- x.melt$SNP.Name[which(!is.na(x.melt$value))]
      
      snp.order.hold <- x.melt$SNP.Name
        
      if(any(lapply(snp.order.list, function (x) all(snp.order.hold == x)) == TRUE)){
        message(k)
        if(k == nrow(flip.res2)) stop("All orders tried previously: check results for best SNP order")

      } else {
        break
      }
      
    }
    
    if(is.null(snp.order.hold)){
      break("All orders tried previously: check results for best SNP order")
    } else {
      
      snp.order <- snp.order.hold
      
      iteration <- iteration + 1
    }
    
  }
}

system("rm *")

setwd("..")

system(paste0("rmdir ", flips.tab$window.name[i]))
