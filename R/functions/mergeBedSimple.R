## mergeBed.R
# Yoichiro Sugimoto
# Usage Rscript mergeBed.R "input directory name" "output bed file name"
options("scipen"=100)

args <- commandArgs(TRUE)

inputDir <- args[1]
outName <- args[2]

# inputDir <- "~/Documents/Staufen/Stau1/Hybridization_iCLIP/2ndTime_with_Transfer/nonHybrid"
# outName <- "~/Desktop/LUs27_nonHybrid.bed" 

runMerge <- function(inputDir, outName, pattern.file = "\\.bed$"){
  bedList <- list.files(inputDir, pattern = pattern.file, full.names = TRUE)
  bedList <- bedList[grep("unmatched", bedList, invert = TRUE)]

  i <- 0
  for(bedFile in bedList){
    
    bed <- read.table(bedFile, sep = "\t", skip = 1)
    
    bed[, 5] <- "+"
    bed[bed[, 4] < 0, 5] <- "-"
    
    colnames(bed) <- c("chr", "start", "end", "score", "strand")	
    
    ## print(head(bed))
    if(i == 0){
      mergedBed <- bed
    } else {
      mergedBed <- (merge(mergedBed, bed, by = c("chr", "start", "end", "strand"), all = TRUE))
	}
    i <- i + 1
  }
  
  mergedBed[is.na(mergedBed)] <- 0
  resultBed <- mergedBed[, 1:3]
  resultBed[, 4] <- rowSums(mergedBed[, 5:ncol(mergedBed)])
  colnames(resultBed)[4] <- "score"
  resultBed  <- resultBed[with(resultBed, order(resultBed$chr, resultBed$start, -1 * resultBed$score)), ]
  resultBed[, 4] <- as.numeric(resultBed[, 4])
  cat("In total, ")
  cat(paste(sum(abs(resultBed[, 4])), " nonHybrid reads were obtained.\n"))
  print(sum(abs(resultBed[, 4])))
  
  resultBed <- resultBed[grep("^ENSG", resultBed[, 1]), ]
  
  print("Among the nonHybrid reads, those mapped to genes are:")
  print(sum(abs(resultBed[, 4])))
  
  createHeader <- function(bedFileName, species = "hg19"){
    if (species == "hg19"){
      header_temp <- 'track type=bedGraph name="bedFileName" description="bedFileName" db=hg19 color="120,101,172" lib_id="bedFileName" maxHeightPixels="100:50:0" visibility="full" mapped_to="hg19" species="Hg" condition="" res_type="G" altColor="200,120,59" priority="20" regions="" replicate="None"\n'
    }
    header <- gsub(pattern = "bedFileName", replacement = bedFileName, header_temp)
    header
  }
  
  zz <- file(outName, "w")

  cat(createHeader(unlist(strsplit(outName, "/"))[length(unlist(strsplit(outName, "/")))]), file=zz)
  write.table(resultBed, file=zz, append=TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
  
  close(zz) 
  
  system(paste("gzip -f ", outName, sep =""))
}

if(!is.na(inputDir)){
  # runMerge(inputDir, outName)
}


## example
# Rscript /Users/yoichiro/Desktop/updated/STAU1_unwinding/mergeBed-2.0.R /Users/yoichiro/Desktop/updated/STAU1_unwinding/wig/raw/wt /Users/yoichiro/Desktop/updated/STAU1_unwinding/wig/merged/wt_merged.wig
# Rscript /Users/yoichiro/Desktop/updated/STAU1_unwinding/mergeBed-2.0.R /Users/yoichiro/Desktop/updated/STAU1_unwinding/wig/raw/KD /Users/yoichiro/Desktop/updated/STAU1_unwinding/wig/merged/KD_merged.wig
# Rscript /Users/yoichiro/Desktop/updated/STAU1_unwinding/mergeBed-2.0.R /Users/yoichiro/Desktop/updated/STAU1_unwinding/wig/raw/RC /Users/yoichiro/Desktop/updated/STAU1_unwinding/wig/merged/RC_merged.wig
