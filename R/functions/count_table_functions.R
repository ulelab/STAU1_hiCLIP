import_nonHybrid <- function(file){
  df <- read.table(file, sep = "\t", skip = 1)
  gr <- GRanges(seqnames = Rle(df[, 1]), 
                ranges = IRanges(df[, 3], width = 1), 
                strand = Rle(ifelse(df[, 4] > 0, "+", "-")),
                score = df[, 4]
                )
  return(gr)
}

createCountTableSub <- function(bed.file, bed, gr, seqlevels.in.use){
  gr <- gr[seqnames(gr) %in% seqlevels.in.use]
  seqlevels(gr) <- seqlevels.in.use

  a.bed <- subsetByOverlaps(bed, gr)
  
  a.count.df <- aggregate(score(a.bed),
                            list(as.character(seqnames(a.bed)))
                        , sum)
  colnames(a.count.df) <- c("id", gsub(".bed", "", strsplit(bed.file, "/")[[1]][length(strsplit(bed.file, "/")[[1]])]) )
  
  return(a.count.df)
}

createCountTable <- function(bed.dir, match.pattern, out.pre, mRNA.df, a.CDS){
  bed.list <- list.files(bed.dir, pattern = match.pattern, full.names = TRUE)
  i <- 0
  for(bed.file in bed.list){
    
    bed <- import_nonHybrid(bed.file)
    bed <- bed[seqnames(bed) %in% mRNA.df$gene_id]
    bed <- bed[strand(bed) == "+"]
    
    seqlevels.in.use <- unique(
      as.character(
        seqnames(bed)
        )
      )
    
    seqlevels(bed) <- seqlevels.in.use
    
    count.df <- aggregate(score(bed),
                          list(as.character(seqnames(bed)))
                          , sum)
    
    colnames(count.df) <- c("id", gsub(".bed", "", strsplit(bed.file, "/")[[1]][length(strsplit(bed.file, "/")[[1]])]) )
    
    
   	a.count.df <- createCountTableSub(bed.file, bed, a.CDS, seqlevels.in.use)
    
    if(i == 0){
      total.count <- count.df
      a.total.count <- a.count.df
      
    } else {
      total.count <- merge(total.count, count.df, by = "id", all = TRUE)
      total.count[is.na(total.count)] <- 0
      
      a.total.count <- merge(a.total.count, a.count.df, by = "id", all = TRUE)
      a.total.count[is.na(a.total.count)] <- 0
    }
    i <- i + 1	
  }
  
  write.table(total.count, paste(out.pre, "_total_count.txt", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  write.table(a.total.count, paste(out.pre, "_trim30.CDS_count.txt", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}


createCountTableUtr3 <- function(bed.dir, match.pattern, out.pre, out.post, mRNA.df, utr3.gr){
  bed.list <- list.files(bed.dir, pattern = match.pattern, full.names = TRUE)
  bed.list <- bed.list[-grep("unmatched", bed.list)]
  i <- 0
  for(bed.file in bed.list){
    
    bed <- import_nonHybrid(bed.file)
    bed <- bed[seqnames(bed) %in% mRNA.df$gene_id]
    bed <- bed[strand(bed) == "+"]
    
    seqlevels.in.use <- unique(
      as.character(
        seqnames(bed)
        )
      )
    
    seqlevels(bed) <- seqlevels.in.use
    
    
    utr3.gr <- utr3.gr[seqnames(utr3.gr) %in% seqlevels.in.use]
    seqlevels(utr3.gr) <- seqlevels.in.use
    
    a.bed <- subsetByOverlaps(bed, utr3.gr)
    
    a.count.df <- aggregate(score(a.bed),
                            list(as.character(seqnames(a.bed)))
                        , sum)
    colnames(a.count.df) <- c("id", gsub(".bed", "", strsplit(bed.file, "/")[[1]][length(strsplit(bed.file, "/")[[1]])]) )
    
    
    if(i == 0){
      a.total.count <- a.count.df
    } else {
      
      a.total.count <- merge(a.total.count, a.count.df, by = "id", all = TRUE)
      a.total.count[is.na(a.total.count)] <- 0
    }
    i <- i + 1	
  }

  write.table(a.total.count, paste(out.pre, out.post, sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

createCountTableNonHybrid <- function(bed.dir, match.pattern, out.pre, mRNA.df){
  bed.list <- list.files(bed.dir, pattern = match.pattern, full.names = TRUE)
  bed.list <- bed.list[-grep("unmatched", bed.list)]
  
  i <- 0
  for(bed.file in bed.list){
    
    bed <- import_nonHybrid(bed.file)
    bed <- bed[seqnames(bed) %in% mRNA.df$gene_id]
    bed <- bed[strand(bed) == "+"]
    
    seqlevels.in.use <- unique(
      as.character(
        seqnames(bed)
        )
      )
    
    seqlevels(bed) <- seqlevels.in.use
    
    count.df <- aggregate(score(bed),
                          list(as.character(seqnames(bed)))
                          , sum)
    
    colnames(count.df) <- c("id", gsub(".bed", "", strsplit(bed.file, "/")[[1]][length(strsplit(bed.file, "/")[[1]])]) )
    
    
    if(i == 0){
      total.count <- count.df
    } else {
      total.count <- merge(total.count, count.df, by = "id", all = TRUE)
      total.count[is.na(total.count)] <- 0
    }
    i <- i + 1	
  }
  
  write.table(total.count, paste(out.pre, "_total_count_mRNA.txt", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}
