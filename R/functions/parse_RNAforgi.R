cleanForgiOut <- function(df){
  maxN <- function(x, N=2){
    ## http://stackoverflow.com/questions/2453326/fastest-way-to-find-second-third-highest-lowest-value-in-vector-or-column
    len <- length(x)
    if(N>len){
      warning('N greater than length(x).  Setting N=length(x)')
      N <- length(x)
    }
    sort(x,partial=len-N+1)[len-N+1]
  }
  
  
  stem.index <- intersect(grep("^s", df[, 2]), grep("^define", df[, 1]))
  stems <- df[stem.index, 2:ncol(df)]
  stems[, 2:5] <- apply(stems[, 2:5], 2, as.integer)
  stems$range <- apply(stems[, 2:5], 1, function(x){maxN(x, N = 2) - maxN(x, N = 3) - 1})
  return(stems)
}

parseRNAforgi <- function(filename){
  forgi.lines <- readLines(filename)
  forgi.start <- grep("^>", forgi.lines)
  forgi.mark <- vector('integer', length(forgi.lines))
  forgi.mark[forgi.start] <- 1
  forgi.mark <- cumsum(forgi.mark)

  df <- lapply(split(forgi.lines, forgi.mark),
               function(.data){
                 .input <- read.table(
                   textConnection(.data),
                   skip = 3,
                   fill = TRUE,
                   stringsAsFactors = FALSE,
                   col.names = c("define", "name", "s1_start", "s1_end", "a1_start", "a1_end")  
                   )
                 attr(.input, 'name') <- gsub(">", "", .data[1])  # save the names
                 attr(.input, 'length') <- as.integer(strsplit(.data[3], " ")[[1]][2])
                 .input
               }
               )

  names(df) <- sapply(df, attr, 'name')
  df <- lapply(df, cleanForgiOut)
  
  return(df)
}

forgiToGR.L <- function(forgi.df){

  forgi.gene.ids <- sapply(strsplit(rownames(forgi.df), "\\."), "[", 1)

  forgi.unsorted.gr.L <- GRanges(
    rname = paste(forgi.gene.ids, 1:nrow(forgi.df), sep = "_"),
    seqnames = Rle(forgi.gene.ids), 
    ranges = IRanges(start = forgi.df$s1_start,
      end = forgi.df$s1_end),
    strand = Rle("+"),
    score = as.integer(1),
    category = "unknown",
    annot = "unknown",
    sequence = "unknown"
    )
  
  return(forgi.unsorted.gr.L)
}

forgiToGR.R <- function(forgi.df){

  forgi.gene.ids <- sapply(strsplit(rownames(forgi.df), "\\."), "[", 1)

  forgi.unsorted.gr.R <- GRanges(
    rname = paste(forgi.gene.ids, 1:nrow(forgi.df), sep = "_"),
    seqnames = Rle(forgi.gene.ids), 
    ranges = IRanges(start = forgi.df$a1_start,
      end = forgi.df$a1_end),
    strand = Rle("+"),
    score = as.integer(1),
    category = "unknown",
    annot = "unknown",
    sequence = "unknown"
    )
  
  return(forgi.unsorted.gr.R)
}

forgiToHGRL <- function(forgi.out.list){

  forgis <- do.call("rbind", forgi.out.list)	

  forgiGR.L <- forgiToGR.L(forgis)
  forgiGR.R <- forgiToGR.R(forgis)

  if(!all(elementMetadata(forgiGR.L)$rname == elementMetadata(forgiGR.R)$rname)){
    stop("rname for left and righ reads are different")
  }

  forgi.unsorted.grL <- GRangesList()
  seqlevels(forgi.unsorted.grL) <- unique(
    names(forgi.out.list)
    )
  
  forgi.unsorted.grL$L <- forgiGR.L
  forgi.unsorted.grL$R <- forgiGR.R

  forgi.HGRL <- new("HybridGRL", forgi.unsorted.grL)

  return(forgi.HGRL)
}
