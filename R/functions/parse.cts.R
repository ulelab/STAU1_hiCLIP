parse.rRNA.ct.file <- function(file.name, rRNA.name = "rRNA_18S"){
  if(rRNA.name == "rRNA_18S"){
    rRNA.length <- 1869
  } else if(rRNA.name == "rRNA_28S"){
    rRNA.length <- 5070
  } else {
    stop("Currently only applicable to parse 18S and 28S rRNA ct file")
  }


  ct <- read.table(file.name, sep = " ", skip = 1)

  indexConsecVector <- function(x, incremental = TRUE){
    
    if(incremental){
      Breaks <- c(0, which(diff(x) != 1), length(x))
    } else {
      Breaks <- c(0, which(diff(x) != -1), length(x))
    }
    consc.list <- sapply(seq(length(Breaks) - 1), function(i) x[(Breaks[i] + 1):Breaks[i+1]])
    
    if(is.list(consc.list)){
      consc.ind <- sapply(consc.list, length)
      consc.indexes <- rep(1:length(consc.ind), consc.ind)
    } else {
      consc.ind <- sapply(consc.list, length)
      consc.indexes <- consc.ind
    }
    
    return(consc.indexes)
  }


  ## Ignored buldges
  ct <- ct[ct[, 5] != 0, ]
  
  ct$l.grp <- indexConsecVector(ct[, 6], incremental = TRUE)
  ct$r.grp <- indexConsecVector(ct[, 5], incremental = FALSE)
  ct$index <-as.integer(factor(paste(ct$l.grp, ct$r.grp, sep = "_")))

  ct.list <- split(ct, f=ct[, "index"])

  extract.duplex.position <- function(ct.part){
    left.strat <- ct.part[1, 6]
    left.end <- ct.part[nrow(ct.part), 6]
    right.strat <- ct.part[nrow(ct.part), 5]
    right.end <- ct.part[1, 5]
    interval.range <- right.strat - left.end - 1
    
    vec <- c(left.strat,
             left.end,
             right.strat,
             right.end,
             interval.range
             )
    
    return(vec)
  }

  rRNA.duplexes <- as.data.frame(t(sapply(ct.list, extract.duplex.position)))
  colnames(rRNA.duplexes) <- c("s1_start", "s1_end", "a1_start", "a1_end", "range")
  rownames(rRNA.duplexes) <- paste(rRNA.name, 1:nrow(rRNA.duplexes), sep = ".")
  rRNA.duplexes$range.ratio <- rRNA.duplexes$range / rRNA.length

  rRNA.duplexes <- rRNA.duplexes[order(rRNA.duplexes$s1_start), ]
  rRNA.duplexes <- rRNA.duplexes[rRNA.duplexes$range > 0, ]

  result.df <- rRNA.duplexes

  rRNA.GR.L <- forgiToGR.L(rRNA.duplexes)
  rRNA.GR.R <- forgiToGR.R(rRNA.duplexes)

  rRNA.grL <- GRangesList()
  seqlevels(rRNA.grL) <- unique(
    c(as.character(seqnames(rRNA.GR.L)),
    as.character(seqnames(rRNA.GR.R))
    )
    )
  
  rRNA.grL$L <- rRNA.GR.L
  rRNA.grL$R <- rRNA.GR.R

  rRNA.HGRL <- new("HybridGRL", rRNA.grL)
  
  result.list <- list(rRNA.HGRL = rRNA.HGRL, rRNA.df = result.df)

  return(result.list)
}
