##############################################################
#' addColumnHGRL
#'
#' addColumnHGRL extend regions of HybridGRL object 
#' @param \code{hgrl}. HybridGRL object to be analyzed. 
#' @param \code{fasta}. ShortRead object imported from fasta file of transcript
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' addColumnHGRL(hgrl, fasta)
setGeneric(
  name = "addColumnHGRL",
  def = function(object, c.name, default.value = NA){standardGeneric("addColumnHGRL")}
  )

setMethod(
  f = "addColumnHGRL",
  signature = "HybridGRL",
  definition = function(object, c.name, default.value = NA){
      unlisted <- unlist(object, use.names=FALSE)
      values(unlisted)[, c.name] <- default.value
      r.grl <- relist(unlisted, object)
      hgrl <- new( "HybridGRL", r.grl)
      return(hgrl)
  }
  )


##############################################################
#' confidentIsland
#'
#' confidentIsland identified RNA duplexed identified by more than \code{n} hybrid reads.
#' @param \code{hgrL}. HybridGRL object to be analyzed. 
#' @param \code{low.val}. The number of independent hybrid reads which identified the duplex to be defined as confident duplex. 
#" @param \code{width.min}. Minimun overlap width required.
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' confidentIsland(hgrL, low.val = 2, width.min = 8)
setGeneric(
  name = "confidentIsland",
  def = function(object, low.val = 2, width.min = 8){standardGeneric("confidentIsland")}
  )

setMethod(
  f = "confidentIsland",
  signature = "HybridGRL",
  definition = function(object, low.val, width.min){
        
    temp.gr.L <- object$L
    temp.gr.L <- unique(temp.gr.L)
    cov.L <- coverage(temp.gr.L)
    island.L <- slice(cov.L, lower = low.val)
    gr.L <- GRanges(seqnames = Rle(rep(names(sapply(island.L, length)), as.integer(sapply(island.L, length)))),
                    ranges = IRanges(
                      start = unlist(start(island.L)),
                      end = unlist(end(island.L))
                      ),
                    strand = Rle(rep("+", length(unlist(start(island.L))))),
                    score = unlist(viewMaxs(island.L))
                    )
    gr.L <- gr.L[width(gr.L) > width.min]
    rep.hyb <- GRangesList(L = GRanges(), R = GRanges())
    seqlevels(rep.hyb) <-  unique(as.character(seqnames(gr.L)))
    
    for(i in 1:length(gr.L)){
      temp.hybrid.gr.L <- object$L[seqnames(object$L) == as.character(seqnames(gr.L[i]))]
      seqlevels(temp.hybrid.gr.L) <- as.character(seqnames(gr.L[i]))
      
      temp.hybrid.gr.L <- subsetByOverlaps(temp.hybrid.gr.L, gr.L[i])
      temp.hybrid.gr.R <- object$R[elementMetadata(object$R)$rname %in% elementMetadata(temp.hybrid.gr.L)$rname]
      seqlevels(temp.hybrid.gr.R) <- as.character(seqnames(gr.L[i]))
      
      temp.gr.R <- temp.hybrid.gr.R
      temp.gr.R <- unique(temp.gr.R)
      temp.cov.R <- coverage(temp.gr.R)
      island.R <- slice(temp.cov.R, lower = low.val)
      if(length(unlist(start(island.R))) > 0){
        gr.R <- GRanges(seqnames = Rle(rep(names(sapply(island.R, length)), as.integer(sapply(island.R, length)))),
                        ranges = IRanges(
                          start = unlist(start(island.R)),
                          end = unlist(end(island.R))
                          ),
                        strand = Rle(rep("+", length(unlist(start(island.R))))),
                        score = unlist(viewMaxs(island.R))
                        )
        
        gr.R <- gr.R[width(gr.R) > width.min]
        if(length(gr.R) > 0){
          gr.R <- gr.R[score(gr.R) == max(score(gr.R))]
          gr.R <- sort(gr.R)[1]
          seqlevels(gr.L) <- unique(as.character(seqnames(gr.L)))
          seqlevels(gr.R) <- unique(as.character(seqnames(gr.L)))
          elementMetadata(gr.L)$score[i] <- min(score(gr.L[i]), score(gr.R))
          elementMetadata(gr.R)$score <- min(score(gr.L[i]), score(gr.R))
          
          rep.hyb$L <- c(rep.hyb$L, gr.L[i])
          rep.hyb$R <- c(rep.hyb$R, gr.R)
        }
      }
    }
    
    rep.hyb.hgrL <- new("HybridGRL", rep.hyb)
    
    rep.hyb.hgrL <- addColumnHGRL(rep.hyb.hgrL, "rname", default.value = 1)
    elementMetadata(rep.hyb.hgrL$L)$rname <- 1:length(rep.hyb.hgrL$L)
    elementMetadata(rep.hyb.hgrL$R)$rname <- 1:length(rep.hyb.hgrL$R)
    
    rep.hyb.hgrL <- addColumnHGRL(rep.hyb.hgrL, "category", "unknown")
    rep.hyb.hgrL <- addColumnHGRL(rep.hyb.hgrL, "annot", "unknown")
    rep.hyb.hgrL <- addColumnHGRL(rep.hyb.hgrL, "sequence", "unknown")
    
    return(rep.hyb.hgrL)
  }
  )

##############################################################
#' findDuplex2
#'
#' findDuplex2 find the longest RNA duplex using RNAhybrid unlike findDuplex, which uses guugle.
#' 
#' @param \code{hgrL}. HybridGRL object to be analyzed. 
#' @param \code{fasta}. ShortRead object imported from fasta file
#' @param \code{t.gtf}. Transcriptome coodrinate gene annotation as GenomicRange object.
#' @param \code{out.file}. Directory for output and temporary file.
#' @param \code{hiCdir}. Directory of package.  
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' findDuplex2(hgrL, fasta, out.file, hiCdir)
setGeneric(
  name = "findDuplex2",
  def = function(object, fasta, out.file, hiCdir){standardGeneric("findDuplex2")}
  )

setMethod(
  f = "findDuplex2",
  signature = "HybridGRL",
  definition = function(object, fasta, out.file, hiCdir){

    ## A function to parse the RNAhybrid compact output is loaded.
    parseRNAhybridVec <- function(RNAhybrid.out.vec){
      
      RNAhybrid.out.vec <- as.character(RNAhybrid.out.vec)

      if(!all(c(nchar(RNAhybrid.out.vec[8]) == nchar(RNAhybrid.out.vec[9]), nchar(RNAhybrid.out.vec[9]) == nchar(RNAhybrid.out.vec[10]), nchar(RNAhybrid.out.vec[10]) == nchar(RNAhybrid.out.vec[11])))){
	print(RNAhybrid.out.vec)
	stop("Unexpected RNAhybrid output")
      } 

      pos.index <- 1:nchar(RNAhybrid.out.vec[8])

      target.start.pos <- as.integer(RNAhybrid.out.vec[7])
      shown.seq.length <- nchar(RNAhybrid.out.vec[8])

      pos.df <- data.frame(
	p.i = pos.index, 
	t.u = strsplit(RNAhybrid.out.vec[8], NULL)[[1]],
	t.p = strsplit(RNAhybrid.out.vec[9], NULL)[[1]],
	q.p = strsplit(RNAhybrid.out.vec[10], NULL)[[1]],
	q.u = strsplit(RNAhybrid.out.vec[11], NULL)[[1]]
	)

      ## positional index for the duplex
      dup.index <- pos.df$p.i
      dup.index[pos.df$t.p == " "] <- 0
      pos.df$dup.index <- dup.index

      ## positional index for the target
      target.bools <- c(pos.df$t.u != " ") | c(pos.df$t.p != " ")
      pos.df$target.index <- 0
      pos.df$target.index[target.bools] <- 1:sum(target.bools)

      ## positional index for the query
      query.bools <- c(pos.df$q.u != " ") | c(pos.df$q.p != " ")
      pos.df$query.index <- 0
      pos.df$query.index[query.bools] <- sum(query.bools):1
      query.seq.len <- sum(query.bools)

      ##
      pos.df <- pos.df[pos.df$dup.index != 0, ]

      indexConsecVector <- function(x, incremental = TRUE){
        
        if(incremental){
          Breaks <- c(0, which(diff(x) != 1), length(x))
        } else {
          Breaks <- c(0, which(diff(x) != -1), length(x))
        }
        consc.list <- lapply(seq(length(Breaks) - 1), function(i) x[(Breaks[i] + 1):Breaks[i+1]])
        
        
        consc.ind <- sapply(consc.list, length)
        consc.indexes <- rep(1:length(consc.ind), consc.ind)
        
        return(consc.indexes)
      }
      
      pos.df$dup.id <- indexConsecVector(pos.df$dup.index, incremental = TRUE)

      index.longest <- names(table(pos.df$dup.id))[table(pos.df$dup.id) == max(table(pos.df$dup.id))][1]
      pos.df <- pos.df[pos.df$dup.id == as.integer(index.longest), ]

      longest.dup.pos.target <- target.start.pos - 1 + pos.df$target.index[1]
      longest.dup.pos.query <- pos.df$query.index[nrow(pos.df)]
      length.longest.dup <- nrow(pos.df)

      return.vec <- c(longest.dup.pos.target, longest.dup.pos.query, length.longest.dup)

      return(return.vec)

    }


    ## A function to run the RNAhybrid program.
    RunRNAhybrid <- function(rep.hyb, RNAhybridOutName, RNAhybridTempOut, hiCdir){
      if(file.exists(RNAhybridOutName)){
  	system(paste("rm", RNAhybridOutName))
      }
      
      for(i in 1:length(rep.hyb$L)){
        fasta.L <- ShortRead(
          sread = DNAStringSet(elementMetadata(rep.hyb$L)$sequence[i]),
          id = BStringSet(paste(i, "L", sep = "_"))
          )
        
        fasta.R <- ShortRead(
          sread = DNAStringSet(elementMetadata(rep.hyb$R)$sequence)[i],
          id = BStringSet(paste(i, "R", sep = "_"))
          )
        
        L.temp <- paste(RNAhybridTempOut, "/L.temp.fa", sep = "")
        R.temp <- paste(RNAhybridTempOut, "/R.temp.fa", sep = "")
        
        writeFasta(fasta.L, L.temp)
        writeFasta(fasta.R, R.temp)
        
        RNAhybrid.command <- paste(hiCdir, "/inst/bin/RNAhybrid/bin/RNAhybrid -s 3utr_human -m 1000 -n 1000 -c -t ", L.temp, " -q ", R.temp, " >> ", RNAhybridOutName, sep = "")
        ## print(RNAhybrid.command)
        system(RNAhybrid.command)	
      }
      return(0)
    }


    object <- addSeqHgrl(object, fasta)

    RNAhybrid.dir <- paste(out.file, "RNAHybridDuplex", sep = "/")
    dir.create(file.path(RNAhybrid.dir), showWarnings = FALSE, recursive = TRUE)
    
    RNAhybridOutName <- paste(RNAhybrid.dir, "hybrid.transcripts.RNAhybrid", sep = "/")
    RNAhybridTempOut <- paste(RNAhybrid.dir, "temp_RNAhybrid", sep = "/")
    
    dir.create(file.path(RNAhybridTempOut), showWarnings = FALSE, recursive = TRUE)




    RunRNAhybrid(object, RNAhybridOutName, RNAhybridTempOut, hiCdir)

    hybrid.df <- read.table(RNAhybridOutName, sep = ":", header = FALSE, stringsAsFactors = FALSE)
    
    if(!all(sapply(strsplit(hybrid.df$V1, "_"), "[", 1) == sapply(strsplit(hybrid.df$V1, "_"), "[", 1))){
      stop("RNAhybrid did not identify a duplex for all the pairs: 1")
    } else if (nrow(hybrid.df) != as.integer(sapply(strsplit(hybrid.df$V1, "_"), "[", 1)[nrow(hybrid.df)])){
      stop("RNAhybrid did not identify a duplex for all the pairs: 2")
    }
    
    
    dup.summary <- apply(hybrid.df, 1, parseRNAhybridVec)
    

    left.stem <- GRanges(seqnames = seqnames(object$L),
                         ranges = IRanges(
                           start = start(object$L) + dup.summary[1, ] - 1,
                           width = dup.summary[3, ]
                           ),
                         strand = strand(object$L),
                         score = score(object$L),
                         rname = elementMetadata(object$L)$rname,
                         category = "unknown",
                         annot = "unknown",
                         sequence = "unknown"
                         )
    
    right.stem <- GRanges(seqnames = seqnames(object$R),
                          ranges = IRanges(
                            start = start(object$R) + dup.summary[2, ] - 1,
                            width = dup.summary[3, ]
                            ),
                          strand = strand(object$R),
                          score = score(object$R),
                          rname = elementMetadata(object$R)$rname,
                          category = "unknown",
                          annot = "unknown",
                          sequence = "unknown"
                          )
    
    
    left.stem <- left.stem[dup.summary[3, ] > 0]
    right.stem <- right.stem[dup.summary[3, ] > 0]
    
    stem.grL <- GRangesList(L = left.stem,
                            R = right.stem)
    
    
    
    stem.hgrL <- new("HybridGRL", stem.grL)

    
    return(stem.hgrL)
  }
  )

##############################################################
#' addSeqHgrl
#'
#' addSeqHgrl add sequence information to HybridGRL object 
#' @param \code{hgrl}. HybridGRL object to be analyzed. 
#' @param \code{fasta}. ShortRead object imported from fasta file of transcript
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' addSeqHgrl(hgrl, fasta)
setGeneric(
  name = "addSeqHgrl",
  def = function(object, fasta){standardGeneric("addSeqHgrl")}
  )

setMethod(
  f = "addSeqHgrl",
  signature = "HybridGRL",
  definition = function(object, fasta){
   
    getSeq.transcript <- function(fasta, gene_id, start, end){          
      require("ShortRead")
      if(!is.character(gene_id)){
        stop("gene_id must be character object")
      }
      fasta.vec <- 1:length(fasta)
      names(fasta.vec) <- as.character(id(fasta))
      return.seq <- substr(sread(fasta[fasta.vec[gene_id]]), start, end)
      return(return.seq)
    }
    
    if(!("sequence" %in% names(elementMetadata(object$L)))){
    	object <- addColumnHGRL(object, "sequence", default.value = "NA")
    }
  
    elementMetadata(object$L)$sequence <- getSeq.transcript(
      fasta,
      as.character(seqnames(object$L)),
      start(object$L),
      end(object$L)
      )
    
    elementMetadata(object$R)$sequence <- getSeq.transcript(
      fasta,
      as.character(seqnames(object$R)),
      start(object$R),
      end(object$R)
      )
  
    return(object)
  }
  )
