##############################################################
#' addColumnGR
#'
#'  
#' @param \code{object}. Add new column for GenomicRange object. 
#' 
#' @export
#' @docType methods
#' @rdname ShortRead-methods
#'
#' @examples
#' addColumnGR(object, c.name, default.value = NA)
setGeneric(
  name = "addColumnGR",
  def = function(object, c.name, default.value = NA){standardGeneric("addColumnGR")}
  )

setMethod(
  f = "addColumnGR",
  signature = "GenomicRanges",
  definition = function(object, c.name, default.value = NA){
	values(object)[, c.name] <- default.value
    return(object)
  }
  )


##############################################################
#' addSeqGR
#'
#'  
#' @param \code{object}. ShortRead object to be exported. 
#' 
#' @export
#' @docType methods
#' @rdname ShortRead-methods
#'
#' @examples
#' addSeqGR(object)
setGeneric(
  name = "addSeqGR",
  def = function(object, fasta){standardGeneric("addSeqGR")}
  )

setMethod(
  f = "addSeqGR",
  signature = "GenomicRanges",
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
  
    rep.hyb <- addColumnGR(object, "sequence")
  elementMetadata(object)$sequence <- getSeq.transcript(
    fasta,
    as.character(seqnames(object)),
    start(object),
    end(object)
    )  
  return(object)
  
  }
  )


##############################################################
#' sortFasta
#'
#'  
#' @param \code{object}. ShortRead object to be exported. 
#' 
#' @export
#' @docType methods
#' @rdname ShortRead-methods
#'
#' @examples
#' sortFasta(object)
setGeneric(
  name = "sortFasta",
  def = function(object){standardGeneric("sortFasta")}
  )

setMethod(
  f = "sortFasta",
  signature = "ShortRead",
  definition = function(object){
      object <- object[order(as.character(id(object)))]
      return(object)
  }
  )

##############################################################
#' writeFastaGR
#'
#' writeFasraGR export Fastafile from GRanges object 
#' @param \code{object}. GRanges object to be exported. 
#' 
#' @export
#' @docType methods
#' @rdname GRanges-methods
#'
#' @examples
#' writeFastaGR(object)
setGeneric(
  name = "writeFastaGR",
  def = function(object, filename){standardGeneric("writeFastaGR")}
  )

setMethod(
  f = "writeFastaGR",
  signature = "GRanges",
  definition = function(object, filename){
    if(!("sequence" %in% colnames(elementMetadata(object)))){
      stop("sequence have to be appended")
    } else {}
    
    fasta <- ShortRead(sread = DNAStringSet(elementMetadata(object)$sequence), 
								id = BStringSet(as.character(seqnames(object)))
								)
    fasta <- sortFasta(fasta)
    
    writeFasta(fasta, filename)
    return(fasta)
  }
  )
  
##############################################################
#' extendHGRL
#'
#'  
#' @param \code{object}. ShortRead object to be exported. 
#' 
#' @export
#' @docType methods
#' @rdname ShortRead-methods
#'
#' @examples
#' extendHGRL(object)
setGeneric(
  name = "extendHGRL",
  def = function(object, ext.range){standardGeneric("extendHGRL")}
  )

setMethod(
  f = "extendHGRL",
  signature = "HybridGRL",
  definition = function(object, ext.range){
  	start(object$L) <- start(object$L) - ext.range
    end(object$L) <- end(object$L) + ext.range
    start(object$R) <- start(object$R) - ext.range
    end(object$R) <- end(object$R) + ext.range

	return(object)
  }
  )


##############################################################
#' extendHGRL2
#'
#'  
#' @param \code{object}. HybridGRL object to be extenede 
#' 
#' @export
#' @docType methods
#' @rdname ShortRead-methods
#'
#' @examples
#' extendHGRL2(object)
setGeneric(
  name = "extendHGRL2",
  def = function(object, l.ext.range, r.ext.range){standardGeneric("extendHGRL2")}
  )

setMethod(
  f = "extendHGRL2",
  signature = "HybridGRL",
  definition = function(object, l.ext.range, r.ext.range){
  	start(object$L) <- start(object$L) - l.ext.range
    end(object$L) <- end(object$L) + r.ext.range
    start(object$R) <- start(object$R) - l.ext.range
    end(object$R) <- end(object$R) + r.ext.range

	return(object)
  }
  )
  



##############################################################
#' concatenateHGRL
#'
#'  
#' @param \code{object}. ShortRead object to be exported. 
#' 
#' @export
#' @docType methods
#' @rdname ShortRead-methods
#'
#' @examples
#' concatenateHGRL(object)
setGeneric(
  name = "concatenateHGRL",
  def = function(object, fasta, filename, sep.seq = "AAAAAAAAAA", prefix.name = "hybrid"){standardGeneric("concatenateHGRL")}
  )

setMethod(
  f = "concatenateHGRL",
  signature = "HybridGRL",
  definition = function(object, fasta, filename, sep.seq = "AAAAAAAAAA", prefix.name = "hybrid"){
  	
  	# Duplicated, remove
  	addColumnGRL <- function(grL, c.name, default.value = NA){
      unlisted <- unlist(grL, use.names=FALSE)
      values(unlisted)[, c.name] <- default.value
      r.grl <- relist(unlisted, grL)
      return(r.grl)
    }
    
    object <- addColumnGRL(object, "sequence")
    object <- as(object, "HybridGRL")
       
    object <- addSeqHgrl(object, fasta)
    
    prop.df <- data.frame(id = rep("", length(object$L)), hybrid.width = rep(0, length(object$L)), stringsAsFactors = FALSE)
    
    prop.df.file <- paste(filename, ".tab", sep = "")
    fa.filename <- paste(filename, ".fa", sep = "")
    
    
    fa.df <- data.frame(
      seq.names = paste(as.character(seqnames(object$L)), 1:length(object$L), sep = "_"), 
      conc.seq = paste(elementMetadata(object$L)$sequence, elementMetadata(object$R)$sequence, sep = "&")
      )
      
    fa.df$seq.names <- paste(">", fa.df$seq.names, sep = "")
    
    write.table(fa.df, fa.filename, sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
    
	write.table(prop.df, prop.df.file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	
	return(prop.df)
  }
  )
  
  
##############################################################
#' randomPermutatons
#'
#' randomPermutatons perfoms the randoomization of the position of hybrid reads, but faster than the old version of randomPermutatons for multiply.factor > 1
#' Replaced on 21/June/2014
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' exportHybrid(hgrl)
setGeneric(
  name = "randomPermutatons",
  def = function(object, gtf, multiply.factor = 1, seed.n = 1){standardGeneric("randomPermutatons")}
  )

setMethod(
  f = "randomPermutatons",
  signature = "HybridGRL",
  definition = function(object, gtf, multiply.factor = 1, seed.n = 1){
	randomRepositioningDf <- function(gr, gtf, ol, multiply.factor){
      df <- data.frame(
        gene_id = rep(as.character(seqnames(gtf[ol[, 2]])), each = multiply.factor),
        start = 0,
        width = rep(width(gr[ol[, 1]]), each = multiply.factor)
        )
      
      for(i in 1:nrow(ol)){
      	df.start.index <- seq(from = 1, to = multiply.factor * nrow(ol), by = multiply.factor)[i]
      	
      	range.start <- start(gtf[ol[i, 2]])
      	range.end <- end(gtf[ol[i, 2]]) - width(gr[ol[i, 1]]) + 1
      	
      	if(range.start >= range.end){
      		## This is for a case if hybrid reads overlapping with the boundary of segments (5' UTR and CDS) and the reads are longer than 5' UTR.
      		range.end <- end(gtf[ol[i, 2]])
      	}
      	
        df$start[df.start.index:(df.start.index + multiply.factor - 1)] <- sample(
          range.start:range.end,
          size = multiply.factor,
          replace = TRUE
          )
      }
      
      r.gr <- GRanges(seqnames = Rle(df$gene_id), 
                      ranges = IRanges(df$start,
                        width = df$width), 
                      strand = Rle(rep("+", nrow(df))),
                      rname = 1:nrow(df),
                      score = as.integer(1),
                      category = "unknown",
                      annot = "unknown",
                      sequence = "unknown"
                      )
      return(r.gr)
    }
    
    
      set.seed(seed.n)
      object.seqlevels <- unique(c(
        as.character(seqnames(object$L)),
        as.character(seqnames(object$R))
        ))
      
      gtf <- gtf[seqnames(gtf) %in% object.seqlevels]
      seqlevels(gtf) <- object.seqlevels
      
      seqlevels(object) <- object.seqlevels
      
      ol.L <- as.matrix(findOverlaps(object$L, gtf))
      ol.R <- as.matrix(findOverlaps(object$R, gtf))
      
      ol.L <- ol.L[!duplicated(ol.L[, 1]), ]
      ol.R <- ol.R[!duplicated(ol.R[, 1]), ]
      
      r.object <- GRangesList(L = GRanges(), R = GRanges())
      seqlevels(r.object) <- object.seqlevels
      
      temp.L <- randomRepositioningDf(object$L, gtf, ol.L, multiply.factor)
      seqlevels(temp.L) <- object.seqlevels
      r.object$L <- temp.L
      
      temp.R <- randomRepositioningDf(object$R, gtf, ol.R, multiply.factor)
      seqlevels(temp.R) <- object.seqlevels
      r.object$R <- temp.R

      r.hobject <- new("HybridGRL", r.object)
      return(r.hobject)
      
  }
  )