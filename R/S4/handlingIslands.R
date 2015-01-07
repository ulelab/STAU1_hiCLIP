##############################################################
#' Extendhybrid
#'
#' extendHybrid merge two HybridGRL objects 
#' @param \code{hgrl}. HybridGRL object to be extended. 
#' @param \code{min.width}. Minimum width of the arms. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' extendHybrid(hgrl, min.width)
setGeneric(
  name = "extendHybrid",
  def = function(object, min.width = 17){standardGeneric("extendHybrid")}
  )

setMethod(
  f = "extendHybrid",
  signature = "HybridGRL",
  definition = function(object, min.width = 17){

    cor.df <- data.frame(
      left.start = rep(0, length(object$L)),
      left.end = rep(0, length(object$L)),
      right.start = rep(0, length(object$L)),
      right.end = rep(0, length(object$L))
      )


    cor.df$left.width <- width(object$L)
    cor.df$right.width <- width(object$R)

    cor.df$left.b <- cor.df$left.width < min.width
    cor.df$right.b <- cor.df$right.width < min.width

    cor.df$left.cor.fac <- min.width - cor.df$left.width
    cor.df$right.cor.fac <- min.width - cor.df$right.width

    cor.df$left.start <- ceiling(cor.df$left.cor.fac/2)
    cor.df$left.end <- floor(cor.df$left.cor.fac/2)

    cor.df$right.start <- ceiling(cor.df$right.cor.fac/2)
    cor.df$right.end <- floor(cor.df$right.cor.fac/2)

    cor.df$left.start[!cor.df$left.b] <- 0
    cor.df$left.end[!cor.df$left.b] <- 0

    cor.df$right.start[!cor.df$right.b] <- 0
    cor.df$right.end[!cor.df$right.b] <- 0

    left.stem <- GRanges(seqnames = seqnames(object$L),
                         ranges = IRanges(
                           start = start(object$L) -  cor.df$left.start,
                           end = end(object$L) +  cor.df$left.end
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
                            start = start(object$R) - cor.df$right.start,
                            end = end(object$R) +  cor.df$right.end
                            ),
                          strand = strand(object$R),
                          score = score(object$R),
                          rname = elementMetadata(object$R)$rname,
                          category = "unknown",
                          annot = "unknown",
                          sequence = "unknown"
                          )
    
    
    
    stem.grL <- GRangesList(L = left.stem,
                            R = right.stem)
    
    
    
    stem.hgrL <- new("HybridGRL", stem.grL)

    return(stem.hgrL)
  }
  )

##############################################################
#' findOverlapHGRL
#'
#' findOverlapHGRL finds the overlapping hybrids and returns them, and can also be used to subtract overlapping hybrid reads.
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' exportHybrid(hgrl)
setGeneric(
  name = "findOverlapHGRL",
  def = function(object, subject, type = c("intersect", "subtract")){standardGeneric("findOverlapHGRL")}
  )

setMethod(
  f = "findOverlapHGRL",
  signature = "HybridGRL",
  definition = function(object, subject, type = c("intersect", "subtract")){
    object.seqlevels <- unique(c(
      as.character(seqnames(object$L)),
      as.character(seqnames(object$R)),
      as.character(seqnames(subject$L)),
      as.character(seqnames(subject$R))
      ))
    seqlevels(object) <- object.seqlevels
    seqlevels(subject) <- object.seqlevels
    
    ol.L <- as.matrix(findOverlaps(object$L, subject$L, type = "any"))
    ol.R <- as.matrix(findOverlaps(object$R, subject$R, type = "any"))
    
    ol.L.df <- as.data.frame(ol.L)
    ol.R.df <- as.data.frame(ol.R)
    
    ol.L.df$indicator.L <- TRUE
    ol.R.df$indicator.R <- TRUE
    
    res <- merge(ol.L.df, ol.R.df, all=TRUE)
    res$indicator.both <- apply(res[, c("indicator.L", "indicator.R")], 1, all)
    res$indicator.both[is.na(res$indicator.both)] <- FALSE
    
    res.both <- res[res$indicator.both, ]
    index.both <- sort(unique(res.both$queryHits))
    
    if(type == "intersect"){
      object$L <- object$L[index.both]
      object$R <- object$R[index.both]
    } else if(type == "subtract"){
      object$L <- object$L[!(1:length(object$L) %in% index.both)]
      object$R <- object$R[!(1:length(object$R) %in% index.both)]
    }
    
    return(object)

  }
  )

##############################################################
#' collapsedIdenticalHGRL
#'
#' collapsedIdenticalHGRL collapsed identical hybrid reads 
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' exportHybrid(hgrl)
setGeneric(
  name = "collapsedIdenticalHGRL",
  def = function(object){standardGeneric("collapsedIdenticalHGRL")}
  )

setMethod(
  f = "collapsedIdenticalHGRL",
  signature = "HybridGRL",
  definition = function(object){
    left_cord <- paste(seqnames(object$L), start(object$L), end(object$L), sep = ".")
    right_cord <- paste(seqnames(object$R), start(object$R), end(object$R), sep = ".")
    
    total_cord <- paste(left_cord, right_cord, sep = "_")
    
    identical.index <- duplicated(total_cord)
    
    elementMetadata(object$L)$score <- table(total_cord)[total_cord]
    elementMetadata(object$R)$score <- table(total_cord)[total_cord]
    
    object$L <- object$L[!identical.index]
    object$R <- object$R[!identical.index]
    
    return(object)
  }
  )

