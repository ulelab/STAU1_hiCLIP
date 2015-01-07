##############################################################
#' mergeHybrid
#'
#' mergeHybrid merge two HybridGRL objects 
#' @param \code{hgrl1}. HybridGRL object to be merged. 
#' @param \code{hgrl2}. HybridGRL object to be merged. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' mergeHybrid(hgrl1, hgrl2)
setGeneric(
  name = "mergeHybrid",
  def = function(object, hgrl2){standardGeneric("mergeHybrid")}
  )

setMethod(
  f = "mergeHybrid",
  signature = "HybridGRL",
  definition = function(object, hgrl2){
    all.seq <- unique(c(as.character(seqlevels(object$L)),
                        as.character(seqlevels(hgrl2$L))
                        )
                      )
    seqlevels(object) <- all.seq
    seqlevels(hgrl2) <- all.seq
    
    merged.grL <- GRangesList(L = c(object$L, hgrl2$L),
                              R = c(object$R, hgrl2$R)
                              )
    
    elementMetadata(merged.grL$L)$rname <- 1:length(merged.grL$L)
    elementMetadata(merged.grL$R)$rname <- 1:length(merged.grL$R)
    
    merged.grL <- new( "HybridGRL", merged.grL)
    
    return(merged.grL) 
  }
  )
  
##############################################################
#' selectHybridByGeneName
#'
#' selectHybridByGeneName returns valid number of hybrid reads 
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' exportHybrid(hgrl)
setGeneric(
  name = "selectHybridByGeneName",
  def = function(object, gene.vec){standardGeneric("selectHybridByGeneName")}
  )

setMethod(
  f = "selectHybridByGeneName",
  signature = "HybridGRL",
  definition = function(object, gene.vec){
    if(length(object$L) != length(object$R)){
      stop("")
    } else {}
    
    if(!all(seqnames(object$L) == seqnames(object$R))){
    	stop("Only functional for intra-gene hybrid")
    }
    
    object$L <- object$L[as.character(seqnames(object$L)) %in% gene.vec]
    object$R <- object$R[as.character(seqnames(object$R)) %in% gene.vec]
    
    return(object)
  }
  )

##############################################################
#' selectHybridByIndex
#'
#' selectHybridByIndex returns valid number of hybrid reads 
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' exportHybrid(hgrl)
setGeneric(
  name = "selectHybridByIndex",
  def = function(object, indexes){standardGeneric("selectHybridByIndex")}
  )

setMethod(
  f = "selectHybridByIndex",
  signature = "HybridGRL",
  definition = function(object, indexes){
    if(length(object$L) != length(object$R)){
      stop("")
    } else {}
    
    object$L <- object$L[indexes]
    object$R <- object$R[indexes]
    
    return(object)
  }
  )

##############################################################
#' sortHybrid
#'
#' sortHybrid sort HybridGRL files. Only menaning full for intra gene hybridds 
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' exportHybrid(hgrl)
setGeneric(
  name = "sortHybrid",
  def = function(object){standardGeneric("sortHybrid")}
  )

setMethod(
  f = "sortHybrid",
  signature = "HybridGRL",
  definition = function(object){
  	hybrid.unsorted.grL <- object
    # 1. sort hybrid.unsorted.grL$L and $R by the position with sort() function
    # 2. If necessary, wwap left and right reads so that left reads located to the 5' side of transcripts
    # 3. Return sorted GRangesList object (hybrid.grL) with $L and $R

    if(any(start(hybrid.unsorted.grL$L) > start(hybrid.unsorted.grL$R))){
      temp.ids <- (1:length(hybrid.unsorted.grL$L))[start(hybrid.unsorted.grL$L) > start(hybrid.unsorted.grL$R)]
      
      hybrid.grL <- GRangesList()
      hybrid.grL$L <- hybrid.unsorted.grL$L
      hybrid.grL$R <- hybrid.unsorted.grL$R
      
      hybrid.grL$R[temp.ids] <- hybrid.unsorted.grL$L[temp.ids]
      hybrid.grL$L[temp.ids] <- hybrid.unsorted.grL$R[temp.ids]
    } else {
      hybrid.grL <- hybrid.unsorted.grL
    }
    
    hybrid.grL <- new( "HybridGRL", hybrid.grL)
    return(hybrid.grL)

  }
  )

##############################################################
#' expandHGRL
#'
#' expandHGRL expand HybridGRL object by the score.
#' @param \code{hgrl}. HybridGRL object to be expanded. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' expandHGRL(hgrl)
setGeneric(
  name = "expandHGRL",
  def = function(object){standardGeneric("expandHGRL")}
  )

setMethod(
  f = "expandHGRL",
  signature = "HybridGRL",
  definition = function(object){
    
    expandGr <- function(gr){
	  gr <- gr[rep(1:length(gr), score(gr), each = TRUE)]
      # elementMetadata(gr)[["score"]] <- 1
      return(gr)
	}
	
	object$L <- expandGr(object$L)
	object$R <- expandGr(object$R)
    
    return(object) 
  }
  )
 
