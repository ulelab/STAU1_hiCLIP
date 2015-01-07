##############################################################
#' nonhybridCountOnDuplex
#'
#'  
#' @param \code{object}. HybridGRL object to be analyzed 
#' 
#' @export
#' @docType methods
#' @rdname ShortRead-methods
#'
#' @examples
#' nonhybridCountOnDuplex(object, nonHybrid.gr, ext.left, ext.right)
setGeneric(
  name = "nonhybridCountOnDuplex",
  def = function(object, nonHybrid.gr, ext.left, ext.right){standardGeneric("nonhybridCountOnDuplex")}
  )

setMethod(
  f = "nonhybridCountOnDuplex",
  signature = "HybridGRL",
  definition = function(object, nonHybrid.gr, ext.left, ext.right){
  	object <- extendHGRL2(object, ext.left, ext.right)
  	
  	seqlevel.in.use <- unique(c(as.character(seqlevels(object$L)),
  	as.character(seqlevels(object$R)),
  	as.character(seqlevels(nonHybrid.gr))
  	))

  	gr <- GRanges()
  	seqlevels(gr) <- seqlevel.in.use
  	seqlevels(object$L) <- seqlevel.in.use
  	seqlevels(object$R) <- seqlevel.in.use

  	for(i in 1:length(object$L)){
  		gr.temp <- c(object$L[i], object$R[i])
  		gr <- c(gr, reduce(gr.temp))
  	}
  	
	seqlevels(nonHybrid.gr) <- seqlevel.in.use

  	total.count <- countOverlaps(gr, nonHybrid.gr)
  	
  	return(sum(total.count))  
  
  }
  )
  
