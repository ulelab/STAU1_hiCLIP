##############################################################
#' readNumber
#'
#' readNumber returns valid number of hybrid reads 
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' exportHybrid(hgrl)
setGeneric(
  name = "readNumber",
  def = function(object){standardGeneric("readNumber")}
  )

setMethod(
  f = "readNumber",
  signature = "HybridGRL",
  definition = function(object){
    if(length(object$L) != length(object$R)){
      stop("")
    } else {}
    return(length(object$L))
  }
  )

