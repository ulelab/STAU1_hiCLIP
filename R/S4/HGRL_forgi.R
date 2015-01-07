##############################################################
#' findOverlapsHGRL
#'
#' findOverlapsHGRL identified query hybrid reads overlapping with the subjct reads. 
#' @param \code{query}. HybridGRL object to be examined. 
#' @param \code{subjct}. HybridGRL object to be used as the reference.
#' @param \code{overlap.type}. See the manual for findOverlaps function in GRanges.
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' findOverlapsHGRL(query, subjct, overlap.type = "within")
setGeneric(
  name = "findOverlapsHGRL",
  def = function(object, subject.HGRL, overlap.type = "within"){standardGeneric("findOverlapsHGRL")}
  )

setMethod(
  f = "findOverlapsHGRL",
  signature = "HybridGRL",
  definition = function(object, subject.HGRL, overlap.type = "within"){
    object <- addColumnHGRL(object, c.name = "RNAfold.predicted", default.value = FALSE)
    
    seqlevels(object) <- seqlevels(subject.HGRL)
    left.match <- as.data.frame(findOverlaps(object$L, subject.HGRL$L, type= overlap.type))
    right.match <- as.data.frame(findOverlaps(object$R, subject.HGRL$R, type= overlap.type))
    
    left.match$merge <- with(left.match, paste(queryHits, subjectHits, sep = "_"))
    
    right.match$merge <- with(right.match, paste(queryHits, subjectHits, sep = "_"))
    
    common_elements <- intersect(left.match$merge, right.match$merge)
    
    overlapped_ids <- sapply(strsplit(common_elements, "_"), "[", 1)
    bin.results <- 1:length(object$L) %in% as.integer(overlapped_ids)
    
    elementMetadata(object$L)$RNAfold.predicted <- bin.results
    elementMetadata(object$R)$RNAfold.predicted <- bin.results
    
    return(object)
  }
  )

##############################################################
#' compareForgiRNAfold
#'
#' compareForgiRNAfold compared the hybrid identified structures and RNAfold predicted structures.
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' exportHybrid(hgrl)
setGeneric(
  name = "compareForgiRNAfold",
  def = function(object, forgi.HGRL){standardGeneric("compareForgiRNAfold")}
  )

setMethod(
  f = "compareForgiRNAfold",
  signature = "HybridGRL",
  definition = function(object, forgi.HGRL){
    object <- findOverlapsHGRL(object, forgi.HGRL, overlap.type = "within")

    RNAfold.df <- data.frame(range = start(object$R) - end(object$L) - 1, is.predicted = elementMetadata(object$L)$RNAfold.predicted)
    
    print(paste("In total ", readNumber(object), " structures are examined.", sep = ""))
    print(paste(sum(elementMetadata(object$L)$RNAfold.predicted), " structures identified by the hybrid reads are predicted by RNAfold."))
    print(paste(sum(!elementMetadata(object$L)$RNAfold.predicted), " structures identified by the hybrid reads are not predicted by RNAfold."))
    
    p1 <- ggplot(data = RNAfold.df) + geom_density(aes(x = log(range, 10), group = is.predicted, fill = is.predicted), alpha = 0.3)

    p2 <- ggplot(data = RNAfold.df) + geom_boxplot(aes(y = log(range, 10), x = is.predicted, fill = is.predicted), alpha = 0.3)

    wil.out <- wilcox.test(RNAfold.df$range[RNAfold.df$is.predicted], RNAfold.df$range[!RNAfold.df$is.predicted])
    print(wil.out)
    print(paste("p-value is:",
                  wil.out$p.value
                  ))

    out.plot <- list(a = p1, b = p2)
    return(out.plot)
  }
  )


