##############################################################
#' interGenePartner
#'
#' interGenePartner returns the interacting partner of specified gene type. 
#' @param \code{hgrl}. HybridGRL object to be examined.
#' @param \code{gene_type}. gene type where the interacting partner will be examinded
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' interGenePartner(hgrl, "lncRNA")
setGeneric(
  name = "interGenePartner",
  def = function(object, gene_type){standardGeneric("interGenePartner")}
  )

setMethod(
  f = "interGenePartner",
  signature = "HybridGRL",
  definition = function(object, gene_type){

    concatenateFlag <- function(vec){
      vec <- paste(sort(vec), collapse = "-")
      vec
    }
    
    inter.annot.df <- data.frame(
      left.gene = as.character(seqnames(object$L)),
      right.gene = as.character(seqnames(object$R)),
      left = elementMetadata(object$L)$annot,
      right = elementMetadata(object$R)$annot
      )
    
    inter.annot.df$inter.annot <- apply(inter.annot.df[, 3:4], 1, concatenateFlag)
    inter.annot.df$inter.gene <- apply(inter.annot.df[, 1:2], 1, concatenateFlag)
    
    inter.annot.df$inter.annot <- factor(inter.annot.df$inter.annot, 
                                         levels=names(sort(table(inter.annot.df$inter.annot), 
                                           decreasing=TRUE))
                                         )
    
    selected.inter.df <- inter.annot.df[grep(gene_type, inter.annot.df$inter.annot), ]
    
    cat(paste("The number of inter gene interactions, where ", gene_type, " involves \n", sep = ""))
    print(
      round(
        table(as.character(selected.inter.df$inter.annot))
        , digit = 1)
      )
    
    print(
      ggplot(selected.inter.df) + geom_bar(aes(x = inter.annot, col = inter.annot))
      )
  }
  )

