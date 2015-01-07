##############################################################
#' selectGeneMappedHybrid
#'
#' Select hybrid reads both of arms were uniquly mapped to genes after excluding reads from rRNAs, tRNAs, rDNAs, and mitochondria transcript.
#' @param \code{hgrl}. HybridGRL object to be analyzed.
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' selectGeneMappedHybrid(hgrl)
setGeneric(
  name = "selectGeneMappedHybrid",
  def = function(object){standardGeneric("selectGeneMappedHybrid")}
  )

setMethod(
  f = "selectGeneMappedHybrid",
  signature = "HybridGRL",
  definition = function(object){
      intersected.id <- intersect(grep("^ENSG", seqnames(object$L)),
                                  grep("^ENSG", seqnames(object$R))
                                  )
      object$L <- object$L[intersected.id]
      object$R <- object$R[intersected.id]
  
      strand.flag <- (as.character(strand(object$L)) == "+") &
        (as.character(strand(object$R)) == "+")

      object$L <- object$L[strand.flag]
      object$R <- object$R[strand.flag]
  
      return(object)
    }
  )

##############################################################
#' isGeneMapped
#'
#' Test whether HbridGRL object contains those after selectGeneMappedHybrid().
#' @param \code{hgrl}. HybridGRL object to be analyzed.
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' isGeneMapped(hgrl)
setGeneric(
  name = "isGeneMapped",
  def = function(object){standardGeneric("isGeneMapped")}
  )

setMethod(
  f = "isGeneMapped",
  signature = "HybridGRL",
  definition = function(object){
    if(any(strand(object$L) == "-") | any(strand(object$R) == "-")){
      stop("This function is only meaningful for the analysis of gene mapped hybrid reads.")
    } else {}

    if(
      (length(
        grep("^ENSG", seqnames(object$L)) 
        ) != length(object$L)
       ) |
      (length(
        grep("^ENSG", seqnames(object$R))
        ) != length(object$R)
       )
      ){
      stop("This function is only meaningful for the analysis of gene mapped hybrid reads.")
    } else {}
  }
  )


##############################################################
#' plotInterVsIntra
#'
#' plotInterVsIntra plot the ratio of hybrid reads from inter-molecular or intra-molecular RNA structures.
#' @param \code{hgrl}: HybridGRL object to be analyzed.
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' plotInterVsIntra(hgrl)
setGeneric(
  name = "plotInterVsIntra",
  def = function(object){standardGeneric("plotInterVsIntra")}
  )

setMethod(
  f = "plotInterVsIntra",
  signature = "HybridGRL",
  definition = function(object){

    isGeneMapped(object)
    
    inter.intra <- data.frame(id = 1:length(object$L),
                              inter.flag = factor(ifelse(as.character(seqnames(object$L)) == as.character(seqnames(object$R)), "intra", "inter"), levels = c("intra", "inter"))
                              )

    cat("The ratio of intra-molecular vs inter-molecular RNA structures [%].\n")
    print(
      round(
        prop.table(table(as.character(inter.intra$inter.flag))) * 100
        , digit = 1)
      )

    plotGG <- ggplot(inter.intra) +
      geom_bar(aes(x = factor(1), fill = inter.flag), width = 1) +
        coord_polar(theta = "y") +
          theme(axis.ticks = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y =  element_blank(),
                axis.text.x =  element_blank(),
                axis.title.x = element_blank(),
                panel.background = element_blank(),
                legend.title=element_blank()
                ) 
    print(plotGG)
  }
  )

##############################################################
#' selectHybridType
#'
#' Select hybrid reads from either inter- or intra- molecular RNA structures.
#' @param \code{hgrl} HybridGRL object to be analyzed.
#' @param \code{hybridType} hybridType, either "inter" or "intra".
#'
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' selectHybridType(hgrl, "inter")
setGeneric(
  name = "selectHybridType",
  def = function(object, hybridType){standardGeneric("selectHybridType")}
  )

setMethod(
  f = "selectHybridType",
  signature = "HybridGRL",
  definition = function(object, hybridType){
    isGeneMapped(object)

    if(!(hybridType %in% c("inter", "intra"))){
      stop("hybridType has to be 'inter' or 'intra'")
    }

    if(hybridType == "inter"){
      inter.gene.flag <- c(as.character(seqnames(object$L)) !=
                           as.character(seqnames(object$R))
                           )
      object$L <- object$L[inter.gene.flag]
      object$R <- object$R[inter.gene.flag]
  
      return(object)
    } else if(hybridType == "intra"){
      intra.gene.flag <- as.character(seqnames(object$L)) ==
        as.character(seqnames(object$R))
      
      object$L <- object$L[intra.gene.flag]
      object$R <- object$R[intra.gene.flag]
  
      return(object)
    }
  }
  )

##############################################################
#' selectHybridAnnot
#'
#' selectHybridAnnot returns valid number of hybrid reads 
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' @param \code{RNAannot} RNA annotation of hybrid reads, typically "utr3", "CDS", and "utr5".
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' selectHybridAnnot(hgrl, "utr3")
setGeneric(
  name = "selectHybridAnnot",
  def = function(object, t.gtf, annot){standardGeneric("selectHybridAnnot")}
  )

setMethod(
  f = "selectHybridAnnot",
  signature = "HybridGRL",
  definition = function(object, t.gtf, annot){

    ## Method specific function...
    
    subAnnotHybrid <- function(hybrid.grL, t.gtf, type = "start"){
      
      if(all(seqnames(hybrid.grL$L) != seqnames(hybrid.grL$R))){
        stop("L and R are not matched")
      }
      
      temp.hyb <- GRanges(seqnames = seqnames(hybrid.grL$L),
                          ranges = IRanges(
                            start = start(hybrid.grL$L),
                            width = 1
                            ),
                          strand = strand(hybrid.grL$L),
                          score = score(hybrid.grL$L),
                          rname = elementMetadata(hybrid.grL$L)$rname
                          )
      if(type == "all"){
        end(temp.hyb) <- end(hybrid.grL$R)
      } else if (type == "start"){
        end(temp.hyb) <- start(hybrid.grL$R)
      } else {
        stop("type must be start or all")
      }
      
      seqlevels.in.use <- unique(as.character(seqnames(temp.hyb)))
      sub.t.gtf <- t.gtf[seqnames(t.gtf) %in% seqlevels.in.use]
      seqlevels(sub.t.gtf) <- seqlevels.in.use
      seqlevels(temp.hyb) <- seqlevels.in.use
  
      mapDf <- as.matrix(findOverlaps(temp.hyb, sub.t.gtf, type = "within"))
      
      if(nrow(mapDf[duplicated(mapDf[, 1]), ]) != 0){
        stop("Unexpected duplicate")
      }
      
      elementMetadata(temp.hyb)$annot <- "unknown"
      elementMetadata(temp.hyb)$annot[mapDf[, 1]] <- as.character(elementMetadata(sub.t.gtf)$annot[mapDf[, 2]])
      temp.hyb <- temp.hyb[mapDf[, 1]]
      return(temp.hyb)
    }

    ## Till here, method specific function...

    temp.hyb <- subAnnotHybrid(object, t.gtf, type = "start")
    
    temp.hyb <- temp.hyb[!is.na(elementMetadata(temp.hyb)$annot)]
    selected.hyb <- temp.hyb[elementMetadata(temp.hyb)$annot == annot]
    selected.rname <- elementMetadata(selected.hyb)$rname
    
    selected.grL <- GRangesList(
      L = object$L[elementMetadata(object$L)$rname %in% selected.rname],
      R = object$R[elementMetadata(object$R)$rname %in% selected.rname]
      )
    
    selected.hgrL <- new("HybridGRL", selected.grL)
    
    return(selected.hgrL)
    
  }
  )


##############################################################
#' select.rRNAsHybrid
#'
#' Select hybrid reads both of arms were mapped to rRNAs.
#' @param \code{hgrl}. HybridGRL object to be analyzed.
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' select.rRNAsHybrid(hgrl)
setGeneric(
  name = "select.rRNAsHybrid",
  def = function(object){standardGeneric("select.rRNAsHybrid")}
  )

setMethod(
  f = "select.rRNAsHybrid",
  signature = "HybridGRL",
  definition = function(object){  
      rRNA.flag <- intersect(grep("rRNA", seqnames(object$L)),
                                  grep("rRNA", seqnames(object$R))
                                  )
      
      object$L <- object$L[rRNA.flag]
      object$R <- object$R[rRNA.flag]
      
      strand.flag <- (as.character(strand(object$L)) == "+") &
        (as.character(strand(object$R)) == "+")

      object$L <- object$L[strand.flag]
      object$R <- object$R[strand.flag]
  
      return(object)
    }
  )
  
  
##############################################################
#' grepHybrid
#'
#' Regular expression for hybrid reads.
#' @param \code{hgrl}. HybridGRL object to be analyzed.
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' grepHybrid(hgrl)
setGeneric(
  name = "grepHybrid",
  def = function(object, pattern, type = c("all", "intra", "inter")){standardGeneric("grepHybrid")}
  )

setMethod(
  f = "grepHybrid",
  signature = "HybridGRL",
  definition = function(object, pattern, type = c("all", "intra", "inter")){  
      grep.flag <- intersect(grep(pattern, seqnames(object$L)),
                                  grep(pattern, seqnames(object$R))
                                  )
      
      object$L <- object$L[grep.flag]
      object$R <- object$R[grep.flag]
      
      strand.flag <- (as.character(strand(object$L)) == "+") &
        (as.character(strand(object$R)) == "+")

      object$L <- object$L[strand.flag]
      object$R <- object$R[strand.flag]
      
      if(type == "intra"){
      	intra.index <- seqnames(object$L) == seqnames(object$R)
      	object$L <- object$L[intra.index]
      	object$R <- object$R[intra.index]
      } else if(type == "inter"){
      	inter.index <- !(seqnames(object$L) == seqnames(object$R))
      	object$L <- object$L[inter.index]
      	object$R <- object$R[inter.index]
      } else {
      	print("Both intra and inter gene hybrid reads are returned.")
      }
  
      return(object)
    }
  )
