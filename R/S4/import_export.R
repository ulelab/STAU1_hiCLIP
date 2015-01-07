##############################################################
#' exportHybrid
#'
#' exportHybrid export HybridGRL objects as IGV compatible sam file. 
#' This function is specific for HybridGRL objects for intra-molecular RNA strctures. See below.
#'
#' @param \code{hgrl}. HybridGRL object to be exported. HybridGRL objects should be those represnting intra-molecular RNA structures. 
#' @param \code{filename}: Output filename.
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' exportHybrid(hgrl, "~/test.sam")
setGeneric(
  name = "exportHybrid",
  def = function(object, filename){standardGeneric("exportHybrid")}
  )

setMethod(
  f = "exportHybrid",
  signature = "HybridGRL",
  definition = function(object, filename){
    # Export hybrid.grL objects as sam file inorder to vesualize on genome browser such as IGV
    # Check validity of input
    # 1. Are they only from intra-molecular RNA strctures?
    if(!all(seqnames(object$L) == seqnames(object$R))){
    	stop("This method is designed to export intra-moecular hybrid reads.")
    } else {}
    # 2. Are the HybridGRL objects are sorted?
    if(any(start(object$L) > start(object$R))){
    	stop("HybridGRL objects must be sorted.")
    } else {}
    
    # If left arm and right arm of hybrid reads are overlapped, ignore them.
    gap.width.init <- start(object$R) - end(object$L) - 1
    
    print(sum(gap.width.init < 0))
    print(table(gap.width.init[gap.width.init < 0]))

    object$L <- object$L[gap.width.init >= 0]
    object$R <- object$R[gap.width.init >= 0]
  
    sam.df <- data.frame(
      read.ids = elementMetadata(object$L)$rname,
      strand = 0,
      seqnames = as.character(seqnames(object$L)),
      start = start(object$L),
      quality = 255,
      length = paste(
        width(object$L),
        "M",
        gap.width.init[gap.width.init >=0],
        "N",
        width(object$R),
        "M",
        sep = ""
        ),
      V7 = "*",
      V8 = 0,
      V9 = 0,
      V10 = "*",
      V11 = "*"
      )

    sam.df <- sam.df[order(sam.df$seqnames, sam.df$start), ]
    write.table(sam.df, filename, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)  
    }
)

#' importHybridSam: Constructor for HybridGRL class.
#'
#' This function cleate HybridGRL class from bed like and sam like file. 
#'
#' @param hybrid.samfile: exported file by exportHybrid
#' @examples
importHybridSam <- function(hybrid.samfile){
  # Import ".hybrid.sam" file created by exportHybrid function and generate hybrid.grL object
    # read sam file and extract the start site as GRanges object
  require(GenomicRanges)
  sam <- read.table(hybrid.samfile, sep = "\t")
  sam$left.length <- as.integer(sapply(strsplit(as.character(sam[, 6]), "[NM]"), "[", 1))
  sam$all.length <- sapply(strsplit(as.character(sam[, 6]), "[NM]"), function(x){sum(as.integer(x))})
  sam$right.length <- as.integer(sapply(strsplit(as.character(sam[, 6]), "[NM]"), "[", 3))
  
  sam$left.start <- sam[,4]
  sam$right.end <- sam[, 4] + sam$all.length - 1
  
  left.gr <- GRanges(seqnames = Rle(sam[, 3]),
                     ranges = IRanges(start = sam$left.start,
                           width = sam$left.length
                           ),
                     strand = Rle(strand(rep("+", nrow(sam)))),
                     score = rep(1, nrow(sam)),
                     rname = 1:nrow(sam)
                         )
                         
  right.gr <- GRanges(seqnames = Rle(sam[, 3]),
                      ranges = IRanges(end = sam$right.end,
                        width = sam$right.length
                        ),
                      strand = Rle(strand(rep("+", nrow(sam)))),
                      score = rep(1, nrow(sam)),
                      rname = 1:nrow(sam)
                      )
                         
  hybrid.grL <- GRangesList(L = left.gr, R = right.gr)
  hgrL <- new( "HybridGRL", hybrid.grL)
  return(hgrL)
}
