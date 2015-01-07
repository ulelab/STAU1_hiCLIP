##############################################################
#' analyzeRNADuplexRange
#'
#' analyzeRNADuplexRange perform statistical test for the range of RNA strctures in CDS and 3' UTR as well as produce corresponding plot.
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' @param \code{t.gtf}. Annotation of transcripts as GenomicRanges object
#' @param \code{out.dir}. Location of output directory where the result table will be created.
#' @param \code{hiCdir}. Directory of package
#' @param \code{print.flag}. If TRUE, perform statistical anlysis and areate plot, otherwise only output summary of RNA structure range.
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' analyzeRNADuplexRange(hgrl, t.gtf, out.dir, hiCdir)
setGeneric(
  name = "analyzeRNADuplexRange",
  def = function(object, t.gtf, out.dir, hiCdir, print.flag = FALSE){standardGeneric("analyzeRNADuplexRange")}
  )

setMethod(
  f = "analyzeRNADuplexRange",
  signature = "HybridGRL",
  definition = function(object, t.gtf, out.dir, hiCdir, print.flag){
    ### Method specific function..
    
    rangeDuplexDf <- function(hybrid.grL, t.gtf, annot){

      if(all(seqnames(hybrid.grL$L) != seqnames(hybrid.grL$R))
         &
         all(score(hybrid.grL$L) != score(hybrid.grL$R))
         ){
        stop("L and R have to have the smae seqnames")
      }

	selected.grL <- selectHybridAnnot(hybrid.grL, t.gtf, annot)
	## In order to consider the score.
	selected.grL <- expandHGRL(selected.grL)

      distance.df <- data.frame(
        gene_id = as.character(seqnames(selected.grL$L)),
        distance = start(selected.grL$R) - end(selected.grL$L) - 1,
        score = score(selected.grL$L)
        )
      distance.df <- distance.df[distance.df$distance > 0, ]

      ## anlyze width of gene with specific annoyation such as CDS and utr3
      t.gtf.annot <- t.gtf[elementMetadata(t.gtf)$annot == annot]
      width.vec <- width(t.gtf.annot)
      names(width.vec) <- as.character(seqnames(t.gtf.annot))
      
      distance.df$gene.width <- width.vec[as.character(distance.df$gene_id)]

      if(!all(complete.cases(distance.df))){
        stop("Error, annotation is likely to be wrong...")
      }
      
      distance.df$distance.ratio <- distance.df$distance / distance.df$gene.width
      distance.df$annot <- annot
      
      return(distance.df)
    }

    ### Till here, method specific function.

    utr3.range.df <- rangeDuplexDf(object, t.gtf, "utr3")
    CDS.range.df <- rangeDuplexDf(object, t.gtf, "CDS")

    if(print.flag){
      cat(paste("N for CDS is:", nrow(CDS.range.df), "\n"))
      cat(paste("N for 3' UTR is:", nrow(utr3.range.df), "\n"))
    }
    
    merged.df <- rbind(CDS.range.df, utr3.range.df)
    
    merged.df <- merged.df[order(merged.df$annot, merged.df$distance, decreasing = TRUE), ]
    
    dir.create(out.dir, showWarnings = FALSE)

    table.file <- paste(out.dir, "range_summary.tab", sep = "/")
    write.table(merged.df, table.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

    if(print.flag){
      cat("Mann-Whitney U test for the range of RNA structures \n in 3' UTR and CDS.\n")
      raw.distance.wilcox <- wilcox.test(merged.df$distance[merged.df$annot == "CDS"], merged.df$distance[merged.df$annot == "utr3"])
      # print(raw.distance.wilcox)
      print(paste("p-value is:",
                  raw.distance.wilcox$p.value
                  ))

      plot1 <- ggplot(merged.df) + geom_density(aes(x = log(distance, 10), fill = annot), alpha = 0.3) + theme(legend.position="top")
      
    
      cat("Mann-Whitney U test for the normalized range \n(normalized by the width of gene segment) \nof RNA structures in 3' UTR and CDS.\n")
      ratio.distance.wicox <- wilcox.test(merged.df$distance.ratio[merged.df$annot == "CDS"], merged.df$distance.ratio[merged.df$annot == "utr3"])
      # print(ratio.distance.wicox)
      print(paste("p-value for ratio is: ",
                ratio.distance.wicox$p.value
                  ))
      
      
      plot2 <- ggplot(merged.df) + geom_density(aes(x = log(distance.ratio, 10), fill = annot), alpha = 0.3) + theme(legend.position="top")
      
      multiplot(plot1, plot2, cols = 2)
    }

    
  }
  )

##############################################################
#' analyzeRNADuplexRangeUnexpanded
#'
#' This is specific function for the analysis of the comparison of RNAfold and hiCLIP
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' @param \code{t.gtf}. Annotation of transcripts as GenomicRanges object
#' @param \code{out.dir}. Location of output directory where the result table will be created.
#' @param \code{hiCdir}. Directory of package
#' @param \code{print.flag}. If TRUE, perform statistical anlysis and areate plot, otherwise only output summary of RNA structure range.
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' analyzeRNADuplexRangeUnexpanded(hgrl, t.gtf, out.dir, hiCdir)
setGeneric(
  name = "analyzeRNADuplexRangeUnexpanded",
  def = function(object, t.gtf, out.dir, hiCdir, print.flag = FALSE){standardGeneric("analyzeRNADuplexRangeUnexpanded")}
  )

setMethod(
  f = "analyzeRNADuplexRangeUnexpanded",
  signature = "HybridGRL",
  definition = function(object, t.gtf, out.dir, hiCdir, print.flag){
    ### Method specific function..
    
    rangeDuplexDf <- function(hybrid.grL, t.gtf, annot){

      if(all(seqnames(hybrid.grL$L) != seqnames(hybrid.grL$R))
         &
         all(score(hybrid.grL$L) != score(hybrid.grL$R))
         ){
        stop("L and R have to have the smae seqnames")
      }

	selected.grL <- selectHybridAnnot(hybrid.grL, t.gtf, annot)
	## In order to consider the score.
	## selected.grL <- expandHGRL(selected.grL)

      distance.df <- data.frame(
        gene_id = as.character(seqnames(selected.grL$L)),
        distance = start(selected.grL$R) - end(selected.grL$L) - 1,
        score = score(selected.grL$L)
        )
      distance.df <- distance.df[distance.df$distance > 0, ]

      ## anlyze width of gene with specific annoyation such as CDS and utr3
      t.gtf.annot <- t.gtf[elementMetadata(t.gtf)$annot == annot]
      width.vec <- width(t.gtf.annot)
      names(width.vec) <- as.character(seqnames(t.gtf.annot))
      
      distance.df$gene.width <- width.vec[as.character(distance.df$gene_id)]

      if(!all(complete.cases(distance.df))){
        stop("Error, annotation is likely to be wrong...")
      }
      
      distance.df$distance.ratio <- distance.df$distance / distance.df$gene.width
      distance.df$annot <- annot
      
      return(distance.df)
    }

    ### Till here, method specific function.

    utr3.range.df <- rangeDuplexDf(object, t.gtf, "utr3")
    CDS.range.df <- rangeDuplexDf(object, t.gtf, "CDS")

    if(print.flag){
      cat(paste("N for CDS is:", nrow(CDS.range.df), "\n"))
      cat(paste("N for 3' UTR is:", nrow(utr3.range.df), "\n"))
    }
    
    merged.df <- rbind(CDS.range.df, utr3.range.df)
    
    merged.df <- merged.df[order(merged.df$annot, merged.df$distance, decreasing = TRUE), ]
    
    dir.create(out.dir, showWarnings = FALSE)

    table.file <- paste(out.dir, "unexpanded_range_summary.tab", sep = "/")
    write.table(merged.df, table.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

    if(print.flag){
      cat("Mann-Whitney U test for the range of RNA structures \n in 3' UTR and CDS.\n")
      raw.distance.wilcox <- wilcox.test(merged.df$distance[merged.df$annot == "CDS"], merged.df$distance[merged.df$annot == "utr3"])
      # print(raw.distance.wilcox)
      print(paste("p-value is:",
                  raw.distance.wilcox$p.value
                  ))

      plot1 <- ggplot(merged.df) + geom_density(aes(x = log(distance, 10), fill = annot), alpha = 0.3) + theme(legend.position="top")
      
    
      cat("Mann-Whitney U test for the normalized range \n(normalized by the width of gene segment) \nof RNA structures in 3' UTR and CDS.\n")
      ratio.distance.wicox <- wilcox.test(merged.df$distance.ratio[merged.df$annot == "CDS"], merged.df$distance.ratio[merged.df$annot == "utr3"])
      # print(ratio.distance.wicox)
      print(paste("p-value for ratio is: ",
                ratio.distance.wicox$p.value
                  ))
      
      
      plot2 <- ggplot(merged.df) + geom_density(aes(x = log(distance.ratio, 10), fill = annot), alpha = 0.3) + theme(legend.position="top")
      
      multiplot(plot1, plot2, cols = 2)
    }

    
  }
  )


