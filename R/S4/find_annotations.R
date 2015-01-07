##############################################################
#' annotateHybrid
#'
#' annotateHybrid annotates hybrid reads.
#' @param \code{hgrl}. HybridGRL object to be examined.
#' @param \code{a.df}. Annotation data.frame contains catergory and annotation of each genes.
#' @param \code{t.gtf.gr}. Transcriptome coordinate of protein coding gene as GRanges object
#' @param \code{mRNA.with.intron.grL}. Genomic coordinate of genes with the annotation of intron as GRanges object. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' annotateHybrid(hgrl, a.df, t.gtf.gr, mRNA.with.intron.grL)
setGeneric(
  name = "annotateHybrid",
  def = function(object, a.df, t.gtf.gr, mRNA.with.intron.grL){standardGeneric("annotateHybrid")}
  )

setMethod(
  f = "annotateHybrid",
  signature = "HybridGRL",
  definition = function(object, a.df, t.gtf.gr, mRNA.with.intron.grL){
    ## Define Method specific functions
    annotateProteinCoding <- function(gr.input, t.gtf){
      gr <- gr.input
      sub.t.gtf <- t.gtf[seqnames(t.gtf) %in% unique(as.character(seqnames(gr)))]
      seqlevels(sub.t.gtf) <- unique(as.character(seqnames(sub.t.gtf)))
      seqlevels(gr) <- as.character(seqlevels(sub.t.gtf))
      gr.temp <- gr
      ## The next line was added to define the hybrid reads by the position of the first base.
      end(gr.temp) <- start(gr.temp)
      ol <- findOverlaps(gr.temp, sub.t.gtf, type = "any")
      ol <- as.matrix(ol)
      ol <- ol[!duplicated(ol[, 1]), ]
      elementMetadata(gr)$annot[ol[, 1]] <- as.character(elementMetadata(sub.t.gtf)$annot[ol[, 2]])
      return(gr)
    }

    addCategoryAnnot <- function(gr, a.df, t.gtf.gr, mRNA.with.intron.grL){
      
      annotation.df <- a.df
      t.gtf <- t.gtf.gr
      mRNA.with.intron <- mRNA.with.intron.grL
      
      if(!all(c("category", "annot") %in% names(elementMetadata(gr)))){
        stop("column have to be pre-prepared")
      }
      
      if(!all(
        c(
          is.character(annotation.df[, 1]),
          is.character(annotation.df[, 2]),
          is.character(annotation.df[, 3])
          )
        )
         )
        {
          stop("annotation data frame only accept character")
        }

      # 1st: annotate rRNA and tRNA
      elementMetadata(gr)$category[grep("^rRNA", seqnames(gr))] <- "rRNA"
      elementMetadata(gr)$annot[grep("^rRNA", seqnames(gr))] <- as.character(seqnames(gr))[grep("^rRNA", seqnames(gr))]
  
      elementMetadata(gr)$category[grep("^trna", seqnames(gr))] <- "tRNA"
      elementMetadata(gr)$annot[grep("^trna", seqnames(gr))] <- "tRNA"

      # 2nd: annotate protein_coding and ncRNAs
      ensg.category <- annotation.df$category
      names(ensg.category) <- annotation.df$gene_id
      
      ensg.annot <- annotation.df$annot
      names(ensg.annot) <- annotation.df$gene_id
      
      elementMetadata(gr)$category[grep("^ENSG", seqnames(gr))] <-  ensg.category[as.character(seqnames(gr[grep("^ENSG", seqnames(gr))]))]
      elementMetadata(gr)$annot[grep("^ENSG", seqnames(gr))] <-  ensg.annot[as.character(seqnames(gr[grep("^ENSG", seqnames(gr))]))]

      # elementMetadata(gr)$annot[elementMetadata(gr)$category == "protein_coding"]

      gr[elementMetadata(gr)$category == "protein_coding"] <- annotateProteinCoding(gr[elementMetadata(gr)$category == "protein_coding"], t.gtf)
      
      # if mapped to minus strand of gene => intergenic
      elementMetadata(gr)$category[(as.character(strand(gr)) == "-") & (elementMetadata(gr)$annot %in% c("protein_coding", "lncRNA", "other_ncRNAs", "miRNA"))] <- "intergenic"
      elementMetadata(gr)$annot[(as.character(strand(gr)) == "-") & (elementMetadata(gr)$annot %in% c("protein_coding", "lncRNA", "other_ncRNAs", "miRNA"))] <- "intergenic"

      # 3rd: genomic regions as intron or intergenic
      elementMetadata(gr)$annot[grep("^chr", seqnames(gr))] <- "intergenic"
      
      if(length(gr[grep("^chr", seqnames(gr))]) != 0){
      	gr[grep("^chr", seqnames(gr))] <- annotateProteinCoding(gr[grep("^chr", seqnames(gr))], mRNA.with.intron)
      }
      
      elementMetadata(gr)$category[grep("^chr", seqnames(gr))] <- elementMetadata(gr)$annot[grep("^chr", seqnames(gr))]
      
      return(gr)
    }
    ## Define Method specofic functions: Up to here

    object$L <- addCategoryAnnot(object$L, a.df, t.gtf.gr, mRNA.with.intron.grL)
    object$R <- addCategoryAnnot(object$R, a.df, t.gtf.gr, mRNA.with.intron.grL)

    validObject(object)
    return(object)
  }
  )

##############################################################
#' plotHybridRNASource
#'
#' plotHybridRNASource plot RNA source of hybrid reads 
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' plotHybridRNASource(hgrl)
setGeneric(
  name = "plotHybridRNASource",
  def = function(object){standardGeneric("plotHybridRNASource")}
  )

setMethod(
  f = "plotHybridRNASource",
  signature = "HybridGRL",
  definition = function(object){
    annotation.each.reads.df <- data.frame(
      category = c(elementMetadata(object$L)$category, elementMetadata(object$R)$category),
      annot = c(elementMetadata(object$L)$annot, elementMetadata(object$R)$annot)
      )
    
    annot.collapase.rRNA <- annotation.each.reads.df
    annot.collapase.rRNA$annot <- as.character(annot.collapase.rRNA$annot)
    annot.collapase.rRNA$annot[grep("^rRNA", annot.collapase.rRNA$annot)] <- "rRNA"

    annot.collapase.rRNA$annot <- factor(as.character(annot.collapase.rRNA$annot), levels = c("utr5", "CDS", "utr3", "protein_coding", "lncRNA", "miRNA", "other_ncRNAs", "rRNA", "tRNA", "intron", "intergenic"))
    
    cat("RNA source of hybrid reads [total]\n")
    print(table(as.character(annot.collapase.rRNA$annot)))
	cat("\n")

    cat("RNA source of hybrid reads [%]\n")
    print(round(
      prop.table(table(as.character(annot.collapase.rRNA$annot))) * 100
      , digits = 1)
          )
    cat("\n")


    print(ggplot(annot.collapase.rRNA) +
          geom_bar(aes(x = factor(1), fill = annot), width = 1) +
          coord_polar(theta = "y") +
          theme(axis.ticks = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y =  element_blank(),
                axis.text.x =  element_blank(),
                axis.title.x = element_blank(),
                panel.background = element_blank(),
                legend.title=element_blank()
                )
          )
          
  }
  )


##############################################################
#' plotHybridPairRNASource
#'
#' plotHybridPairRNASource is simlar to plotHybridRNASource, but the results are calcualted from the pair of RNA sources.
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' plotHybridPairRNASource(hgrl)
setGeneric(
  name = "plotHybridPairRNASource",
  def = function(object){standardGeneric("plotHybridPairRNASource")}
  )

setMethod(
  f = "plotHybridPairRNASource",
  signature = "HybridGRL",
  definition = function(object){
  	if(!all(seqnames(object$L) == seqnames(object$R))){
  	  stop("This function is only applicable to those mapped to intra_molecular duplexes")
  	}
  	
    annotation.each.reads.df <- data.frame(
      category = paste(elementMetadata(object$L)$category, elementMetadata(object$R)$category, sep = "_"),
      annot = paste(elementMetadata(object$L)$annot, elementMetadata(object$R)$annot, sep = "_")
      )
    
    plot.annotation.df <- as.data.frame(table(annotation.each.reads.df$annot))
    plot.annotation.df$Var2 <- hotfactor(plot.annotation.df, n = 3)
    
    name.dict <- as.character(plot.annotation.df$Var2)
    names(name.dict) <- as.character(plot.annotation.df$Var1)

    
    cat("RNA source of hybrid reads [total]\n")
    print(table(as.character(annotation.each.reads.df$annot)))
	cat("\n")

    cat("RNA source of hybrid reads [%]\n")
    print(round(
      prop.table(table(as.character(annotation.each.reads.df$annot))) * 100
      , digits = 1)
          )
    cat("\n")


	annotation.each.reads.df$annot_collapsed <- name.dict[as.character(annotation.each.reads.df$annot)]

    print(ggplot(annotation.each.reads.df) +
          geom_bar(aes(x = factor(1), fill = annot_collapsed), width = 1) +
          coord_polar(theta = "y") +
          theme(axis.ticks = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y =  element_blank(),
                axis.text.x =  element_blank(),
                axis.title.x = element_blank(),
                panel.background = element_blank(),
                legend.title=element_blank()
                ) 
          )
          
  }
  )


