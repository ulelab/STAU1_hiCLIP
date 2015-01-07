
##############################################################
#' shiftCD
#'
#' shiftCD was used to adjust the coordinate to the "3' UTR base" from "transcript base" 
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' exportHybrid(hgrl)
setGeneric(
  name = "shiftCD",
  def = function(object, gr.utr3){standardGeneric("shiftCD")}
  )

setMethod(
  f = "shiftCD",
  signature = "HybridGRL",
  definition = function(object, gr.utr3){

    start.vec <- start(gr.utr3)
    names(start.vec) <- as.character(seqnames(gr.utr3))
    
    shift.gr <- function(gr, start.vec){
      start(gr) <- start(gr) - start.vec[as.character(seqnames(gr))] + 1
      end(gr) <- end(gr) - start.vec[as.character(seqnames(gr))] + 1
      return(gr)
    }
    
    object$L <- shift.gr(object$L, start.vec)
    object$R <- shift.gr(object$R, start.vec)
    
    
    return(object)
  }
  )



##############################################################
#' selectRNAfoldPredictableGenes
#'
#' selectRNAfoldPredictableGenes was used to find whether a gene structure can be predicted by RNAfold with a constaint using hybrid data. 
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' exportHybrid(hgrl)
setGeneric(
  name = "selectRNAfoldPredictableGenes",
  def = function(object, selected = TRUE){standardGeneric("selectRNAfoldPredictableGenes")}
  )

setMethod(
  f = "selectRNAfoldPredictableGenes",
  signature = "HybridGRL",
  definition = function(object, selected = TRUE){

    ## Functions specific for this method
    sub.is.RNAfoldPredictableGenes <- function(temp.object){
      b = end(temp.object$L[1])
      c = start(temp.object$R[1])
      d = end(temp.object$R[1])
      
      o = end(temp.object$L[2])
      r = end(temp.object$R[2])
      
      Rbf <- FALSE
      
      if(d < o){
        Rbf <- TRUE
      } else {
        if((r < c) & (b < o)){
          Rbf <- TRUE
        } else {
          Rbf <- FALSE
        }
      }
      
      return(Rbf)
    }

    is.RNAfoldPredictableGenes <- function(object){
      if(length(object$L) == 1){
        b.res <- TRUE
      } else {
        
        n.dup <- length(object$L)
        n.comb <- combn(1:n.dup, 2)
        
        b.results <- c()
        for(i in 1:ncol(n.comb)){
          temp.object <- selectHybridByIndex(object, indexes = n.comb[, i])
          Rbf <- sub.is.RNAfoldPredictableGenes(temp.object)
          b.results <- c(b.results, Rbf)
        }
        
        b.res <- all(b.results)
      }
      
      return(b.res)
    }

    ## Caution: HGRL object should be sorted before running this function
    gene.vec <- unique(as.character(seqnames(object$L)))
    bf.genes <- c()
    
    for(i.g in gene.vec){
      tmp.object <- selectHybridByGeneName(object, i.g)
      tmp.bf.genes <- is.RNAfoldPredictableGenes(tmp.object)
      
      bf.genes <- c(bf.genes, tmp.bf.genes)
    }
    
    if(selected){
      selected.gene.vec <- gene.vec[bf.genes]
    } else {
      selected.gene.vec <- gene.vec[!bf.genes]
    }
    return(selected.gene.vec)

  }
  )


##############################################################
#' createDB 
#'
#' createDB returns structure constraint by hybrid
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' exportHybrid(hgrl)
setGeneric(
  name = "createDB",
  def = function(object, gr.utr3, filename){standardGeneric("createDB")}
  )

setMethod(
  f = "createDB",
  signature = "HybridGRL",
  definition = function(object, gr.utr3, filename){


    mergeDB <- function(vec1, vec2){
      if(length(vec1) != length(vec2)){
        stop("Vector length should be the same.")
      }
      
      left_elements <- lapply(strsplit(vec1, "\\_"), as.integer)
      right_elements <- lapply(strsplit(vec2, "\\_"), as.integer)
      
      sum_elements <- mapply("+", left_elements, right_elements, SIMPLIFY = FALSE)
      elements_conct <- sapply(sum_elements, function(x){paste(x, collapse = "_")})
      
      return(elements_conct)  
    }
    
    
    createDB <- function(gr, gr.utr3, bracket = 1){
      for(i in 1:length(gr)){
      	
      	utr3.logical <- as.character(seqnames(gr.utr3)) == as.character(seqnames(gr[i]))
      	utr3.length <- end(gr.utr3[utr3.logical]) - start(gr.utr3[utr3.logical]) + 1
      	
      	left <-paste(rep(0, (start(gr[i]) - 1)), collapse = "_")
      	mid <- paste(rep(bracket, (end(gr[i]) - start(gr[i]) + 1)), collapse = "_")      	
      	right <- paste(rep(0, (utr3.length - end(gr[i]))), collapse = "_")

        all.elements <- paste(c(left, mid, right), collapse = "_")
        all.elements <- gsub("^\\_", "", all.elements)
        all.elements <- gsub("\\_$", "", all.elements)
        
      	elementMetadata(gr)$DB[i] <- all.elements
      }
      
      return(gr)
    }
    
    createDB.df <- function(object){
      DB.df <- data.frame(
        gene_id = as.character(seqnames(object$L)),
        DB = "NA",
        stringsAsFactors = FALSE
        )
      
      elements_conct <- mergeDB(elementMetadata(object$L)$DB, elementMetadata(object$R)$DB)
      

      
      if(length(grep("3", elements_conct)) != 0){
    	stop("Conflicting duplexes exist")
      }

      DB.df$DB <- elements_conct
      return(DB.df)
    }
    
    compressDBdf <- function(DB.df){
      duplicated_id <- unique(
        DB.df$gene_id[duplicated(DB.df$gene_id)]
        )

      unique.df <- DB.df[!(DB.df$gene_id %in% duplicated_id), ]
      duplicated.df <- DB.df[DB.df$gene_id %in% duplicated_id, ]

      compressed.df <- data.frame(
        gene_id = unique(duplicated.df$gene_id),
        DB = "NA",
        stringsAsFactors = FALSE
        )
      
      for(gene in compressed.df$gene_id){
        temp.df <- duplicated.df[duplicated.df$gene_id %in% gene, ]
        merged.DB <- temp.df$DB[1]
        
        for(i in 1:(nrow(temp.df) - 1)){
          merged.DB <- mergeDB(merged.DB, temp.df$DB[i + 1])
        }
        
        compressed.df$DB[compressed.df$gene_id == gene] <- merged.DB
      }

      result.df <- rbind(unique.df, compressed.df)
      result.df <- result.df[order(result.df$gene_id), ]
      
      return(result.df)
    }
    
    
    convertIntoDB <- function(vec){
      temp_dp <- gsub("_", "", vec)
      temp_dp_1 <- gsub("0", ".", temp_dp)
      temp_dp_2 <- gsub("1", "(", temp_dp_1)
      dp_vec <- gsub("2", ")", temp_dp_2)
      
      return(dp_vec)
    } 
    
    object <- addColumnHGRL(object, "DB", default.value = "NA")
    object$L <- createDB(object$L, gr.utr3, 1)
    object$R <- createDB(object$R, gr.utr3, 2)  
    
    DB.df <- createDB.df(object)
    compressed.DB.df <- compressDBdf(DB.df)
    compressed.DB.df$DB <- convertIntoDB(compressed.DB.df$DB)
    
    filename.faconst <- paste(filename, "faconst", sep = ".")
    filename.const <- paste(filename, "const", sep = ".")

    sink(filename.faconst)
    for(i in 1:nrow(compressed.DB.df)){
      line.id <- paste(">", compressed.DB.df[i, 1], "\n")
      line.constrain <- paste(compressed.DB.df[i, 2], "\n")
      
      cat(line.id)
      cat(elementMetadata(utr3.selected)$sequence[as.character(seqnames(utr3.selected)) == compressed.DB.df[i, 1]])
      cat("\n")
      cat(line.constrain)
      
    }
    sink()
    
    
    sink(filename.const)
    for(i in 1:nrow(compressed.DB.df)){
      line.id <- paste(">", compressed.DB.df[i, 1], "\n")
      line.constrain <- paste(compressed.DB.df[i, 2], "\n")
      
      cat(line.id)
      cat(line.constrain)
      
    }
    sink()

    return(DB.df)
  }
  )
