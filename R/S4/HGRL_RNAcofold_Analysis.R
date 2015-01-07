##############################################################
#' RNAcoFoldAnalysis
#'
#' RNAcoFoldAnalysis performs folding energy analysis of hybrid reads
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' @param \code{all.gtf}. Transcripts annotation including ncRNAs, as transcriptome coordinate GenomicRanges object.
#' @param \code{fasta}. ShortRead objet of transcript sequence from fasta file.
#' @param \code{out.dir}. The directory for output file.
#' @param \code{hiCdir}. The directory of the package.
#' @param \code{pdir}. Python directory.
#' @param \code{rdir}. R directory.
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' RNAcoFoldAnalysis(hgrl, all.gtf, fasta, out.dir, hiCdir, pdir, rdir))
setGeneric(
  name = "RNAcoFoldAnalysis",
  def = function(object, all.gtf, fasta, out.dir, hiCdir, pdir, rdir){standardGeneric("RNAcoFoldAnalysis")}
  )

setMethod(
  f = "RNAcoFoldAnalysis",
  signature = "HybridGRL",
  definition = function(object, all.gtf, fasta, out.dir, hiCdir, pdir, rdir){
  	
    ### Method specific functions
    randomRepositioning <- function(grL, gtf){
      set.seed(1)
      grL.seqlevels <- unique(c(
        as.character(seqnames(grL$L)),
        as.character(seqnames(grL$R))
        ))
      
      gtf <- gtf[seqnames(gtf) %in% grL.seqlevels]
      seqlevels(gtf) <- grL.seqlevels
      
      seqlevels(grL) <- grL.seqlevels
      
      ol.L <- as.matrix(findOverlaps(grL$L, gtf))
      ol.R <- as.matrix(findOverlaps(grL$R, gtf))
      
      ol.L <- ol.L[!duplicated(ol.L[, 1]), ]
      ol.R <- ol.R[!duplicated(ol.R[, 1]), ]
      
      r.grL <- GRangesList(L = GRanges(), R = GRanges())
      seqlevels(r.grL) <- grL.seqlevels
      
      temp.L <- randomRepositioningDf(grL$L, gtf, ol.L)
      seqlevels(temp.L) <- grL.seqlevels
      r.grL$L <- temp.L
      
      temp.R <- randomRepositioningDf(grL$R, gtf, ol.R)
      seqlevels(temp.R) <- grL.seqlevels
      r.grL$R <- temp.R

      r.hgrL <- new("HybridGRL", r.grL)
      return(r.hgrL)
    }
    
    randomRepositioningDf <- function(gr, gtf, ol){
      df <- data.frame(
        gene_id = as.character(seqnames(gtf[ol[, 2]])),
        start = 0,
        width = width(gr[ol[, 1]])
        )
      
      for(i in 1:nrow(ol)){
        df$start[i] <- sample(
          start(gtf[ol[i, 2]]):(end(gtf[ol[i, 2]]) - width(gr[ol[i, 1]]) + 1),
          1
          )
      }
      
      r.gr <- GRanges(seqnames = Rle(df$gene_id), 
                      ranges = IRanges(df$start,
                        width = df$width), 
                      strand = Rle(rep("+", nrow(df))),
                      rname = 1:nrow(df),
                      score = as.integer(1),
                      category = "unknown",
                      annot = "unknown",
                      sequence = "unknown"
                      )
      return(r.gr)
    }
    
    generateFastaFromHybrid <- function(grL, fasta.file){
      
      fasta.both <- ShortRead(
        sread = DNAStringSet(c(elementMetadata(grL$L)$sequence,
          elementMetadata(grL$R)$sequence)),
        id = BStringSet(c(paste(elementMetadata(grL$L)$rname, "L", sep = "_"),
          paste(elementMetadata(grL$R)$rname, "R", sep = "_")
          )
          )
        )
      
      fasta.sorted <- sortFastaByReadname(fasta.both)
      writeFasta(fasta.sorted, fasta.file, width = 200)
    }
    
    sortFastaByReadname <- function(fasta){
      fasta <- fasta[order(sapply(strsplit(as.character(id(fasta)), "_"), "[[", 2))]
      fasta <- fasta[order(as.integer(sapply(strsplit(as.character(id(fasta)), "_"), "[[", 1)))]
      return(fasta)
	}


  ### Method specific functions, till here
  
  random.grL <- randomRepositioning(object, all.gtf)
  random.grL <- addSeqHgrl(random.grL, fasta)
  object <- addSeqHgrl(object, fasta)

  temp.file <- paste(out.dir, 'temp', sep = "/")
  dir.create(file.path(out.dir), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(temp.file), showWarnings = FALSE, recursive = TRUE)
  
  hybrid.fasta.file <- paste(temp.file, "grL.fa", sep = "/")
  random.fasta.file <- paste(temp.file, "random.grL.fa", sep = "/")
  coFoldOutDir <- out.dir

  generateFastaFromHybrid(object, hybrid.fasta.file)
  generateFastaFromHybrid(random.grL, random.fasta.file)

  cofoldCommand <- paste(pdir, "/python ", hiCdir, "/inst/Python/RNAcofoldAnalysis.py -i ", hybrid.fasta.file, " -n ", random.fasta.file, " -o ",  coFoldOutDir, " -d ", hiCdir, " -r ", rdir, sep = "")
  # print(cofoldCommand)
  system(cofoldCommand)


  # Start plot from here
  hybridCofoldOut <- gsub("temp/grL.fa", "RNAcofold/cofoldOut/Hybrid_grL.cofoldOut", hybrid.fasta.file)
  nonHybridCofoldOut <- gsub("/temp/random.grL.fa", "/RNAcofold/cofoldOut/nonHybrid_random.grL.cofoldOut", random.fasta.file)
  
  require(ggplot2)
  
  extractMFE <- function(RNAcofoldOut){
	table <- read.table(RNAcofoldOut, fill = TRUE, header = FALSE, sep = "\n")
	table <- table[seq(3, nrow(table), 3), ]
	table <- data.frame(do.call('rbind', strsplit(as.character(table), " (", fixed = TRUE)))
	MFE <- table[, 2]
	MFE <- gsub(")", "", MFE)
	MFE <- as.numeric(MFE)
	MFE
  }
  
  hybridMFE <- extractMFE(hybridCofoldOut)
  nonHybridMFE <- extractMFE(nonHybridCofoldOut)
  
  wtR <- wilcox.test(hybridMFE, nonHybridMFE, alternative = "two.sided")
  print(wtR)
  print(wtR[3])
  pVal <- as.numeric(wtR[3])

  print(paste("N (hybrid)", length(hybridMFE)))
  print(paste("N (non-hybrid)", length(nonHybridMFE)))

  df <- data.frame(variable = factor(c(rep("hybrid", length(hybridMFE)), rep("non-hybrid", length(nonHybridMFE)))), value = c(hybridMFE, nonHybridMFE))

  print(ggplot(df, aes(value, fill = variable)) + geom_density(alpha = 0.2) + labs(title = paste("two-sample Wilcoxon tests, p-val: ", pVal, sep ="")) + xlab("Ensemble free energy [kcal/mol]") + coord_cartesian(xlim = c(-80, 5), ylim = c(-0.005, 0.08)))

  }
  )

##############################################################
#' RNAcoFoldAnalysis2
#'
#' RNAcoFoldAnalysis2 performs folding energy analysis of hybrid reads
#' Unlike RNAcoFoldAnalysis, shuffled sequences were used as the control
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' @param \code{all.gtf}. Transcripts annotation including ncRNAs, as transcriptome coordinate GenomicRanges object.
#' @param \code{fasta}. ShortRead objet of transcript sequence from fasta file.
#' @param \code{out.dir}. The directory for output file.
#' @param \code{hiCdir}. The directory of the package.
#' @param \code{pdir}. Python directory.
#' @param \code{rdir}. R directory.
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' RNAcoFoldAnalysis2(hgrl, all.gtf, fasta, out.dir, hiCdir, pdir, rdir))
setGeneric(
  name = "RNAcoFoldAnalysis2",
  def = function(object, all.gtf, fasta, out.dir, hiCdir, pdir, rdir){standardGeneric("RNAcoFoldAnalysis2")}
  )

setMethod(
  f = "RNAcoFoldAnalysis2",
  signature = "HybridGRL",
  definition = function(object, all.gtf, fasta, out.dir, hiCdir, pdir, rdir){
    generateFastaFromHybrid <- function(grL, fasta.file){
      
      fasta.both <- ShortRead(
        sread = DNAStringSet(c(elementMetadata(grL$L)$sequence,
          elementMetadata(grL$R)$sequence)),
        id = BStringSet(c(paste(elementMetadata(grL$L)$rname, "L", sep = "_"),
          paste(elementMetadata(grL$R)$rname, "R", sep = "_")
          )
          )
        )
      
      fasta.sorted <- sortFastaByReadname(fasta.both)
      writeFasta(fasta.sorted, fasta.file, width = 200)
    }
    
    sortFastaByReadname <- function(fasta){
      fasta <- fasta[order(sapply(strsplit(as.character(id(fasta)), "_"), "[[", 2))]
      fasta <- fasta[order(as.integer(sapply(strsplit(as.character(id(fasta)), "_"), "[[", 1)))]
      return(fasta)
	}


  ### Method specific functions, till here
  object <- addSeqHgrl(object, fasta)

  temp.file <- paste(out.dir, 'temp', sep = "/")
  dir.create(file.path(out.dir), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(temp.file), showWarnings = FALSE, recursive = TRUE)
  
  hybrid.fasta.file <- paste(temp.file, "grL.fa", sep = "/")
  random.fasta.file <- paste(temp.file, "random.grL.fa", sep = "/")
  
  generateFastaFromHybrid(object, hybrid.fasta.file)
  
  randomizationCommand <- paste(pdir, "/python ", hiCdir, "/inst/Python/shuffle_RNAcoFold.py ", hybrid.fasta.file, " ", random.fasta.file, sep = "")
  system(randomizationCommand)
  
  random.fasta <- readFasta(random.fasta.file)
  random.sorted.fasta <- sortFastaByReadname(random.fasta)
  unlink(random.fasta.file)
  writeFasta(random.sorted.fasta, random.fasta.file, width = 200)
 
  coFoldOutDir <- out.dir

  cofoldCommand <- paste(pdir, "/python ", hiCdir, "/inst/Python/RNAcofoldAnalysis.py -i ", hybrid.fasta.file, " -n ", random.fasta.file, " -o ",  coFoldOutDir, " -d ", hiCdir, " -r ", rdir, sep = "")
  # print(cofoldCommand)
  system(cofoldCommand)


  # Start plot from here
  hybridCofoldOut <- gsub("temp/grL.fa", "RNAcofold/cofoldOut/Hybrid_grL.cofoldOut", hybrid.fasta.file)
  nonHybridCofoldOut <- gsub("/temp/random.grL.fa", "/RNAcofold/cofoldOut/nonHybrid_random.grL.cofoldOut", random.fasta.file)
  
  require(ggplot2)
  
  extractMFE <- function(RNAcofoldOut){
	table <- read.table(RNAcofoldOut, fill = TRUE, header = FALSE, sep = "\n")
	table <- table[seq(3, nrow(table), 3), ]
	table <- data.frame(do.call('rbind', strsplit(as.character(table), " (", fixed = TRUE)))
	MFE <- table[, 2]
	MFE <- gsub(")", "", MFE)
	MFE <- as.numeric(MFE)
	MFE
  }
  
  hybridMFE <- extractMFE(hybridCofoldOut)
  nonHybridMFE <- extractMFE(nonHybridCofoldOut)
  
  wtR <- wilcox.test(hybridMFE, nonHybridMFE, alternative = "two.sided")
  print(wtR)
  print(wtR[3])
  pVal <- as.numeric(wtR[3])

  print(paste("N (hybrid)", length(hybridMFE)))
  print(paste("N (non-hybrid)", length(nonHybridMFE)))

  df <- data.frame(variable = factor(c(rep("hybrid", length(hybridMFE)), rep("non-hybrid", length(nonHybridMFE)))), value = c(hybridMFE, nonHybridMFE))

  print(ggplot(df, aes(value, fill = variable)) + geom_density(alpha = 0.2) + labs(title = paste("two-sample Wilcoxon tests, p-val: ", pVal, sep ="")) + xlab("Ensemble free energy [kcal/mol]") + coord_cartesian(xlim = c(-80, 5), ylim = c(-0.005, 0.08)))

  }
  )
