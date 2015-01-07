##############################################################
#' RNAhybridAnalysis
#'
#' RNAhybridAnalysis performs folding energy analysis of hybrid reads using RNAhybrid program.
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' @param \code{all.gr}. Transcripts annotation including ncRNAs, as transcriptome coordinate GenomicRanges object.
#' @param \code{fasta}. ShortRead objet of transcript sequence from fasta file.
#' @param \code{out.dir}. The directory for output file.
#' @param \code{hiCdir}. The directory of the package.
#' @param \code{pdir}. Python directory.
#' 
#' @export
#' @docType methods
#' @rdname hybridGRL-methods
#'
#' @examples
#' RNAhybridAnalysis(hgrl, all.gr, fasta, out.dir, hiCdir, pdir))
setGeneric(
  name = "RNAhybridAnalysis",
  def = function(object, all.gr, fasta, out.dir, hiCdir, pdir){standardGeneric("RNAhybridAnalysis")}
  )

setMethod(
  f = "RNAhybridAnalysis",
  signature = "HybridGRL",
  definition = function(object, all.gr, fasta, out.dir, hiCdir, pdir){
  
  random.grL <- randomPermutatons(object, all.gr, 1)
  random.grL <- addSeqHgrl(random.grL, fasta)
  object <- addSeqHgrl(object, fasta)

  temp.file <- paste(out.dir, 'temp', sep = "/")
  dir.create(file.path(out.dir), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(temp.file), showWarnings = FALSE, recursive = TRUE)
  
  hybrid.file <- paste(temp.file, "grL", sep = "/")
  random.file <- paste(temp.file, "random.grL", sep = "/")
  
  hybrid.fasta.file <- paste(hybrid.file, ".fa", sep = "")
  random.fasta.file <- paste(random.file, ".fa", sep = "")
  
  
  coFoldOutDir <- out.dir

  df1 <- concatenateHGRL(object, 
                      	fasta, 
                      	hybrid.file, 
                      	sep.seq = "&")
                      	
  df2 <- concatenateHGRL(random.grL, 
                      	fasta, 
                      	random.file, 
                      	sep.seq = "&")
  
                     	
  hybrid.RNAhybrid.out <- paste(temp.file, "/hybrid.RNAhybrid.out", sep = "")                  	
  random_control.RNAhybrid.out <- paste(temp.file, "/random_control.RNAhybrid.out", sep = "")
   

  RNAhybridCommand <- paste(pdir, "/python ", hiCdir, "/inst/Python/runRNAhybrid_parallel.py -i ", hybrid.fasta.file, " -o ", hybrid.RNAhybrid.out, " -d ", hiCdir, sep = "")
  # print(RNAhybridCommand)
  system(RNAhybridCommand)
  
  RNAhybridCommand2 <- paste(pdir, "/python ", hiCdir, "/inst/Python/runRNAhybrid_parallel.py -i ", random.fasta.file, " -o ", random_control.RNAhybrid.out, " -d ", hiCdir, sep = "")
  # print(RNAhybridCommand2)
  system(RNAhybridCommand2)


  # Start plot from here
  require(ggplot2)
  
  hybridMFE <- readRNAhybrid(hybrid.RNAhybrid.out)
  nonHybridMFE <- readRNAhybrid(random_control.RNAhybrid.out)
  
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
#' RNAhybridAnalysis2
#'
#' RNAhybridAnalysis2 performs free energy analysis of hybrid reads using RNAhybrid program
#' Unlike RNAcoFoldAnalysis, shuffled sequences were used as the control
#' @param \code{hgrl}. HybridGRL object to be examined. 
#' @param \code{all.gr}. Transcripts annotation including ncRNAs, as transcriptome coordinate GenomicRanges object.
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
#' RNAhybridAnalysis2(hgrl, all.gr, fasta, out.dir, hiCdir, pdir))
setGeneric(
  name = "RNAhybridAnalysis2",
  def = function(object, all.gr, fasta, out.dir, hiCdir, pdir){standardGeneric("RNAhybridAnalysis2")}
  )

setMethod(
  f = "RNAhybridAnalysis2",
  signature = "HybridGRL",
  definition = function(object, all.gr, fasta, out.dir, hiCdir, pdir){
  
  random.grL <- randomPermutatons(object, all.gr, 1)
  random.grL <- addSeqHgrl(random.grL, fasta)
  object <- addSeqHgrl(object, fasta)

  temp.file <- paste(out.dir, 'temp', sep = "/")
  dir.create(file.path(out.dir), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(temp.file), showWarnings = FALSE, recursive = TRUE)
  
  hybrid.file <- paste(temp.file, "grL", sep = "/")
  random.file <- paste(temp.file, "random.grL", sep = "/")
  
  hybrid.fasta.file <- paste(hybrid.file, ".fa", sep = "")
  random.fasta.file <- paste(random.file, ".fa", sep = "")
  
  
  coFoldOutDir <- out.dir

  df1 <- concatenateHGRL(object, 
                      	fasta, 
                      	hybrid.file, 
                      	sep.seq = "&")
                      	

  randomizationCommand <- paste(pdir, "/python ", hiCdir, "/inst/Python/shuffle_RNAcoFold.py ", hybrid.fasta.file, " ", random.fasta.file, sep = "")
  system(randomizationCommand)
                     	
  hybrid.RNAhybrid.out <- paste(temp.file, "/hybrid.RNAhybrid.out", sep = "")                  
  random_control.RNAhybrid.out <- paste(temp.file, "/random_control.RNAhybrid.out", sep = "")

  RNAhybridCommand <- paste(pdir, "/python ", hiCdir, "/inst/Python/runRNAhybrid_parallel.py -i ", hybrid.fasta.file, " -o ", hybrid.RNAhybrid.out, " -d ", hiCdir, sep = "")
  # print(RNAhybridCommand)
  system(RNAhybridCommand)
  
  RNAhybridCommand2 <- paste(pdir, "/python ", hiCdir, "/inst/Python/runRNAhybrid_parallel.py -i ", random.fasta.file, " -o ", random_control.RNAhybrid.out, " -d ", hiCdir, sep = "")
  # print(RNAhybridCommand2)
  system(RNAhybridCommand2)


  # Start plot from here
  require(ggplot2)
  
  hybridMFE <- readRNAhybrid(hybrid.RNAhybrid.out)
  nonHybridMFE <- readRNAhybrid(random_control.RNAhybrid.out)
  
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
