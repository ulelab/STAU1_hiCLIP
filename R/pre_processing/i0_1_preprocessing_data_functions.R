# Modified from Quantify RNA-Seq data. version 0.01 {https://r-forge.r-project.org/scm/viewvc.php/pkg/R/?root=qrnaseq}

loadGTF <- function(file, feature=c("exon", "CDS", "intron", "utr"), attributes=c("gene_id", "transcript_id", "gene_name", "transcript_name", "protein_id")){
  require(GenomicRanges)
  
  cat("load GTF file ... \n")
  GTF <- read.table(file, sep="\t", header=FALSE, stringsAsFactors=FALSE, colClasses=c(rep("character", 3), rep("numeric", 2), rep("character", 4)))
  cat("parse attributes ... \n")
  attrs <- GTF[,9]
  attrs <- strsplit(attrs, split=" |;")
  for(i in 1:length(attrs)){
    ai <- match(attributes, attrs[[i]])
    if(sum(is.na(ai))==0){
      break
    }
  }
  Attrs <- do.call("cbind", lapply(attrs, function(x)x[ai+1]))
  Attrs <- t(Attrs)
  colnames(Attrs) <- attributes
  GR <- GRanges(seqnames=GTF[,1], IRanges(start=GTF[,4], end=GTF[,5]), strand=GTF[,7], source=GTF[,2], feature=GTF[,3], score=GTF[,6], frame=GTF[,8], Attrs)
  
  ## add transcript.gene_id
  elementMetadata(GR)$transcript.gene_id <- paste(elementMetadata(GR)$transcript_id, elementMetadata(GR)$gene_id, sep = "_")

  return(GR)
}

write.GTF <- function(GR, file=""){
  GRtable <- cbind(seqname=as.character(seqnames(GR)), start=as.character(start(GR)), end=as.character(end(GR)), strand=as.character(strand(GR)), as.data.frame(GR@elementMetadata))
  idx <- match(c("seqname", "source", "feature", "start", "end", "score", "strand", "frame"), colnames(GRtable))
  idx1 <- na.omit(idx)
  attributes <- as.matrix(GRtable[,-idx1])
  atn <- colnames(attributes)

  atts <- c()
  for(i in 1:ncol(attributes)){
    atts1 <- paste(atn[i], attributes[,i])
    atts <- paste(atts, atts1, sep="; ")
  }
  atts <- sub("^;", "", atts)
  atts <- sub("\\s*\\w*\\sNA;?", "", atts)
  
  if(length(grep("gene_id", atts)) == length(atts)){
  	cat("GTF has full attributes...\n")
  	atts <- substr(atts, 2, nchar(atts))
  } else {
  	cat("GTF has limited attributes...\n")
  	atts <- paste('gene_id ', substr(atts, 3, nchar(atts)), '', sep = "")
  }

  GRTable <- matrix(".", nrow=nrow(GRtable), ncol=9)
  if(is.null(na.action(idx1))){
    GRTable[,1:8] <- as.matrix(GRtable[,idx1])
  }else{
    GRTable[,c(1:8)[-na.action(idx1)]] <- as.matrix(GRtable[,idx1])
  }
  GRTable[,9] <- atts

  # Added on 2013/Jan/13th
  if(!("feature" %in% names(elementMetadata(GR)))){
    GRTable[,3] <- "exon"
  }
  # Till here
  
  colnames(GRTable) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
  write.table(GRTable, file=file, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
}


diffGR <- function(gr1, gr2, by="gene_id"){
  idx <- match(by, colnames(gr1@elementMetadata))
  GR1 <- GRanges(paste(seqnames(gr1), gr1@elementMetadata[,idx], sep="|"), ranges(gr1), strand(gr1))
  GR2 <- GRanges(paste(seqnames(gr2), gr2@elementMetadata[,idx], sep="|"), ranges(gr2), strand(gr2))
  GRd <- suppressWarnings(setdiff(GR1, GR2))
  chrid <- do.call("rbind", strsplit(as.character(seqnames(GRd)), split="\\|"))
  GRD <- GRanges(seqnames=chrid[,1], ranges=ranges(GRd), strand(GRd), chrid[,2])
  colnames(GRD@elementMetadata) <- by
  return(GRD)
}

reduceGR.transcript.gene_id <- function(GR){
  GR.ids.temp <- GRanges(seqnames=paste(seqnames(GR), sapply(strsplit(GR@elementMetadata$transcript.gene_id, "_"), "[", 2), sep="|"), ranges=ranges(GR), strand=strand(GR))
  GR.ids <- reduce(GR.ids.temp)
  chrid <- do.call("rbind", strsplit(as.character(seqnames(GR.ids)), split="\\|"))
  GR.ids <- GRanges(seqnames=chrid[,1], ranges=ranges(GR.ids), strand=strand(GR.ids), gene_id=chrid[,2])
  return(GR.ids)
}

selectLongestGenes <- function(file){
  # file <- "~/Desktop/gtf_test/chr22.gtf"
  # feature=c("exon", "CDS", "intron", "utr")
  # attributes=c("gene_id", "transcript_id", "gene_name", "transcript_name", "protein_id")

  require(GenomicRanges)

  GR <- loadGTF(file)
  GR.exon <- GR[elementMetadata(GR)$feature == "exon"]
  
  
  # Calculate the length of each transcript
  transcript.width <- aggregate(width(GR.exon), by = list(elementMetadata(GR.exon)$transcript_id), FUN = sum)
  transcript.width$Group.1 <- as.character(transcript.width$Group.1)

  # Add gene_id
  t_g.id.dict <- createDict(key = elementMetadata(GR)$transcript_id[!duplicated(elementMetadata(GR)$transcript_id)], value = elementMetadata(GR)$gene_id[!duplicated(elementMetadata(GR)$transcript_id)])

  transcript.width$gene_id <- t_g.id.dict[transcript.width$Group.1]

  # Find longest transcript among the same gene_id
  transcript.width <- transcript.width[order(transcript.width$gene_id, -1 * transcript.width$x), ]
  longest.transcript_id <- transcript.width$Group.1[!duplicated(transcript.width$gene_id)]

  GR.longest <- GR[elementMetadata(GR)$transcript_id %in% longest.transcript_id]
  
  longest.transcript_id.df <- data.frame(gene_id = elementMetadata(GR.longest)$gene_id, transcript_id = elementMetadata(GR.longest)$transcript_id)

  longest.transcript_id.df <- longest.transcript_id.df[!duplicated(longest.transcript_id.df$gene_id), ]
 
  GR.longest.list <- list(GR = GR.longest, id.df = longest.transcript_id.df)

  return(GR.longest.list)
}

createDict <- function(key = key, value = value){
  dict <- value
  names(dict) <- key
  if(any(duplicated(key))){stop("key is duplicated.")}
  return(dict)
}
