lenGr <- function(gr){
  gr.list <- split(gr, elementMetadata(gr)$gene_id)
  gr.len <- sapply(gr.list, FUN = function(x){sum(width(x))})
  return(gr.len)
}

write.IGV.GTF <- function(GR, file=""){
  # Modified from Quantify RNA-Seq data. version 0.01.
  GRtable <- cbind(seqname=as.character(seqnames(GR)), start=as.character(start(GR)), end=as.character(end(GR)), strand=as.character(strand(GR)), as.data.frame(GR@elementMetadata))

  idx <- match(c("seqname", "source", "feature", "start", "end", "score", "strand", "frame"), colnames(GRtable))
  idx1 <- na.omit(idx)
  attributes <- as.matrix(GRtable[,-idx1])

  atts <- paste('gene_id "', as.character(seqnames(GR)), '"; transcript_id "', as.character(seqnames(GR)), '"', sep = "")
  
  GRTable <- matrix(".", nrow=nrow(GRtable), ncol=9)
  if(is.null(na.action(idx1))){
    GRTable[,1:8] <- as.matrix(GRtable[,idx1])
  }else{
    GRTable[,c(1:8)[-na.action(idx1)]] <- as.matrix(GRtable[,idx1])
  }
  GRTable[,9] <- atts
  GRTable[,8] <- ifelse(elementMetadata(GR)$annot != "exon", "0", "")
  GRTable[,3] <- elementMetadata(GR)$annot
  GRTable[,2] <- "protein_coding"
  
  colnames(GRTable) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

  write.table(GRTable, file=file, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
}