pasteC <- function(vec){
  vec <- paste(vec, collapse = "_")
}


selectLongest <- function(hybrid, prefix){
  hybrid$gi_annot <- apply(hybrid[, c("gene_id", "annot")], 1, pasteC)
  hybrid <- hybrid[order(hybrid$gi_annot, hybrid$distance, decreasing = TRUE), ]
  hybrid <- hybrid[!duplicated(hybrid$gi_annot),]
  hybrid <- hybrid[, c("gene_id", "distance", "score", "gene.width", "distance.ratio", "annot")]

  CDS <- hybrid[hybrid$annot == "CDS", ]
  colnames(CDS) <- paste(prefix, "CDS", colnames(CDS), sep = ".")
  colnames(CDS)[1] <- "gene_id"
  
  utr3 <- hybrid[hybrid$annot == "utr3", ]
  colnames(utr3) <- paste(prefix, "utr3", colnames(utr3), sep = ".")
  colnames(utr3)[1] <- "gene_id"

  merged.hybrid <- merge(CDS, utr3, by = "gene_id", all = TRUE)
  return(merged.hybrid)
}
