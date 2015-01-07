source("~/R_program/exec/transcriptome_mapping/functions/annotation_functions.R")

library(GenomicRanges)
library(rtracklayer)
library(ggplot2)


drawRNAmaps <- function(RPFile, mRNASeqFile, startCodons, stopCodons, gtfGenes, file.name, searchRange = 50){
  
  print("load RP data")
  RP <- import_nonHybrid(RPFile)
  print("load mRNA-Seq data")
  mRNASeq <- import_nonHybrid(mRNASeqFile)
  
  print("draw RNAmap for start codon")
  p1 <- drawRNAmap(RP, mRNASeq, gtfGenes, startCodons, c(2, 3, 1)[3 - searchRange %% 3], searchRange)
  print("draw RNAmap for stop codon")
  p2 <- drawRNAmap(RP, mRNASeq, gtfGenes, stopCodons, c(1, 2, 3)[3 - searchRange %% 3], searchRange)
  print("plotting")

  pdf(paste(file.name, "/", unlist(strsplit(RPFile, "/"))[length(unlist(strsplit(RPFile, "/")))], "_start.pdf", sep = ""), width = 7, height = 7, useDingbats = FALSE)
  print(p1 + labs(title = unlist(strsplit(RPFile, "/"))[length(unlist(strsplit(RPFile, "/")))]))
  dev.off()

  pdf(paste(file.name, "/", unlist(strsplit(RPFile, "/"))[length(unlist(strsplit(RPFile, "/")))], "_stop.pdf", sep = ""), width = 7, height = 7, useDingbats = FALSE)
  print(p2 + labs(title = unlist(strsplit(RPFile, "/"))[length(unlist(strsplit(RPFile, "/")))]))
  dev.off()
}

RNAmapDf <- function(RP, gtfGenes, startCodons, searchRange = 50){
  # RP is GRange object of risbosome profiling.
  # gtfGenes is GRangeList object of gtf files
  # in order to increase speed for iteration, these files should be GRange files and pre-procesed
  # extract position of start codons
  normalisationFactor <- sum(elementMetadata(RP)$score)

  RP.seqlevels <- as.character(seqlevels(RP))
  seqlevels.in.use <- RP.seqlevels[RP.seqlevels %in% as.character(seqlevels(gtfGenes))]
  
  RP <- RP[seqnames(RP) %in% seqlevels.in.use]
  seqlevels(RP) <- seqlevels.in.use

  gtfGenes <- gtfGenes[seqnames(gtfGenes) %in% seqlevels.in.use]
  seqlevels(gtfGenes) <- seqlevels.in.use

  startCodons <- startCodons[start(startCodons) - searchRange > 0, ]
  startCodons <- startCodons[seqnames(startCodons) %in% seqlevels.in.use]
  seqlevels(startCodons) <- seqlevels.in.use
                               
  start(startCodons) <- start(startCodons) - searchRange
  end(startCodons) <- end(startCodons) + searchRange
  
  overlaps=as.matrix(findOverlaps(startCodons, RP))
  flank.size <- searchRange
  overlaps=as.data.frame(overlaps)
  overlaps$clust.start=start(RP)[overlaps[,2]]
  overlaps$strand=as.vector(strand(RP))[overlaps[,2]]
  overlaps$score=elementMetadata(RP)$score[overlaps[,2]]
  overlaps$start_codon=start(startCodons)[overlaps[,1]]+flank.size
  
  # get distances to exon and separate according to strand
  distances=mapply(function(x, y, z){ rep(x-z, y) }, overlaps$clust.start, overlaps$score, overlaps$start_codon)
  #distances=overlaps$clust.start-overlaps$exon.start
  distances=split(distances, overlaps$strand)
  
  # split by strands
  #distances=c(distances[["+"]], (-1)*distances[["-"]])
  distances=c(unlist(distances[["+"]]), (-1)*unlist(distances[["-"]]))
  
  # keep only those that are within the range (for clusters that might extent over the region end)
  leftMost <- searchRange
  rightMost <- searchRange
  distances=distances[distances %in% c(-leftMost:(rightMost))]
  
  # determine frequencies
  distances=as.data.frame( table(distances) )
  names(distances)=c("rel.pos","xnts")
  
  # get plotting positions (ensure that all positions are present)
  pos=c(-leftMost:(rightMost-1))
  distances=data.frame(rel.pos=pos, 
                       xnts=distances$xnts[match(pos,distances$rel.pos)])
  distances$xnts[is.na(distances$xnts)]=0
  
  # get relative counts
  distances$rel.count = 1000000 * distances$xnts / normalisationFactor
  
  distances
}

drawRNAmap <- function(RP, mRNASeq, gtfGenes, startCodons, codonPos, searchRange = 50){
  RPDistances <- RNAmapDf(RP, gtfGenes, startCodons, searchRange)
  colnames(RPDistances) <- c("rel.pos", "RP.xnts", "RP.rel.count")
  mRNASeqDistances <- RNAmapDf(mRNASeq, gtfGenes, startCodons, searchRange)
  colnames(mRNASeqDistances) <- c("rel.pos", "mRNASeq.xnts", "mRNASeq.rel.count")
  
  codonPat <- c(1, 2, 3, 1, 2, 3)
  
  mergedDistances <- merge(RPDistances, mRNASeqDistances)
  
  mergedDistances$color <- factor(rep(codonPat[codonPos:(codonPos + 2)], nrow(mergedDistances)/3 + 1)[1:nrow(mergedDistances)], levels = 1:3)
  
  plotRNAmap <- ggplot(mergedDistances) +
    geom_line(aes(x = rel.pos, y = RP.rel.count)) +
      geom_point(aes(x = rel.pos, y = RP.rel.count, color = color)) +
        geom_line(aes(x = rel.pos, y = mRNASeq.rel.count)) +
          coord_cartesian(ylim = c(-50, 7000))
}

extractStart <- function(gtf){
	gtf <- sort(gtf)
	cds <- gtf[elementMetadata(gtf)$annot == "CDS"]
	startCodons <- cds
	end(startCodons) <- start(startCodons)
	
        return(startCodons)
}

extractStop <- function(gtf){
	gtf <- sort(gtf)
	cds <- gtf[elementMetadata(gtf)$annot == "CDS"]
	stopCodons <- cds
	start(stopCodons) <- end(stopCodons)
        
	return(stopCodons)
}

extractSeqnames <- function(GrangeGene){
  seqnames <- unique(as.character(seqnames(GrangeGene)))
  seqnames
}

extractStrand <- function(GrangeGene){
  strand <- unique(as.character(strand(GrangeGene)))
  strand
}

import_iCLIP <- function(file.name, genome="hg19"){
  # from Kathi
  xnts=import.bed(file.name,genome=genome,asRangedData=FALSE)
  elementMetadata(xnts)$score=as.numeric(elementMetadata(xnts)$name)
  strand(xnts)=ifelse(elementMetadata(xnts)$score >=0, "+", "-")
  elementMetadata(xnts)$score=abs(elementMetadata(xnts)$score)
  elementMetadata(xnts)$name=NULL
  return(xnts)
}

RNAmapDfAadapt <- function(RP, startCodons, normalisationFactor, searchRange = 50){
  # RP is GRange object of risbosome profiling.
  # gtfGenes is GRangeList object of gtf files
  # in order to increase speed for iteration, these files should be GRange files and pre-procesed
  # extract position of start codons
  
  start(startCodons) <- start(startCodons) - searchRange
  end(startCodons) <- end(startCodons) + searchRange
  
  overlaps=as.matrix(findOverlaps(startCodons, RP))
  flank.size <- searchRange
  overlaps=as.data.frame(overlaps)
  overlaps$clust.start=start(RP)[overlaps[,2]]
  overlaps$strand=as.vector(strand(RP))[overlaps[,2]]
  overlaps$score=elementMetadata(RP)$score[overlaps[,2]]
  overlaps$start_codon=start(startCodons)[overlaps[,1]]+flank.size
  
  # get distances to exon and separate according to strand
  distances=mapply(function(x, y, z){ rep(x-z, y) }, overlaps$clust.start, overlaps$score, overlaps$start_codon)
  #distances=overlaps$clust.start-overlaps$exon.start
  distances=split(distances, overlaps$strand)
  
  # split by strands
  #distances=c(distances[["+"]], (-1)*distances[["-"]])
  distances=c(unlist(distances[["+"]]), (-1)*unlist(distances[["-"]]))
  
  # keep only those that are within the range (for clusters that might extent over the region end)
  leftMost <- searchRange
  rightMost <- searchRange
  distances=distances[distances %in% c(-leftMost:(rightMost))]
  
  # determine frequencies
  distances=as.data.frame( table(distances) )
  names(distances)=c("rel.pos","xnts")
  
  # get plotting positions (ensure that all positions are present)
  pos=c(-leftMost:(rightMost-1))
  distances=data.frame(rel.pos=pos, 
                       xnts=distances$xnts[match(pos,distances$rel.pos)])
  distances$xnts[is.na(distances$xnts)]=0
  
  # get relative counts
  distances$rel.count = distances$xnts / normalisationFactor
  
  distances
}

