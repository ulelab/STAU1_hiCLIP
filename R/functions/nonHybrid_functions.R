import_iCLIP <- function(file.name, genome="hg19"){
  xnts=import.bed(file.name,genome=genome,asRangedData=FALSE)
  elementMetadata(xnts)$score=as.numeric(elementMetadata(xnts)$name)
  strand(xnts)=ifelse(elementMetadata(xnts)$score >=0, "+", "-")
  elementMetadata(xnts)$score=abs(elementMetadata(xnts)$score)
  elementMetadata(xnts)$name=NULL
  return(xnts)
}


randomRepositioningiCLIP <- function(gr, gtf, rep.num = 10){
  set.seed(1)
  gr.seqlevels <- unique(
    as.character(seqnames(gr))
    )
      
  gtf <- gtf[seqnames(gtf) %in% gr.seqlevels]
  seqlevels(gtf) <- gr.seqlevels
  
  seqlevels(gr) <- gr.seqlevels
  
  ol <- as.matrix(findOverlaps(gr, gtf))
  
  ol <- ol[!duplicated(ol[, 1]), ]
      
  r.gr <- randomRepositioningDf(gr, gtf, ol, rep.num)
  
  return(r.gr)
}
    
randomRepositioningDf <- function(gr, gtf, ol, rep.num = 20){
  require(multicore)
	
  temp.df <- data.frame(matrix(0, nrow = length(gtf[ol[, 2]]), ncol = rep.num + 2))
  colnames(temp.df)[1:2] <- c("gene_id", "width")
  temp.df$gene_id <- as.character(seqnames(gtf[ol[, 2]]))
  temp.df$width <- width(gr[ol[, 1]])
  # The column for width is designed for hybrid which have ranges
  
  df <- temp.df
  
  df[, 3:(rep.num + 2)] <- do.call("rbind", mclapply(1:nrow(df), function(x) sample(start(gtf[ol[x, 2]]):(end(gtf[ol[x, 2]]) - width(gr[ol[x, 1]]) + 1), rep.num, replace = TRUE)))
  
  # TO DO, add loop for creating GRangesLists
  
  rL.gr <- list()
  
  for(i in 1:rep.num){

	  	  r.gr <- GRanges(seqnames = Rle(df$gene_id), 
	                  ranges = IRanges(df[, 2 + i],
	                    width = df$width), 
	                  strand = Rle(rep("+", nrow(df))),
	                  rname = 1:nrow(df),
	                  score = as.integer(1),
	                  category = "unknown",
	                  annot = "unknown",
	                  sequence = "unknown"
	                  )
	rL.gr <- c(rL.gr, r.gr)
  }
  
  # CAUTION, rL.gr is just list
  return(rL.gr)
}


compareHybridvsiCLIP <- function(iCLIP, duplex, gtf, rep.num = 10, y.axis.lim = 130){
  iCLIP.ol.hybrid <- iCLIP[seqnames(iCLIP) %in% as.character(seqnames(duplex$L))]
  seqlevels(iCLIP.ol.hybrid) <- unique(
    as.character(seqnames(iCLIP.ol.hybrid))
    )
  
  r.iCLIP.ol.hybrid <- randomRepositioningiCLIP(iCLIP.ol.hybrid, gtf, rep.num)
  
  plot.out.list <- list()
  
  for(i in 1:4){
    sam.flag <- i
    plot.out.list[[i]] <- compareSites(iCLIP, r.iCLIP.ol.hybrid, duplex, sam.flag, y.axis.lim)
  }

  print(multiplot(plot.out.list[[1]], plot.out.list[[2]], plot.out.list[[3]], plot.out.list[[4]], cols = 2))
}



## Following functions are used inside the compareHybridvsiCLIP() function.
compareSitesDf <- function(iCLIP, CLIP, searchRange = 100, XS.flag = TRUE){
  
  if(XS.flag){
    elementMetadata(iCLIP)[["score"]] <- 1
  }
  
  rnaMap <- CLIP
  
  start(rnaMap) <- start(CLIP) - searchRange
  end(rnaMap) <- end(CLIP) + searchRange
  
  overlaps=as.matrix(findOverlaps(rnaMap, iCLIP))
  flank.size <- searchRange
  
  overlaps=as.data.frame(overlaps)
  overlaps$iCLIP.tr=start(iCLIP)[overlaps[,2]]
  overlaps$strand=as.vector(strand(iCLIP))[overlaps[,2]]
  overlaps$score=elementMetadata(iCLIP)$score[overlaps[,2]]
  overlaps$CLIP.del=start(rnaMap)[overlaps[,1]]+flank.size

  ## get distances to exon and separate according to strand
  distances=mapply(function(x, y, z){ rep(x-z, y) }, overlaps$iCLIP.tr, overlaps$score, overlaps$CLIP.del)
  distances=split(distances, overlaps$strand)

  ## split by strands
  distances=c(unlist(distances[["+"]]), (-1)*unlist(distances[["-"]]))
  ## keep only those that are within the range (for clusters that might extent over the region end)
  leftMost <- searchRange
  rightMost <- searchRange
  distances=distances[distances %in% c(-leftMost:(rightMost))]

                                        # determine frequencies
  distances=as.data.frame( table(distances) )
  names(distances)=c("rel.pos","xnts")

  ## get plotting positions (ensure that all positions are present)
  pos=c(-leftMost:(rightMost-1))
  distances=data.frame(rel.pos=pos, 
    xnts=distances$xnts[match(pos,distances$rel.pos)])
  distances$xnts[is.na(distances$xnts)]=0

  ## get relative counts
  distances$rel.count=distances$xnts/length(rnaMap)
  
  return(distances)
}

compareSites <- function(iCLIP, random.iCLIP, hybrid.grL, sam.flag, y.axis.lim){
  
  require(ggplot2)
  
  print(sam.flag)
  
  if(sam.flag == 1){
    sam.label <- "L, start"
    CLIP <- hybrid.grL$L
    end(CLIP) <- start(CLIP)	
  } else if(sam.flag == 2){
    sam.label <- "L, end"
    CLIP <- hybrid.grL$L
    start(CLIP) <- end(CLIP)	
  } else if(sam.flag == 3){
    sam.label <- "R, start"
    CLIP <- hybrid.grL$R
    end(CLIP) <- start(CLIP)	
  } else if(sam.flag == 4){
    sam.label <- "R, end"
    CLIP <- hybrid.grL$R
    start(CLIP) <- end(CLIP)	
  } else {
    stop("Define sam.flag")
  }
  
  iCLIP.df <- compareSitesDf(iCLIP, CLIP, searchRange = 100)
  colnames(iCLIP.df)[2] <- "iCLIP.xnts"
  
  for(i in 1:length(random.iCLIP)){
    random.iCLIP.df <- compareSitesDf(random.iCLIP[[i]], CLIP, searchRange = 100)
    colnames(random.iCLIP.df)[2] <- paste("random.iCLIP.xnts", i, sep = "")
    iCLIP.df <- cbind(iCLIP.df, random.iCLIP.df[, 2])
  }
  
  iCLIP.df$random.mean <- rowMeans(iCLIP.df[, 4:ncol(iCLIP.df)])
  iCLIP.df$random.sd <- apply(iCLIP.df[, 4:ncol(iCLIP.df)], 1, sd) 
  
  
  if((sam.flag %% 2) == 1){
    plot.out <- ggplot(iCLIP.df[1:111, ]) + 
      geom_area(aes(x = rel.pos, y = iCLIP.xnts), fill = "blue", alpha = 0.3) + 
        geom_line(aes(x = rel.pos, y = iCLIP.xnts), col = "blue") + 
          geom_line(aes(x = rel.pos, y = random.mean)) + 
            geom_ribbon(aes(x = rel.pos, ymin = random.mean - random.sd, ymax = random.mean + random.sd), alpha = 0.3) + 
              coord_cartesian(ylim=c(0, y.axis.lim))
  } else if ((sam.flag %% 2) == 0){
    plot.out <- ggplot(iCLIP.df[91:nrow(iCLIP.df), ]) + 
      geom_area(aes(x = rel.pos, y = iCLIP.xnts), fill = "blue", alpha = 0.3) + 
        geom_line(aes(x = rel.pos, y = iCLIP.xnts), col = "blue") + 
  	  geom_line(aes(x = rel.pos, y = random.mean)) + 
            geom_ribbon(aes(x = rel.pos, ymin = random.mean - random.sd, ymax = random.mean + random.sd), alpha = 0.3) + 
              coord_cartesian(ylim=c(0, y.axis.lim))
  }
  
  return(plot.out)
}
