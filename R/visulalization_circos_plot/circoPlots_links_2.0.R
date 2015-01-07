main.dir <- "~/Documents/hi-ClIP/results/processed/"
setwd(main.dir)


load("all.collapsed.duplex.hgrL.Rdata")
hybrids <- all.collapsed.duplex.hgrL@unlistData
library(GenomicRanges)

## all mRNA
cds <- read.table("gr.tc.longest.mRNA_df_cds.tab")
utr5 <- read.table("gr.tc.longest.mRNA_df_5UTR.tab")
utr3 <- read.table("gr.tc.longest.mRNA_df_3UTR.tab")

## targeted mRNAs
targets <- unique(as.vector(seqnames(hybrids)))
cds <- cds[cds[,1]%in%targets,]
utr5 <- utr5[utr5[,1]%in%targets,]
utr3 <- utr3[utr3[,1]%in%targets,]

## make a GRanges object
mRNAs <- rbind(utr5,cds,utr3)
mRNAs.GR <- GRanges(seqnames=Rle(mRNAs$seqnames),ranges=IRanges(start=mRNAs$start,end=mRNAs$end))
values(mRNAs.GR) <- data.frame(annot=mRNAs$annot)


## find overlap 
ovl <- findOverlaps(hybrids,mRNAs.GR,select="first")
hybrids <- hybrids[!is.na(ovl)]
ovl <- ovl[!is.na(ovl)]
values(hybrids) <- cbind(rnames=as.data.frame(values(hybrids))$rname,annot=as.vector(as.data.frame(values(mRNAs.GR))$annot[ovl]),loop=(start(hybrids)[((length(hybrids)/2)+1):length(hybrids)]-start(hybrids[1:(length(hybrids)/2)])))
mRNAs.GR <- mRNAs.GR[ovl]



## transform coordinates (normalized)
end.feature <- end(mRNAs.GR) - start(mRNAs.GR)
center.arm <- as.integer((end(hybrids)+start(hybrids))/2)
positions.arms <- center.arm - start(mRNAs.GR)
relative.pos <- as.integer((positions.arms*1000)/end.feature)
relative.GR <- GRanges(seqnames=seqnames(hybrids),ranges=IRanges(start=relative.pos,end=relative.pos+1))
values(relative.GR) <- cbind(as.data.frame(values(hybrids)),end.feature = end.feature)
annot <- as.vector(as.data.frame(values(relative.GR))$annot)
end(relative.GR)[annot=="utr3"] <- end(relative.GR)[annot=="utr3"]*6
start(relative.GR)[annot=="utr3"] <- start(relative.GR)[annot=="utr3"]*6
end(relative.GR)[annot=="CDS"] <- end(relative.GR)[annot=="CDS"]*2
start(relative.GR)[annot=="CDS"] <- start(relative.GR)[annot=="CDS"]*2


## transform coordinates (extend)

end.feature <- end(mRNAs.GR) - start(mRNAs.GR)
center.arm <- as.integer((end(hybrids)+start(hybrids))/2)
positions.arms <- center.arm - start(mRNAs.GR)
relative.GR <- GRanges(seqnames=seqnames(hybrids),ranges=IRanges(start=positions.arms,end=positions.arms+1))
values(relative.GR) <- cbind(as.data.frame(values(hybrids)),end.feature = end.feature)

utrs.GR <- relative.GR[as.vector(as.data.frame(values(relative.GR))$annot)!="CDS"]
cds.GR <- relative.GR[as.vector(as.data.frame(values(relative.GR))$annot)=="CDS"]
features <- as.vector(as.data.frame(values(relative.GR))$annot)
max.CDS <- max(end.feature[as.vector(as.data.frame(values(mRNAs.GR))$annot)=="CDS"])
end.cds<- end(cds.GR) + (max.CDS - end.feature[features=="CDS"])
start.cds <- start(cds.GR) + (max.CDS - end.feature[features=="CDS"])

ranges(cds.GR) <- IRanges(start=start.cds,end=end.cds)

relative.GR <- c(cds.GR,utrs.GR)


## build links
duplex.u <- unique(as.vector(as.data.frame(values(relative.GR))$rnames))
duplex <- as.vector(as.data.frame(values(relative.GR))$rnames)
colors <- rep(NA,length(relative.GR))
l <- rep(NA,length(relative.GR))
gloub <- 0
for (i in duplex.u){
	temp <- as.vector(as.data.frame(values(relative.GR))$annot)[duplex==i]
	coords <- start(relative.GR)[duplex==i]
	temp <- temp[order(temp)]
	if (temp[2] == "utr5" & temp[1] != "utr5"){
		gloub <- gloub + 1
	}	
	temp <- paste(temp[1],temp[2],sep=".")
	loop <- as.integer(as.vector(as.data.frame(values(relative.GR))$loop)[duplex==i][1])
	#print(relative.GR[duplex==i])
	#print(temp)
	if (temp %in% c("utr5.utr5","CDS.CDS","utr3.utr3")){
		if (loop <= 200){
			col <- "intra10"
		}
		else if (loop > 200 & loop <= 500){
			col <- "intra20"
		}
		else if (loop > 500 & loop <= 1000){
			col <- "intra50"
		}
		else if (loop > 1000 & loop <= 2000){
			col <- "intra100"
		}
		else if (loop > 2000 & loop <= 2500){
			col <- "intra200"
		}
		else if (loop > 2500 & loop <= 3000){
			col <- "intra500"
		}
		else if (loop > 3000 & loop <= 4000){
			col <- "intra1000"
		}
		else if (loop > 4000 & loop <= 5000){
			col <- "intra2000"
		}
		else if (loop > 5000 & loop <= 6000){
			col <- "intra4000"
		}
		else if (loop > 6000 & loop <= 8000){
			col <- "intra8000"
		}
		else if (loop > 8000){
			col <- "intra10000"
		}
	}
	else{
		print("inter")
		if (loop <= 100 ){
			col <- "inter100"
		}
		else if (loop > 100 & loop <= 200){
			col <- "inter200"
		}
		else if (loop > 200 & loop <= 300){
			col <- "inter300"
		}
		else if (loop > 300 & loop <= 500){
			col <- "inter500"
		}
		else if (loop > 500 & loop <= 1000){
			col <- "inter1000"
		}
		else if (loop > 1000 & loop <= 2000){
			col <- "inter2000"
		}
		else if (loop > 2000 & loop <= 3000){
			col <- "inter3000"
		}
		else if (loop > 3000 & loop <= 4000){
			col <- "inter4000"
		}
		else if (loop > 4000 & loop <= 5000){
			col <- "inter5000"
		}
		else if (loop > 5000 & loop <= 6000){
			col <- "inter6000"
		}
		else if (loop > 6000 & loop <= 8000){
			col <- "inter8000"
		}
		else if (loop > 8000){
			col <- "inter10000"
		}
		#print(relative.GR[duplex==i])
		print(hybrids[as.data.frame(values(hybrids))$rnames==i])
		print(loop)
		print(col)
		
	}
	
	colors[duplex==i] <- col
	l[duplex==i] <- loop
}


values(relative.GR) <- cbind(as.data.frame(values(relative.GR))[,-4],colors=colors)
relative.GR <- relative.GR[order(l)]


## write interactions
## /Users/darbo01/Documents/hi-ClIP/results/circoPlots/data/HiCLIP_all_duplexes_interactions.txt
vals <- as.data.frame(values(relative.GR))
int.names <- paste("hybrid",as.vector(vals$rnames),sep="")
feature <- as.vector(vals$annot)
colors <- paste("color=",as.vector(vals$colors),sep="")

to.write <- data.frame(int.names,feature,start(relative.GR),end(relative.GR),colors)

write.table(to.write,"/Users/darbo01/Documents/hi-ClIP/results/circoPlots/data/HiCLIP_all_duplexes_interactions_colors_normalized.txt",sep="\t")


## create heatmap with number of sequences
end.feature <- end(mRNAs.GR) - start(mRNAs.GR)
features <- as.vector(as.data.frame(mRNAs.GR)$annot)
max.CDS <- max(end.feature[features=="CDS"])
max.utr3 <- max(end.feature[features=="utr3"])
max.utr5 <- 4422

heatmap.CDS <- matrix(0,ncol=max.CDS,nrow=length(which(features=="CDS")))
heatmap.utr3 <- matrix(0,ncol=max.utr3,nrow=length(which(features=="utr3")))
heatmap.utr5 <- matrix(0,ncol=max.utr5,nrow=length(which(features=="utr5")))


feats <- c("utr5","CDS","utr3")
feats <- "utr5"
for (f in feats){
	mat <- get(paste("heatmap",f,sep="."))
	for (i in 1:nrow(mat)){
		print(i)
		l.feature <- end.feature[features==f][i]
		if (f == "CDS"){
			temp <- c(rep(0,(ncol(mat)-l.feature)),rep(1,l.feature))
			mat[i,] <- temp
			print(mat[i,ncol(mat)])
		}
		else{
			mat[i,1:l.feature] <- 1
		}
	}
	counts <- round(apply(mat,2,sum) / nrow(mat) * 100,digit=1)
	colors <- rep("vvvlr",length(counts))
	colors[counts>=12.5] <- "vvlr"
	colors[counts>=25] <- "vlr"
	colors[counts>=37.5] <- "lr"
	colors[counts>=50] <- "dr"
	colors[counts>=62.5] <- "vdr"
	colors[counts>=75] <- "vvdr"
	colors[counts>=87.5] <- "vvvdr"
	to.write <- data.frame(chr=rep(f,ncol(mat)),start=1:(ncol(mat)),end=2:(ncol(mat)+1),value=counts,colors=paste("color=",colors,sep=""))
	write.table(to.write,paste("/Users/darbo01/Documents/hi-ClIP/results/circoPlots/data/",f,"_HiCLIP_all_duplexes_heatmap.txt",sep=""))
}




