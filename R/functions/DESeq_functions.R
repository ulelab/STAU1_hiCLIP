normaliseDf <- function(df, min.val = 30){
	QCdf <- df
	QCdf <- QCdf[apply(QCdf, 1, min) > min.val, ]
	df <- df[rownames(df) %in% rownames(QCdf), ]
	df
}

RunDESeq <- function(file.name, out.pre){
	require("DESeq")
	
	df <- read.table(file.name, row.names = 1, header = 1, sep = "\t")
	df <- normaliseDf(df)
	conds <- factor(c("STAU1", "STAU1", "KD", "KD", "STAU1", "STAU1"), levels = c("STAU1", "KD"))
	
	## 1 Input data and preparation
	cds <- newCountDataSet(df, conds)
	cds <- estimateSizeFactors(cds)
	write.table(sizeFactors(cds), paste(out.pre, "_size_factor.tab", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
	
	write.table(counts(cds, normalized=TRUE), paste(out.pre, "_DESeqNormalisedCounts.tab", sep = ""), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}

RunDESeqNonHybrid <- function(file.name, out.pre){
	require("DESeq")
	
	df <- read.table(file.name, row.names = 1, header = 1, sep = "\t")
	# df <- normaliseDf(df)
	conds <- factor(c("STAU1", "STAU1", "STAU1"))
	
	## 1 Input data and preparation
	cds <- newCountDataSet(df, conds)
	cds <- estimateSizeFactors(cds)
	write.table(sizeFactors(cds), paste(out.pre, "_size_factor.tab", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
	
	write.table(counts(cds, normalized=TRUE), paste(out.pre, "_DESeqNormalisedCounts.tab", sep = ""), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}

