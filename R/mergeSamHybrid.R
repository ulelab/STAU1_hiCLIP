args<-commandArgs(TRUE)

P1.file <- args[1]
P3.file <- args[2]
P4.file <- args[3]
out.file <- args[4]

# Except P4, all unmapped reads should be removed

P1 <- read.table(P1.file, sep = "\t", stringsAsFactors = FALSE, comment.char = "@", fill = TRUE, col.names=paste("V", 1:14, sep=""))
P3 <- read.table(P3.file, sep = "\t", stringsAsFactors = FALSE, comment.char = "@", fill = TRUE, col.names=paste("V", 1:14, sep=""))
P4 <- read.table(P4.file, sep = "\t", stringsAsFactors = FALSE, comment.char = "@", fill = TRUE, col.names=paste("V", 1:14, sep=""))

P1 <- P1[P1$V2 != 4, ]
P3 <- P3[P3$V2 != 4, ]

P1 <- P1[, 1:11]
P3 <- P3[, 1:11]
P4 <- P4[, 1:11]

print(nrow(P1))
print(nrow(P3))
print(nrow(P4))

P.merged <- rbind(P1, P3, P4)

P.merged <- P.merged[order(as.integer(sapply(strsplit(P.merged$V1, "_"), "[[", 1)), 
	as.character(sapply(strsplit(P.merged$V1, "_"), "[[", 2))
	), ]
	
	
zz <- file(out.file, "w")
cat("@HD header was removed\n", file=zz)
write.table(P.merged, file=zz, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
close(zz)