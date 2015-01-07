# QC of Ribosome Profiling data
.libPaths("~/R-2.15.1/library/")
source("~/R_program/exec/transcriptome_mapping/functions/annotation_functions.R")
source("~/R_program/exec/transcriptome_mapping/functions/RNAmapFunctions_transcripts.R")

load("~/Mapping/biodata/Rdata/annotations/transcript.coordinate.gtf.Rdata")

RP.file <- "/netscr/yoichiro/RibosomeProfiling/bed/total/id1_wt_1_RP_total.bed"
mRNASeq.file <- "/netscr/yoichiro/mRNASeq/bed/total/id1_wt_1_mRNASeq_total.bed"

start.codon <- extractStart(t.gtf)
stop.codon <- extractStop(t.gtf)

drawRNAmaps(RP.file, mRNASeq.file, start.codon, stop.codon, t.gtf, "~/STAU1_project/QC", searchRange = 50)

sessionInfo()
