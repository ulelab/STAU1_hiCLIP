hiCLIP.dir <- "/netscr/yoichiro/STAU1_hiCLIP"
R.dir <- "~/R-2.15.1"
python.dir <- "~/Python-2.7.1/"
lib.path <- paste(R.dir, "/library/", sep = "")
.libPaths(lib.path)

library(GGally)
library(reshape2)
library(plyr)
library(GenomicRanges)
library(ShortRead)
library(ggplot2)

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}

S4dir <- paste(hiCLIP.dir, "/R/S4", sep = "")
sourceDir(S4dir)

functions.dir <- paste(hiCLIP.dir, "/R/functions", sep = "")
sourceDir(functions.dir)