library(TFBSTools)
library(JASPAR2014)
library(CNEr)

pfm <- getMatrixByID(JASPAR2014, ID="MA0139")
pwm <- toPWM(pfm)

### Timing for Axt alignments
axtFns <- list.files("/export/downloads/ucsc/axtNet/hg19",
                     pattern=".*hg19.mm10.net.axt.gz",
                     full.names=TRUE)
axt <- readAxt(axtFns)

cputime <- proc.time()
sitePairSet <-  searchAln(pwm, axt, min.score="80%",
                          windowSize=51L, cutoff=0.7, strand="+",
                          type="any", conservation=NULL, mc.cores=8)
cputime <- proc.time() - cputime
## > cputime
##   user     system    elapsed
## 181650.713    206.973  30450.453
## 25Mbp/CPUHour

GRangesTFBS <- toGRangesList(sitePairSet, axt)


### Timing for BSGenome
library(rtracklayer)
library(JASPAR2014)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)
chain <- import.chain("~/Downloads/hg19ToMm10.over.chain")
cputime <- proc.time()
sitePairSet <- searchPairBSgenome(pwm, BSgenome.Hsapiens.UCSC.hg19,
                                 BSgenome.Mmusculus.UCSC.mm10,
                                 chr1="chr1", chr2="chr1",
                                 min.score="80%", strand="+", chain=chain)
cputime <- proc.time() - cputime
## > cputime
##   user  system elapsed
##  26.696   0.706  27.514

### Timing for Sequence
cputime <- proc.time()
siteset <- searchSeq(pwm, BSgenome.Hsapiens.UCSC.hg19[["chr1"]],
                     min.score="80%", strand="+")
cputime <- proc.time() - cputime
## > cputime
##   user  system elapsed
##   12.818   0.009  12.828
## 74Gbp/CPU hour


