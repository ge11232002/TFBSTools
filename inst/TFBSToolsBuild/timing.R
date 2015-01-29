library(TFBSTools)
library(JASPAR2014)
library(CNEr)

pfm <- getMatrixByID(JASPAR2014, ID="MA0139")
pwm <- toPWM(pfm)

### Timing for Axt alignments
axtFns <- list.files("/export/downloads/ucsc/axtNet/hg19",
                     pattern=".*hg19.mm10.net.axt.gz",
                     full.names=TRUE)
axt <- readAxt(axtFns[1])
axt <- axt[1:2000] ## 1124711bp

cputime <- proc.time()
sitePairSet <-  searchAln(pwm, axt, min.score="80%",
                          windowSize=51L, cutoff=0.7, strand="+",
                          type="any", conservation=NULL, mc.cores=8)
cputime <- proc.time() - cputime
## > cputime
##   user  system elapsed
## 138.703   1.692  24.353
## 30Mb per CPU hour
GRangesTFBS <- toGRangesList(sitePairSet, axt)


### Timing for BSGenome
library(rtracklayer)
library(JASPAR2014)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)
chain <- import.chain("/mnt/biggley/home/gtan/work2/projects/hg19ToMm10.over.chain")
cputime <- proc.time()
sitePairSet <- searchPairBSgenome(pwm, BSgenome.Hsapiens.UCSC.hg19,
                                 BSgenome.Mmusculus.UCSC.mm10,
                                 chr1="chr1", chr2="chr1",
                                 min.score="80%", strand="+", chain=chain)
cputime <- proc.time() - cputime

