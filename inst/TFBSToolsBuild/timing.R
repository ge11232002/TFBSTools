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
axt <- axt[1:2000] ## 1124711bp

cputime <- proc.time()
sitePairSet <-  searchAln(pwm, axt, min.score="80%",
                          windowSize=51L, cutoff=0.7, strand="+",
                          type="any", conservation=NULL, mc.cores=4)
cputime <- proc.time() - cputime


