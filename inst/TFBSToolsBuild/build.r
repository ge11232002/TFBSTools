############ create package ##################
rm *.pdf
R CMD Rd2pdf .

remove.packages("TFBSTools")
rm TFBSTools_0.99.3.tar.gz
R CMD check TFBSTools_1.1.3.tar.gz
R CMD build TFBSTools
R CMD INSTALL TFBSTools_1.1.3.tar.gz

library(TFBSTools)

library(JASPAR2014)

pfm = getMatrixByID(JASPAR2014, ID="MA0003")

###########################debug ##########################
pfm = PFMatrix(ID="MA0004.1", name="Arnt", matrixClass="Zipper-Type", strand="+",
				bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
				tags=list(family="Helix-Loop-Helix", species="10090", tax_group="vertebrates",
				medline="7592839", type="SELEX", ACC="P53762", pazar_tf_id="TF0000003",
				TFBSshape_ID="11", TFencyclopedia_ID="580"),
				matrix=matrix(c(4L,  19L, 0L,  0L,  0L,  0L,
								16L, 0L,  20L, 0L,  0L,  0L,
								0L,  1L,  0L,  20L, 0L,  20L,
								0L,  0L,  0L,  0L,  20L, 0L), 
								byrow=TRUE, nrow=4, dimnames=list(c("A", "C", "G", "T")))
				)

pfmList = PFMatrixList(pfm1=pfm, pfm2=pfm)

## to ICM, PWM
pwm = toPWM(pfm, pseudocounts=sqrt(ncol(Matrix(pfm))))
icm = toICM(pfm)
seqLogo(icm)

## searchAln
s1 <- DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
s2 <-    DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
globalAlign <-
    pairwiseAlignment(s1, s2, substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)
x = DNAStringSet(c(as.character(pattern(globalAlign)), as.character(subject(globalAlign))))
min.score = "50%"
searchAln(pwm, x, min.score="50%", windowSize=11L, cutoff=0.5)

## MatrixList
pfmList = PFMatrixList(a=pfm, b=pfm)
icmList = ICMatrixList(a=icm, b=icm)
pwmList = PWMatrixList(a=pwm, b=pwm)
searchSeq(pwmList, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", min.score="60%")
searchAln(pwmList, x, min.score="50%", windowSize=11L, cutoff=0.5)

## searchSeq
site = searchSeq(pwm, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", min.score="60%")

## SiteList
siteList = SiteList(site, site)
siteList = searchSeq(pwmList, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", min.score="60%")

## Alignments
pwm = toPWM(pfm)
aln = PairwiseAlignmentTFBS(pattern="ACTTCACCAGCTCCCTGGCGGTAAGTTGATC---AAAGG---AAACGCAAAGTTTTCAAG", subject="GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC", windowSize=13L)

sitePair = searchAln(pwm, aln1=aln, min.score="50%", windowSize=13L, cutoff=0.4)
sitePairList = searchAln(pwmList, aln1=aln, min.score="50%", windowSize=13L, cutoff=0.4)


## DB
pfm = getMatrixByID(JASPAR2014, ID="MA0003")
pwm = toPWM(pfm)
icm = toICM(pfm)
pfm2 = getMatrixByName(JASPAR2014, name="TFAP2A")

opts = list()
opts[["species"]] = 9606
opts[["name"]] = "RUNX1"
#opts[["class"]] = "Ig-fold"
opts[["type"]] = "SELEX"
opts[["all_versions"]] = TRUE
pfmList = getMatrixSet(JASPAR2014, opts)


## MEME wrappers
/usr/local/Cellar/meme/4.9.0-p4/bin/meme /Users/gtan/src/meme_4.9.1/tests/crp0.s -text -dna -mod oops -pal 2>/dev/null

foo = runMEME("/Users/gtan/src/meme_4.9.1/tests/crp0.s", binary="/usr/local/Cellar/meme/4.9.0-p4/bin/meme", arguments="-nmotifs 3")

## test on long chromosome
library(BSgenome.Hsapiens.UCSC.hg19)
library(TFBSTools)

library(JASPAR2014)
chr2 = unmasked(BSgenome.Hsapiens.UCSC.hg19[[2]])
#db = "/Users/gtan/CSC/JASPAR/JASPAR2014/inst/extdata/JASPAR2014.sqlite"
pfm = getMatrixByName(JASPAR2014, name="CTCF")
pwm = toPWM(pfm)
icm = toICM(pfm)

system.time((ans = searchSeq(pwm, chr2, min.score="90%")))
  user  system elapsed
 21.100   0.006  21.110

## searchMatrix
dyn.load("matrixAlignerDynamic.so")
.Call("matrixAligner", matrixQuery, matrixSubject, 3, 0.01)

library(TFBSTools)
library(JASPAR2014)
pfmSubject = getMatrixByID(JASPAR2014, ID="MA0055")
pfmQuery = getMatrixByID(JASPAR2014, ID="MA0048")
searchMatrix(pfmSubject, pfmQuery)	

