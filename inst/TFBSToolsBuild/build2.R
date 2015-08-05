
#R CMD build --no-build-vignettes --no-manual TFBSTools
#R_dev CMD build TFBSTools
#R_dev CMD check TFBSTools_1.7.1.tar.gz
#R_dev CMD INSTALL TFBSTools_1.7.1.tar.gz

###########################debug ##########################
pfm = PFMatrix(ID="M0001", name="MyProfile", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), matrix=matrix(as.integer(c(12, 3, 0, 0, 4, 0, 0, 0, 0, 11, 7, 0, 0, 9, 12, 0, 0, 0, 0, 0, 0, 1, 1, 12)), byrow=TRUE, nrow=4, dimnames=list(c("A", "C", "G", "T"))))

## to ICM, PWM
pwm = toPWM(pfm, pseudocounts=sqrt(12))
icm = toICM(pfm)
plotLogo(icm)

## searchSeq
searchSeq(pwm, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", min.score="60%")

## searchAln
s1 <- 
    DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
s2 <-
    DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
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

## Site
site = searchSeq(pwm, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", min.score="60%")

## SiteList
siteList = SiteList(site, site)
searchSeq(pwmList, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", min.score="60%")

## Alignments
foo = PairwiseAlignments(as.character(pattern(globalAlign)), as.character(subject(globalAlign)))

foo = PairwiseAlignmentTFBS(pattern="ACTTCACCAGCTCCCTGGCGGTAAGTTGATC---AAAGG---AAACGCAAAGTTTTCAAG", subject="GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC", windowSize=13L)
doSiteSearch(pwm, aln1=foo, min.score="50%")
searchAln(pwm, aln1=foo,min.score="50%", windowSize=13L, cutoff=0.7)
searchAln(pwmList, aln1=foo,min.score="50%", windowSize=13L, cutoff=0.7)


## DB
db = "/Users/gtan/CSC/JASPAR/JASPARData/inst/extdata/JASPAR_2010.sqlite"
pfm = get_Matrix_by_ID("/Users/gtan/CSC/JASPAR/JASPARData/inst/extdata/JASPAR_2010.sqlite", ID="MA0003", type="PFM")
pwm = toPWM(pfm, pseudocounts=sqrt(185))
icm = get_Matrix_by_ID(db, ID="MA0003", type="ICM")

pfm2 = get_Matrix_by_name(db, name="TFAP2A", type="PFM")
icm2 = get_Matrix_by_name(db, name="TFAP2A", type="ICM")

## DB .get_IDlist_by_query
opts = list()
opts[["species"]] = 9606
opts[["name"]] = "RUNX1"
#opts[["class"]] = "Ig-fold"
opts[["type"]] = "SELEX"
opts[["all_versions"]] = TRUE
foo = get_MatrixSet(con, opts)


## MEME wrappers
/usr/local/Cellar/meme/4.9.0-p4/bin/meme crp0.s -text -dna -mod oops -pal 2>/dev/null
foo=runMEME("/Users/gtan/src/meme_4.9.1/tests/crp0.s", binary="/usr/local/Cellar/meme/4.9.0-p4/bin/meme")


