
#R CMD build --no-build-vignettes --no-manual TFBSTools
#R CMD build TFBSTools
#R CMD check TFBSTools_1.7.7.tar.gz
#R CMD BiocCheck TFBSTools_1.7.7.tar.gz
#R CMD INSTALL TFBSTools_1.7.7.tar.gz


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
# /usr/local/Cellar/meme/4.9.0-p4/bin/meme crp0.s -text -dna -mod oops -pal 2>/dev/null
foo=runMEME("/Users/gtan/src/meme_4.9.1/tests/crp0.s", 
            binary="/usr/local/Cellar/meme/4.9.0-p4/bin/meme")


