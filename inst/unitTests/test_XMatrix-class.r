
test_Accessor <- function(){
  pfm <- PFMatrix(ID="MA0004.1", name="Arnt", matrixClass="Zipper-Type", 
                  strand="+",
                  bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                  tags=list(family="Helix-Loop-Helix", 
                            species="10090", tax_group="vertebrates",
                            medline="7592839", type="SELEX", ACC="P53762", 
                            pazar_tf_id="TF0000003",
                            TFBSshape_ID="11", TFencyclopedia_ID="580"),
                  matrix=matrix(c(4L,  19L, 0L,  0L,  0L,  0L,
                                  16L, 0L,  20L, 0L,  0L,  0L,
                                  0L,  1L,  0L,  20L, 0L,  20L,
                                  0L,  0L,  0L,  0L,  20L, 0L),
                                byrow=TRUE, nrow=4, 
                                dimnames=list(c("A", "C", "G", "T")))
                  )
  checkIdentical(ID(pfm), "MA0004.1")
  checkIdentical(name(pfm), "Arnt")
  checkIdentical(matrixClass(pfm), "Zipper-Type")

}

