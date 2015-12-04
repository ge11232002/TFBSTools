### -----------------------------------------------------------------
### make matrix_list.txt file from JASPAR database.
### The matrix_list.txt has columns: ID, InformationContent, Name, Class, Tags
### Exported!
makeMatrixList <- function(JASPAR){
  pfms <- getMatrixSet(JASPAR, opts=list(all=TRUE))
  pfmIDs <- ID(pfms)
  pfmICs <- sapply(pfms, function(x){sum(totalIC(toICM(x)))})
  pfmNames <- name(pfms)
  ## There are cases that have more than two classes
  pfmClasses <- sapply(matrixClass(pfms), paste, collapse="/")
  #tags <- c("acc", "collection", "family", "medline", "species",
  #          "tax_group", "type")
  pfmTags <- lapply(pfms, function(x){sapply(x@tags, paste, 
                                              collapse="/")})
  pfmTags <- sapply(pfmTags, function(x){x[order(names(x))]})
  pfmTags <- sapply(pfmTags, function(x){paste(";", paste0(names(x), " \"", x, 
                                                "\" ", collapse="; "))})
  ans <- cbind(data.frame(IDs=pfmIDs,
                          ICs=pfmICs,
                          Names=pfmNames,
                          Classes=pfmClasses,
                          Tags=pfmTags
                          ))
  write.table(ans, file="matrix_list.txt", quote=FALSE, sep="\t",
              col.names = FALSE, row.names = FALSE)
}