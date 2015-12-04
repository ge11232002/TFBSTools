### -----------------------------------------------------------------
### make FlatFileDir from JASPAR database: including *.pfm and matrix_list.txt
### The matrix_list.txt has columns: ID, InformationContent, Name, Class, Tags
### Exported!
makeFlatFileDir <- function(JASPAR){
  outputDir <- "FlatFileDir"
  dir.create(outputDir)
  
  # Make the matrix_list.txt
  pfms <- getMatrixSet(JASPAR, opts=list(all=TRUE))
  pfmIDs <- ID(pfms)
  pfmICs <- sapply(pfms, function(x){sum(totalIC(toICM(x)))})
  pfmNames <- name(pfms)
  ## There are cases that have more than two classes
  pfmClasses <- sapply(matrixClass(pfms), paste, collapse="/")
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
  write.table(ans, file=file.path(outputDir, "matrix_list.txt"),
              quote=FALSE, sep="\t",
              col.names = FALSE, row.names = FALSE)
  
  # Make the *.pfm
  lapply(pfms, function(x){write.table(x@profileMatrix,
                                       file=file.path(outputDir,
                                                      paste0(x@ID, ".pfm")),
                                       quote=FALSE, sep="\t",
                                       col.names = FALSE, row.names = FALSE)})
  return("success")
}

