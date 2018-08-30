# TFBSTools
Software Package for Transcription Factor Binding Site (TFBS) Analysis

## Installation of the stable version of `TFBSTools` from Bioconductor

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("TFBSTools")
```

## Installation of the development version of `TFBSTools` from github
**Prerequsite**:

  * Mac: Install "Command Line Tools" via `gcc` on terminal
  * Linux: Install a compiler and various development libraries 
    (details vary across different flavors of Linux).
  * Windows: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

```R
devtools::install_github("ge11232002/TFBSTools")
```
