# agalmar: A Collection of Tools for Analyzing Agalma Expression Projects

[![Build Status](https://travis-ci.org/caseywdunn/agalmar.svg?branch=master)](https://travis-ci.org/caseywdunn/agalmar)

This package contains tools for analyzing phylogenetic and expression projects 
from [Agalma](https://bitbucket.org/caseywdunn/agalma). The interface between 
Agalma and agalmar is the json file 
[exported](https://bitbucket.org/caseywdunn/agalma/src/master/scripts/agalma-export-expression) 
from Agalma.


## Documentation

See the manual, 
[agalmar-manual.pdf](https://github.com/caseywdunn/agalmar/raw/master/agalmar-manual.pdf)

## Example usage

This package comes with a test dataset called `janedoe`. To process 
it, you can run: 

    e <- lapply(janedoe$expression, Expression)

## Citing

To find out how to cite hutan, run the following in R:

    citation("agalmar")

## Installing

### From the git repository

First, install the [devtools](https://github.com/hadley/devtools) package. Then, run the following in R:

    source("https://bioconductor.org/biocLite.R")
    biocLite("edgeR")
    biocLite("DESeq2")
    library(devtools)
    install_github('caseywdunn/hutan')
    install_github('caseywdunn/agalmar')


## Development

This package is built with the excellent devtools 
(https://github.com/hadley/devtools). Extensive explanations on using devtools 
to create R packages is available in the book 
[R Packages](http://r-pkgs.had.co.nz/).

Development typically involves `cd`ing to the package directory, launching R, 
and running some combination of the following: 
	
	options(error=traceback) # Get line numbers for errors
    library(devtools)
    load_all()
    test()
    document()
    check() # A wrapper for R CMD check, see http://r-pkgs.had.co.nz/check.html#check
    build() # Create package bundle, including executed vignettes

To regenerate the pdf manual, run the following shell command in the package directory:

    R CMD Rd2pdf . --force --output=agalmar-manual.pdf
    
