% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\docType{class}
\name{Expression-class}
\alias{Expression-class}
\title{An S4 class to represent gene expression data for multiple samples
for a given species. Assumes that all expession data are derived
from mapping to the same reference. Applies to data for g genes
across s samples (ie, sequenced libraries).}
\description{
The fields that apply to the s samples correspond to those
}
\section{Slots}{

\describe{
\item{\code{species}}{The species}

\item{\code{edgeR}}{An edgeR DGEList object holding data for all samples. Expression
matrix of dimension g,s.}

\item{\code{lengths}}{Summary of length, in bp, of transctipts for each gene. Length g.}

\item{\code{individual}}{Factors indicating which individual each sample is from. Length s.}

\item{\code{treatment}}{Factors indicating which treatment applies to each sample. Length s.}

\item{\code{id}}{Factors indicating the unique sequencing run id of each sample, eg HWI-ST625-75-D0PBDACXX-6-ATCACG. Length s.}

\item{\code{library_id}}{Factors indicating the unique library id of each sample, eg FEG365. Length s.}

\item{\code{sample_prep}}{Sample prep strategy. Length s.}

\item{\code{genome_type}}{Character indicating genome type, eg nuclear. See agalma documentation. Length g.}

\item{\code{molecule_type}}{Character indicating encoded molecule type, eg protein. See agalma documentation. Length g.}

\item{\code{blast_hit}}{Blast hit. Length g.}

\item{\code{rRNA}}{Fraction of reads in sample that are rRNA. Length s.}

\item{\code{protein}}{Fraction of reads in sample that are protein coding. Length s.}

\item{\code{x}}{The counts matrix}
}}

