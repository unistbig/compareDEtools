% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edgeR.GLM.Robust.createRmd.R
\name{edgeR.GLM.Robust.createRmd}
\alias{edgeR.GLM.Robust.createRmd}
\title{Generate a \code{.Rmd} file containing code to perform differential expression analysis with the edgeR GLM Robust dispersion estimation}
\usage{
edgeR.GLM.Robust.createRmd(data.path, result.path, codefile,
  norm.method = "TMM", disp.estimation = "CoxReid")
}
\arguments{
\item{data.path}{The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.}

\item{result.path}{The path to the file where the result object will be saved.}

\item{codefile}{The path to the file where the code will be written.}

\item{norm.method}{The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}.}

\item{disp.estimation}{method for estimating the dispersion. Possible values are "CoxReid", "Pearson" or "deviance".}
}
\value{
The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
}
\description{
A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the edgeR package.
}
\details{
For more information about the methods and the interpretation of the parameters, see the \code{edgeR} package and the corresponding publications.
}
\examples{
tmpdir <- normalizePath(tempdir(), winslash = "/")
SyntheticDataSimulation(simul.data= 'KIRC', dataset = file.path(tmpdir, "mydata.rds"), n.var = 500,
                                    samples.per.cond = 3, n.diffexp = 100, dispType='same', mode='D',
                                    fraction.upregulated=0.5, dataset.parameters=generateDatasetParameter())
compcodeR::runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "edgeR.GLM.Robust",
           Rmdfunction = "edgeR.GLM.Robust.createRmd",
           output.directory = tmpdir)
}
\references{
Robinson MD, McCarthy DJ and Smyth GK (2010): edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140

Zhou X, Lindsay H, Robinson MD (2014). Robustly detecting differential expression in RNA sequencing data using observation weights. Nucleic Acids Research, 42(11), e91.
}
\author{
BuKyung Baik
}
