% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/voom.qn.limma.createRmd.R
\name{voom.qn.limma.createRmd}
\alias{voom.qn.limma.createRmd}
\title{Generate a \code{.Rmd} file containing code to perform differential expression analysis with the quantile-normalized voom+limma}
\usage{
voom.qn.limma.createRmd(data.path, result.path, codefile,
  norm.method = "quantile")
}
\arguments{
\item{data.path}{The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.}

\item{result.path}{The path to the file where the result object will be saved.}

\item{codefile}{The path to the file where the code will be written.}

\item{norm.method}{The microarray-style normalization method to be applied to the logCPM values. Choices are as for the method argument of voom function.}
}
\value{
The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
}
\description{
A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) by applying the voom transformation with quantile normalization (from the limma package) followed by differential expression analysis with limma.
}
\details{
For more information about the methods and the interpretation of the parameters, see the \code{limma} package and the corresponding publications.
}
\examples{
tmpdir <- normalizePath(tempdir(), winslash = "/")
SyntheticDataSimulation(simul.data= 'KIRC', dataset = file.path(tmpdir, "mydata.rds"), n.var = 500,
                                    samples.per.cond = 3, n.diffexp = 100, dispType='same', mode='D',
                                    fraction.upregulated=0.5, dataset.parameters=generateDatasetParameter())
compcodeR::runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "voom.qn.limma",
           Rmdfunction = "voom.qn.limma.createRmd",
           output.directory = tmpdir,)
}
\references{
Smyth GK (2005): Limma: linear models for microarray data. In: 'Bioinformatics and Computational Biology Solutions using R and Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds), Springer, New York, pages 397-420

 Law CW, Chen Y, Shi W and Smyth GK (2014): voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology 15, R29
}
\author{
BuKyung Baik
}
