% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlationheatmap.R
\name{correlation_heatmap}
\alias{correlation_heatmap}
\title{make correlation heatmap function}
\usage{
correlation_heatmap(working.dir, figure.dir, simul.data, nsample, topgenes,
  AnalysisMethods, rep)
}
\arguments{
\item{working.dir}{Input file location.}

\item{figure.dir}{Figure save location.}

\item{simul.data}{Type of dataset(e.g. KIRC and Bottomly)}

\item{nsample}{Number of samples. Input as a number.}

\item{topgenes}{Number of top significant genes obtained from each methods to estimate rank correlation.}

\item{AnalysisMethods}{DEmethods used for figures. Input as a character vector. (e.g. 'edgeR','edgeR.ql','edgeR.rb','DESeq.pd','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')}

\item{rep}{Number of replication each test contain.}
}
\description{
make correlation heatmap function
}
