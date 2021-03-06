% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performanceplot.R
\name{performance_realdata_plot}
\alias{performance_realdata_plot}
\title{make realdata plot function}
\usage{
performance_realdata_plot(working.dir, figure.dir, simul.data, rep,
  nsample, rowType, AnalysisMethods)
}
\arguments{
\item{working.dir}{Input file location}

\item{figure.dir}{Figure save location}

\item{simul.data}{Type of dataset (e.g. KIRC ,Bottomly and SEQC)}

\item{rep}{Number of replication each test contain.}

\item{nsample}{Number of samples. Input as a numeric vector.}

\item{rowType}{Type of measures. Combination of DetectedDE and FPC. (e.g. c('DetectedDE','FP.count'))}

\item{AnalysisMethods}{DEmethods used for figures. Input as character vectors
(e.g. 'edgeR','edgeR.ql','edgeR.rb','DESeq.pd','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')}
}
\description{
make realdata plot function
}
