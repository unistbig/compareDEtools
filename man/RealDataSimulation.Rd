% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimulationDataGeneration.R
\name{RealDataSimulation}
\alias{RealDataSimulation}
\title{Generate real data for analysis}
\usage{
RealDataSimulation(simul.data, dataset, samples.per.cond, fpc = FALSE,
  dataset.parameters)
}
\arguments{
\item{simul.data}{Characters indicating which dataset will be used for simulation data generation
"KIRC" for KIRC dataset.
"Bottomly" for Bottomly dataset.
"SEQC" for SEQC ERCC Spike-In dataset.}

\item{dataset}{Characters specifying the file name of simulation dataset.}

\item{samples.per.cond}{An integer indicating number of samples for each condition (e.g. 3).}

\item{fpc}{A logical indicating whether simulation data is made from single sample condition(e.g. normal) to calculate False positive counts. Default is FALSE.}

\item{dataset.parameters}{A list containing estimated mean and dispersion parameters and filtered count from original count dataset.}
}
\description{
Generate real data for analysis
}
