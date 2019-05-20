#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with baySeq
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the \code{baySeq} package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{baySeq} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. Possible values are \code{"quantile"}, \code{"total"} and \code{"edgeR"}.
#' @param equaldisp Logical parameter indicating whether or not equal dispersion should be assumed across all conditions.
#' @param sample.size The size of the sample used to estimate the priors (default 5000).
#' @param estimation The approach used to estimate the priors. Possible values are \code{"QL"} (default), \code{"ML"} and \code{"edgeR"}.
#' @param pET The method used to re-estimate the priors. Possible values are \code{"BIC"} (default), \code{"none"} and \code{"iteratively"}.
#' @export
#' @author Charlotte Soneson
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references
#' Hardcastle TJ (2012): baySeq: Empirical Bayesian analysis of patterns of differential expression in count data. R package
#'
#' Hardcastle TJ and Kelly KA (2010): baySeq: Empirical Bayesian methods for identifying differential expression in sequence count data. BMC Bioinformatics 11:422
#' @examples
#' try(
#' if (require(baySeq)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' SyntheticDataSimulation(simul.data= 'KIRC', dataset = file.path(tmpdir, "mydata.rds"), n.var = 500,
#'                                     samples.per.cond = 3, n.diffexp = 100, dispType='same', mode='D',
#'                                     fraction.upregulated=0.5, dataset.parameters=generateDatasetParameter())
#' ## Note! In the interest of speed, we set sample.size=10 in this example.
#' ## In a real analysis, much larger sample sizes are recommended (the default is 5000).
#' compcodeR::runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "baySeq",
#'            Rmdfunction = "baySeq.createRmd",
#'            output.directory = tmpdir, norm.method = "edgeR",
#'            equaldisp = TRUE, sample.size = 10)
#' })
baySeq.createRmd <- function(data.path, result.path, codefile,
                             norm.method, equaldisp, sample.size = 5000,
                             estimation = "QL", pET = "BIC") {
  codefile <- file(codefile, open = 'w')
  writeLines("### baySeq", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = TRUE, error = TRUE, warning = TRUE}",
               "require(baySeq)",
               "require(edgeR)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               "baySeq.cd <- new('countData', data = count.matrix(cdata), replicates = sample.annotations(cdata)$condition, groups = list(NDE = rep(1, length(sample.annotations(cdata)$condition)), DE = sample.annotations(cdata)$condition))",
               paste("libsizes(baySeq.cd) <- baySeq::getLibsizes(baySeq.cd, estimationType = '", norm.method, "')", sep = '')), codefile)
  writeLines(c(paste("baySeq.cd <- baySeq::getPriors.NB(baySeq.cd, samplesize =", sample.size, ", equalDispersions = ", equaldisp, ", estimation = '", estimation, "', cl = NULL)", sep = ''),
               paste("baySeq.cd <- baySeq::getLikelihoods(baySeq.cd, prs = c(0.5, 0.5), pET = '", pET, "', cl = NULL)", sep = ''),
               "baySeq.cd@annotation <- data.frame(rowID=rownames(baySeq.cd@data),row.names=rownames(baySeq.cd@data))"), codefile)
  writeLines(c("baySeq.posteriors.DE <- exp(baySeq.cd@posteriors)[, 'DE']",
               "baySeq.FDR <- baySeq::topCounts(baySeq.cd, group = 'DE', FDR = 1)$FDR.DE[match(rownames(count.matrix(cdata)), rownames(baySeq::topCounts(baySeq.cd, group = 'DE', FDR = 1)))]",
               "baySeq.score <- 1 - baySeq.FDR",
               "result.table <- data.frame('FDR' = baySeq.FDR, 'score' = baySeq.score, 'posterior.DE' = baySeq.posteriors.DE)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('baySeq,', packageVersion('baySeq'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'baySeq', 'full.name' = '", paste('baySeq.', packageVersion('baySeq'), '.', norm.method, '.', ifelse(equaldisp == TRUE, 'equaldisp', 'nonequaldisp'), '.samplesize', sample.size, '.', estimation, '.', pET, sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}
