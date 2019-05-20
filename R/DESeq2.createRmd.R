#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with DESeq2
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the DESeq2 package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{DESeq2} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param fit.type The fitting method used to get the dispersion-mean relationship. Possible values are \code{"parametric"}, \code{"local"} and \code{"mean"}.
#' @param test The test to use. Possible values are \code{"Wald"} and \code{"LRT"}.
#' @param beta.prior Whether or not to put a zero-mean normal prior on the non-intercept coefficients. Default is \code{TRUE}.
#' @param independent.filtering Whether or not to perform independent filtering of the data. With independent filtering=TRUE, the adjusted p-values for genes not passing the filter threshold are set to NA.
#' @param cooks.cutoff The cutoff value for the Cook's distance to consider a value to be an outlier. Set to Inf or FALSE to disable outlier detection. For genes with detected outliers, the p-value and adjusted p-value will be set to NA.
#' @param impute.outliers Whether or not the outliers should be replaced by a trimmed mean and the analysis rerun.
#' @export
#' @author Charlotte Soneson
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references
#' Anders S and Huber W (2010): Differential expression analysis for sequence count data. Genome Biology 11:R106
#' @examples
#' try(
#' if (require(DESeq2)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' SyntheticDataSimulation(simul.data= 'KIRC', dataset = file.path(tmpdir, "mydata.rds"), n.var = 500,
#'                                     samples.per.cond = 3, n.diffexp = 100, dispType='same', mode='D',
#'                                     fraction.upregulated=0.5, dataset.parameters=generateDatasetParameter())
#' compcodeR::runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "DESeq2",
#'            Rmdfunction = "DESeq2.createRmd",
#'            output.directory = tmpdir, fit.type = "parametric",
#'            test = "Wald")
#' })
DESeq2.createRmd <- function(data.path, result.path, codefile,
                             fit.type, test, beta.prior = TRUE,
                             independent.filtering = TRUE, cooks.cutoff = TRUE,
                             impute.outliers = TRUE) {
  codefile <- file(codefile, open = 'w')
  writeLines("### DESeq2", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = TRUE, error = TRUE, warning = TRUE}",
               "require(DESeq2)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               "DESeq2.ds <- DESeq2::DESeqDataSetFromMatrix(countData = count.matrix(cdata), colData = data.frame(condition = factor(sample.annotations(cdata)$condition)), design = ~ condition)",
               paste("DESeq2.ds <- DESeq2::DESeq(DESeq2.ds, fitType = '", fit.type, "', test = '", test, "', betaPrior = ", beta.prior, ")", sep = "")), codefile)
  if (impute.outliers == TRUE) {
    writeLines(c("DESeq2.ds.clean <- DESeq2::replaceOutliersWithTrimmedMean(DESeq2.ds)",
                 paste("DESeq2.ds.clean <- DESeq2::DESeq(DESeq2.ds.clean, fitType = '", fit.type, "', test = '", test, "', betaPrior = ", beta.prior, ")", sep = ""),
                 "DESeq2.ds <- DESeq2.ds.clean"), codefile)
  }
  writeLines(c(paste("DESeq2.results <- DESeq2::results(DESeq2.ds, independentFiltering = ", independent.filtering, ", cooksCutoff = ", cooks.cutoff, ")", sep = ""),
               "DESeq2.pvalues <- DESeq2.results$pvalue",
               "DESeq2.adjpvalues <- DESeq2.results$padj",
               "DESeq2.logFC <- DESeq2.results$log2FoldChange",
               "DESeq2.score <- 1 - DESeq2.pvalues",
               "result.table <- data.frame('pvalue' = DESeq2.pvalues, 'adjpvalue' = DESeq2.adjpvalues, 'logFC' = DESeq2.logFC, 'score' = DESeq2.score)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('DESeq2,', packageVersion('DESeq2'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'DESeq2', 'full.name' = '", paste('DESeq2.', packageVersion('DESeq2'), '.', fit.type, '.', test, '.', ifelse(beta.prior == TRUE, 'bp', 'nobp'), '.', ifelse(independent.filtering == TRUE, 'indf', 'noindf'), paste(".cook_", cooks.cutoff, sep = ""), ifelse(impute.outliers, ".imp", ".noimp"), sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}
