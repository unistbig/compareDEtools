#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with the DESeq nbinom approach
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the nbinom test from the DESeq package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{DESeq} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param sharing.mode The method used to select between the individually estimated dispersion and the dispersion estimate obtained by fitting a dispersion-mean relationship to the estimated values for all genes. Possible values are \code{"fit-only"} (use the fitted value), \code{"maximum"} (take the maximum of the fitted and the estimated value) and \code{"gene-est-only"} (use the estimated value).
#' @param disp.method The method used to estimate the dispersion. Possible values are \code{"pooled"}, \code{"per-condition"} and \code{"blind"}.
#' @param fit.type The fitting method used to get the dispersion-mean relationship. Possible values are \code{"parametric"} and \code{"local"}.
#' @export
#' @author Charlotte Soneson
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references
#' Anders S and Huber W (2010): Differential expression analysis for sequence count data. Genome Biology 11:R106
#' @examples
#' try(
#' if (require(DESeq)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' SyntheticDataSimulation(simul.data= 'KIRC', dataset = file.path(tmpdir, "mydata.rds"), n.var = 500,
#'                                     samples.per.cond = 3, n.diffexp = 100, dispType='same', mode='D',
#'                                     fraction.upregulated=0.5, dataset.parameters=generateDatasetParameter())
#' compcodeR::runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "DESeq.nbinom",
#'            Rmdfunction = "DESeq.nbinom.createRmd",
#'            output.directory = tmpdir, sharing.mode = "maximum",
#'            disp.method = "pooled", fit.type = "parametric")
#' })
DESeq.nbinom.createRmd <- function(data.path, result.path, codefile,
                                   sharing.mode, disp.method, fit.type) {
  codefile <- file(codefile, open = 'w')
  cdata <- readRDS(data.path)
  writeLines("### DESeq", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = TRUE, error = TRUE, warning = TRUE}",
               "require(DESeq)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(cdata)) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               "DESeq.cds <- DESeq::newCountDataSet(countData = count.matrix(cdata), conditions = factor(sample.annotations(cdata)$condition))",
               "DESeq.cds <- DESeq::estimateSizeFactors(DESeq.cds)",
               paste("DESeq.cds <- DESeq::estimateDispersions(DESeq.cds, sharingMode = '", sharing.mode, "', method = '",
                     disp.method, "', fitType = '", fit.type, "')", sep = '')), codefile)
  if (disp.method == 'pooled') {
    writeLines(c("dispersion.S1 <- fData(DESeq.cds)$disp_pooled",
                 "dispersion.S2 <- fData(DESeq.cds)$disp_pooled"), codefile)
  } else {
    if (disp.method == 'blind') {
      writeLines(c("dispersion.S1 <- fData(DESeq.cds)$disp_blind",
                   "dispersion.S2 <- fData(DESeq.cds)$disp_blind"), codefile)
    } else {
      writeLines(c(paste("dispersion.S1 <- fData(DESeq.cds)['disp_", levels(factor(cdata@sample.annotations$condition))[1], "']", sep = ""),
                   paste("dispersion.S2 <- fData(DESeq.cds)['disp_", levels(factor(cdata@sample.annotations$condition))[2], "']", sep = "")), codefile)
    }
  }
  writeLines(c("DESeq.test <- DESeq::nbinomTest(DESeq.cds, levels(factor(sample.annotations(cdata)$condition))[1], levels(factor(sample.annotations(cdata)$condition))[2])",
               "DESeq.pvalues <- DESeq.test$pval",
               "DESeq.adjpvalues <- p.adjust(DESeq.pvalues, method = 'BH')",
               "DESeq.logFC <- DESeq.test$log2FoldChange",
               "DESeq.score <- 1 - DESeq.pvalues",
               "result.table <- data.frame('pvalue' = DESeq.pvalues, 'adjpvalue' = DESeq.adjpvalues, 'logFC' = DESeq.logFC, 'score' = DESeq.score, 'dispersion.S1' = dispersion.S1, 'dispersion.S2' = dispersion.S2)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('DESeq,', packageVersion('DESeq'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'DESeq', 'full.name' = '", paste('DESeq.', packageVersion('DESeq'), '.nbinom.', disp.method, '.', sharing.mode, '.', fit.type, sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}
