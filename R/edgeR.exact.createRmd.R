#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with the edgeR exact test
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the exact test functionality from the edgeR package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{edgeR} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}.
#' @param trend.method The method used to estimate the trend in the mean-dispersion relationship. Possible values are \code{"none"}, \code{"movingave"} and \code{"loess"}
#' @param disp.type The type of dispersion estimate used. Possible values are \code{"common"}, \code{"trended"} and \code{"tagwise"}.
#' @references
#' Robinson MD, McCarthy DJ and Smyth GK (2010): edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#' @export
#' @author Charlotte Soneson
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' SyntheticDataSimulation(simul.data= 'KIRC', dataset = file.path(tmpdir, "mydata.rds"), n.var = 500,
#'                                     samples.per.cond = 3, n.diffexp = 100, dispType='same', mode='D',
#'                                     fraction.upregulated=0.5, dataset.parameters=generateDatasetParameter())
#' compcodeR::runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "edgeR.exact",
#'            Rmdfunction = "edgeR.exact.createRmd",
#'            output.directory = tmpdir, norm.method = "TMM",
#'            trend.method = "movingave", disp.type = "tagwise")
edgeR.exact.createRmd <- function(data.path, result.path, codefile,
                                  norm.method, trend.method, disp.type) {
  ## Write the code for applying edgeR exact test to an Rmd file
  codefile <- file(codefile, open = "w")
  writeLines("### edgeR", codefile)
  writeLines(paste("Data file: ", data.path, sep = ""), codefile)

  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(edgeR)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               "edgeR.dgelist <- edgeR::DGEList(counts = count.matrix(cdata), group = factor(sample.annotations(cdata)$condition))",
               paste("edgeR.dgelist <- edgeR::calcNormFactors(edgeR.dgelist, method = '", norm.method, "')", sep = ''),
               "edgeR.dgelist <- edgeR::estimateCommonDisp(edgeR.dgelist)"), codefile)
  if (disp.type == 'tagwise') {
    writeLines(c(paste("edgeR.dgelist <- edgeR::estimateTagwiseDisp(edgeR.dgelist, trend = '", trend.method, "')", sep = ''),
                 "dispersion.S1 <- edgeR.dgelist$tagwise.dispersion",
                 "dispersion.S2 <- edgeR.dgelist$tagwise.dispersion"), codefile)
  } else {
    if (disp.type == 'common') {
      writeLines(c("dispersion.S1 <- rep(edgeR.dgelist$common.dispersion, nrow(count.matrix(cdata)))",
                   "dispersion.S2 <- rep(edgeR.dgelist$common.dispersion, nrow(count.matrix(cdata)))"), codefile)
    } else if (disp.type == 'trended') {
      writeLines(c("edgeR.dgelist <- edgeR::estimateTrendedDisp(edgeR.dgelist)",
                   "dispersion.S1 <- edgeR.dgelist$trended.dispersion",
                   "dispersion.S2 <- edgeR.dgelist$trended.dispersion"), codefile)
    }
  }
  writeLines(c("edgeR.exacttest <- edgeR::exactTest(edgeR.dgelist)",
               "edgeR.pvalues <- edgeR.exacttest$table$PValue",
               "edgeR.adjpvalues <- p.adjust(edgeR.pvalues, method = 'BH')",
               "edgeR.logFC <- edgeR.exacttest$table$logFC",
               "edgeR.score <- 1 - edgeR.pvalues",
               "result.table <- data.frame('pvalue' = edgeR.pvalues, 'adjpvalue' = edgeR.adjpvalues, 'logFC' = edgeR.logFC, 'score' = edgeR.score, 'dispersion.S1' = dispersion.S1, 'dispersion.S2' = dispersion.S2)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('edgeR,', packageVersion('edgeR'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'edgeR', 'full.name' = '", paste('edgeR.', packageVersion('edgeR'), '.exact.', norm.method, '.', trend.method, '.', disp.type, sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}
