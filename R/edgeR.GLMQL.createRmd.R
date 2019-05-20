#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with the edgeR GLM Quasi-likelihood Tests
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the edgeR package.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{edgeR} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method A Character parameter for calcNormFactors function. Possible values are "TMM", "TMMwzp", "RLE", "upperquartile", or "none".
#' @param disp.estimation A Character parameter for estimateGLMCommonDisp function. Possible values are "CoxReid", "Pearson" or "deviance".
#' @references
#' Robinson MD, McCarthy DJ and Smyth GK (2010): edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#'
#' Lun, ATL, Chen, Y, and Smyth, GK (2016). It’s DE-licious: a recipe for differential expression analyses of RNA-seq experiments using quasi-likelihood methods in edgeR. Methods in Molecular Biology 1418, 391–416.
#' @export
#' @author BuKyung Baik
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' SyntheticDataSimulation(simul.data= 'KIRC', dataset = file.path(tmpdir, "mydata.rds"), n.var = 500,
#'                                     samples.per.cond = 3, n.diffexp = 100, dispType='same', mode='D',
#'                                     fraction.upregulated=0.5, dataset.parameters=generateDatasetParameter())
#' compcodeR::runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "edgeR.GLMQL",
#'            Rmdfunction = "edgeR.GLMQL.createRmd",
#'            output.directory = tmpdir)
edgeR.GLMQL.createRmd <- function(data.path, result.path, codefile,
                                  norm.method='TMM', disp.estimation='CoxReid') {


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
               paste("edgeR.dgelist <- edgeR::estimateGLMCommonDisp(edgeR.dgelist, design = model.matrix(~factor(sample.annotations(cdata)$condition)), method = '", disp.estimation, "')", sep = '')), codefile)
  writeLines("edgeR.dgelist <- edgeR::estimateGLMTrendedDisp(edgeR.dgelist, design = model.matrix(~factor(sample.annotations(cdata)$condition)), method = 'auto')", codefile)



  writeLines(c("edgeR.glmqlfit <- edgeR::glmQLFit(edgeR.dgelist, design = model.matrix(~factor(sample.annotations(cdata)$condition)), robust=FALSE)",
               "edgeR.glmqlft <- edgeR::glmQLFTest(edgeR.glmqlfit, coef = 2)",
               "edgeR.pvalues <- edgeR.glmqlft$table$PValue",
               "edgeR.adjpvalues <- p.adjust(edgeR.pvalues, method = 'BH')",
               "edgeR.logFC <- edgeR.glmqlft$table$logFC",
               "edgeR.score <- 1 - edgeR.pvalues",
               "result.table <- data.frame('pvalue' = edgeR.pvalues, 'adjpvalue' = edgeR.adjpvalues, 'logFC' = edgeR.logFC, 'score' = edgeR.score, 'trended_dispersion' = edgeR.dgelist$trended.dispersion )",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('edgeR,', packageVersion('edgeR'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'edgeR', 'full.name' = '", paste('edgeR.', packageVersion('edgeR'), '.GLM.', norm.method, '.trended.', disp.estimation, '.', sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}
