#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with the edgeR GLM Robust dispersion estimation
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the edgeR package.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{edgeR} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}.
#' @param disp.estimation method for estimating the dispersion. Possible values are "CoxReid", "Pearson" or "deviance".
#' @references
#' Robinson MD, McCarthy DJ and Smyth GK (2010): edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#'
#' Zhou X, Lindsay H, Robinson MD (2014). Robustly detecting differential expression in RNA sequencing data using observation weights. Nucleic Acids Research, 42(11), e91.
#' @export
#' @author BuKyung Baik
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' SyntheticDataSimulation(simul.data= 'KIRC', dataset = file.path(tmpdir, "mydata.rds"), n.var = 500,
#'                                     samples.per.cond = 3, n.diffexp = 100, dispType='same', mode='D',
#'                                     fraction.upregulated=0.5, dataset.parameters=generateDatasetParameter())
#' compcodeR::runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "edgeR.GLM.Robust",
#'            Rmdfunction = "edgeR.GLM.Robust.createRmd",
#'            output.directory = tmpdir)
edgeR.GLM.Robust.createRmd <- function(data.path, result.path, codefile,
                                       norm.method = 'TMM', disp.estimation = 'CoxReid') {

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

  writeLines(paste("edgeR.dgelist <- edgeR::estimateGLMRobustDisp(edgeR.dgelist, design = model.matrix(~factor(sample.annotations(cdata)$condition)), trend.method = 'auto')", sep = ''), codefile)

  writeLines(c("dispersion.S1 <- edgeR.dgelist$tagwise.dispersion",
               "dispersion.S2 <- edgeR.dgelist$tagwise.dispersion"), codefile)


  writeLines(c("edgeR.glmfit <- edgeR::glmFit(edgeR.dgelist, design = model.matrix(~factor(sample.annotations(cdata)$condition)))",
               "edgeR.glmlrt <- edgeR::glmLRT(edgeR.glmfit, coef = 2)",
               "edgeR.pvalues <- edgeR.glmlrt$table$PValue",
               "edgeR.adjpvalues <- p.adjust(edgeR.pvalues, method = 'BH')",
               "edgeR.logFC <- edgeR.glmlrt$table$logFC",
               "edgeR.score <- 1 - edgeR.pvalues",
               "result.table <- data.frame('pvalue' = edgeR.pvalues, 'adjpvalue' = edgeR.adjpvalues, 'logFC' = edgeR.logFC, 'score' = edgeR.score, 'dispersion.S1' = dispersion.S1, 'dispersion.S2' = dispersion.S2)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('edgeR,', packageVersion('edgeR'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'edgeR', 'full.name' = '", paste('edgeR.', packageVersion('edgeR'), '.GLM.', norm.method, '.trend.', disp.estimation, '.', 'tagwise', sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}
