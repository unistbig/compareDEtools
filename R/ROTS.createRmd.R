#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with Reproducibility-Optimized Test Statistic(ROTS)
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the ROTS package.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{ROTS} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param normalize A logical indicating whether to run ROTS with TMM normalized count or raw count data. Default value is True (TMM normalized).
#' @param transformation A logical indicating whether to run ROTS with voom transformation. Default value is True.
#' @param B An integer specifying the number of bootstrap and permutation resamplings (default 1000).
#' @param K An integer indicating the largest top list size considered. If no value is given, 1/4 of the features are used.
#' @param log A logical parameter indicating whether input data is log2 scaled. Choices are as for the method argument of ROTS function. This information is only used to calculate log fold change.
#' @references
#' L. L. Elo, S. Filen, R. Lahesmaa and T. Aittokallio: Reproducibility-optimized test statistic for ranking genes in microarray studies. IEEE/ACM Transactions on Computational Biology and Bioinformatics 5: 423–431, 2008.
#' @export
#' @author BuKyung Baik
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' SyntheticDataSimulation(simul.data= 'KIRC', dataset = file.path(tmpdir, "mydata.rds"), n.var = 500,
#'                                     samples.per.cond = 3, n.diffexp = 100, dispType='same', mode='D',
#'                                     fraction.upregulated=0.5, dataset.parameters=generateDatasetParameter())
#' compcodeR::runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "ROTS",
#'            Rmdfunction = "ROTS.createRmd", transformation=TRUE, normalize = TRUE,
#'            output.directory = tmpdir, B = 1000,
#'            K = 500, log = FALSE)
ROTS.createRmd <- function(data.path, result.path, codefile, normalize=TRUE, transformation=TRUE,
                           B = 1000, K = NULL, log = FALSE){

  ## Write the code for applying edgeR exact test to an Rmd file
  codefile <- file(codefile, open = "w")
  writeLines("### ROTS", codefile)
  writeLines(paste("Data file: ", data.path, sep = ""), codefile)

  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(edgeR)",
               "require(ROTS)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')"), codefile)
  if(normalize&&transformation){
    if(transformation){
      log=TRUE
      writeLines(c("nf <- edgeR::calcNormFactors(count.matrix(cdata), method = 'TMM' )",
                   "voom.data <- limma::voom(count.matrix(cdata), design = model.matrix(~factor(sample.annotations(cdata)$condition)), lib.size = colSums(count.matrix(cdata)) * nf)",
                   "Exp <- voom.data$E"), codefile)
    }else{
      writeLines(c("edgeR.dgelist <- edgeR::DGEList(counts = count.matrix(cdata), group = factor(sample.annotations(cdata)$condition))",
                   "edgeR.dgelist <- edgeR::calcNormFactors(edgeR.dgelist, method = 'TMM')",
                   "Factors <- edgeR.dgelist$samples$lib.size * edgeR.dgelist$samples$norm.factors",
                   "Exp <- t(t(edgeR.dgelist$counts)/Factors) * mean(Factors)"), codefile)
    }
  }else{
    writeLines('Exp <- count.matrix(cdata)', codefile)
  }


  writeLines(c(paste("results <- ROTS::ROTS(data = Exp, groups = as.numeric(as.character(factor(sample.annotations(cdata)$condition))), K = " ,K , ", B = " ,B , ", log = ",log , ")",sep=''),
               "rots.pvalues <- results$pvalue",
               "rots.logFC <- results$logfc",
               "rots.FDR <- results$FDR",
               "result.table <- data.frame('pvalue' = rots.pvalues, 'logFC' = rots.logFC, 'FDR' = rots.FDR)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('ROTS,', packageVersion('ROTS'), ';', 'edgeR,', packageVersion('edgeR'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'rots', 'full.name' = '",
                     paste('rots.', packageVersion('ROTS'), if(normalize){'TMM'}, '.', K, '.', B, '.', if(log){'log2 scaled'}, sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}
