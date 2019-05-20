#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with a Poisson log linear model
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the PoissonSeq package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{PoissonSeq} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @references
#' Li J, Witten DM, Johnstone I, Tibshirani R (2012). Normalization, testing, and false discovery rate estimation for RNA-sequencing data. Biostatistics 13(3): 523-38.
#' @export
#' @author BuKyung Baik
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' SyntheticDataSimulation(simul.data= 'KIRC', dataset = file.path(tmpdir, "mydata.rds"), n.var = 500,
#'                                     samples.per.cond = 3, n.diffexp = 100, dispType='same', mode='D',
#'                                     fraction.upregulated=0.5, dataset.parameters=generateDatasetParameter())
#' compcodeR::runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "PoissonSeq",
#'            Rmdfunction = "PoissonSeq.createRmd",
#'            output.directory = tmpdir)
PoissonSeq.createRmd = function(data.path, result.path, codefile)
{
  codefile = file(codefile, open="w")
  cdata = readRDS(data.path)
  writeLines("### PoissonSeq", codefile)
  writeLines(paste("Data file: ", data.path, sep=""), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = TRUE, error = TRUE, warning = TRUE}",
               "require(PoissonSeq)", paste("cdata <- readRDS('", data.path,"')", sep = "")), codefile)
  if(is.list(cdata)) {
    writeLines("cdata = convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)", "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               "dat = list(n = count.matrix(cdata), y=sample.annotations(cdata)$condition, pair = FALSE, type = 'twoclass', gname=rownames(count.matrix(cdata)))",
               "PS.test = PoissonSeq::PS.Main(dat)"),codefile)

  writeLines(c("PS.nc = PS.test$nc",
               "PS.gname = rownames(PS.test)",
               "PS.tt = PS.test$tt",
               "PS.pvalue = PS.test$pval",
               "PS.FDR = PS.test$fdr",
               "PS.logfc =  PS.test$log.fc",
               "result.table = data.frame('nc' = PS.nc, 'Genename' = PS.gname, 'tt' = PS.tt, 'pvalue' = PS.pvalue, 'FDR'=PS.FDR, 'log.fc' = PS.logfc)"), codefile)

  writeLines(c("result.table(cdata) = result.table",
               "package.version(cdata) = paste('PoissonSeq', packageVersion('PoissonSeq'))",
               "analysis.date(cdata) = date()",
               paste("method.names(cdata) <- list('short.name' = 'PoissonSeq', 'full.name' = '", paste('PoissonSeq.', packageVersion('PoissonSeq'), sep = ''), "')", sep = ''),
               "is.valid = check_compData_results(cdata)",
               "if(!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep="")), codefile)

  writeLines("print(paste('Unique data set ID:', info.parameter(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}
