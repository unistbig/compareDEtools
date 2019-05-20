#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with SAMseq
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using SAMseq. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{SAMseq} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @export
#' @author Charlotte Soneson
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references
#' Li J and Tibshirani R (2011): Finding consistent patterns: a nonparametric approach for identifying differential expression in RNA-Seq data. Statistical Methods in Medical Research
#' @examples
#' try(
#' if (require(samr)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' SyntheticDataSimulation(simul.data= 'KIRC', dataset = file.path(tmpdir, "mydata.rds"), n.var = 500,
#'                                     samples.per.cond = 3, n.diffexp = 100, dispType='same', mode='D',
#'                                     fraction.upregulated=0.5, dataset.parameters=generateDatasetParameter())
#' compcodeR::runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "SAMseq",
#'            Rmdfunction = "SAMseq.createRmd",
#'            output.directory = tmpdir)
#' })
SAMseq.createRmd <- function(data.path, result.path, codefile) {
  codefile <- file(codefile, open = 'w')
  writeLines("### SAMseq (samr)", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(samr)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               "condition.12 <- rep(1, length(sample.annotations(cdata)$condition))",
               "condition.12[which(sample.annotations(cdata)$condition == levels(factor(sample.annotations(cdata)$condition))[2])] <- 2",
               "SAMseq.test <- samr::SAMseq(count.matrix(cdata), condition.12, resp.type = 'Two class unpaired', geneid = rownames(count.matrix(cdata)), genenames = rownames(count.matrix(cdata)), nperms = 100, nresamp = 20, fdr.output = 1)",
               "SAMseq.result <- rbind(SAMseq.test$siggenes.table$genes.up, SAMseq.test$siggenes.table$genes.lo)",
               "SAMseq.statistic <- rep(0, nrow(count.matrix(cdata)))",
               "SAMseq.statistic[match(SAMseq.result[, 1], rownames(count.matrix(cdata)))] <- as.numeric(SAMseq.result[, 3])",
               "SAMseq.FDR <- rep(1, nrow(count.matrix(cdata)))",
               "SAMseq.FDR[match(SAMseq.result[, 1], rownames(count.matrix(cdata)))] <- as.numeric(SAMseq.result[, 5])/100",
               "SAMseq.score <- 1 - SAMseq.FDR",
               "result.table <- data.frame('statistic' = SAMseq.statistic, 'FDR' = SAMseq.FDR, 'score' = SAMseq.score)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('SAMseq,', packageVersion('samr'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'SAMseq', 'full.name' = '",
                     paste('SAMseq.', packageVersion('samr'), sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}
