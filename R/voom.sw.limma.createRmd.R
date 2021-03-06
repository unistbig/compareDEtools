#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with the quantile-normalized voom+limma
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) by applying the voom transformation with quantile normalization (from the limma package) followed by differential expression analysis with limma.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{limma} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. The normalization factors are calculated using the \code{calcNormFactors} of the \code{edgeR} package. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}
#' @references
#' Smyth GK (2005): Limma: linear models for microarray data. In: 'Bioinformatics and Computational Biology Solutions using R and Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds), Springer, New York, pages 397-420
#'
#' Law CW, Chen Y, Shi W and Smyth GK (2014): voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology 15, R29
#'
#' Liu, R., Holik, A. Z., Su, S., Jansz, N., Chen, K., Leong, H. S., Blewitt, M. E., Asselin-Labat, M.-L., Smyth, G. K., Ritchie, M. E. (2015). Why weight? Combining voom with estimates of sample quality improves power in RNA-seq analyses. Nucleic Acids Research 43, e97.
#'
#' @export
#' @author BuKyung Baik
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' SyntheticDataSimulation(simul.data= 'KIRC', dataset = file.path(tmpdir, "mydata.rds"), n.var = 500,
#'                                     samples.per.cond = 3, n.diffexp = 100, dispType='same', mode='D',
#'                                     fraction.upregulated=0.5, dataset.parameters=generateDatasetParameter())
#' compcodeR::runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "voom.sw.limma",
#'            Rmdfunction = "voom.sw.limma.createRmd",
#'            output.directory = tmpdir,)
voom.sw.limma.createRmd <- function(data.path, result.path, codefile, norm.method = 'TMM') {
  codefile <- file(codefile, open = 'w')
  writeLines("### voom + limma", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(limma)",
               "require(edgeR)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               "pt<-proc.time()",
               paste("nf <- edgeR::calcNormFactors(count.matrix(cdata), method = '", norm.method, "')", sep = ''),
               "voom.data <- limma::voom(count.matrix(cdata), design = model.matrix(~factor(sample.annotations(cdata)$condition)), lib.size = colSums(count.matrix(cdata)) * nf)",
               "voom.data$genes <- rownames(count.matrix(cdata))",
               "voom.sw<-limma::voomWithQualityWeights(count.matrix(cdata),design = model.matrix(~factor(sample.annotations(cdata)$condition)), lib.size = colSums(count.matrix(cdata)) * nf)",
               "voom.fitlimma <- limma::lmFit(voom.sw, design = model.matrix(~factor(sample.annotations(cdata)$condition)))",
               "voom.fitbayes <- limma::eBayes(voom.fitlimma)",
               "voom.pvalues <- voom.fitbayes$p.value[, 2]",
               "voom.adjpvalues <- p.adjust(voom.pvalues, method = 'BH')",
               "cpu.time<-as.numeric(proc.time()-pt)[3]",
               "voom.logFC <- voom.fitbayes$coefficients[, 2]",
               "voom.score <- 1 - voom.pvalues",
               "result.table <- data.frame('pvalue' = voom.pvalues, 'adjpvalue' = voom.adjpvalues, 'logFC' = voom.logFC, 'score' = voom.score, 'cpu.time'=cpu.time)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('limma,', packageVersion('limma'), ';', 'edgeR,', packageVersion('edgeR'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'voom', 'full.name' = '",
                     paste('voom.', packageVersion('limma'), '.limma.', norm.method , sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}
