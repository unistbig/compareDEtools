#' make correlation heatmap function
#' @param working.dir Input file location.
#' @param figure.dir Figure save location.
#' @param simul.data Type of dataset(e.g. KIRC and Bottomly)
#' @param nsample Number of samples. Input as a number.
#' @param topgenes Number of top significant genes obtained from each methods to estimate rank correlation.
#' @param AnalysisMethods DEmethods used for figures. Input as a character vector. (e.g. 'edgeR','edgeR.ql','edgeR.rb','DESeq.pd','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')
#' @param rep Number of replication each test contain.
#' @export
correlation_heatmap<-function(working.dir, figure.dir ,simul.data, nsample, topgenes, AnalysisMethods, rep){
  # working.dir='/home/node01/test/real/'
  # figure.dir='/home/node01/test/fig/'
  # simul.data='KIRC'
  # nsample=5
  # AnalysisMethods=AnalysisMethods
  # rep=1
  if(simul.data=='KIRC'){
    nvar=16621
  }else if(simul.data=='Bottomly'){
    nvar=8550
  }
  DEmethods<-unlist(sapply(AnalysisMethods,DEanalysis::select_tool),use.names = FALSE)
  whole_corr_matrix=array(rep(0,length(rep*length(AnalysisMethods)*length(AnalysisMethods))) ,dim=c(rep,length(AnalysisMethods),length(AnalysisMethods)))

  for(k in 1:rep){
    result_matrix=matrix(0,ncol=length(AnalysisMethods),nrow=nvar)
    rank_matrix=matrix(0,ncol=length(AnalysisMethods),nrow=nvar)

    for(i in 1:length(AnalysisMethods)){
      fileName=paste(working.dir, simul.data, '_', nsample,'spc_rep_',k,'_',DEmethods[i],'.rds',sep='')
      result = try(readRDS(fileName), silent=T)
      if(class(result)=='try-error'){next}
      result = result@result.table
      if(nrow(result)==0){next}
      if(DEmethods[i] == "PoissonSeq"){rownames(result)=as.character(result$Genename)}

      if(!is.null(result$FDR)){
        FDR=result$FDR
      }else if(!is.null(result$adjpvalue)){
        FDR=result$adjpvalue
      }

      if(DEmethods[i]=='PoissonSeq'){
        real_FDR=rep(0,nvar)
        for(j in 1:nvar){
          if(paste('g',j,sep='') %in% result$Genename){
            real_FDR[j]=FDR[match(paste('g',j,sep=''), result$Genename)]
          }else{
            real_FDR[j]=1
          }
        }
        FDR<-real_FDR
      }
      rank_matrix[,i]=order(order(FDR))
      result_matrix[,i]=(FDR)
    }
    inter=vector()
    for(i in 1:length(AnalysisMethods)){
      temp=which(rank_matrix[,i]<=topgenes)
      inter=union(inter,temp)
    }

    result_matrix<-result_matrix[inter,]

    colnames(result_matrix)=AnalysisMethods
    rownames(result_matrix)=paste("g",1:dim(result_matrix)[1],sep="")
    cor_matrix=matrix(0,ncol=length(AnalysisMethods),nrow=length(AnalysisMethods))
    rownames(cor_matrix)=AnalysisMethods
    colnames(cor_matrix)=AnalysisMethods
    for(i in 1:length(AnalysisMethods)){
      for(j in 1:length(AnalysisMethods)){
        cor_matrix[i,j]=(cor.test( x=result_matrix[,i],y=result_matrix[,j],
                                   method = "spearman",
                                   continuity = FALSE,
                                   conf.level = 0.95)$estimate)
        if(cor_matrix[i,j]<0){
          cor_matrix[i,j]=0
        }
      }
    }
    whole_corr_matrix[k,,]<-cor_matrix
  }
  correlation=matrix(0,ncol=length(AnalysisMethods),nrow=length(AnalysisMethods))
  for(i in 1:length(AnalysisMethods)){
    for(j in 1:length(AnalysisMethods)){
      correlation[i,j]=mean(whole_corr_matrix[,i,j])
    }
  }
  rownames(correlation)<-AnalysisMethods
  colnames(correlation)<-AnalysisMethods

  pdf(paste(figure.dir,simul.data,'_top ',topgenes,'_clustering heatmap.pdf',sep=''))
  colors = c(seq(0,1,length=50))
  my_palette <- colorRampPalette(c("white",brewer.pal(9,"Blues")))(n = 49)
  heatmap.2(correlation,distfun = function(x) as.dist(1-x),  symm=TRUE,keysize=1,density.info=c('none'),trace='none',revC=TRUE, breaks = colors, col=my_palette)
  dev.off()
}


