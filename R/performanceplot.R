#' tool selector
#'
#' Take input methods and return extensions of analysis rds files.
#' @param method_in A character parameter indicating the object method to be switched into its file extension.
#' @export
select_tool = function(method_in)
{
  answer = switch(method_in,
                  "DESeq.pc"="DESeq.pc",
                  "DESeq2"="DESeq2",
                  "edgeR"="edgeR.exact",
                  "edgeR.ql"="edgeR.GLMQL",
                  "edgeR.rb"="edgeR.GLM.Robust",
                  "voom.tmm"="voom.limma",
                  "voom.qn"="voom.qn.limma",
                  "voom.sw"="voom.sw.limma",
                  "ROTS"="ROTS",
                  "BaySeq"="BaySeq",
                  "BaySeq.qn"="BaySeq.qn",
                  "PoissonSeq"="PoissonSeq",
                  "SAMseq"="SAMseq"
  )
  return(answer)
}

#' color selector
#'
#' Take input methods and return pre-assigned color for figure.
#' @param method_color A character parameter indicating the object method to be switched into its corresponding color in the performance plot.
#' @export
select_color = function(method_color)
{
  answer = switch(method_color,
                  "DESeq.pc"="#801279",
                  "DESeq2"="#1B0FFF",
                  "edgeR"="#1BAFFF",
                  "edgeR.glm"="#0FC4FF",
                  "edgeR.ql"="#022604",
                  "edgeR.rb"="#1E7E23",
                  "voom.tmm"="#FCFF0F",
                  "voom.qn"="#FF9933",
                  "voom.sw"= "#E094D9",
                  "ROTS"="#CE00B3",
                  "BaySeq"="#6D700D",
                  "BaySeq.qn"="#28D1B0",
                  "PoissonSeq"="#0FFF63",
                  "SAMseq"="#FF0000"
  )
  return(answer)
}




#' make synthetic data plot function
#' @param working.dir Input file location
#' @param figure.dir Figure save location
#' @param fixedfold A logical indicating whether analyzed dataset used random fold changes following exponential distribution or fixed fold change values. Possible values are TRUE or FALSE. If fixedfold is TRUE, fraction.upregulated is automatically fixed to 0.67.
#' @param simul.data Type of dataset (e.g. KIRC, Bottomly, mBdK and mKdB)
#' @param rep An integer specifying iterations DE analysis methods run for each condition.
#' @param nsample An integer vector indicating number of samples in each sample group.
#' @param nvar An integer indicating how many genes are in the analyzed dataset.
#' @param nDE An integer vector indicating how many DE genes are in the analyzed dataset.
#' @param fraction.upregulated A numeric vector specifying proportions of upregulated DE genes among total DE genes in the analzyed dataset. (e.g. 0.5, 0.7 and 0.9)
#' @param disp.Type A vector indicating how is the dispersion parameter assumed to be for each sample group in the analzyed dataset. Possible values are 'same' and 'different'.
#' @param mode A character specifying a test condition used for analyzed dataset.
#' "D" for basic simulation (not adding outliers).
#' "R" for adding 5% of random outlier.
#' "OS" for adding outlier sample to each sample group.
#' "DL" for decreasing KIRC simulation dispersion 22.5 times (similar to SEQC data dispersion) to compare with SEQC data.
#' @param rowType A character vector indicating which results are shown in performance plot. Combination of AUC, TPR and trueFDR. (e.g. c('AUC','TPR'))
#' @param AnalysisMethods A character vector specifying DE methods used for the analysis.
#' (e.g. 'edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','BaySeq.qn','PoissonSeq','SAMseq')
#' @export
performance_plot = function(working.dir, figure.dir, fixedfold=FALSE, simul.data, rep, nsample, nvar, nDE, fraction.upregulated, disp.Type, mode, rowType, AnalysisMethods){


  if(length(simul.data)!=1){stop('simul.data must have one element.')}
  if(length(mode)!=1){stop('mode must have one element.')}
  if(length(disp.Type)!=1){stop('disp.Type must have one element.')}
  if(fixedfold){fraction.upregulated=0.67}

  if(disp.Type=='same'){
    test.cond=mode
  }else if(disp.Type=='different'){
    test.cond=paste('DiffDisp_', mode, sep='')
  }
  type = switch(test.cond, 'D'='No Random outlier test / same dispersion ',
                'DiffDisp_D'= 'No Random outlier test / different dispersion ',
                'R'='Random outlier test / same dispersion ',
                'DiffDisp_R' = 'Random outlier test / different dispersion ',
                'OS' = 'Outlier dispersion sample test / same dispersion ',
                'DiffDisp_OS' = 'Outlier dispersion sample test / different dispersion ')


  tpr = tfdr = auc = ts = tc = NDE = NSAMPLE = UPPROP = METHOD = COLOR = REPEAT = nDE_Factor =NULL

  for(nde in nDE)
  {
    for(s in nsample)
    {
      for(prop in fraction.upregulated)
      {
        for(tools in  AnalysisMethods)
        {
          tools2=select_tool((tools))
          tpr_temp=c()
          tfdr_temp=c()
          auc_temp=c()
          for(i in 1:rep)
          {
            if(fixedfold){
              fileName = paste(working.dir,simul.data,'_', test.cond,'_',nde,'DE_',s,'spc_fixedfold_upFrac_',prop,'_rep_',i,"_",tools2,'.rds',sep='')
            }else{
              fileName = paste(working.dir,simul.data,'_', test.cond,'_',nde,'DE_',s,'spc_upFrac_',prop,'_rep_',i,"_",tools2,'.rds',sep='')
            }


            result = try(readRDS(fileName), silent=T)
            if(class(result)=='try-error'){next}
            result = result@result.table
            if(nrow(result)==0){next}
            if(tools == "PoissonSeq"){rownames(result)=as.character(result$Genename)}


            ts = append(ts, setdiff(tools2, ts))
            mColor = select_color(tools)
            tc = append(tc, setdiff(mColor, tc))


            if(!is.null(result$FDR)){
              FDR=result$FDR
            }else if(!is.null(result$adjpvalue)){
              FDR=result$adjpvalue
            }

            GeneName = rownames(result)
            if(nde > 0){
              if(nde < nrow(result)){
                TrueGene = paste('g',1:nde,sep="")
                FalseGene = paste('g', (nde+1):(nrow(result)),sep="")
              }else if(nde == nrow(result)){
                TrueGene = paste('g',1:nrow(result),sep="")
                FalseGene = ''
              }else{
                stop("nde cannot exceed number of total genes")
              }
            }else if(nde==0){
              TrueGene = ''
              FalseGene = paste('g', (nde+1):(nrow(result)),sep="")
            }else{
              stop("nde cannot have values below 0")
            }

            indexTrue = which(GeneName%in%TrueGene)
            indexFalse = which(GeneName%in%FalseGene)
            tpr_temp = append(tpr_temp, length(which(FDR[indexTrue]<0.1))/length(FDR[indexTrue]))
            tfdr_temp = append(tfdr_temp, length(setdiff(which(FDR<0.1),indexTrue))/length(which(FDR<0.1)))
            label = rep(0, nrow(result))
            label[indexTrue] = 1
            pred = prediction(predictions = 1-FDR, labels = label)
            auc_temp = append(auc_temp, performance(pred, 'auc')@y.values[[1]][1])
            REPEAT = append(REPEAT,i)
            NDE = append(NDE, paste("pDE = ",round(nde*100/nvar,2),"%",sep=""))
            NSAMPLE = append(NSAMPLE,s)
            UPPROP = append(UPPROP, paste("upDE = ",round(prop*100,2),"%",sep=""))
            METHOD = append(METHOD, tools)
            COLOR = append(COLOR, mColor)
          }
          tpr = append(tpr, tpr_temp)
          tfdr = append(tfdr, tfdr_temp)
          auc = append(auc, auc_temp)
        }
      }
    }
  }
  res = data.frame(Methods = METHOD, nSample=NSAMPLE, Repeat=REPEAT , nDE=NDE, upDE = UPPROP, TPR = tpr, trueFDR = tfdr, AUC = auc, Color = COLOR)
  res$Color=factor(res$Color)
  res$nDE = paste(res$nDE, sep="")
  res2 = melt(res, measure.vars=c("AUC","TPR","trueFDR"))
  for(nde in NDE){
    nDE_Factor = append(nDE_Factor, setdiff(paste("pDE = ",round(nDE*100/nvar,2),"%",sep=""),nDE_Factor))
  }

  res2$nDE = factor(res$nDE, levels=nDE_Factor)
  default_order=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','BaySeq.qn','PoissonSeq','SAMseq')
  axis_order = intersect(default_order, AnalysisMethods)
  for(size in nsample)
  {
    for(up in fraction.upregulated)
    {
      up=paste("upDE = ",round(up*100,2),"%",sep="")
      sub.res.temp = res2[res2$nSample == size,]
      sub.res = sub.res.temp[sub.res.temp$upDE == up,]
      rowtypes.index<-which(sub.res$variable %in% rowType)
      sub.res = sub.res[rowtypes.index,]
      miss = which(is.na(sub.res$value))
      if(length(miss)>0){sub.res = sub.res[-miss,]}
      pd = position_dodge(width=0.0)
      gbase = ggplot(sub.res, aes(y=value, x=Methods, color=Methods))+geom_boxplot(position=pd,outlier.shape=NA)+
        facet_grid(variable~nDE, scales='free')+
        scale_x_discrete(limits = axis_order)+
        theme(axis.text.x=element_text(angle=90, hjust=1))+
        scale_colour_manual(name = "Methods",
                            labels = ts[order(ts)],
                            values = tc[order(ts)])

      gline = gbase
      if(fixedfold){
        tt = paste(simul.data,' / ',disp.Type,' / SS = ',size,' / ',up,' / fixedfold',sep="")
      }else{
        tt = paste(simul.data,' / ',disp.Type,' / SS = ',size,' / ',up,sep="")
      }

      tt=gsub(pattern = 'upDE', replacement = 'Bal', x = tt)
      print(gline+aes(x=Methods)+labs(x='Methods', y=up)+ggtitle(tt))
      figurename=gsub(pattern = " / ", replacement = "_", x = tt)
      figurename=gsub(pattern=" = ",replacement="_",x=figurename,fixed=T)
      figurename=gsub(pattern="%",replacement="percent",x=figurename,fixed=T)
      figurename=paste(figurename,".pdf",sep = "")
      ggsave(file=paste(figure.dir,"/",figurename,sep=""),width = 10,height = 8)
      dev.off()

    }
  }
}


#' make synthetic data False Positive count plot function
#' @param working.dir Input file location
#' @param figure.dir Figure save location
#' @param simul.data A character parameter indicating which given dataset analyzed dataset is based on. ‘KIRC’, ‘Bottomly’, ‘mBdK’ and ‘mKdB’ are available.
#' @param rep An integer specifying iterations DE analysis methods run for each condition.
#' @param nsample An integer vector indicating number of samples in each sample group.
#' @param disp.Type A vector indicating how is the dispersion parameter assumed to be for each sample group in the analzyed dataset. Possible values are 'same' and 'different'.
#' @param modes A character specifying a test condition used for analyzed dataset.
#' "D" for basic simulation (not adding outliers).
#' "R" for adding 5% of random outlier.
#' "OS" for adding outlier sample to each sample group.
#' "DL" for decreasing KIRC simulation dispersion 22.5 times (similar to SEQC data dispersion) to compare with SEQC data.
#' @param AnalysisMethods A character vector specifying DE methods used for the analysis.
#' (e.g. 'edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')
#' @export
fpc_performance_plot = function(working.dir, figure.dir, simul.data, rep, nsample, disp.Type, modes, AnalysisMethods){

  if(length(simul.data)!=1){stop('simul.data must have one element.')}
  if(length(disp.Type)!=1){stop('disp.Type must have one element.')}



  fpc = ts = tc = NSAMPLE = METHOD = COLOR = COND = REPEAT = NULL

  for(mode in modes)
  {
    for(s in nsample)
    {
      for(tools in  AnalysisMethods)
      {
        tools2=select_tool((tools))
        fpc_temp=c()
        for(i in 1:rep)
        {

          if(disp.Type=='same'){
            test.cond=mode
          }else if(disp.Type=='different'){
            test.cond=paste('DiffDisp_', mode, sep='')
          }

          fileName = paste(working.dir,simul.data,'_', test.cond,'_0DE_',s,'spc_rep_',i,"_",tools2,'.rds',sep='')
          result = try(readRDS(fileName), silent=T)
          if(class(result)=='try-error'){next}
          result = result@result.table
          if(nrow(result)==0){next}
          if(tools == "PoissonSeq"){rownames(result)=as.character(result$Genename)}



          ts = append(ts, setdiff(tools2, ts))
          mColor = select_color(tools)
          tc = append(tc, setdiff(mColor, tc))


          if(!is.null(result$FDR)){
            FDR=result$FDR
          }else if(!is.null(result$adjpvalue)){
            FDR=result$adjpvalue
          }

          GeneName = rownames(result)
          FalseGene = paste('g', 1:(nrow(result)),sep="")
          indexFalse = which(GeneName%in%FalseGene)
          fpc_temp = append(fpc_temp, (length(which(FDR[indexFalse]<0.1))))

          REPEAT = append(REPEAT,i)
          COND = append(COND, mode)
          NSAMPLE = append(NSAMPLE,s)
          METHOD = append(METHOD, tools)
          COLOR = append(COLOR, mColor)
        }
        fpc = append(fpc, fpc_temp)
      }
    }
  }
  res = data.frame(Methods = METHOD, nSample=NSAMPLE, Repeat=REPEAT , Condition=COND, FPC = fpc, Color = COLOR)
  res$Color=factor(res$Color)
  res2 = melt(res, measure.vars=c("FPC"))
  res2<-res2[,-which(names(res2)=='Condition')]
  res2$variable=res$Condition
  default_order=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','BaySeq.qn','PoissonSeq','SAMseq')
  axis_order = intersect(default_order, AnalysisMethods)

  sub.res = res2
  miss = which(is.na(sub.res$value))
  if(length(miss)>0){sub.res = sub.res[-miss,]}
  pd = position_dodge(width=0.0)
  gbase = ggplot(sub.res, aes(y=value, x=Methods, color=Methods))+geom_boxplot(position=pd,outlier.shape=NA)+
    facet_grid(variable~nSample, scales='free')+
    scale_x_discrete(limits = axis_order)+
    theme(axis.text.x=element_text(angle=90, hjust=1))+
    scale_colour_manual(name = "Methods",
                        labels = ts[order(ts)],
                        values = tc[order(ts)])

  gline = gbase
  tt = paste(simul.data,' / False Positive Counts / ',disp.Type, ' dispersion ' ,sep="")
  print(gline+aes(x=Methods)+labs(x='Methods', y='False Positive counts')+ggtitle(tt))
  figurename=gsub(pattern = " / ", replacement = "_", x = tt)
  figurename=paste(figurename,".pdf",sep = "")
  ggsave(file=paste(figure.dir,"/",figurename,sep=""),width = 10,height = 8)
  dev.off()
}







#' make realdata plot function
#' @param working.dir Input file location
#' @param figure.dir Figure save location
#' @param simul.data A character parameter indicating which given dataset analyzed dataset is based on. ‘KIRC’, ‘Bottomly’ and ‘SEQC’ are available.
#' @param rep An integer specifying iterations DE analysis methods run for each condition.
#' @param nsample An integer vector indicating number of samples in each sample group.
#' @param rowType A character vector indicating which results are shown in performance plot. Combination of DetectedDE and FPC. (e.g. c('DetectedDE','FP.count')) If simul.data is 'SEQC', only combination of 'AUC', 'TPR' and 'trueFDR' is available.
#' @param AnalysisMethods A character vector specifying DE methods used for the analysis.
#' (e.g. 'edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')
#' @export
performance_realdata_plot = function(working.dir, figure.dir, simul.data, rep, nsample, rowType, AnalysisMethods){

  timer_start = Sys.time()
  if(length(simul.data)!=1){stop('simul.data must have one element.')}

  if(simul.data=='SEQC'){
    spikeinloca <-system.file("extdata", "ERCC_Controls_Analysis.txt",package = "compareDEtools")
    spikein <- read.table(spikeinloca, stringsAsFactors=FALSE, header=T, sep="\t", row.names=2)

    tpr = tfdr = auc = ts = tc = NDE = NSAMPLE = METHOD = COLOR = REPEAT = nDE_Factor =NULL

    nde=48
    nvar=17961
    for(s in nsample)
    {
      for(tools in  AnalysisMethods)
      {
        tools2=select_tool((tools))
        tpr_temp=c()
        tfdr_temp=c()
        auc_temp=c()
        for(i in 1:rep)
        {

          fileName = paste(working.dir,simul.data,'_',s,'spc_rep_',i,'_',tools2,'.rds',sep="")


          result = try(readRDS(fileName), silent=T)
          if(class(result)=='try-error'){next}
          result = result@result.table
          if(nrow(result)==0){next}
          if(tools == "PoissonSeq"){rownames(result)=as.character(result$Genename)}


          ts = append(ts, setdiff(tools2, ts))
          mColor = select_color(tools)
          tc = append(tc, setdiff(mColor, tc))


          if(!is.null(result$FDR)){
            FDR=result$FDR
          }else if(!is.null(result$adjpvalue)){
            FDR=result$adjpvalue
          }

          GeneName = rownames(result)
          TrueGene = rownames(spikein[which(spikein$log2.Mix.1.Mix.2.!=0),])
          FalseGene = rownames(spikein[which(spikein$log2.Mix.1.Mix.2.==0),])
          indexTrue = which(GeneName%in%TrueGene)
          indexFalse = which(GeneName%in%FalseGene)
          FDR<-FDR[c(indexTrue,indexFalse)]
          label = rep(0, nrow(result))
          label[indexTrue] = 1
          label<-label[c(indexTrue,indexFalse)]
          pred = prediction(predictions = 1-FDR, labels = label)


          tpr_temp = append(tpr_temp, length(which(FDR[indexTrue]<0.1))/length(FDR[indexTrue]))
          tfdr_temp = append(tfdr_temp, length(setdiff(which(FDR<0.1),indexTrue))/length(which(FDR<0.1)))

          auc_temp = append(auc_temp, performance(pred, 'auc')@y.values[[1]][1])
          REPEAT = append(REPEAT,i)
          NDE = append(NDE, paste("pDE = ",round(nde*100/nvar,2),"%",sep=""))
          NSAMPLE = append(NSAMPLE,s)
          METHOD = append(METHOD, tools)
          COLOR = append(COLOR, mColor)
        }
        tpr = append(tpr, tpr_temp)
        tfdr = append(tfdr, tfdr_temp)
        auc = append(auc, auc_temp)
      }
    }

    res = data.frame(Methods = METHOD, nSample=NSAMPLE, Repeat=REPEAT , nDE=NDE , TPR = tpr, trueFDR = tfdr, AUC = auc, Color = COLOR)
    res$Color=factor(res$Color)
    res$nDE = paste(res$nDE, sep="")
    res2 = melt(res, measure.vars=c("AUC","TPR","trueFDR"))
    for(nde in NDE){
      nDE_Factor = append(nDE_Factor, setdiff(paste("pDE = ",round(nDE*100/nvar,2),"%",sep=""),nDE_Factor))
    }

    res2$nDE = factor(res$nDE, levels=nDE_Factor)
    default_order=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','BaySeq.qn','PoissonSeq','SAMseq')
    axis_order = intersect(default_order, AnalysisMethods)
    for(size in nsample)
    {
      sub.res = res2[res2$nSample == size,]
      rowtypes.index<-which(sub.res$variable %in% rowType)
      sub.res = sub.res[rowtypes.index,]
      miss = which(is.na(sub.res$value))
      if(length(miss)>0){sub.res = sub.res[-miss,]}
      pd = position_dodge(width=0.0)
      gbase = ggplot(sub.res, aes(y=value, x=Methods, color=Methods))+geom_boxplot(position=pd,outlier.shape=NA)+
        facet_grid(variable~nDE, scales='free')+
        scale_x_discrete(limits = axis_order)+
        theme(axis.text.x=element_text(angle=90, hjust=1))+
        scale_colour_manual(name = "Methods",
                            labels = ts[order(ts)],
                            values = tc[order(ts)])

      gline = gbase
      tt = paste(simul.data,'_Real',sep="")
      print(gline+aes(x=Methods)+labs(x='Methods')+ggtitle(tt))
      figurename=gsub(pattern = " / ", replacement = "_", x = tt)
      figurename=gsub(pattern = " (",replacement = '_',x = figurename, fixed=T)
      figurename=gsub(pattern = ") ",replacement = '',x = figurename,fixed=T)
      figurename=gsub(pattern=" = ",replacement="_",x=figurename,fixed=T)
      figurename=paste(figurename,".pdf",sep = "")
      ggsave(file=paste(figure.dir,"/",figurename,sep=""),width = 10,height = 8)
      dev.off()
    }
  }else{
    gene_num = fpc = ts = tc = REPEAT = NSAMPLE = METHOD= COLOR = nsample_Factor =NULL

    for(s in nsample)
    {
      for(tools in  AnalysisMethods)
      {
        tools2=select_tool((tools))
        gene_num_temp=c()
        fpc_temp=c()
        for(i in 1:rep)
        {
          fileName = paste(working.dir,simul.data,'_',s,'spc_rep_',i,'_',tools2,'.rds',sep="")
          result = try(readRDS(fileName), silent=T)
          if(class(result)=='try-error'){next}
          result = result@result.table
          if(nrow(result)==0){next}
          if(tools == "PoissonSeq"){rownames(result)=as.character(result$Genename)}


          ts = append(ts, setdiff(tools, ts))
          mColor = select_color(tools)
          tc = append(tc, setdiff(mColor, tc))



          if(!is.null(result$FDR)){
            FDR=result$FDR
          }else if(!is.null(result$adjpvalue)){
            FDR=result$adjpvalue
          }


          GeneName = rownames(result)

          gene_num_temp=append(gene_num_temp, length(which(FDR<0.1)))

          REPEAT = append(REPEAT,i)
          NSAMPLE = append(NSAMPLE,paste(s," Sample",sep=""))
          METHOD = append(METHOD, tools)
          COLOR = append(COLOR, mColor)
        }
        if(!(simul.data == 'Bottomly' && s > 5)){
          for(i in 1:rep)
          {
            fileName = paste(working.dir,simul.data,'_',s,'spc_rep_',i,'_fpc_',tools2,'.rds',sep="")
            result = try(readRDS(fileName), silent=T)
            if(class(result)=='try-error'){next}
            result = result@result.table
            if(nrow(result)==0){next}
            if(tools == "PoissonSeq"){rownames(result)=as.character(result$Genename)}


            if(!is.null(result$FDR)){
              FDR=result$FDR
            }else if(!is.null(result$adjpvalue)){
              FDR=result$adjpvalue
            }

            GeneName = rownames(result)
            FalseGene = paste('g', 1:(nrow(result)),sep="")
            indexFalse = which(GeneName%in%FalseGene)
            fpc_temp = append(fpc_temp, (length(which(FDR[indexFalse]<0.1))))
          }
        }
        gene_num = append(gene_num, gene_num_temp)
        fpc = append(fpc, fpc_temp)
      }
    }

    res_DE = data.frame(Methods = METHOD, nSample=NSAMPLE, Repeat=REPEAT, Color = COLOR, DetectedDE=gene_num)
    res_FP = data.frame(Methods = METHOD[1:length(fpc)], nSample=NSAMPLE[1:length(fpc)], Repeat=REPEAT[1:length(fpc)], FP.count = fpc, Color = COLOR[1:length(fpc)])
    res_DE = melt(res_DE, measure.vars=c("DetectedDE"))
    res_FP = melt(res_FP, measure.vars=c("FP.count"))
    res2 = rbind(res_DE, res_FP)
    res2$Color=factor(res2$Color)

    for(s in nsample){
      nsample_Factor = append(nsample_Factor, setdiff(paste(s," Sample",sep=""),nsample_Factor))
    }

    res2$nSample = factor(res2$nSample, levels=nsample_Factor)

    default_order=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','BaySeq.qn','PoissonSeq','SAMseq')
    axis_order = intersect(default_order, AnalysisMethods)


    sub.res=res2
    rowtypes.index<-which(sub.res$variable %in% rowType)
    sub.res = sub.res[rowtypes.index,]

    miss = which(is.na(sub.res$value))
    if(length(miss)>0){sub.res = sub.res[-miss,]}
    pd = position_dodge(width=0.0)
    gbase = ggplot(sub.res, aes(y=value, x=Methods, color=Methods))+geom_boxplot(position=pd,outlier.shape=NA)+
      facet_grid(variable~nSample, scales='free')+
      scale_x_discrete(limits = axis_order)+
      theme(axis.text.x=element_text(angle=90, hjust=1))+
      scale_colour_manual(name = "Methods",
                          labels = ts[order(ts)],
                          values = tc[order(ts)])


    gline = gbase
    plot(gline)
    tt = paste(simul.data,'_Real',sep="")
    print(gline+aes(x=Methods)+labs(x='Methods')+ggtitle(tt))
    figurename=gsub(pattern = " / ", replacement = "_", x = tt)
    figurename=gsub(pattern = " (",replacement = '_',x = figurename, fixed=T)
    figurename=gsub(pattern = ") ",replacement = '',x = figurename,fixed=T)
    figurename=gsub(pattern=" = ",replacement="_",x=figurename,fixed=T)
    figurename=paste(figurename,".pdf",sep = "")
    ggsave(file=paste(figure.dir,"/",figurename,sep=""),width = 20,height = 18)
    dev.off()

    {
      print(Sys.time() - timer_start)
    }
  }
}


#' run PCA
#'
#' Run PCA on two sample datasets, KIRC and Bottomly.
#' @param datatypes A character parameter indicating the object dataset to be analyzed with PCA.
#' @export
run_PCA = function(datatypes)
{
  if(datatypes=='KIRC'){
    data(kidney, package='SimSeq') # Load TCGA KIRC RNA-seq data. 144 samples (72 cancer and matched normal, respectively)
    count = kidney$counts # RNA-seq read count data
    index.cancer=(1:72)*2 # cancer sample index
    index.normal=index.cancer-1 # normal sample index
    count= count[,c(index.cancer, index.normal)] # Arrange samples for convinience
    mean.total = apply(count,1,mean)
    index.filter = which(mean.total < 10)
    count<-count[-index.filter,]

    DESeq.cds <- DESeq::newCountDataSet(countData = count, conditions=c(rep(1,72),rep(2,72))) # DESeq normalization
    dds <- estimateSizeFactors(DESeq.cds)
    mat <- counts(dds, normalized=TRUE)
    rownames(mat)<-paste('g',1:16621,sep='')
    colnames(mat)<-1:144
    require(graphics)
    model <- prcomp(t(mat), scale = TRUE)
    plot(model$x[,1:2], pch=19, col=c(rep('red',72),rep('blue',72)))
    legend("bottom", pch=19, col=c('red','blue'), legend=c('Cancer','Normal'))
  }else if(datatypes=='Bottomly'){
    load(system.file("extdata", "bottomly_eset.RData",package = "compareDEtools"))# Load Bottomly mouse RNA-seq data. 21 samples
    count<-exprs(bottomly.eset) # RNA-seq read count data
    strain<-pData(bottomly.eset)[,'strain']
    index.C=which(strain=='C57BL/6J') # C57BL/6J sample index
    index.D=which(strain=='DBA/2J') # DBA/2J sample index
    count= count[,c(index.C, index.D)] # Arrange samples for convinience
    mean.total = apply(count,1,mean)
    index.filter = which(mean.total < 10)
    count<-count[-index.filter,]
    DESeq.cds <- DESeq::newCountDataSet(countData = count, conditions=c(rep(1,10),rep(2,11))) # DESeq normalization
    dds <- estimateSizeFactors(DESeq.cds)
    mat <- counts(dds, normalized=TRUE)
    rownames(mat)<-paste('g',1:8550,sep='')
    colnames(mat)<-1:21

    require(graphics)
    plot(model$x[,1:2],col=c(rep('red',10),rep('blue',11)),pch=19)
    legend("bottom", pch=19, col=c('red','blue'), legend=c('C57BL/6J','DBA/2J'))
  }else
    stop("Given name is not one of two sample datasets.")
}


