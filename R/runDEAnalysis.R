#' Generation of Synthetic simulation files
#' @param working.dir Input file location
#' @param data.types A vector parameter indicating types of dataset (e.g. data.type = c(KIRC, Bottomly, mBdK and mKdB))
#' @param fixedfold A logical indicating whether simulation data is made from fixed fold to imitate SEQC counts data. If fixedfold is TRUE, fraction of upregulated genes is automatically fixed to 0.67.
#' @param rep.start An integer specifying start number of replication. Default is 1.
#' @param rep.end An integer specifying how many datasets will be generated from \code{rep.start} for each condition to run DE analysis methods.
#' @param nsample An integer vector indicating how many samples are in each sample group.
#' @param nvar An integer indicating the number of total gene in the synthetic data.
#' @param nDE A vector indicating number of generated DE genes in the synthetic data.
#' @param fraction.upregulated A vector specifying proportion of upregulated DE genes in the synthetic data. Default value is 0.5. (e.g. fraction.upregulated = c(0.5, 0.7 and 0.9)) This vector is available only if fixed fold is FALSE.
#' @param disp.Types A vector indicating how is the dispersion parameter assumed to be for each condition to make a synthetic data. Possible values are 'same' and 'different'.
#' @param modes A vector specifying test conditions we used for simulation data generation.
#' "D" for basic simulation (not adding outliers).
#' "R" for adding 5% of random outlier.
#' "OS"for adding outlier sample to each sample group.
#' "DL" for decreasing KIRC simulation dispersion 22.5 times (similar to SEQC data dispersion) to compare with SEQC data.
#' @param RO.prop An integer specifying random outlier proportion percentage that we generate dataset with.
#' @export
GenerateSyntheticSimulation<-function(working.dir, data.types, fixedfold=FALSE, rep.start=1, rep.end, nsample, nvar, nDE, fraction.upregulated = 0.5 , disp.Types, modes, RO.prop=5, random_sampling=FALSE){
  dataset.parameters<-generateDatasetParameter()
  for(simul.data in data.types){
    for(mode in modes){
      for(disp.Type in disp.Types){
        if(disp.Type=='same'){
          datahead=paste(simul.data,'_',mode, '_', sep='')
        }else if(disp.Type=='different'){
          datahead=paste(simul.data,'_','DiffDisp_',mode ,'_', sep='')
        }
        for(i in rep.start:rep.end){
          for(s in nsample){
            for(nde in nDE){
              if(nde==0){
                SyntheticDataSimulation(simul.data = simul.data, dataset = paste(working.dir,datahead,"0DE_",s,"spc_rep_",i,'.rds',sep=""),
                                        fixedfold = fixedfold, n.var = nvar, samples.per.cond = s, n.diffexp = 0,
                                        fraction.upregulated = 0, dispType=disp.Type, mode = mode, RO.prop=RO.prop, dataset.parameters = dataset.parameters, random_sampling=random_sampling)
              }else if(fixedfold){
                SyntheticDataSimulation(simul.data = simul.data, dataset = paste(working.dir,datahead, nde, "DE_",s,"spc_fixedfold_upFrac_0.67_rep_",i,'.rds',sep=""),
                                        fixedfold = fixedfold, n.var = nvar, samples.per.cond = s, n.diffexp = nde,
                                        fraction.upregulated = 0.67, dispType = disp.Type, mode = mode, RO.prop=RO.prop, dataset.parameters = dataset.parameters, random_sampling=random_sampling)
              }else{
                for(frac in fraction.upregulated){
                  SyntheticDataSimulation(simul.data=simul.data, dataset = paste(working.dir,datahead, nde, "DE_",s,"spc_upFrac_",frac,"_rep_",i,'.rds',sep=""),
                                          fixedfold=fixedfold, n.var=nvar, samples.per.cond=s, n.diffexp=nde,
                                          fraction.upregulated=frac, dispType=disp.Type, mode=mode, RO.prop=RO.prop, dataset.parameters=dataset.parameters, random_sampling=random_sampling)
                }
              }
            }
          }
        }
      }
    }
  }
}

#' Generation of Real simulation files
#' @param working.dir Input file location
#' @param fpc A logical indicating whether simulation data is generated with samples from single sample group(e.g. normal) to calculate False positive counts. Possible values are TRUE and FALSE.
#' @param data.types A vector indicating types of dataset (e.g. data.types = c(KIRC, Bottomly))
#' @param rep.start An integer specifying start number of replication. Default is 1.
#' @param rep.end An integer specifying how many datasets will be generated from \code{rep.start} for each condition to run DE analysis methods.
#' @param nsample A vector indicating number of samples in each sample group. Possible values can be 1~72 for KIRC, 1~10 for Bottomly and 1~5 for SEQC, maximum values decreased to 36 for KIRC and 5 for Bottomly if fpc is true.
#' @export
GenerateRealSimulation<-function(working.dir, fpc=FALSE, data.types, rep.start=1, rep.end, nsample){
  dataset.parameters<-generateDatasetParameter()
  for(simul.data in data.types){
    datahead=simul.data
    for(i in rep.start:rep.end){
      for(s in nsample){
        if(simul.data=='SEQC'){
          RealDataSimulation(simul.data=simul.data, dataset=paste(working.dir,datahead,'_',s,'spc_rep_',i,'.rds',sep=""), samples.per.cond=s, fpc=fpc, dataset.parameters=dataset.parameters)
        }else{
          if(fpc){
            RealDataSimulation(simul.data=simul.data, dataset=paste(working.dir,datahead,'_',s,'spc_rep_',i,'_fpc.rds',sep=""), samples.per.cond=s, fpc=fpc, dataset.parameters=dataset.parameters)
          }else{
            RealDataSimulation(simul.data=simul.data, dataset=paste(working.dir,datahead,'_',s,'spc_rep_',i,'.rds',sep=""), samples.per.cond=s, fpc=fpc, dataset.parameters=dataset.parameters)
          }
        }
      }
    }
  }
}

#' Run Analysismethods
#' @param working.dir Input file location
#' @param output.dir Result file location
#' @param real A logical indicating whether this analysis is for synthetic data or real data.
#' @param fpc A logical indicating whether simulation data is made from a single sample group (e.g. normal) to calculate False positive counts. Only used for real data analysis.
#' @param data.types A character vector indicating which given dataset our target dataset is based on. ‘KIRC’, ‘Bottomly’, ‘mBdK’ and ‘mKdB’ are available for synthetic datasets and ‘KIRC’, ‘Bottomly’ and ‘SEQC’ are available for real datasets.
#' @param fixedfold A logical indicating whether simulation data is made from fixed fold to imitate SEQC counts data. Only applied if this simulation is synthetic data simulation.
#' @param rep.start An integer specifying start number of replication. Default is 1.
#' @param rep.end An integer specifying how many datasets will be generated from \code{rep.start} for each condition to run DE analysis methods.
#' @param nsample An integer vector indicating number of samples in each sample group.
#' @param nDE An integer vector indicating number of generated DE genes in the synthetic data.
#' @param fraction.upregulated A numeric vector specifying proportions of upregulated DE genes among total DE genes in the generated dataset. Default value is 0.5. (e.g. fraction.upregulated = c(0.5, 0.7 and 0.9)) This vector is available only if fixed fold is FALSE.
#' @param disp.Types A vector indicating how is the dispersion parameter assumed to be for each sample group to make a synthetic data. Possible values are 'same' and 'different'.
#' @param modes A vector specifying test conditions we used for simulation data generation.
#' "D" for basic simulation (not adding outliers).
#' "R" for adding 5% of random outlier.
#' "OS"for adding outlier sample to each sample group.
#' "DL" for decreasing KIRC simulation dispersion 22.5 times (similar to SEQC data dispersion) to compare with SEQC data.
#' @param AnalysisMethods A character vector specifying DE methods used for the analysis.
#' (e.g. 'edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','BaySeq.qn,'PoissonSeq','SAMseq')
#' @param para A list parameter indicating the parameters to run each DE analysis methods. It contains lists corresponding each method and each list contain the parameters for each DE analysis methods. The analysis methods not in the para list will be run with default parameters. (e.g. para=list(ROTS=list(transformation=FALSE, normalize=FALSE)))
#' @export
runSimulationAnalysis<-function(working.dir, output.dir, real=FALSE, fpc=FALSE, data.types, fixedfold=FALSE, rep.start=1, rep.end, nsample, nDE, fraction.upregulated, disp.Types, modes, AnalysisMethods, para=list()){

  for(i in setdiff(AnalysisMethods, ls(para))){
    para[[i]]=list()
  }

  if(real){
    for(simul.data in data.types){
      for(i in rep.start:rep.end)
      {
        for(s in nsample)
        {
          if(fpc){
            obj = paste(working.dir,simul.data,'_',s,"spc_rep_",i,'_fpc.rds',sep="")
          }else{
            obj = paste(working.dir,simul.data,'_',s,"spc_rep_",i,'.rds',sep="")
          }
          dat = readRDS(obj)

          temp=0
          while(TRUE)
          {

            if(temp==1){break
            }else{
              test = try(simul_methods(obj,output.dir, AnalysisMethods, para=para), silent = T)
              if(class(test)=="try-error"){temp=0}else{temp=1}
            }
          }
        }
      }
    }
  }else{
    for(simul.data in data.types){
      for(mode in modes)
      {
        for(disp.Type in disp.Types){
          out=switch(paste(disp.Type,mode,sep=''), 'differentD'='DiffDisp_D_',
                     'differentR'='DiffDisp_R_',
                     'differentOS'='DiffDisp_OS_',
                     'differentDL'='DiffDisp_DL_',
                     'sameD'='D_',
                     'sameR'='R_',
                     'sameOS'='OS_',
                     'sameDL'='DL_'
          )
          datahead=paste(simul.data,'_',out ,sep='')
          for(i in rep.start:rep.end)
          {
            for(s in nsample)
            {
              if(fpc){
                obj = paste(working.dir,datahead,"0DE_",s,"spc_rep_",i,'.rds',sep="")
                dat = readRDS(obj)

                temp=0
                while(TRUE)
                {
                  if(temp==1){break
                  }else{
                    test = try(simul_methods(obj,output.dir, AnalysisMethods, para=para), silent = T)
                    if(class(test)=="try-error"){temp=0}else{temp=1}
                    print(test)
                  }
                }
              }else{
                for(nde in nDE)
                {
                  if(fixedfold){
                    obj = paste(working.dir,datahead,nde,"DE_",s,"spc_fixedfold_upFrac_0.67_rep_",i,'.rds',sep="")
                    dat = readRDS(obj)

                    temp=0
                    while(TRUE)
                    {
                      if(temp==1){break
                      }else{
                        test = try(simul_methods(obj,output.dir, AnalysisMethods, para=para), silent = T)
                        if(class(test)=="try-error"){temp=0}else{temp=1}
                        print(test)
                      }
                    }
                  }else{
                    for(frac in fraction.upregulated)
                    {
                      obj = paste(working.dir,datahead,nde,"DE_",s,"spc_upFrac_",frac,"_rep_",i,'.rds',sep="")
                      dat = readRDS(obj)

                      temp=0
                      while(TRUE)
                      {
                        if(temp==1){break
                        }else{
                          test = try(simul_methods(obj,output.dir, AnalysisMethods, para=para), silent = T)
                          if(class(test)=="try-error"){temp=0}else{temp=1}
                          print(test)
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}


#' Run Analysismethods
#' @param obj Name and place of input files.
#' @param output.dir Result file save location.
#' @param AnalysisMethods DE methods used for analysis. Input as character vectors
#' (e.g. 'edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','BaySeq.qn,'PoissonSeq','SAMseq')
#' @param para A list parameter indicating the parameters to run each DE analysis methods. It contains lists corresponding each method and each list contain the parameters for each DE analysis methods. The analysis methods not in the para list will be run with default parameters. (e.g. para=list(ROTS=list(transformation=FALSE, normalize=FALSE)))
#' @export
simul_methods=function(obj,output.dir, AnalysisMethods, para=list()){
  if('edgeR' %in% AnalysisMethods){
    edgeR.parameters=list(result.extent = "edgeR.exact", norm.method = 'TMM', trend.method = 'movingave', disp.type = 'tagwise')
    for(i in setdiff(ls(edgeR.parameters),ls(para$edgeR))){
      para$edgeR[[i]]=edgeR.parameters[[i]]
    }
    runDiffExp(data.file = obj, result.extent = para$edgeR$result.extent,
               Rmdfunction = "compareDEtools::edgeR.exact.createRmd", output.directory = output.dir,
               norm.method = para$edgeR$norm.method, trend.method = para$edgeR$trend.method, disp.type = para$edgeR$disp.type) # run edgeR exact test
  }
  if('edgeR.ql' %in% AnalysisMethods){
    edgeR.ql.parameters=list(result.extent = 'edgeR.GLMQL', norm.method = 'TMM', disp.estimation = 'CoxReid')
    for(i in setdiff(ls(edgeR.ql.parameters),ls(para$edgeR.ql))){
      para$edgeR.ql[[i]]=edgeR.ql.parameters[[i]]
    }
    runDiffExp(data.file = obj, result.extent = para$edgeR.ql$result.extent,
               Rmdfunction = 'compareDEtools::edgeR.GLMQL.createRmd', norm.method=para$edgeR.ql$norm.method,
               disp.estimation=para$edgeR.ql$disp.estimation, output.directory = output.dir) #run edgeR quasi-likelihood test
  }
  if('edgeR.rb' %in% AnalysisMethods){
    edgeR.rb.parameters=list(result.extent = 'edgeR.GLM.Robust', norm.method = 'TMM', disp.estimation = 'CoxReid')
    for(i in setdiff(ls(edgeR.rb.parameters),ls(para$edgeR.rb))){
      para$edgeR.rb[[i]]=edgeR.rb.parameters[[i]]
    }
    runDiffExp(data.file = obj, result.extent = para$edgeR.rb$result.extent,
               Rmdfunction = 'compareDEtools::edgeR.GLM.Robust.createRmd', norm.method=para$edgeR.rb$norm.method,
               disp.estimation=para$edgeR.rb$disp.estimation, output.directory = output.dir) #run edgeR GLM Robust
  }
  if('DESeq.pc' %in% AnalysisMethods){
    DESeq.pc.parameters=list(result.extent = "DESeq.pc", sharing.mode = "maximum", disp.method = "per-condition", fit.type = "local")
    for(i in setdiff(ls(DESeq.pc.parameters),ls(para$DESeq.pc))){
      para$DESeq.pc[[i]]=DESeq.pc.parameters[[i]]
    }
    runDiffExp(data.file = obj, result.extent = para$DESeq.pc$result.extent,
               Rmdfunction = "compareDEtools::DESeq.nbinom.createRmd", output.directory = output.dir,
               sharing.mode = para$DESeq.pc$sharing.mode, disp.method = para$DESeq.pc$disp.method,
               fit.type = para$DESeq.pc$fit.type)  # run DESeq.pc
  }
  if('DESeq2' %in% AnalysisMethods){
    DESeq2.parameters=list(result.extent = "DESeq2", fit.type = "parametric", test="Wald")
    for(i in setdiff(ls(DESeq2.parameters),ls(para$DESeq2))){
      para$DESeq2[[i]]=DESeq2.parameters[[i]]
    }
    runDiffExp(data.file=obj, result.extent = para$DESeq2$result.extent,
               Rmdfunction = "compareDEtools::DESeq2.createRmd", output.directory = output.dir,
               fit.type = para$DESeq2$fit.type, test=para$DESeq2$test) # run DESeq2
  }
  if('voom.tmm' %in% AnalysisMethods){
    voom.tmm.parameters=list(result.extent = "voom.limma", norm.method = "TMM")
    for(i in setdiff(ls(voom.tmm.parameters),ls(para$voom.tmm))){
      para$voom.tmm[[i]]=voom.tmm.parameters[[i]]
    }
    runDiffExp(data.file=obj, result.extent = para$voom.tmm$result.extent,
               Rmdfunction = "compareDEtools::voom.limma.createRmd", output.directory = output.dir,
               norm.method = para$voom.tmm$norm.method) # run voom limma (TMM normalization)
  }
  if('voom.qn' %in% AnalysisMethods){
    voom.qn.parameters=list(result.extent = "voom.qn.limma", norm.method = 'quantile')
    for(i in setdiff(ls(voom.qn.parameters),ls(para$voom.qn))){
      para$voom.qn[[i]]=voom.qn.parameters[[i]]
    }
    runDiffExp(data.file = obj, result.extent = para$voom.qn$result.extent,
               Rmdfunction = "compareDEtools::voom.qn.limma.createRmd", output.directory = output.dir,
               norm.method = para$voom.qn$norm.method) # run voom limma (quantile normalization)
  }
  if('voom.sw' %in% AnalysisMethods){
    voom.sw.parameters=list(result.extent = "voom.sw.limma", norm.method = 'TMM')
    for(i in setdiff(ls(voom.sw.parameters),ls(para$voom.sw))){
      para$voom.sw[[i]]=voom.sw.parameters[[i]]
    }
    runDiffExp(data.file=obj, result.extent = para$voom.sw$result.extent,
               Rmdfunction = 'compareDEtools::voom.sw.limma.createRmd', output.directory = output.dir,
               norm.method = para$voom.sw$norm.method) # run voom limma(quality weighted)
  }
  if('ROTS' %in% AnalysisMethods){
    ROTS.parameters=list(result.extent = "ROTS", transformation=TRUE, normalize=TRUE, B = 1000, K = NULL, log = FALSE)
    for(i in setdiff(ls(ROTS.parameters),ls(para$ROTS))){
      para$ROTS[[i]]=ROTS.parameters[[i]]
    }
    runDiffExp(data.file=obj, result.extent = para$ROTS$result.extent,
               Rmdfunction = "compareDEtools::ROTS.createRmd", output.directory = output.dir,
               normalize = para$ROTS$normalize, B = para$ROTS$B,
               K = para$ROTS$K, log = para$ROTS$log)  # run ROTS
  }
  if('BaySeq' %in% AnalysisMethods){
    BaySeq.parameters=list(result.extent = "BaySeq", norm.method = 'edgeR', equaldisp = TRUE)
    for(i in setdiff(ls(BaySeq.parameters),ls(para$BaySeq))){
      para$BaySeq[[i]]=BaySeq.parameters[[i]]
    }
    runDiffExp(data.file = obj, result.extent = para$BaySeq$result.extent,
               Rmdfunction = "compareDEtools::baySeq.createRmd", output.directory = output.dir,
               norm.method = para$BaySeq$norm.method, equaldisp = para$BaySeq$equaldisp) # run Bayseq
  }
  if('BaySeq.qn' %in% AnalysisMethods){
    BaySeq.qn.parameters=list(result.extent = "BaySeq.qn", norm.method = 'quantile', equaldisp = TRUE)
    for(i in setdiff(ls(BaySeq.qn.parameters),ls(para$BaySeq))){
      para$BaySeq.qn[[i]]=BaySeq.qn.parameters[[i]]
    }
    runDiffExp(data.file = obj, result.extent = para$BaySeq.qn$result.extent,
               Rmdfunction = "compareDEtools::baySeq.createRmd", output.directory = output.dir,
               norm.method = para$BaySeq.qn$norm.method, equaldisp = para$BaySeq.qn$equaldisp) # run Bayseq
  }
  if('PoissonSeq' %in% AnalysisMethods){
    PoissonSeq.parameters=list(result.extent = "PoissonSeq")
    for(i in setdiff(ls(PoissonSeq.parameters),ls(para$PoissonSeq))){
      para$PoissonSeq[[i]]=PoissonSeq.parameters[[i]]
    }
    runDiffExp(data.file = obj, result.extent = para$PoissonSeq$result.extent,
               Rmdfunction = "compareDEtools::PoissonSeq.createRmd", output.directory = output.dir) # run PoissionSeq
  }
  if('SAMseq' %in% AnalysisMethods){
    SAMseq.parameters=list(result.extent = "SAMseq")
    for(i in setdiff(ls(SAMseq.parameters),ls(para$SAMseq))){
      para$SAMseq[[i]]=SAMseq.parameters[[i]]
    }
    runDiffExp(data.file = obj, result.extent = para$SAMseq$result.extent,
               Rmdfunction = "compareDEtools::SAMseq.createRmd", output.directory = output.dir) # run SAMseq
  }
}
