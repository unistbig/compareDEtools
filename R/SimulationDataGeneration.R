#' Get random dispersion value based on the read count size (mean).
#' @param mean read count size
#' @param mean.condition read count vector of a certain sample condition
#' @param disp.condition dispersion vector of a certain sample condition
getDisp = function(mean, mean.condition, disp.condition)
{
  pool = disp.condition[which(mean.condition>(mean-20) & mean.condition<(mean+20))]
  if(length(pool)==0){value = disp.condition[which.min(abs(mean.condition-mean))]
  }else{value = sample(pool,1)}
  return(value)
}

#' Generate List containing estimated mean and dispersion parameters and filtered count from original count dataset.
#' @export
generateDatasetParameter = function(){
  #  Effect size 1.2/1.5~
  data(kidney, package='SimSeq') # Load TCGA KIRC RNA-seq data. 144 samples (72 cancer and matched normal, respectively)
  k_count = kidney$counts # RNA-seq read count data
  index.cancer=(1:72)*2 # cancer sample index
  index.normal=index.cancer-1 # normal sample index

  k_count= k_count[,c(index.cancer, index.normal)] # Arrange samples for convinience

  # Get count mean and dispersion using edgeR package

  # Mean and dispersion values are obtained separately from cancer and normal samples when different dispersion is assummed between two sample types.
  # Mean and dispersion values from normal samples
  dge.normal=DGEList(counts=k_count[,73:144], group = factor(rep(2,72)))
  dge.normal=calcNormFactors(dge.normal)
  dge.normal=estimateCommonDisp(dge.normal)
  dge.normal=estimateTagwiseDisp(dge.normal)
  disp.normal = dge.normal$tagwise.dispersion # Dispersion
  mean.normal = apply(k_count[,73:144],1,mean)


  # Mean and dispersion values from cancer samples
  dge.cancer=DGEList(counts=k_count[,1:72], group=factor(rep(1,72)))
  dge.cancer=calcNormFactors(dge.cancer)
  dge.cancer=estimateCommonDisp(dge.cancer)
  dge.cancer=estimateTagwiseDisp(dge.cancer)
  disp.cancer = dge.cancer$tagwise.dispersion
  mean.cancer = apply(k_count[,1:72],1,mean)


  # Gene filtering: Genes having small read count (<10) are filtered
  k_mean.total = apply(k_count,1,mean)
  k_index.filter = which(k_mean.total < 10)
  k_mean.total = k_mean.total[-k_index.filter]
  disp.normal = disp.normal[-k_index.filter]
  disp.cancer = disp.cancer[-k_index.filter]
  mean.normal = mean.normal[-k_index.filter]
  mean.cancer = mean.cancer[-k_index.filter]


  # Mean and dispersion values obtained using all samples when same dispersion is assummed between two sample types..
  k_dge.total = DGEList(counts = k_count, group = factor(c(rep(1,72),rep(2,72))))
  k_dge.total = calcNormFactors(k_dge.total)
  k_dge.total = estimateCommonDisp(k_dge.total)
  k_dge.total = estimateTagwiseDisp(k_dge.total)
  k_disp.total = k_dge.total$tagwise.dispersion
  k_disp.total = k_disp.total[-k_index.filter]



  ###############################################
  #############bottomly###########################
  ################################################
  load(system.file("extdata", "bottomly_eset.RData",package = "compareDEtools"))# Load Bottomly mouse RNA-seq data. 21 samples
  b_count<-exprs(bottomly.eset) # RNA-seq read count data
  strain<-pData(bottomly.eset)[,'strain']
  index.C=which(strain=='C57BL/6J') # C57BL/6J sample index
  index.D=which(strain=='DBA/2J') # DBA/2J sample index


  b_count= b_count[,c(index.C, index.D)] # Arrange samples for convinience


  # Get count mean and dispersion using edgeR package


  # Mean and dispersion values are obtained separately from C57BL/6J and DBA/2J samples when different dispersion is assummed between two sample types.
  # Mean and dispersion values from C57BL/6J samples
  dge.C=DGEList(counts=b_count[,1:10], group=factor(rep(1,10)))
  dge.C=calcNormFactors(dge.C)
  dge.C=estimateCommonDisp(dge.C)
  dge.C=estimateTagwiseDisp(dge.C)
  disp.C = dge.C$tagwise.dispersion
  mean.C = apply(b_count[,1:10],1,mean)


  # Mean and dispersion values from DBA/2J samples
  dge.D=DGEList(counts=b_count[,11:21], group = factor(rep(2,11)))
  dge.D=calcNormFactors(dge.D)
  dge.D=estimateCommonDisp(dge.D)
  dge.D=estimateTagwiseDisp(dge.D)
  disp.D = dge.D$tagwise.dispersion # Dispersion
  mean.D = apply(b_count[,11:21],1,mean)


  # Gene filtering: Genes having small read count (<10) are filtered
  b_mean.total = apply(b_count,1,mean)
  b_index.filter = which(b_mean.total < 10)
  b_mean.total = b_mean.total[-b_index.filter]
  disp.C = disp.C[-b_index.filter]
  disp.D = disp.D[-b_index.filter]
  mean.C = mean.C[-b_index.filter]
  mean.D = mean.D[-b_index.filter]

  # Mean and dispersion values obtained using all samples when same dispersion is assummed between two sample types..
  b_dge.total = DGEList(counts = b_count, group = factor(c(rep(1,10),rep(2,11))))
  b_dge.total = calcNormFactors(b_dge.total)
  b_dge.total = estimateCommonDisp(b_dge.total)
  b_dge.total = estimateTagwiseDisp(b_dge.total)
  b_disp.total = b_dge.total$tagwise.dispersion
  b_disp.total = b_disp.total[-b_index.filter]


  ###############################################
  #############   SEQC   ########################
  ###############################################
  SEQC<-system.file("extdata", "GSE49712_HTSeq.txt",package = "compareDEtools")
  s_count<-read.table(SEQC, header = T)
  s_count <- s_count[grep('no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique' ,
                      rownames(s_count), invert=TRUE),]
  s_mean.total = apply(s_count,1,mean)
  s_index.filter = which(s_mean.total < 10)


  dataset.parameters=list(k_count=k_count, disp.normal=disp.normal, mean.normal=mean.normal, disp.cancer=disp.cancer,
                          mean.cancer=mean.cancer, k_mean.total=k_mean.total, k_index.filter=k_index.filter, k_disp.total=k_disp.total,
                          b_count=b_count, disp.D=disp.D, mean.D=mean.D, disp.C=disp.C, mean.C=mean.C, b_mean.total=b_mean.total,
                          b_index.filter=b_index.filter, b_disp.total=b_disp.total, s_count=s_count, s_mean.total=s_mean.total, s_index.filter=s_index.filter)

  return(dataset.parameters)
}

#' Generate synthetic count data for analysis
#' @param simul.data Characters indicating which dataset will be used for simulation data generation
#' "KIRC" for KIRC dataset.
#' "Bottomly" for Bottomly dataset.
#' "mKdB" for hybrid dataset combining mean of KIRC and dispersion of Bottomly datatset.
#' "mBdK" for hybrid dataset combining mean of Bottomly and dispersion of KIRC datatset.
#' @param dataset Characters specifying the file name of simulation dataset.
#' @param fixedfold A logical indicating whether this dataset is generated by fold changes with fixed values or random folds following exponential distribution.
#' @param samples.per.cond An integer indicating number of samples for each sample group (e.g. 3).
#' @param n.var An integer indicating the number of total gene in the synthetic data.
#' @param n.diffexp An integer indicating number of generated DE genes in the synthetic data.
#' @param fraction.upregulated Proportion of upregulated DE genes in the synthetic data. (e.g. 0.5).
#' @param dispType A character parameter indicating how is the dispersion parameter assumed to be for each condition to make a synthetic data. Possible values are 'same' and 'different'.
#' @param mode Characters specifying test conditions used for simulation data generation.
#' "D" for basic simulation (not adding outliers).
#' "R" for adding 5% of random outlier.
#' "OS" for adding outlier sample to each sample group.
#' "DL" for decreasing KIRC simulation dispersion 22.5 times (similar to SEQC data dispersion) to compare with SEQC data.
#' @param dataset.parameters A list containing estimated mean and dispersion parameters and filtered count from original count dataset.
#' @export
SyntheticDataSimulation = function(simul.data, dataset, fixedfold=FALSE, samples.per.cond, n.var, n.diffexp, fraction.upregulated, dispType, mode, RO.prop=5, dataset.parameters)
{
  # Generate simulation data
  if(mode != 'D' && mode!='S' && mode!='R' && mode!='OS' && mode!='DL'){stop('mode must be "D" (DE variation test), "S" (single outlier test) or "R" (random outlier test) or "OS" (dispersion outlier sample test) or "DL" (dispersion lowered test)')}
  datasetName = dataset
  s = samples.per.cond
  k_count=dataset.parameters$k_count
  b_count=dataset.parameters$b_count
  disp.normal=dataset.parameters$disp.normal
  disp.cancer=dataset.parameters$disp.cancer
  disp.D=dataset.parameters$disp.D
  disp.C=dataset.parameters$disp.C
  mean.normal=dataset.parameters$mean.normal
  mean.cancer=dataset.parameters$mean.cancer
  mean.D=dataset.parameters$mean.D
  mean.C=dataset.parameters$mean.C
  k_mean.total=dataset.parameters$k_mean.total
  k_index.filter=dataset.parameters$k_index.filter
  k_disp.total=dataset.parameters$k_disp.total
  b_mean.total=dataset.parameters$b_mean.total
  b_index.filter=dataset.parameters$b_index.filter
  b_disp.total=dataset.parameters$b_disp.total




  if(simul.data=='KIRC'||simul.data=='mKdB'){
    random.index = sample(1:length(k_mean.total), size=n.var)
    sample.mean1 = k_mean.total[random.index]
    if(simul.data=='mKdB'){
      sub.random.index =sample(1:length(b_mean.total), size=n.var)
      sub.sample.mean1 = b_mean.total[sub.random.index]
      sub.sample.mean2 = sub.sample.mean1
    }
  }else if(simul.data=='Bottomly'||simul.data=='mBdK'){
    random.index = sample(1:length(b_mean.total), size=n.var)
    sample.mean1 = b_mean.total[random.index]
    if(simul.data=='mBdK'){
      sub.random.index =sample(1:length(k_mean.total), size=n.var)
      sub.sample.mean1 = k_mean.total[sub.random.index]
      sub.sample.mean2 = sub.sample.mean1
    }
  }
  sample.mean2 = sample.mean1



  if(n.diffexp!=0){
    if(fixedfold){
      if(simul.data!='KIRC'){
        stop('Simulation with fixed fold must be based on KIRC dataset.')
      }
      factor1 = 1.15
      factor2 = 1.3
      factor3 = 1.6

      sample.mean2[1:round(n.diffexp/3)] = sample.mean2[1:round(n.diffexp/3)]*factor1
      sample.mean2[round((n.diffexp/3)+1):round(2*n.diffexp/3)] = sample.mean2[round((n.diffexp/3)+1):round(2*n.diffexp/3)]*factor2
      sample.mean2[round((2*n.diffexp/3)+1):(n.diffexp)] = sample.mean2[round((2*n.diffexp/3)+1):(n.diffexp)]/factor3
    }else{
      upindex = 1:round(n.diffexp*fraction.upregulated)
      dnindex = round(n.diffexp*fraction.upregulated+1):n.diffexp
      if(s<=3){
        factor1 = 1.5+rexp(n = length(upindex), rate=1) #3 sample fold change 1.5~
        factor2 = 1.5+rexp(n = length(dnindex), rate=1)
      }else if(s<=5){
        factor1 = 1.3+rexp(n = length(upindex), rate=1) #5 sample fold change 1.3~
        factor2 = 1.3+rexp(n = length(dnindex), rate=1)
      }
      else{
        factor1 =  1.2+rexp(n = length(upindex), rate=1) #10 sample fold change 1.2~
        factor2 =  1.2+rexp(n = length(dnindex), rate=1)
      }
      sample.mean2[upindex] = sample.mean2[upindex]*factor1
      sample.mean2[dnindex] = sample.mean2[dnindex]/factor2
      if(simul.data=='mBdK'||simul.data=='mKdB'){
        sub.sample.mean2[upindex] = sub.sample.mean2[upindex]*factor1
        sub.sample.mean2[dnindex] = sub.sample.mean2[dnindex]/factor2
      }
    }
  }

  if(simul.data=='KIRC'||simul.data=='mBdK'){
    mean.condition1=mean.normal
    mean.condition2=mean.cancer
    disp.condition1=disp.normal
    disp.condition2=disp.cancer
    disp.total=k_disp.total
  }else if(simul.data=='Bottomly'||simul.data=='mKdB'){
    mean.condition1=mean.D
    mean.condition2=mean.C
    disp.condition1=disp.D
    disp.condition2=disp.C
    disp.total=b_disp.total
  }

  if(dispType == 'different')
  {
    if(simul.data=='KIRC'||simul.data=='Bottomly'){
      sample.disp1 = sapply(sample.mean1, FUN = getDisp, mean.condition = mean.condition1, disp.condition = disp.condition1, simplify = T, USE.NAMES = F)
      sample.disp2 = sapply(sample.mean2, FUN = getDisp, mean.condition = mean.condition2, disp.condition = disp.condition2, simplify = T, USE.NAMES = F)
    }else if(simul.data=='mKdB'||simul.data=='mBdK'){
      sample.disp1 = sapply(sub.sample.mean1[order(sub.sample.mean1)], FUN = getDisp, mean.condition = mean.condition1, disp.condition = disp.condition1, simplify = T, USE.NAMES = F)
      sample.disp1 <-sample.disp1[order(order(sample.mean1))]
      sample.disp2 = sapply(sub.sample.mean2[order(sub.sample.mean2)], FUN = getDisp, mean.condition = mean.condition2, disp.condition = disp.condition2, simplify = T, USE.NAMES = F)
      sample.disp2 <-sample.disp2[order(order(sample.mean2))]
    }
  }else if(dispType == 'same')
  {
    if(simul.data=='KIRC'||simul.data=='Bottomly'){
      sample.disp1 = disp.total[random.index]
    }else if(simul.data=='mKdB'||simul.data=='mBdK'){
      sample.disp1 = disp.total[sub.random.index[order(sub.sample.mean1)]][order(order(sample.mean1))]
    }
    sample.disp2 = sample.disp1
  }

  counts = matrix(nrow=n.var, ncol = 2*s)
  if(mode == "OS"){
    for(i in 1:n.var){
      counts[i,1:round(s/3)]=rnbinom(round(s/3), 1/(5*sample.disp2[i]),mu=sample.mean2[i])
      counts[i,(round(s/3)+1):s] = rnbinom((s-round(s/3)), 1/sample.disp2[i], mu = sample.mean2[i])
      counts[i,(s+1):(s+round(s/3))] = rnbinom(round(s/3), 1/(5*sample.disp1[i]), mu = sample.mean1[i])
      counts[i,(s+round(s/3)+1):(2*s)] = rnbinom((s-round(s/3)), 1/sample.disp1[i], mu = sample.mean1[i])
    }
  }else if(mode == "DL"){
    for(i in 1:n.var)
    {
      counts[i,1:s] = rnbinom(s, 22.5/sample.disp2[i], mu = sample.mean2[i])
      counts[i,(s+1):(2*s)] = rnbinom(s, 22.5/sample.disp1[i], mu = sample.mean1[i])
    }
  }else{
    for(i in 1:n.var)
    {
      counts[i,1:s] = rnbinom(s, 1/sample.disp2[i], mu = sample.mean2[i])
      counts[i,(s+1):(2*s)] = rnbinom(s, 1/sample.disp1[i], mu = sample.mean1[i])
    }
  }


  ### Random Outlier
  if(mode == "R")
  {
    RO = matrix(runif(n.var*2*s , min = 0, max = 100), nrow = n.var, ncol = 2*s)
    index.outlier = which(RO<RO.prop)
    counts[index.outlier] = counts[index.outlier]*runif(n = length(index.outlier), min=5, max=10)
    counts = round(counts)
  }


  rownames(counts) = paste("g",1:n.var,sep="")

  sample.annot = data.frame(condition = c(rep(1,s),rep(2,s)))
  colnames(counts) = rownames(sample.annot)
  info.parameters = list(dataset = datasetName, uID = datasetName)

  # save as compdata
  cpd = compData(count.matrix = counts,
                 sample.annotations = sample.annot,
                 info.parameters = info.parameters)
  saveRDS(cpd, datasetName)
}


#' Generate real data for analysis
#' @param simul.data Characters indicating which dataset will be used for simulation data generation
#' "KIRC" for KIRC dataset.
#' "Bottomly" for Bottomly dataset.
#' "SEQC" for SEQC ERCC Spike-In dataset.
#' @param dataset A character parameter specifying the file name of simulation dataset.
#' @param samples.per.cond An integer indicating number of samples for each sample group (e.g. 3).
#' @param fpc A logical indicating whether simulation data is made from single sample group(e.g. normal) to calculate False positive counts. Default is FALSE.
#' @param dataset.parameters A list containing estimated mean and dispersion parameters and filtered count from original count dataset.
#' @export
RealDataSimulation = function(simul.data ,dataset, samples.per.cond, fpc=FALSE, dataset.parameters)
{
  # Generate simulation data
  datasetName = dataset
  s=samples.per.cond

  k_count=dataset.parameters$k_count
  b_count=dataset.parameters$b_count
  k_index.filter=dataset.parameters$k_index.filter
  b_index.filter=dataset.parameters$b_index.filter
  s_count=dataset.parameters$s_count
  s_index.filter=dataset.parameters$s_index.filter


  if(fpc){
    if(simul.data=='KIRC'){
      random.index_row = sample(1:72, size=2*s)
      counts<-k_count[-k_index.filter,random.index_row]
    }else if(simul.data=='Bottomly'){
      random.index_row = sample(1:10, size=2*s)
      counts<-b_count[-b_index.filter,random.index_row]
    }else if(simul.data=='SEQC'){
      stop('SEQC is not used for FP count calculation')
    }
  }else{
    if(simul.data=='KIRC'){
      random.index_1 = sample(1:72, size=s)
      random.index_2 = sample(73:144, size=s)
      count<-k_count[-k_index.filter,]
      random.index_row<-append(random.index_1,random.index_2)
      counts<-count[,random.index_row]
    }else if(simul.data=='Bottomly'){
      random.index_1 = sample(1:10, size=s)
      random.index_2 = sample(11:21, size=s)
      count<-b_count[-b_index.filter,]
      random.index_row<-append(random.index_1,random.index_2)
      counts<-count[,random.index_row]
    }else if(simul.data=='SEQC'){
      count<-s_count[-s_index.filter,]
      random.index_1 = sample(1:5, size=s)
      random.index_2 = sample(6:10, size=s)
      random.index_row<-append(random.index_1,random.index_2)
      counts<-count[,random.index_row]
    }
  }

  if(simul.data!='SEQC'){
    rownames(counts) = paste("g",1:nrow(counts),sep="")
  }


  sample.annot = data.frame(condition = c(rep(1,s),rep(2,s)))
  colnames(counts) = rownames(sample.annot)
  info.parameters = list(dataset = datasetName, uID = datasetName)

  # save as compdata
  cpd = compData(count.matrix = counts,
                 sample.annotations = sample.annot,
                 info.parameters = info.parameters)
  saveRDS(cpd, datasetName)
}
