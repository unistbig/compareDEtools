#install required packages
install.packages('devtools')
library(devtools)

BiocManager::install(c('baySeq','Biobase','compcodeR','DESeq','DESeq2','edgeR,','impute','limma','ROTS'))
install.packages('gplots','gtools','ggplot2','PoissonSeq','reshape','RColorBrewer','ROCR','samr','SimSeq','statmod','XML')

#install compareDEtools and load package
install_github('unistbig/compareDEtools')
library(compareDEtools)


#Fig1
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','BaySeq.qn','PoissonSeq','SAMseq') #Set DE methods to analyze with
dataset.dir='~/test/example/Fig1/Dataset/' #Set dataset save directory
analysis.dir='~/test/example/Fig1/Analysis/' #Set DE analysis results save directory
figure.dir='~/test/example/Fig1/Figure/' #Set result figure save directory

#A
GenerateRealSimulation(working.dir=dataset.dir, data.types='SEQC', rep.end=1, nsample=5) #Generate SEQC real dataset
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=TRUE, data.types='SEQC', rep.end=1, nsample=c(5), AnalysisMethods = AnalysisMethods, para=list()) #Run DE analysis for preset methods
performance_realdata_plot(working.dir=analysis.dir,figure.dir=figure.dir,simul.data='SEQC', rep.end=1, nsample=c(5), AnalysisMethods=AnalysisMethods, rowType = 'AUC') #Draw real data performance plot for AUC

#B
GenerateSyntheticSimulation(working.dir=dataset.dir, data.types='KIRC', rep.end=10, nsample=c(3,5), nvar=10000, nDE=27, fraction.upregulated = 0.67, disp.Types = 'same', modes='D') #Generate KIRC synthetic data with high proportion of DE genes
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types='KIRC', rep.end=10, nsample=c(3,5), nDE=27, fraction.upregulated=0.67, disp.Types='same', modes='D', AnalysisMethods = AnalysisMethods, para=list()) #Run DE analysis for preset methods
performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,fixedfold=F,simul.data='KIRC', rep.end=10, nsample=c(3,5), nvar=10000, nDE=27, fraction.upregulated = 0.67, disp.Type = 'same', mode='D', AnalysisMethods=AnalysisMethods, rowType = c('AUC')) #Draw syntehtic data performance plot for AUC

#C
GenerateSyntheticSimulation(working.dir=dataset.dir, data.types='KIRC', fixedfold = T, rep.end=10, nsample=c(3,5), nvar=10000, nDE=27, disp.Types = 'same', modes='DL') #Generate KIRC synthetic data with fixed foldchange and lowered dispersion to compare with SEQC dataset
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types='KIRC', fixedfold = T, rep.end=10, nsample=c(3,5), nDE=27, disp.Types='same', modes='DL', AnalysisMethods = AnalysisMethods, para=list()) #Run DE analysis for preset methods
performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,fixedfold=T,simul.data='KIRC', rep.end=10, nsample=c(3,5), nvar=10000, nDE=27, fraction.upregulated = 0.67, disp.Type = 'same', mode='DL', AnalysisMethods=AnalysisMethods, rowType = c('AUC')) #Draw synthetic data performance plot for AUC


#Fig2
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') #Set DE methods to analyze with
dataset.dir='~/test/example/Fig2/Dataset/' #Set dataset save directory
analysis.dir='~/test/example/Fig2/Analysis/' #Set DE analysis results save directory
figure.dir='~/test/example/Fig2/Figure/' #Set result figure save directory

GenerateSyntheticSimulation(working.dir=dataset.dir, data.types='KIRC', rep.end=10, nsample=c(3,10), nvar=10000, nDE=c(500,1000,3000,6000), fraction.upregulated = 0.5, disp.Types = 'same', modes=c('D','R','OS')) #Generate KIRC synthetic dataset
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types='KIRC', rep.end=10, nsample=c(3,10), nDE=c(500,1000,3000,6000), fraction.upregulated=0.5, disp.Types='same', modes=c('D','R','OS'), AnalysisMethods = AnalysisMethods, para=list()) #Run DE analysis for preset methods

for(mode in c('D','R','OS')){
  performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,fixedfold=F,simul.data='KIRC', rep.end=10, nsample=c(3,10), nvar=10000, nDE=c(500,1000,3000,6000), fraction.upregulated = 0.5, disp.Type = 'same', mode=mode, AnalysisMethods=AnalysisMethods, rowType = c('AUC','TPR','trueFDR')) #Draw synthetic data performance plot
}

#Fig3
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') #Set DE methods to analyze with
dataset.dir='~/test/example/Fig3/Dataset/' #Set dataset save directory
analysis.dir='~/test/example/Fig3/Analysis/' #Set DE analysis results save directory
figure.dir='~/test/example/Fig3/Figure/' #Set result figure save directory

GenerateSyntheticSimulation(working.dir=dataset.dir, data.types='KIRC', rep.end=10, nsample=c(3,10), nvar=10000, nDE=0, fraction.upregulated = 0.5, disp.Types = 'same', modes=c('D','R','OS')) #Generate KIRC synthetic dataset without DE genes to calculate false positive counts
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types='KIRC', rep.end=10, nsample=c(3,10), fpc=T, nDE=c(0), disp.Types='same', modes=c('D','R','OS'), AnalysisMethods = AnalysisMethods, para=list()) #Run DE analysis for preset methods

fpc_performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,simul.data='KIRC', rep.end=10, nsample=c(3,10), disp.Type = 'same', modes=c('D','R','OS'), AnalysisMethods=AnalysisMethods) #Draw synthetic data false positive counts performance plot


#Fig4
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') #Set DE methods to analyze with
dataset.dir='~/test/example/Fig4/Dataset/' #Set dataset save directory
analysis.dir='~/test/example/Fig4/Analysis/' #Set DE analysis results save directory
figure.dir='~/test/example/Fig4/Figure/' #Set result figure save directory

#A-B
GenerateRealSimulation(working.dir=dataset.dir, data.types='KIRC', rep.end=10, nsample=c(3,5,10,20)) #Generate KIRC real dataset
GenerateRealSimulation(working.dir=dataset.dir, data.types='KIRC', fpc=TRUE, rep.end=10, nsample=c(3,5,10,20)) #Generate KIRC real dataset with counts from single sample condition to calculate false positive counts
GenerateRealSimulation(working.dir=dataset.dir, data.types='Bottomly', rep.end=10, nsample=c(3,5,10)) #Generate Bottomly real dataset without DE genes
GenerateRealSimulation(working.dir=dataset.dir, data.types='Bottomly', fpc=TRUE, rep.end=10, nsample=c(3,5)) #Generate Bottomly real dataset with counts from single sample condition to calculate false positive counts

runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=TRUE, data.types='KIRC', rep.end=10, nsample=c(3,5,10,20), AnalysisMethods = AnalysisMethods, para=list()) #Run DE analysis for preset methods
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=TRUE, fpc=TRUE, data.types='KIRC', rep.end=10, nsample=c(3,5,10,20), AnalysisMethods = AnalysisMethods, para=list())
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=TRUE, data.types='Bottomly', rep.end=10, nsample=c(3,5,10), AnalysisMethods = AnalysisMethods, para=list())
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=TRUE, fpc=TRUE, data.types='Bottomly', rep.end=10, nsample=c(3,5), AnalysisMethods = AnalysisMethods, para=list())


performance_realdata_plot(working.dir=analysis.dir,figure.dir=figure.dir,simul.data='KIRC', rep.end=10, nsample=c(3,5,10,20), AnalysisMethods=AnalysisMethods, rowType = c("DetectedDE","FP.count")) #Draw KIRC real data performance plot
performance_realdata_plot(working.dir=analysis.dir,figure.dir=figure.dir,simul.data='Bottomly', rep.end=10, nsample=c(3,5,10), AnalysisMethods=AnalysisMethods, rowType = c("DetectedDE","FP.count")) #Draw Bottomly real data performance plot
#C-D

correlation_heatmap(working.dir=analysis.dir, figure.dir=figure.dir,simul.data='KIRC', nsample=5, topgenes=5000, AnalysisMethods=AnalysisMethods, rep.end=10) #Draw correlation heatmap for DE methods run with KIRC real data analysis
correlation_heatmap(working.dir=analysis.dir, figure.dir=figure.dir,simul.data='Bottomly', nsample=5, topgenes=500, AnalysisMethods=AnalysisMethods, rep.end=10) #Draw correlation heatmap for DE methods run with Bottomly real data analysis


#FigS1
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') #Set DE methods to analyze with
dataset.dir='~/test/example/FigS1/Dataset/' #Set dataset save directory
analysis.dir='~/test/example/FigS1/Analysis/' #Set DE analysis results save directory
figure.dir='~/test/example/FigS1/Figure/' #Set result figure save directory

GenerateSyntheticSimulation(working.dir=dataset.dir, data.types='KIRC', rep.end=10, nsample=c(3,10), nvar=10000, nDE=c(500,1000,3000,6000), fraction.upregulated = 0.5, disp.Types = 'different', modes=c('D','R')) #Generate KIRC synthetic dataset, different dispersions to each sample condition are assumed to generate dataset.
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types='KIRC', rep.end=10, nsample=c(3,10), nDE=c(500,1000,3000,6000), fraction.upregulated=0.5, disp.Types='different', modes=c('D','R'), AnalysisMethods = AnalysisMethods, para=list()) #Run DE analysis for preset methods

for(mode in c('D','R')){
  performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,fixedfold=F,simul.data='KIRC', rep.end=10, nsample=c(3,10), nvar=10000, nDE=c(500,1000,3000,6000), fraction.upregulated = 0.5, disp.Type = 'different', mode=mode, AnalysisMethods=AnalysisMethods, rowType = c('AUC','TPR','trueFDR')) #Draw synthetic data performance plot
}

#FigS2
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') #Set DE methods to analyze with
dataset.dir='~/test/example/FigS2/Dataset/' #Set dataset save directory
analysis.dir='~/test/example/FigS2/Analysis/' #Set DE analysis results save directory
figure.dir='~/test/example/FigS2/Figure/' #Set result figure save directory

GenerateSyntheticSimulation(working.dir=dataset.dir, data.types='KIRC', rep.end=10, nsample=c(3,10), nvar=10000, nDE=c(500,1000,3000,6000), fraction.upregulated = c(0.7, 0.9), disp.Types = 'same', modes=c('D')) #Generate KIRC synthetic dataset with high proportion of upregulated DE genes.
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types='KIRC', rep.end=10, nsample=c(3,10), nDE=c(500,1000,3000,6000), fraction.upregulated=c(0.7, 0.9), disp.Types='same', modes=c('D'), AnalysisMethods = AnalysisMethods, para=list(ROTS=list(transformation=FALSE, normalize=FALSE))) #Run DE analysis for preset methods. ROTS parameter is set to unnormalized, without voom transformation.


performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,fixedfold=F,simul.data='KIRC', rep.end=10, nsample=c(3,10), nvar=10000, nDE=c(500,1000,3000,6000), fraction.upregulated = c(0.7, 0.9), disp.Type = 'same', mode='D', AnalysisMethods=AnalysisMethods, rowType = c('AUC','TPR','trueFDR')) #Draw synthetic data performance plot

#FigS3
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') #Set DE methods to analyze with
dataset.dir='~/test/example/FigS3/Dataset/' #Set dataset save directory
analysis.dir='~/test/example/FigS3/Analysis/' #Set DE analysis results save directory
figure.dir='~/test/example/FigS3/Figure/' #Set result figure save directory

GenerateSyntheticSimulation(working.dir=dataset.dir, data.types='Bottomly', rep.end=10, nsample=c(3,10), nvar=5000, nDE=c(250,500,1500,3000), fraction.upregulated = 0.5, disp.Types = 'same', modes=c('D','R','OS')) #Generate Bottomly synthetic dataset
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types='Bottomly', rep.end=10, nsample=c(3,10), nDE=c(250,500,1500,3000), fraction.upregulated=0.5, disp.Types='same', modes=c('D','R','OS'), AnalysisMethods = AnalysisMethods, para=list()) #Run DE analysis for preset methods

for(mode in c('D','R','OS')){
  performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,fixedfold=F,simul.data='KIRC', rep.end=10, nsample=c(3,10), nvar=10000, nDE=c(500,1000,3000,6000), fraction.upregulated = 0.5, disp.Type = 'same', mode=mode, AnalysisMethods=AnalysisMethods, rowType = c('AUC','TPR','trueFDR')) #Draw synthetic data performance plot
}

#FigS4
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') #Set DE methods to analyze with
dataset.dir='~/test/example/FigS4/Dataset/' #Set dataset save directory
analysis.dir='~/test/example/FigS4/Analysis/' #Set DE analysis results save directory
figure.dir='~/test/example/FigS4/Figure/' #Set result figure save directory

GenerateSyntheticSimulation(working.dir=dataset.dir, data.types='Bottomly', rep.end=10, nsample=c(3,10), nvar=5000, nDE=c(250,500,1500,3000), fraction.upregulated = 0.5, disp.Types = 'different', modes=c('D','R')) #Generate Bottomly synthetic dataset, different dispersions to each sample condition are assumed to generate dataset.
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types='Bottomly', rep.end=10, nsample=c(3,10), nDE=c(250,500,1500,3000), fraction.upregulated=0.5, disp.Types='different', modes=c('D','R'), AnalysisMethods = AnalysisMethods, para=list()) #Run DE analysis for preset methods

for(mode in c('D','R')){
  performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,fixedfold=F,simul.data='Bottomly', rep.end=10, nsample=c(3,10), nvar=5000, nDE=c(250,500,1500,3000), fraction.upregulated = 0.5, disp.Type = 'different', mode=mode, AnalysisMethods=AnalysisMethods, rowType = c('AUC','TPR','trueFDR')) #Draw synthetic data performance plot
}


#FigS5
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') #Set DE methods to analyze with
dataset.dir='~/test/example/FigS5/Dataset/' #Set dataset save directory
analysis.dir='~/test/example/FigS5/Analysis/' #Set DE analysis results save directory
figure.dir='~/test/example/FigS5/Figure/' #Set result figure save directory

GenerateSyntheticSimulation(working.dir=dataset.dir, data.types='Bottomly', rep.end=10, nsample=c(3,10), nvar=5000, nDE=c(250,500,1500,3000), fraction.upregulated = c(0.7, 0.9), disp.Types = 'same', modes=c('D')) #Generate Bottomly synthetic dataset with high proportion of upregulated DE genes.
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types='Bottomly', rep.end=10, nsample=c(3,10), nDE=c(250,500,1500,3000), fraction.upregulated=c(0.7, 0.9), disp.Types='same', modes=c('D'), AnalysisMethods = AnalysisMethods, para=list(ROTS=list(transformation=FALSE, normalize=FALSE))) #Run DE analysis for preset methods. ROTS parameter is set to unnormalized, without voom transformation.


performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,fixedfold=F,simul.data='Bottomly', rep.end=10, nsample=c(3,10), nvar=5000, nDE=c(250,500,1500,3000), fraction.upregulated = c(0.7, 0.9), disp.Type = 'same', mode='D', AnalysisMethods=AnalysisMethods, rowType = c('AUC','TPR','trueFDR')) #Draw synthetic data performance plot

#FigS6
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') #Set DE methods to analyze with
dataset.dir='~/test/example/FigS6/Dataset/' #Set dataset save directory
analysis.dir='~/test/example/FigS6/Analysis/' #Set DE analysis results save directory
figure.dir='~/test/example/FigS6/Figure/' #Set result figure save directory

GenerateSyntheticSimulation(working.dir=dataset.dir, data.types='Bottomly', rep.end=10, nsample=c(3,10), nvar=5000, nDE=0, fraction.upregulated = 0.5, disp.Types = 'same', modes=c('D','R','OS')) #Generate Bottomly synthetic dataset without DE genes to calculate false positive counts
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types='Bottomly', rep.end=10, nsample=c(3,10), fpc=T, nDE=0, disp.Types='same', modes=c('D','R','OS'), AnalysisMethods = AnalysisMethods, para=list()) #Run DE analysis for preset methods

fpc_performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,simul.data='Bottomly', rep.end=10, nsample=c(3,10), disp.Type = 'same', modes=c('D','R','OS'), AnalysisMethods=AnalysisMethods) #Draw synthetic data false positive counts performance plot


#FigS7
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') #Set DE methods to analyze with
dataset.dir='~/test/example/FigS7/Dataset/' #Set dataset save directory
analysis.dir='~/test/example/FigS7/Analysis/' #Set DE analysis results save directory
figure.dir='~/test/example/FigS7/Figure/' #Set result figure save directory

GenerateSyntheticSimulation(working.dir=dataset.dir, data.types='mBdK', rep.end=10, nsample=c(3,10), nvar=5000, nDE=c(250,500,1500,3000), fraction.upregulated = 0.5, disp.Types = 'same', modes=c('D','R')) #Generate hybrid dataset with mean counts from Bottomly dataset and dispersion from KIRC dataset
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types='mBdK', rep.end=10, nsample=c(3,10), nDE=c(250,500,1500,3000), fraction.upregulated=0.5, disp.Types='same', modes=c('D','R'), AnalysisMethods = AnalysisMethods, para=list()) #Run DE analysis for preset methods

for(mode in c('D','R')){
  performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,fixedfold=F,simul.data='mBdK', rep.end=10, nsample=c(3,10), nvar=10000, nDE=c(500,1000,3000,6000), fraction.upregulated = 0.5, disp.Type = 'same', mode=mode, AnalysisMethods=AnalysisMethods, rowType = c('AUC','TPR','trueFDR')) #Draw synthetic data performance plot
}


#FigS8
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') #Set DE methods to analyze with
dataset.dir='~/test/example/FigS8/Dataset/' #Set dataset save directory
analysis.dir='~/test/example/FigS8/Analysis/' #Set DE analysis results save directory
figure.dir='~/test/example/FigS8/Figure/' #Set result figure save directory

GenerateSyntheticSimulation(working.dir=dataset.dir, data.types='mKdB', rep.end=10, nsample=c(3,10), nvar=5000, nDE=c(250,500,1500,3000), fraction.upregulated = 0.5, disp.Types = 'same', modes=c('D','R')) #Generate hybrid dataset with mean counts from KIRC dataset and dispersion from Bottomly dataset
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types='mKdB', rep.end=10, nsample=c(3,10), nDE=c(250,500,1500,3000), fraction.upregulated=0.5, disp.Types='same', modes=c('D','R'), AnalysisMethods = AnalysisMethods, para=list()) #Run DE analysis for preset methods

for(mode in c('D','R')){
  performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,fixedfold=F,simul.data='mKdB', rep.end=10, nsample=c(3,10), nvar=10000, nDE=c(500,1000,3000,6000), fraction.upregulated = 0.5, disp.Type = 'same', mode=mode, AnalysisMethods=AnalysisMethods, rowType = c('AUC','TPR','trueFDR')) #Draw synthetic data performance plot
}

#FigS9
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') #Set DE methods to analyze with
analysis.dir='~/test/example/Fig4/Analysis/' #Set DE analysis results save directory
figure.dir='~/test/example/FigS9/Figure/' #Set result figure save directory


correlation_heatmap(working.dir=analysis.dir, figure.dir=figure.dir,simul.data='KIRC', nsample=5, topgenes=3000, AnalysisMethods=AnalysisMethods, rep.end=10) #Draw correlation heatmap for DE methods run with KIRC real data analysis
correlation_heatmap(working.dir=analysis.dir, figure.dir=figure.dir,simul.data='Bottomly', nsample=5, topgenes=300, AnalysisMethods=AnalysisMethods, rep.end=10) #Draw correlation heatmap for DE methods run with Bottomly real data analysis


#FigS10
run_PCA(datatypes = 'KIRC') #run PCA for KIRC reference dataset
run_PCA(datatypes = 'Bottomly') #run PCA for Bottomly reference dataset


