# compareDEtools


compareDEtools is a R package for comparing multiple DE methods with different conditions and datasets.

## Installation

Open R program and type following commands in R console.

```R
install.packages('devtools')
library(devtools)
install_github('unistbig/compareDEtools')
library(compareDEtools)
```

## Example Run


Following codes are to generate synthetic data, run DE analysis tools and show results in boxplot
```R
dataset.dir='~/test/Dataset/' #user defined directory
analysis.dir='~/test/Analysis/' #user defined directory
figure.dir='~/test/Fig/' #user defined directory

#Assign methods for analysis
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')

#Generate dataset for DE analysis
GenerateSyntheticSimulation(working.dir=dataset.dir, data.types='KIRC', rep=10, nsample=c(10), nvar=1000, nDE=c(50), fraction.upregulated = 0.5, disp.Types = 'same', modes=c('D'))

#Run DE methods with generated datasets
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types='KIRC', rep=10, nsample=c(10), nDE=c(50), fraction.upregulated=0.5, disp.Types='same', modes=c('D'), AnalysisMethods = AnalysisMethods, para=list())

#Print boxplots for comparing DE methods performance
performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,fixedfold=F,simul.data='KIRC', rep=10, nsample=c(10), nvar=1000, nDE=c(50), fraction.upregulated = 0.5, disp.Type = 'same', mode='D', AnalysisMethods=AnalysisMethods, rowType = c('AUC','TPR','trueFDR'))

```



## User's Manual


User's manual will be available [here](https://github.com/unistbig/compareDEtools/blob/master/compareDEtools_manual_v.0.0.1.docx)


## Contact: BuKyung Baik (back829@unist.ac.kr)

Any feedback or comments are greatly appreciated!!

## License

[GPL2.0](https://choosealicense.com/licenses/gpl-2.0/)
