################################################################################
#Burn-in

rm(list=ls())
start_time <- Sys.time()
set.seed(333)
library(AlphaSimR); library(writexl)
founderGenomes=quickHaplo(nInd=300, nChr=20, genLen=c(98, 140, 99, 112, 86, 136, 135, 146, 99, 132, 124, 120, 120, 108, 99, 92, 119, 105, 101, 112),  segSites=300, inbred=TRUE)
SP=SimParam$new(founderGenomes)
SP$addTraitAEG (210, mean=3000, var=24000) #Trait in kg/ha
SP$addSnpChip(90) #BARCSoySNP6K Illumina Inc.
Burnincycle1Par=newPop(founderGenomes[1:50]); Burnincycle1Par=setPheno(Burnincycle1Par, h2=0.40, reps=80)
Burnincycle2Par=newPop(founderGenomes[1:50]); Burnincycle2Par=setPheno(Burnincycle2Par, h2=0.40, reps=80)
Burnincycle3Par=newPop(founderGenomes[1:50]); Burnincycle3Par=setPheno(Burnincycle3Par, h2=0.40, reps=80)

################################################################################
#Burnincycle1 (burn-in)

#Parentals, F1s and nurseries
Burnincycle1F1=randCross(Burnincycle1Par, 1225)
Burnincycle1F2=self(Burnincycle1F1, nProgeny=1); Burnincycle1F3=self(Burnincycle1F2, nProgeny=25)
Burnincycle1F4=self(Burnincycle1F3, nProgeny=1); Burnincycle1F5=self(Burnincycle1F4, nProgeny=1)

#Trials 
Burnincycle1Pro=setPheno(Burnincycle1F5, h2=0.10, reps=1); Burnincycle1ProS=selectInd(Burnincycle1Pro, nInd=6000, use='pheno')
Burnincycle1Pre=setPheno(Burnincycle1ProS, h2=0.30, reps=10); Burnincycle1PreS=selectInd(Burnincycle1Pre, nInd=300, use='pheno')
Burnincycle1Adv=setPheno(Burnincycle1PreS, h2=0.40, reps=40); Burnincycle1AdvS=selectInd(Burnincycle1Adv, nInd=50, use='pheno')
Burnincycle1Eli=setPheno(Burnincycle1AdvS, h2=0.40, reps=80); Burnincycle1EliS=selectInd(Burnincycle1Eli, nInd=10, use='pheno')
Burnincycle1Eli1=selectInd(Burnincycle1Eli, nInd=1, use='pheno')

#Estimates and Qtl
calc_stats <- function(data) {corr <- cor(gv(data), pheno(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets <- list(Burnincycle1Pro, Burnincycle1Pre, Burnincycle1EliS)
names(datasets) <- c('Burnincycle1Pro', 'Burnincycle1Pre', 'Burnincycle1EliS')
Burnincycle1 <- do.call(rbind, lapply(datasets, calc_stats)); colnames(Burnincycle1) <- c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Burnincycle1
al <- function(data, simParam) {Qtl <- t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Burnincycle1PreQtl <- al(Burnincycle1Pre, SP); Burnincycle1EliSQtl<- al(Burnincycle1EliS, SP)

################################################################################
#Burnincycle2 (burn-in)

#Parentals, F1s and nurseries
Burnincycle2F1=randCross(Burnincycle2Par, 1225)	
Burnincycle2F2=self(Burnincycle2F1, nProgeny=1); Burnincycle2F3=self(Burnincycle2F2, nProgeny=25)
Burnincycle2F4=self(Burnincycle2F3, nProgeny=1); Burnincycle2F5=self(Burnincycle2F4, nProgeny=1)

#Trials 
Burnincycle2Pro=setPheno(Burnincycle2F5, h2=0.10, reps=1); Burnincycle2ProS=selectInd(Burnincycle2Pro, nInd=6000, use='pheno')
Burnincycle2Pre=setPheno(Burnincycle2ProS, h2=0.30, reps=10); Burnincycle2PreS=selectInd(Burnincycle2Pre, nInd=300, use='pheno')
Burnincycle2Adv=setPheno(Burnincycle2PreS, h2=0.40, reps=40); Burnincycle2AdvS=selectInd(Burnincycle2Adv, nInd=50, use='pheno')
Burnincycle2Eli=setPheno(Burnincycle2AdvS, h2=0.40, reps=80); Burnincycle2EliS=selectInd(Burnincycle2Eli, nInd=10, use='pheno')
Burnincycle2Eli1=selectInd(Burnincycle2Eli, nInd=1, use='pheno')

#Estimates and Qtl
calc_stats <- function(data) {corr <- cor(gv(data), pheno(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets <- list(Burnincycle2Pro, Burnincycle2Pre, Burnincycle2EliS)
names(datasets) <- c('Burnincycle2Pro', 'Burnincycle2Pre', 'Burnincycle2EliS')
Burnincycle2 <- do.call(rbind, lapply(datasets, calc_stats)); colnames(Burnincycle2) <- c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Burnincycle2
al <- function(data, simParam) {Qtl <- t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Burnincycle2PreQtl <- al(Burnincycle2Pre, SP); Burnincycle2EliSQtl<- al(Burnincycle2EliS, SP)

################################################################################
#Burnincycle3 (burn-in)

#Parentals, F1s and nurseries
Burnincycle3F1=randCross(Burnincycle3Par, 1225)
Burnincycle3F2=self(Burnincycle3F1, nProgeny=1); Burnincycle3F3=self(Burnincycle3F2, nProgeny=25)
Burnincycle3F4=self(Burnincycle3F3, nProgeny=1); Burnincycle3F5=self(Burnincycle3F4, nProgeny=1)

#Trials 
Burnincycle3Pro=setPheno(Burnincycle3F5, h2=0.10, reps=1); Burnincycle3ProS=selectInd(Burnincycle3Pro, nInd=6000, use='pheno')
Burnincycle3Pre=setPheno(Burnincycle3ProS, h2=0.30, reps=10); Burnincycle3PreS=selectInd(Burnincycle3Pre, nInd=300, use='pheno')
Burnincycle3Adv=setPheno(Burnincycle3PreS, h2=0.40, reps=40); Burnincycle3AdvS=selectInd(Burnincycle3Adv, nInd=50, use='pheno')
Burnincycle3Eli=setPheno(Burnincycle3AdvS, h2=0.40, reps=80); Burnincycle3EliS=selectInd(Burnincycle3Eli, nInd=10, use='pheno')
Burnincycle3Eli1=selectInd(Burnincycle3Eli, nInd=1, use='pheno')

#Estimates and Qtl
calc_stats <- function(data) {corr <- cor(gv(data), pheno(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets <- list(Burnincycle3Pro, Burnincycle3Pre, Burnincycle3EliS)
names(datasets) <- c('Burnincycle3Pro', 'Burnincycle3Pre', 'Burnincycle3EliS')
Burnincycle3 <- do.call(rbind, lapply(datasets, calc_stats)); colnames(Burnincycle3) <- c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Burnincycle3
al <- function(data, simParam) {Qtl <- t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Burnincycle3PreQtl <- al(Burnincycle3Pre, SP); Burnincycle3EliSQtl<- al(Burnincycle3EliS, SP)

################################################################################
#Burnincycle4 (burn-in)

#Parentals, F1s and nurseries
Burnincycle4Par=Burnincycle1PreS
Burnincycle4F1=randCross(Burnincycle4Par, 1225)
Burnincycle4F2=self(Burnincycle4F1, nProgeny=1); Burnincycle4F3=self(Burnincycle4F2, nProgeny=25)
Burnincycle4F4=self(Burnincycle4F3, nProgeny=1); Burnincycle4F5=self(Burnincycle4F4, nProgeny=1)

#Trials 
Burnincycle4Pro=setPheno(Burnincycle4F5, h2=0.10, reps=1); Burnincycle4ProS=selectInd(Burnincycle4Pro, nInd=6000, use='pheno')
Burnincycle4Pre=setPheno(Burnincycle4ProS, h2=0.30, reps=10); Burnincycle4PreS=selectInd(Burnincycle4Pre, nInd=300, use='pheno')
Burnincycle4Adv=setPheno(Burnincycle4PreS, h2=0.40, reps=40); Burnincycle4AdvS=selectInd(Burnincycle4Adv, nInd=50, use='pheno')
Burnincycle4Eli=setPheno(Burnincycle4AdvS, h2=0.40, reps=80); Burnincycle4EliS=selectInd(Burnincycle4Eli, nInd=10, use='pheno')
Burnincycle4Eli1=selectInd(Burnincycle4Eli, nInd=1, use='pheno')

#Rapid cycling
Par=Burnincycle2PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Burnincycle4CS=selectInd(F3, 150, use='ebv')
Burnincycle4CSF3Acc=cor(F3@gv, F3@ebv); Burnincycle4CSF3Acc
mean(Par@gv); mean(Burnincycle4CS@gv)

ParRound2=Burnincycle4CS
F1=randCross(ParRound2, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12); F2=setPheno(F2, h2=0.10, reps=2)
GSModel=fastRRBLUP(F2, use='pheno'); F3=setEBV(F3, GSModel); Burnincycle4CSRound2=selectInd(F3, 150, use='ebv')
Burnincycle4F3AccRound2=cor(F3@gv, F3@ebv); Burnincycle4F3AccRound2
mean(ParRound2@gv); mean(Burnincycle4CSRound2@gv)

#Estimates and Qtl
calc_stats <- function(data) {corr <- cor(gv(data), pheno(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets <- list(Burnincycle4Pro, Burnincycle4Pre, Burnincycle4EliS)
names(datasets) <- c('Burnincycle4Pro', 'Burnincycle4Pre', 'Burnincycle4EliS')
Burnincycle4 <- do.call(rbind, lapply(datasets, calc_stats)); colnames(Burnincycle4) <- c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Burnincycle4
al <- function(data, simParam) {Qtl <- t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Burnincycle4PreQtl <- al(Burnincycle4Pre, SP); Burnincycle4EliSQtl<- al(Burnincycle4EliS, SP)

################################################################################
#Burnincycle5 (burn-in)

#Parentals, F1s and nurseries
Burnincycle5Par=Burnincycle2PreS
Burnincycle5F1=randCross(Burnincycle5Par, 1225)
Burnincycle5F2=self(Burnincycle5F1, nProgeny=1); Burnincycle5F3=self(Burnincycle5F2, nProgeny=25)
Burnincycle5F4=self(Burnincycle5F3, nProgeny=1); Burnincycle5F5=self(Burnincycle5F4, nProgeny=1)

#Trials 
Burnincycle5Pro=setPheno(Burnincycle5F5, h2=0.10, reps=1); Burnincycle5ProS=selectInd(Burnincycle5Pro, nInd=6000, use='pheno')
Burnincycle5Pre=setPheno(Burnincycle5ProS, h2=0.30, reps=10); Burnincycle5PreS=selectInd(Burnincycle5Pre, nInd=300, use='pheno')
Burnincycle5Adv=setPheno(Burnincycle5PreS, h2=0.40, reps=40); Burnincycle5AdvS=selectInd(Burnincycle5Adv, nInd=50, use='pheno')
Burnincycle5Eli=setPheno(Burnincycle5AdvS, h2=0.40, reps=80); Burnincycle5EliS=selectInd(Burnincycle5Eli, nInd=10, use='pheno')
Burnincycle5Eli1=selectInd(Burnincycle5Eli, nInd=1, use='pheno')

#Rapid cycling
Par=Burnincycle3PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Burnincycle5CS=selectInd(F3, 150, use='ebv')
Burnincycle5CSF3Acc=cor(F3@gv, F3@ebv); Burnincycle5CSF3Acc
mean(Par@gv); mean(Burnincycle5CS@gv)

ParRound2=Burnincycle5CS
F1=randCross(ParRound2, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12); F2=setPheno(F2, h2=0.10, reps=2)
GSModel=fastRRBLUP(F2, use='pheno'); F3=setEBV(F3, GSModel); Burnincycle5CSRound2=selectInd(F3, 150, use='ebv')
Burnincycle5F3AccRound2=cor(F3@gv, F3@ebv); Burnincycle5F3AccRound2
mean(ParRound2@gv); mean(Burnincycle5CSRound2@gv)

#Estimates and Qtl
calc_stats <- function(data) {corr <- cor(gv(data), pheno(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets <- list(Burnincycle5Pro, Burnincycle5Pre, Burnincycle5EliS)
names(datasets) <- c('Burnincycle5Pro', 'Burnincycle5Pre', 'Burnincycle5EliS')
Burnincycle5 <- do.call(rbind, lapply(datasets, calc_stats)); colnames(Burnincycle5) <- c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Burnincycle5
al <- function(data, simParam) {Qtl <- t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Burnincycle5PreQtl <- al(Burnincycle5Pre, SP); Burnincycle5EliSQtl<- al(Burnincycle5EliS, SP)

################################################################################
#Genomic1cycle6
#Training, parentals and F1s
TrPar=Burnincycle3PreS; GSModel=fastRRBLUP(TrPar, use='pheno')
Genomic1cycle6Par=setEBV(TrPar, GSModel); Genomic1cycle6F1=randCross(Genomic1cycle6Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic1cycle6F2=self(Genomic1cycle6F1, nProgeny=1); Genomic1cycle6F3=self(Genomic1cycle6F2, nProgeny=25)
Genomic1cycle6F3=setEBV(Genomic1cycle6F3, GSModel); cor(Genomic1cycle6F3@gv, Genomic1cycle6F3@ebv)
Genomic1cycle6F3S=selectInd(Genomic1cycle6F3, nInd=6000, use='ebv'); mean(Genomic1cycle6F3@gv); mean(Genomic1cycle6F3S@gv)
Genomic1cycle6F4=self(Genomic1cycle6F3S, nProgeny=1); Genomic1cycle6F5=self(Genomic1cycle6F4, nProgeny=1)

#Trials 
Genomic1cycle6Pre=setPheno(Genomic1cycle6F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic1cycle6Pre, traits=1, use='pheno')
Genomic1cycle6Pre=setEBV(Genomic1cycle6Pre, GSModel); Genomic1cycle6PreS=selectInd(Genomic1cycle6Pre, nInd=300, use='ebv')
Genomic1cycle6Adv=setPheno(Genomic1cycle6PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic1cycle6Adv, traits=1, use='pheno')
Genomic1cycle6Adv=setEBV(Genomic1cycle6Adv, GSModel); Genomic1cycle6AdvS=selectInd(Genomic1cycle6Adv, nInd=50, use='ebv')
Genomic1cycle6Eli=setPheno(Genomic1cycle6AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic1cycle6Eli, traits=1, use='pheno')
Genomic1cycle6Eli=setEBV(Genomic1cycle6Eli, GSModel); Genomic1cycle6EliS=selectInd(Genomic1cycle6Eli, nInd=10, use='ebv')
Genomic1cycle6Eli1=selectInd(Genomic1cycle6Eli, nInd=1, use='ebv')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic1cycle6F3, Genomic1cycle6Pre, Genomic1cycle6EliS)
names(datasets)=c('Genomic1cycle6F3', 'Genomic1cycle6Pre', 'Genomic1cycle6EliS')
Genomic1cycle6=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic1cycle6)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic1cycle6
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic1cycle6PreQtl=al(Genomic1cycle6Pre, SP); Genomic1cycle6EliSQtl= al(Genomic1cycle6EliS, SP)

#Genomic1cycle7
#Training, parentals and F1s
TrPar=Burnincycle4PreS; GSModel=fastRRBLUP(TrPar, use='pheno')
Genomic1cycle7Par=setEBV(TrPar, GSModel); Genomic1cycle7F1=randCross(Genomic1cycle7Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic1cycle7F2=self(Genomic1cycle7F1, nProgeny=1); Genomic1cycle7F3=self(Genomic1cycle7F2, nProgeny=25)
Genomic1cycle7F3=setEBV(Genomic1cycle7F3, GSModel); cor(Genomic1cycle7F3@gv, Genomic1cycle7F3@ebv)
Genomic1cycle7F3S=selectInd(Genomic1cycle7F3, nInd=6000, use='ebv'); mean(Genomic1cycle7F3@gv); mean(Genomic1cycle7F3S@gv)
Genomic1cycle7F4=self(Genomic1cycle7F3S, nProgeny=1); Genomic1cycle7F5=self(Genomic1cycle7F4, nProgeny=1)

#Trials 
Genomic1cycle7Pre=setPheno(Genomic1cycle7F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic1cycle7Pre, traits=1, use='pheno')
Genomic1cycle7Pre=setEBV(Genomic1cycle7Pre, GSModel); Genomic1cycle7PreS=selectInd(Genomic1cycle7Pre, nInd=300, use='ebv')
Genomic1cycle7Adv=setPheno(Genomic1cycle7PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic1cycle7Adv, traits=1, use='pheno')
Genomic1cycle7Adv=setEBV(Genomic1cycle7Adv, GSModel); Genomic1cycle7AdvS=selectInd(Genomic1cycle7Adv, nInd=50, use='ebv')
Genomic1cycle7Eli=setPheno(Genomic1cycle7AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic1cycle7Eli, traits=1, use='pheno')
Genomic1cycle7Eli=setEBV(Genomic1cycle7Eli, GSModel); Genomic1cycle7EliS=selectInd(Genomic1cycle7Eli, nInd=10, use='ebv')
Genomic1cycle7Eli1=selectInd(Genomic1cycle7Eli, nInd=1, use='ebv')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic1cycle7F3, Genomic1cycle7Pre, Genomic1cycle7EliS)
names(datasets)=c('Genomic1cycle7F3', 'Genomic1cycle7Pre', 'Genomic1cycle7EliS')
Genomic1cycle7=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic1cycle7)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic1cycle7
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic1cycle7PreQtl=al(Genomic1cycle7Pre, SP); Genomic1cycle7EliSQtl= al(Genomic1cycle7EliS, SP)

#Genomic1cycle8
#Training, parentals and F1s
TrPar=Genomic1cycle6PreS; GSModel=fastRRBLUP(TrPar, use='pheno')
Genomic1cycle8Par=setEBV(TrPar, GSModel); Genomic1cycle8F1=randCross(Genomic1cycle8Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic1cycle8F2=self(Genomic1cycle8F1, nProgeny=1); Genomic1cycle8F3=self(Genomic1cycle8F2, nProgeny=25)
Genomic1cycle8F3=setEBV(Genomic1cycle8F3, GSModel); cor(Genomic1cycle8F3@gv, Genomic1cycle8F3@ebv)
Genomic1cycle8F3S=selectInd(Genomic1cycle8F3, nInd=6000, use='ebv'); mean(Genomic1cycle8F3@gv); mean(Genomic1cycle8F3S@gv)
Genomic1cycle8F4=self(Genomic1cycle8F3S, nProgeny=1); Genomic1cycle8F5=self(Genomic1cycle8F4, nProgeny=1)

#Trials 
Genomic1cycle8Pre=setPheno(Genomic1cycle8F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic1cycle8Pre, traits=1, use='pheno')
Genomic1cycle8Pre=setEBV(Genomic1cycle8Pre, GSModel); Genomic1cycle8PreS=selectInd(Genomic1cycle8Pre, nInd=300, use='ebv')
Genomic1cycle8Adv=setPheno(Genomic1cycle8PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic1cycle8Adv, traits=1, use='pheno')
Genomic1cycle8Adv=setEBV(Genomic1cycle8Adv, GSModel); Genomic1cycle8AdvS=selectInd(Genomic1cycle8Adv, nInd=50, use='ebv')
Genomic1cycle8Eli=setPheno(Genomic1cycle8AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic1cycle8Eli, traits=1, use='pheno')
Genomic1cycle8Eli=setEBV(Genomic1cycle8Eli, GSModel); Genomic1cycle8EliS=selectInd(Genomic1cycle8Eli, nInd=10, use='ebv')
Genomic1cycle8Eli1=selectInd(Genomic1cycle8Eli, nInd=1, use='ebv')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic1cycle8F3, Genomic1cycle8Pre, Genomic1cycle8EliS)
names(datasets)=c('Genomic1cycle8F3', 'Genomic1cycle8Pre', 'Genomic1cycle8EliS')
Genomic1cycle8=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic1cycle8)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic1cycle8
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic1cycle8PreQtl=al(Genomic1cycle8Pre, SP); Genomic1cycle8EliSQtl= al(Genomic1cycle8EliS, SP)

#Genomic1cycle9
#Training, parentals and F1s
TrPar=Genomic1cycle7PreS; GSModel=fastRRBLUP(TrPar, use='pheno')
Genomic1cycle9Par=setEBV(TrPar, GSModel); Genomic1cycle9F1=randCross(Genomic1cycle9Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic1cycle9F2=self(Genomic1cycle9F1, nProgeny=1); Genomic1cycle9F3=self(Genomic1cycle9F2, nProgeny=25)
Genomic1cycle9F3=setEBV(Genomic1cycle9F3, GSModel); cor(Genomic1cycle9F3@gv, Genomic1cycle9F3@ebv)
Genomic1cycle9F3S=selectInd(Genomic1cycle9F3, nInd=6000, use='ebv'); mean(Genomic1cycle9F3@gv); mean(Genomic1cycle9F3S@gv)
Genomic1cycle9F4=self(Genomic1cycle9F3S, nProgeny=1); Genomic1cycle9F5=self(Genomic1cycle9F4, nProgeny=1)

#Trials 
Genomic1cycle9Pre=setPheno(Genomic1cycle9F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic1cycle9Pre, traits=1, use='pheno')
Genomic1cycle9Pre=setEBV(Genomic1cycle9Pre, GSModel); Genomic1cycle9PreS=selectInd(Genomic1cycle9Pre, nInd=300, use='ebv')
Genomic1cycle9Adv=setPheno(Genomic1cycle9PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic1cycle9Adv, traits=1, use='pheno')
Genomic1cycle9Adv=setEBV(Genomic1cycle9Adv, GSModel); Genomic1cycle9AdvS=selectInd(Genomic1cycle9Adv, nInd=50, use='ebv')
Genomic1cycle9Eli=setPheno(Genomic1cycle9AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic1cycle9Eli, traits=1, use='pheno')
Genomic1cycle9Eli=setEBV(Genomic1cycle9Eli, GSModel); Genomic1cycle9EliS=selectInd(Genomic1cycle9Eli, nInd=10, use='ebv')
Genomic1cycle9Eli1=selectInd(Genomic1cycle9Eli, nInd=1, use='ebv')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic1cycle9F3, Genomic1cycle9Pre, Genomic1cycle9EliS)
names(datasets)=c('Genomic1cycle9F3', 'Genomic1cycle9Pre', 'Genomic1cycle9EliS')
Genomic1cycle9=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic1cycle9)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic1cycle9
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic1cycle9PreQtl=al(Genomic1cycle9Pre, SP); Genomic1cycle9EliSQtl= al(Genomic1cycle9EliS, SP)

#Genomic1cycle10
#Training, parentals and F1s
TrPar=Genomic1cycle8PreS; GSModel=fastRRBLUP(TrPar, use='pheno')
Genomic1cycle10Par=setEBV(TrPar, GSModel); Genomic1cycle10F1=randCross(Genomic1cycle10Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic1cycle10F2=self(Genomic1cycle10F1, nProgeny=1); Genomic1cycle10F3=self(Genomic1cycle10F2, nProgeny=25)
Genomic1cycle10F3=setEBV(Genomic1cycle10F3, GSModel); cor(Genomic1cycle10F3@gv, Genomic1cycle10F3@ebv)
Genomic1cycle10F3S=selectInd(Genomic1cycle10F3, nInd=6000, use='ebv'); mean(Genomic1cycle10F3@gv); mean(Genomic1cycle10F3S@gv)
Genomic1cycle10F4=self(Genomic1cycle10F3S, nProgeny=1); Genomic1cycle10F5=self(Genomic1cycle10F4, nProgeny=1)

#Trials 
Genomic1cycle10Pre=setPheno(Genomic1cycle10F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic1cycle10Pre, traits=1, use='pheno')
Genomic1cycle10Pre=setEBV(Genomic1cycle10Pre, GSModel); Genomic1cycle10PreS=selectInd(Genomic1cycle10Pre, nInd=300, use='ebv')
Genomic1cycle10Adv=setPheno(Genomic1cycle10PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic1cycle10Adv, traits=1, use='pheno')
Genomic1cycle10Adv=setEBV(Genomic1cycle10Adv, GSModel); Genomic1cycle10AdvS=selectInd(Genomic1cycle10Adv, nInd=50, use='ebv')
Genomic1cycle10Eli=setPheno(Genomic1cycle10AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic1cycle10Eli, traits=1, use='pheno')
Genomic1cycle10Eli=setEBV(Genomic1cycle10Eli, GSModel); Genomic1cycle10EliS=selectInd(Genomic1cycle10Eli, nInd=10, use='ebv')
Genomic1cycle10Eli1=selectInd(Genomic1cycle10Eli, nInd=1, use='ebv')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic1cycle10F3, Genomic1cycle10Pre, Genomic1cycle10EliS)
names(datasets)=c('Genomic1cycle10F3', 'Genomic1cycle10Pre', 'Genomic1cycle10EliS')
Genomic1cycle10=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic1cycle10)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic1cycle10
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic1cycle10PreQtl=al(Genomic1cycle10Pre, SP); Genomic1cycle10EliSQtl= al(Genomic1cycle10EliS, SP)

#Genomic1cycle11
#Training, parentals and F1s
TrPar=Genomic1cycle9PreS; GSModel=fastRRBLUP(TrPar, use='pheno')
Genomic1cycle11Par=setEBV(TrPar, GSModel); Genomic1cycle11F1=randCross(Genomic1cycle11Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic1cycle11F2=self(Genomic1cycle11F1, nProgeny=1); Genomic1cycle11F3=self(Genomic1cycle11F2, nProgeny=25)
Genomic1cycle11F3=setEBV(Genomic1cycle11F3, GSModel); cor(Genomic1cycle11F3@gv, Genomic1cycle11F3@ebv)
Genomic1cycle11F3S=selectInd(Genomic1cycle11F3, nInd=6000, use='ebv'); mean(Genomic1cycle11F3@gv); mean(Genomic1cycle11F3S@gv)
Genomic1cycle11F4=self(Genomic1cycle11F3S, nProgeny=1); Genomic1cycle11F5=self(Genomic1cycle11F4, nProgeny=1)

#Trials 
Genomic1cycle11Pre=setPheno(Genomic1cycle11F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic1cycle11Pre, traits=1, use='pheno')
Genomic1cycle11Pre=setEBV(Genomic1cycle11Pre, GSModel); Genomic1cycle11PreS=selectInd(Genomic1cycle11Pre, nInd=300, use='ebv')
Genomic1cycle11Adv=setPheno(Genomic1cycle11PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic1cycle11Adv, traits=1, use='pheno')
Genomic1cycle11Adv=setEBV(Genomic1cycle11Adv, GSModel); Genomic1cycle11AdvS=selectInd(Genomic1cycle11Adv, nInd=50, use='ebv')
Genomic1cycle11Eli=setPheno(Genomic1cycle11AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic1cycle11Eli, traits=1, use='pheno')
Genomic1cycle11Eli=setEBV(Genomic1cycle11Eli, GSModel); Genomic1cycle11EliS=selectInd(Genomic1cycle11Eli, nInd=10, use='ebv')
Genomic1cycle11Eli1=selectInd(Genomic1cycle11Eli, nInd=1, use='ebv')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic1cycle11F3, Genomic1cycle11Pre, Genomic1cycle11EliS)
names(datasets)=c('Genomic1cycle11F3', 'Genomic1cycle11Pre', 'Genomic1cycle11EliS')
Genomic1cycle11=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic1cycle11)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic1cycle11
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic1cycle11PreQtl=al(Genomic1cycle11Pre, SP); Genomic1cycle11EliSQtl= al(Genomic1cycle11EliS, SP)

#Genomic1cycle12
#Training, parentals and F1s
TrPar=Genomic1cycle10PreS; GSModel=fastRRBLUP(TrPar, use='pheno')
Genomic1cycle12Par=setEBV(TrPar, GSModel); Genomic1cycle12F1=randCross(Genomic1cycle12Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic1cycle12F2=self(Genomic1cycle12F1, nProgeny=1); Genomic1cycle12F3=self(Genomic1cycle12F2, nProgeny=25)
Genomic1cycle12F3=setEBV(Genomic1cycle12F3, GSModel); cor(Genomic1cycle12F3@gv, Genomic1cycle12F3@ebv)
Genomic1cycle12F3S=selectInd(Genomic1cycle12F3, nInd=6000, use='ebv'); mean(Genomic1cycle12F3@gv); mean(Genomic1cycle12F3S@gv)
Genomic1cycle12F4=self(Genomic1cycle12F3S, nProgeny=1); Genomic1cycle12F5=self(Genomic1cycle12F4, nProgeny=1)

#Trials 
Genomic1cycle12Pre=setPheno(Genomic1cycle12F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic1cycle12Pre, traits=1, use='pheno')
Genomic1cycle12Pre=setEBV(Genomic1cycle12Pre, GSModel); Genomic1cycle12PreS=selectInd(Genomic1cycle12Pre, nInd=300, use='ebv')
Genomic1cycle12Adv=setPheno(Genomic1cycle12PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic1cycle12Adv, traits=1, use='pheno')
Genomic1cycle12Adv=setEBV(Genomic1cycle12Adv, GSModel); Genomic1cycle12AdvS=selectInd(Genomic1cycle12Adv, nInd=50, use='ebv')
Genomic1cycle12Eli=setPheno(Genomic1cycle12AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic1cycle12Eli, traits=1, use='pheno')
Genomic1cycle12Eli=setEBV(Genomic1cycle12Eli, GSModel); Genomic1cycle12EliS=selectInd(Genomic1cycle12Eli, nInd=10, use='ebv')
Genomic1cycle12Eli1=selectInd(Genomic1cycle12Eli, nInd=1, use='ebv')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic1cycle12F3, Genomic1cycle12Pre, Genomic1cycle12EliS)
names(datasets)=c('Genomic1cycle12F3', 'Genomic1cycle12Pre', 'Genomic1cycle12EliS')
Genomic1cycle12=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic1cycle12)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic1cycle12
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic1cycle12PreQtl=al(Genomic1cycle12Pre, SP); Genomic1cycle12EliSQtl= al(Genomic1cycle12EliS, SP)

#Genomic1cycle13
#Training, parentals and F1s
TrPar=Genomic1cycle11PreS; GSModel=fastRRBLUP(TrPar, use='pheno')
Genomic1cycle13Par=setEBV(TrPar, GSModel); Genomic1cycle13F1=randCross(Genomic1cycle13Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic1cycle13F2=self(Genomic1cycle13F1, nProgeny=1); Genomic1cycle13F3=self(Genomic1cycle13F2, nProgeny=25)
Genomic1cycle13F3=setEBV(Genomic1cycle13F3, GSModel); cor(Genomic1cycle13F3@gv, Genomic1cycle13F3@ebv)
Genomic1cycle13F3S=selectInd(Genomic1cycle13F3, nInd=6000, use='ebv'); mean(Genomic1cycle13F3@gv); mean(Genomic1cycle13F3S@gv)
Genomic1cycle13F4=self(Genomic1cycle13F3S, nProgeny=1); Genomic1cycle13F5=self(Genomic1cycle13F4, nProgeny=1)

#Trials 
Genomic1cycle13Pre=setPheno(Genomic1cycle13F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic1cycle13Pre, traits=1, use='pheno')
Genomic1cycle13Pre=setEBV(Genomic1cycle13Pre, GSModel); Genomic1cycle13PreS=selectInd(Genomic1cycle13Pre, nInd=300, use='ebv')
Genomic1cycle13Adv=setPheno(Genomic1cycle13PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic1cycle13Adv, traits=1, use='pheno')
Genomic1cycle13Adv=setEBV(Genomic1cycle13Adv, GSModel); Genomic1cycle13AdvS=selectInd(Genomic1cycle13Adv, nInd=50, use='ebv')
Genomic1cycle13Eli=setPheno(Genomic1cycle13AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic1cycle13Eli, traits=1, use='pheno')
Genomic1cycle13Eli=setEBV(Genomic1cycle13Eli, GSModel); Genomic1cycle13EliS=selectInd(Genomic1cycle13Eli, nInd=10, use='ebv')
Genomic1cycle13Eli1=selectInd(Genomic1cycle13Eli, nInd=1, use='ebv')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic1cycle13F3, Genomic1cycle13Pre, Genomic1cycle13EliS)
names(datasets)=c('Genomic1cycle13F3', 'Genomic1cycle13Pre', 'Genomic1cycle13EliS')
Genomic1cycle13=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic1cycle13)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic1cycle13
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic1cycle13PreQtl=al(Genomic1cycle13Pre, SP); Genomic1cycle13EliSQtl= al(Genomic1cycle13EliS, SP)

#Genomic1cycle14
#Training, parentals and F1s
TrPar=Genomic1cycle12PreS; GSModel=fastRRBLUP(TrPar, use='pheno')
Genomic1cycle14Par=setEBV(TrPar, GSModel); Genomic1cycle14F1=randCross(Genomic1cycle14Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic1cycle14F2=self(Genomic1cycle14F1, nProgeny=1); Genomic1cycle14F3=self(Genomic1cycle14F2, nProgeny=25)
Genomic1cycle14F3=setEBV(Genomic1cycle14F3, GSModel); cor(Genomic1cycle14F3@gv, Genomic1cycle14F3@ebv)
Genomic1cycle14F3S=selectInd(Genomic1cycle14F3, nInd=6000, use='ebv'); mean(Genomic1cycle14F3@gv); mean(Genomic1cycle14F3S@gv)
Genomic1cycle14F4=self(Genomic1cycle14F3S, nProgeny=1); Genomic1cycle14F5=self(Genomic1cycle14F4, nProgeny=1)

#Trials 
Genomic1cycle14Pre=setPheno(Genomic1cycle14F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic1cycle14Pre, traits=1, use='pheno')
Genomic1cycle14Pre=setEBV(Genomic1cycle14Pre, GSModel); Genomic1cycle14PreS=selectInd(Genomic1cycle14Pre, nInd=300, use='ebv')
Genomic1cycle14Adv=setPheno(Genomic1cycle14PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic1cycle14Adv, traits=1, use='pheno')
Genomic1cycle14Adv=setEBV(Genomic1cycle14Adv, GSModel); Genomic1cycle14AdvS=selectInd(Genomic1cycle14Adv, nInd=50, use='ebv')
Genomic1cycle14Eli=setPheno(Genomic1cycle14AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic1cycle14Eli, traits=1, use='pheno')
Genomic1cycle14Eli=setEBV(Genomic1cycle14Eli, GSModel); Genomic1cycle14EliS=selectInd(Genomic1cycle14Eli, nInd=10, use='ebv')
Genomic1cycle14Eli1=selectInd(Genomic1cycle14Eli, nInd=1, use='ebv')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic1cycle14F3, Genomic1cycle14Pre, Genomic1cycle14EliS)
names(datasets)=c('Genomic1cycle14F3', 'Genomic1cycle14Pre', 'Genomic1cycle14EliS')
Genomic1cycle14=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic1cycle14)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic1cycle14
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic1cycle14PreQtl=al(Genomic1cycle14Pre, SP); Genomic1cycle14EliSQtl= al(Genomic1cycle14EliS, SP)

#Genomic1cycle15
#Training, parentals and F1s
TrPar=Genomic1cycle13PreS; GSModel=fastRRBLUP(TrPar, use='pheno')
Genomic1cycle15Par=setEBV(TrPar, GSModel); Genomic1cycle15F1=randCross(Genomic1cycle15Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic1cycle15F2=self(Genomic1cycle15F1, nProgeny=1); Genomic1cycle15F3=self(Genomic1cycle15F2, nProgeny=25)
Genomic1cycle15F3=setEBV(Genomic1cycle15F3, GSModel); cor(Genomic1cycle15F3@gv, Genomic1cycle15F3@ebv)
Genomic1cycle15F3S=selectInd(Genomic1cycle15F3, nInd=6000, use='ebv'); mean(Genomic1cycle15F3@gv); mean(Genomic1cycle15F3S@gv)
Genomic1cycle15F4=self(Genomic1cycle15F3S, nProgeny=1); Genomic1cycle15F5=self(Genomic1cycle15F4, nProgeny=1)

#Trials 
Genomic1cycle15Pre=setPheno(Genomic1cycle15F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic1cycle15Pre, traits=1, use='pheno')
Genomic1cycle15Pre=setEBV(Genomic1cycle15Pre, GSModel); Genomic1cycle15PreS=selectInd(Genomic1cycle15Pre, nInd=300, use='ebv')
Genomic1cycle15Adv=setPheno(Genomic1cycle15PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic1cycle15Adv, traits=1, use='pheno')
Genomic1cycle15Adv=setEBV(Genomic1cycle15Adv, GSModel); Genomic1cycle15AdvS=selectInd(Genomic1cycle15Adv, nInd=50, use='ebv')
Genomic1cycle15Eli=setPheno(Genomic1cycle15AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic1cycle15Eli, traits=1, use='pheno')
Genomic1cycle15Eli=setEBV(Genomic1cycle15Eli, GSModel); Genomic1cycle15EliS=selectInd(Genomic1cycle15Eli, nInd=10, use='ebv')
Genomic1cycle15Eli1=selectInd(Genomic1cycle15Eli, nInd=1, use='ebv')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic1cycle15F3, Genomic1cycle15Pre, Genomic1cycle15EliS)
names(datasets)=c('Genomic1cycle15F3', 'Genomic1cycle15Pre', 'Genomic1cycle15EliS')
Genomic1cycle15=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic1cycle15)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic1cycle15
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic1cycle15PreQtl=al(Genomic1cycle15Pre, SP); Genomic1cycle15EliSQtl= al(Genomic1cycle15EliS, SP)

################################################################################
#Genomic1 exporting estimates and Qtl

Genomic1cycle1to15PreQtl=cbind(Burnincycle1PreQtl, Burnincycle2PreQtl, Burnincycle3PreQtl, Burnincycle4PreQtl, Burnincycle5PreQtl, Genomic1cycle6PreQtl, Genomic1cycle7PreQtl, Genomic1cycle8PreQtl, Genomic1cycle9PreQtl, Genomic1cycle10PreQtl, Genomic1cycle11PreQtl, Genomic1cycle12PreQtl, Genomic1cycle13PreQtl, Genomic1cycle14PreQtl, Genomic1cycle15PreQtl); write.csv(Genomic1cycle1to15PreQtl, 'r1_Genomic1PreQtl.csv')

Genomic1cycle1to15EliSQtl=cbind(Burnincycle1EliSQtl, Burnincycle2EliSQtl, Burnincycle3EliSQtl, Burnincycle4EliSQtl, Burnincycle5EliSQtl, Genomic1cycle6EliSQtl, Genomic1cycle7EliSQtl, Genomic1cycle8EliSQtl, Genomic1cycle9EliSQtl, Genomic1cycle10EliSQtl, Genomic1cycle11EliSQtl, Genomic1cycle12EliSQtl, Genomic1cycle13EliSQtl, Genomic1cycle14EliSQtl, Genomic1cycle15EliSQtl); write.csv(Genomic1cycle1to15EliSQtl, 'r1_Genomic1EliSQtl.csv')

Genomic1cycle1to15=rbind(Burnincycle1, Burnincycle2, Burnincycle3, Burnincycle4, Burnincycle5, Genomic1cycle6, Genomic1cycle7, Genomic1cycle8, Genomic1cycle9, Genomic1cycle10, Genomic1cycle11, Genomic1cycle12, Genomic1cycle13, Genomic1cycle14, Genomic1cycle15); write.csv(Genomic1cycle1to15, 'r1_Genomic1.csv')

################################################################################
#Genomic2cycle6
#Training, parentals and F1s
Tr=list(Burnincycle3PreS, Burnincycle2PreS); Tr=mergePops(Tr)
Par1=Burnincycle3PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Burnincycle4CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic2cycle6Par=setEBV(Par, GSModel); Genomic2cycle6F1=randCross(Genomic2cycle6Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic2cycle6F2=self(Genomic2cycle6F1, nProgeny=1); Genomic2cycle6F3=self(Genomic2cycle6F2, nProgeny=25)
Genomic2cycle6F3=setEBV(Genomic2cycle6F3, GSModel); cor(Genomic2cycle6F3@gv, Genomic2cycle6F3@ebv)
Genomic2cycle6F3S=selectInd(Genomic2cycle6F3, nInd=6000, use='ebv'); mean(Genomic2cycle6F3@gv); mean(Genomic2cycle6F3S@gv)
Genomic2cycle6F4=self(Genomic2cycle6F3S, nProgeny=1); Genomic2cycle6F5=self(Genomic2cycle6F4, nProgeny=1)

#Trials 
Genomic2cycle6Pre=setPheno(Genomic2cycle6F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic2cycle6Pre, traits=1, use='pheno')
Genomic2cycle6Pre=setEBV(Genomic2cycle6Pre, GSModel); Genomic2cycle6PreS=selectInd(Genomic2cycle6Pre, nInd=300, use='ebv')
Genomic2cycle6Adv=setPheno(Genomic2cycle6PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic2cycle6Adv, traits=1, use='pheno')
Genomic2cycle6Adv=setEBV(Genomic2cycle6Adv, GSModel); Genomic2cycle6AdvS=selectInd(Genomic2cycle6Adv, nInd=50, use='ebv')
Genomic2cycle6Eli=setPheno(Genomic2cycle6AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic2cycle6Eli, traits=1, use='pheno')
Genomic2cycle6Eli=setEBV(Genomic2cycle6Eli, GSModel); Genomic2cycle6EliS=selectInd(Genomic2cycle6Eli, nInd=10, use='ebv')
Genomic2cycle6Eli1=selectInd(Genomic2cycle6Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Burnincycle4PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic2cycle6CS=selectInd(F3, 150, use='ebv')
Genomic2cycle6F3Acc=cor(F3@gv, F3@ebv); Genomic2cycle6F3Acc
mean(Par@gv); mean(Genomic2cycle6CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic2cycle6F3, Genomic2cycle6Pre, Genomic2cycle6EliS)
names(datasets)=c('Genomic2cycle6F3', 'Genomic2cycle6Pre', 'Genomic2cycle6EliS')
Genomic2cycle6=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic2cycle6)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic2cycle6
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic2cycle6PreQtl=al(Genomic2cycle6Pre, SP); Genomic2cycle6EliSQtl= al(Genomic2cycle6EliS, SP)

#Genomic2cycle7
#Training, parentals and F1s
Tr=list(Burnincycle4PreS, Burnincycle3PreS); Tr=mergePops(Tr)
Par1=Burnincycle4PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Burnincycle5CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic2cycle7Par=setEBV(Par, GSModel); Genomic2cycle7F1=randCross(Genomic2cycle7Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic2cycle7F2=self(Genomic2cycle7F1, nProgeny=1); Genomic2cycle7F3=self(Genomic2cycle7F2, nProgeny=25)
Genomic2cycle7F3=setEBV(Genomic2cycle7F3, GSModel); cor(Genomic2cycle7F3@gv, Genomic2cycle7F3@ebv)
Genomic2cycle7F3S=selectInd(Genomic2cycle7F3, nInd=6000, use='ebv'); mean(Genomic2cycle7F3@gv); mean(Genomic2cycle7F3S@gv)
Genomic2cycle7F4=self(Genomic2cycle7F3S, nProgeny=1); Genomic2cycle7F5=self(Genomic2cycle7F4, nProgeny=1)

#Trials 
Genomic2cycle7Pre=setPheno(Genomic2cycle7F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic2cycle7Pre, traits=1, use='pheno')
Genomic2cycle7Pre=setEBV(Genomic2cycle7Pre, GSModel); Genomic2cycle7PreS=selectInd(Genomic2cycle7Pre, nInd=300, use='ebv')
Genomic2cycle7Adv=setPheno(Genomic2cycle7PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic2cycle7Adv, traits=1, use='pheno')
Genomic2cycle7Adv=setEBV(Genomic2cycle7Adv, GSModel); Genomic2cycle7AdvS=selectInd(Genomic2cycle7Adv, nInd=50, use='ebv')
Genomic2cycle7Eli=setPheno(Genomic2cycle7AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic2cycle7Eli, traits=1, use='pheno')
Genomic2cycle7Eli=setEBV(Genomic2cycle7Eli, GSModel); Genomic2cycle7EliS=selectInd(Genomic2cycle7Eli, nInd=10, use='ebv')
Genomic2cycle7Eli1=selectInd(Genomic2cycle7Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic2cycle6PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic2cycle7CS=selectInd(F3, 150, use='ebv')
Genomic2cycle7F3Acc=cor(F3@gv, F3@ebv); Genomic2cycle7F3Acc
mean(Par@gv); mean(Genomic2cycle7CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic2cycle7F3, Genomic2cycle7Pre, Genomic2cycle7EliS)
names(datasets)=c('Genomic2cycle7F3', 'Genomic2cycle7Pre', 'Genomic2cycle7EliS')
Genomic2cycle7=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic2cycle7)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic2cycle7
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic2cycle7PreQtl=al(Genomic2cycle7Pre, SP); Genomic2cycle7EliSQtl= al(Genomic2cycle7EliS, SP)

#Genomic2cycle8
#Training, parentals and F1s
Tr=list(Genomic2cycle6PreS, Burnincycle4PreS); Tr=mergePops(Tr)
Par1=Genomic2cycle6PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic2cycle6CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic2cycle8Par=setEBV(Par, GSModel); Genomic2cycle8F1=randCross(Genomic2cycle8Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic2cycle8F2=self(Genomic2cycle8F1, nProgeny=1); Genomic2cycle8F3=self(Genomic2cycle8F2, nProgeny=25)
Genomic2cycle8F3=setEBV(Genomic2cycle8F3, GSModel); cor(Genomic2cycle8F3@gv, Genomic2cycle8F3@ebv)
Genomic2cycle8F3S=selectInd(Genomic2cycle8F3, nInd=6000, use='ebv'); mean(Genomic2cycle8F3@gv); mean(Genomic2cycle8F3S@gv)
Genomic2cycle8F4=self(Genomic2cycle8F3S, nProgeny=1); Genomic2cycle8F5=self(Genomic2cycle8F4, nProgeny=1)

#Trials 
Genomic2cycle8Pre=setPheno(Genomic2cycle8F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic2cycle8Pre, traits=1, use='pheno')
Genomic2cycle8Pre=setEBV(Genomic2cycle8Pre, GSModel); Genomic2cycle8PreS=selectInd(Genomic2cycle8Pre, nInd=300, use='ebv')
Genomic2cycle8Adv=setPheno(Genomic2cycle8PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic2cycle8Adv, traits=1, use='pheno')
Genomic2cycle8Adv=setEBV(Genomic2cycle8Adv, GSModel); Genomic2cycle8AdvS=selectInd(Genomic2cycle8Adv, nInd=50, use='ebv')
Genomic2cycle8Eli=setPheno(Genomic2cycle8AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic2cycle8Eli, traits=1, use='pheno')
Genomic2cycle8Eli=setEBV(Genomic2cycle8Eli, GSModel); Genomic2cycle8EliS=selectInd(Genomic2cycle8Eli, nInd=10, use='ebv')
Genomic2cycle8Eli1=selectInd(Genomic2cycle8Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic2cycle7PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic2cycle8CS=selectInd(F3, 150, use='ebv')
Genomic2cycle8F3Acc=cor(F3@gv, F3@ebv); Genomic2cycle8F3Acc
mean(Par@gv); mean(Genomic2cycle8CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic2cycle8F3, Genomic2cycle8Pre, Genomic2cycle8EliS)
names(datasets)=c('Genomic2cycle8F3', 'Genomic2cycle8Pre', 'Genomic2cycle8EliS')
Genomic2cycle8=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic2cycle8)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic2cycle8
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic2cycle8PreQtl=al(Genomic2cycle8Pre, SP); Genomic2cycle8EliSQtl= al(Genomic2cycle8EliS, SP)

#Genomic2cycle9
#Training, parentals and F1s
Tr=list(Genomic2cycle7PreS, Genomic2cycle6PreS); Tr=mergePops(Tr)
Par1=Genomic2cycle7PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic2cycle7CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic2cycle9Par=setEBV(Par, GSModel); Genomic2cycle9F1=randCross(Genomic2cycle9Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic2cycle9F2=self(Genomic2cycle9F1, nProgeny=1); Genomic2cycle9F3=self(Genomic2cycle9F2, nProgeny=25)
Genomic2cycle9F3=setEBV(Genomic2cycle9F3, GSModel); cor(Genomic2cycle9F3@gv, Genomic2cycle9F3@ebv)
Genomic2cycle9F3S=selectInd(Genomic2cycle9F3, nInd=6000, use='ebv'); mean(Genomic2cycle9F3@gv); mean(Genomic2cycle9F3S@gv)
Genomic2cycle9F4=self(Genomic2cycle9F3S, nProgeny=1); Genomic2cycle9F5=self(Genomic2cycle9F4, nProgeny=1)

#Trials 
Genomic2cycle9Pre=setPheno(Genomic2cycle9F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic2cycle9Pre, traits=1, use='pheno')
Genomic2cycle9Pre=setEBV(Genomic2cycle9Pre, GSModel); Genomic2cycle9PreS=selectInd(Genomic2cycle9Pre, nInd=300, use='ebv')
Genomic2cycle9Adv=setPheno(Genomic2cycle9PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic2cycle9Adv, traits=1, use='pheno')
Genomic2cycle9Adv=setEBV(Genomic2cycle9Adv, GSModel); Genomic2cycle9AdvS=selectInd(Genomic2cycle9Adv, nInd=50, use='ebv')
Genomic2cycle9Eli=setPheno(Genomic2cycle9AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic2cycle9Eli, traits=1, use='pheno')
Genomic2cycle9Eli=setEBV(Genomic2cycle9Eli, GSModel); Genomic2cycle9EliS=selectInd(Genomic2cycle9Eli, nInd=10, use='ebv')
Genomic2cycle9Eli1=selectInd(Genomic2cycle9Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic2cycle8PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic2cycle9CS=selectInd(F3, 150, use='ebv')
Genomic2cycle9F3Acc=cor(F3@gv, F3@ebv); Genomic2cycle9F3Acc
mean(Par@gv); mean(Genomic2cycle9CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic2cycle9F3, Genomic2cycle9Pre, Genomic2cycle9EliS)
names(datasets)=c('Genomic2cycle9F3', 'Genomic2cycle9Pre', 'Genomic2cycle9EliS')
Genomic2cycle9=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic2cycle9)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic2cycle9
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic2cycle9PreQtl=al(Genomic2cycle9Pre, SP); Genomic2cycle9EliSQtl= al(Genomic2cycle9EliS, SP)

#Genomic2cycle10
#Training, parentals and F1s
Tr=list(Genomic2cycle8PreS, Genomic2cycle7PreS); Tr=mergePops(Tr)
Par1=Genomic2cycle8PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic2cycle8CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic2cycle10Par=setEBV(Par, GSModel); Genomic2cycle10F1=randCross(Genomic2cycle10Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic2cycle10F2=self(Genomic2cycle10F1, nProgeny=1); Genomic2cycle10F3=self(Genomic2cycle10F2, nProgeny=25)
Genomic2cycle10F3=setEBV(Genomic2cycle10F3, GSModel); cor(Genomic2cycle10F3@gv, Genomic2cycle10F3@ebv)
Genomic2cycle10F3S=selectInd(Genomic2cycle10F3, nInd=6000, use='ebv'); mean(Genomic2cycle10F3@gv); mean(Genomic2cycle10F3S@gv)
Genomic2cycle10F4=self(Genomic2cycle10F3S, nProgeny=1); Genomic2cycle10F5=self(Genomic2cycle10F4, nProgeny=1)

#Trials 
Genomic2cycle10Pre=setPheno(Genomic2cycle10F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic2cycle10Pre, traits=1, use='pheno')
Genomic2cycle10Pre=setEBV(Genomic2cycle10Pre, GSModel); Genomic2cycle10PreS=selectInd(Genomic2cycle10Pre, nInd=300, use='ebv')
Genomic2cycle10Adv=setPheno(Genomic2cycle10PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic2cycle10Adv, traits=1, use='pheno')
Genomic2cycle10Adv=setEBV(Genomic2cycle10Adv, GSModel); Genomic2cycle10AdvS=selectInd(Genomic2cycle10Adv, nInd=50, use='ebv')
Genomic2cycle10Eli=setPheno(Genomic2cycle10AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic2cycle10Eli, traits=1, use='pheno')
Genomic2cycle10Eli=setEBV(Genomic2cycle10Eli, GSModel); Genomic2cycle10EliS=selectInd(Genomic2cycle10Eli, nInd=10, use='ebv')
Genomic2cycle10Eli1=selectInd(Genomic2cycle10Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic2cycle9PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic2cycle10CS=selectInd(F3, 150, use='ebv')
Genomic2cycle10F3Acc=cor(F3@gv, F3@ebv); Genomic2cycle10F3Acc
mean(Par@gv); mean(Genomic2cycle10CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic2cycle10F3, Genomic2cycle10Pre, Genomic2cycle10EliS)
names(datasets)=c('Genomic2cycle10F3', 'Genomic2cycle10Pre', 'Genomic2cycle10EliS')
Genomic2cycle10=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic2cycle10)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic2cycle10
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic2cycle10PreQtl=al(Genomic2cycle10Pre, SP); Genomic2cycle10EliSQtl= al(Genomic2cycle10EliS, SP)

#Genomic2cycle11
#Training, parentals and F1s
Tr=list(Genomic2cycle9PreS, Genomic2cycle8PreS); Tr=mergePops(Tr)
Par1=Genomic2cycle9PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic2cycle9CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic2cycle11Par=setEBV(Par, GSModel); Genomic2cycle11F1=randCross(Genomic2cycle11Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic2cycle11F2=self(Genomic2cycle11F1, nProgeny=1); Genomic2cycle11F3=self(Genomic2cycle11F2, nProgeny=25)
Genomic2cycle11F3=setEBV(Genomic2cycle11F3, GSModel); cor(Genomic2cycle11F3@gv, Genomic2cycle11F3@ebv)
Genomic2cycle11F3S=selectInd(Genomic2cycle11F3, nInd=6000, use='ebv'); mean(Genomic2cycle11F3@gv); mean(Genomic2cycle11F3S@gv)
Genomic2cycle11F4=self(Genomic2cycle11F3S, nProgeny=1); Genomic2cycle11F5=self(Genomic2cycle11F4, nProgeny=1)

#Trials 
Genomic2cycle11Pre=setPheno(Genomic2cycle11F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic2cycle11Pre, traits=1, use='pheno')
Genomic2cycle11Pre=setEBV(Genomic2cycle11Pre, GSModel); Genomic2cycle11PreS=selectInd(Genomic2cycle11Pre, nInd=300, use='ebv')
Genomic2cycle11Adv=setPheno(Genomic2cycle11PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic2cycle11Adv, traits=1, use='pheno')
Genomic2cycle11Adv=setEBV(Genomic2cycle11Adv, GSModel); Genomic2cycle11AdvS=selectInd(Genomic2cycle11Adv, nInd=50, use='ebv')
Genomic2cycle11Eli=setPheno(Genomic2cycle11AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic2cycle11Eli, traits=1, use='pheno')
Genomic2cycle11Eli=setEBV(Genomic2cycle11Eli, GSModel); Genomic2cycle11EliS=selectInd(Genomic2cycle11Eli, nInd=10, use='ebv')
Genomic2cycle11Eli1=selectInd(Genomic2cycle11Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic2cycle10PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic2cycle11CS=selectInd(F3, 150, use='ebv')
Genomic2cycle11F3Acc=cor(F3@gv, F3@ebv); Genomic2cycle11F3Acc
mean(Par@gv); mean(Genomic2cycle11CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic2cycle11F3, Genomic2cycle11Pre, Genomic2cycle11EliS)
names(datasets)=c('Genomic2cycle11F3', 'Genomic2cycle11Pre', 'Genomic2cycle11EliS')
Genomic2cycle11=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic2cycle11)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic2cycle11
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic2cycle11PreQtl=al(Genomic2cycle11Pre, SP); Genomic2cycle11EliSQtl= al(Genomic2cycle11EliS, SP)

#Genomic2cycle12
#Training, parentals and F1s
Tr=list(Genomic2cycle10PreS, Genomic2cycle9PreS); Tr=mergePops(Tr)
Par1=Genomic2cycle10PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic2cycle10CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic2cycle12Par=setEBV(Par, GSModel); Genomic2cycle12F1=randCross(Genomic2cycle12Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic2cycle12F2=self(Genomic2cycle12F1, nProgeny=1); Genomic2cycle12F3=self(Genomic2cycle12F2, nProgeny=25)
Genomic2cycle12F3=setEBV(Genomic2cycle12F3, GSModel); cor(Genomic2cycle12F3@gv, Genomic2cycle12F3@ebv)
Genomic2cycle12F3S=selectInd(Genomic2cycle12F3, nInd=6000, use='ebv'); mean(Genomic2cycle12F3@gv); mean(Genomic2cycle12F3S@gv)
Genomic2cycle12F4=self(Genomic2cycle12F3S, nProgeny=1); Genomic2cycle12F5=self(Genomic2cycle12F4, nProgeny=1)

#Trials 
Genomic2cycle12Pre=setPheno(Genomic2cycle12F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic2cycle12Pre, traits=1, use='pheno')
Genomic2cycle12Pre=setEBV(Genomic2cycle12Pre, GSModel); Genomic2cycle12PreS=selectInd(Genomic2cycle12Pre, nInd=300, use='ebv')
Genomic2cycle12Adv=setPheno(Genomic2cycle12PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic2cycle12Adv, traits=1, use='pheno')
Genomic2cycle12Adv=setEBV(Genomic2cycle12Adv, GSModel); Genomic2cycle12AdvS=selectInd(Genomic2cycle12Adv, nInd=50, use='ebv')
Genomic2cycle12Eli=setPheno(Genomic2cycle12AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic2cycle12Eli, traits=1, use='pheno')
Genomic2cycle12Eli=setEBV(Genomic2cycle12Eli, GSModel); Genomic2cycle12EliS=selectInd(Genomic2cycle12Eli, nInd=10, use='ebv')
Genomic2cycle12Eli1=selectInd(Genomic2cycle12Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic2cycle11PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic2cycle12CS=selectInd(F3, 150, use='ebv')
Genomic2cycle12F3Acc=cor(F3@gv, F3@ebv); Genomic2cycle12F3Acc
mean(Par@gv); mean(Genomic2cycle12CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic2cycle12F3, Genomic2cycle12Pre, Genomic2cycle12EliS)
names(datasets)=c('Genomic2cycle12F3', 'Genomic2cycle12Pre', 'Genomic2cycle12EliS')
Genomic2cycle12=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic2cycle12)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic2cycle12
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic2cycle12PreQtl=al(Genomic2cycle12Pre, SP); Genomic2cycle12EliSQtl= al(Genomic2cycle12EliS, SP)

#Genomic2cycle13
#Training, parentals and F1s
Tr=list(Genomic2cycle11PreS, Genomic2cycle10PreS); Tr=mergePops(Tr)
Par1=Genomic2cycle11PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic2cycle11CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic2cycle13Par=setEBV(Par, GSModel); Genomic2cycle13F1=randCross(Genomic2cycle13Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic2cycle13F2=self(Genomic2cycle13F1, nProgeny=1); Genomic2cycle13F3=self(Genomic2cycle13F2, nProgeny=25)
Genomic2cycle13F3=setEBV(Genomic2cycle13F3, GSModel); cor(Genomic2cycle13F3@gv, Genomic2cycle13F3@ebv)
Genomic2cycle13F3S=selectInd(Genomic2cycle13F3, nInd=6000, use='ebv'); mean(Genomic2cycle13F3@gv); mean(Genomic2cycle13F3S@gv)
Genomic2cycle13F4=self(Genomic2cycle13F3S, nProgeny=1); Genomic2cycle13F5=self(Genomic2cycle13F4, nProgeny=1)

#Trials 
Genomic2cycle13Pre=setPheno(Genomic2cycle13F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic2cycle13Pre, traits=1, use='pheno')
Genomic2cycle13Pre=setEBV(Genomic2cycle13Pre, GSModel); Genomic2cycle13PreS=selectInd(Genomic2cycle13Pre, nInd=300, use='ebv')
Genomic2cycle13Adv=setPheno(Genomic2cycle13PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic2cycle13Adv, traits=1, use='pheno')
Genomic2cycle13Adv=setEBV(Genomic2cycle13Adv, GSModel); Genomic2cycle13AdvS=selectInd(Genomic2cycle13Adv, nInd=50, use='ebv')
Genomic2cycle13Eli=setPheno(Genomic2cycle13AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic2cycle13Eli, traits=1, use='pheno')
Genomic2cycle13Eli=setEBV(Genomic2cycle13Eli, GSModel); Genomic2cycle13EliS=selectInd(Genomic2cycle13Eli, nInd=10, use='ebv')
Genomic2cycle13Eli1=selectInd(Genomic2cycle13Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic2cycle12PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic2cycle13CS=selectInd(F3, 150, use='ebv')
Genomic2cycle13F3Acc=cor(F3@gv, F3@ebv); Genomic2cycle13F3Acc
mean(Par@gv); mean(Genomic2cycle13CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic2cycle13F3, Genomic2cycle13Pre, Genomic2cycle13EliS)
names(datasets)=c('Genomic2cycle13F3', 'Genomic2cycle13Pre', 'Genomic2cycle13EliS')
Genomic2cycle13=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic2cycle13)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic2cycle13
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic2cycle13PreQtl=al(Genomic2cycle13Pre, SP); Genomic2cycle13EliSQtl= al(Genomic2cycle13EliS, SP)

#Genomic2cycle14
#Training, parentals and F1s
Tr=list(Genomic2cycle12PreS, Genomic2cycle11PreS); Tr=mergePops(Tr)
Par1=Genomic2cycle12PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic2cycle12CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic2cycle14Par=setEBV(Par, GSModel); Genomic2cycle14F1=randCross(Genomic2cycle14Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic2cycle14F2=self(Genomic2cycle14F1, nProgeny=1); Genomic2cycle14F3=self(Genomic2cycle14F2, nProgeny=25)
Genomic2cycle14F3=setEBV(Genomic2cycle14F3, GSModel); cor(Genomic2cycle14F3@gv, Genomic2cycle14F3@ebv)
Genomic2cycle14F3S=selectInd(Genomic2cycle14F3, nInd=6000, use='ebv'); mean(Genomic2cycle14F3@gv); mean(Genomic2cycle14F3S@gv)
Genomic2cycle14F4=self(Genomic2cycle14F3S, nProgeny=1); Genomic2cycle14F5=self(Genomic2cycle14F4, nProgeny=1)

#Trials 
Genomic2cycle14Pre=setPheno(Genomic2cycle14F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic2cycle14Pre, traits=1, use='pheno')
Genomic2cycle14Pre=setEBV(Genomic2cycle14Pre, GSModel); Genomic2cycle14PreS=selectInd(Genomic2cycle14Pre, nInd=300, use='ebv')
Genomic2cycle14Adv=setPheno(Genomic2cycle14PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic2cycle14Adv, traits=1, use='pheno')
Genomic2cycle14Adv=setEBV(Genomic2cycle14Adv, GSModel); Genomic2cycle14AdvS=selectInd(Genomic2cycle14Adv, nInd=50, use='ebv')
Genomic2cycle14Eli=setPheno(Genomic2cycle14AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic2cycle14Eli, traits=1, use='pheno')
Genomic2cycle14Eli=setEBV(Genomic2cycle14Eli, GSModel); Genomic2cycle14EliS=selectInd(Genomic2cycle14Eli, nInd=10, use='ebv')
Genomic2cycle14Eli1=selectInd(Genomic2cycle14Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic2cycle13PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic2cycle14CS=selectInd(F3, 150, use='ebv')
Genomic2cycle14F3Acc=cor(F3@gv, F3@ebv); Genomic2cycle14F3Acc
mean(Par@gv); mean(Genomic2cycle14CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic2cycle14F3, Genomic2cycle14Pre, Genomic2cycle14EliS)
names(datasets)=c('Genomic2cycle14F3', 'Genomic2cycle14Pre', 'Genomic2cycle14EliS')
Genomic2cycle14=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic2cycle14)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic2cycle14
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic2cycle14PreQtl=al(Genomic2cycle14Pre, SP); Genomic2cycle14EliSQtl= al(Genomic2cycle14EliS, SP)

#Genomic2cycle15
#Training, parentals and F1s
Tr=list(Genomic2cycle13PreS, Genomic2cycle12PreS); Tr=mergePops(Tr)
Par1=Genomic2cycle13PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic2cycle13CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic2cycle15Par=setEBV(Par, GSModel); Genomic2cycle15F1=randCross(Genomic2cycle15Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic2cycle15F2=self(Genomic2cycle15F1, nProgeny=1); Genomic2cycle15F3=self(Genomic2cycle15F2, nProgeny=25)
Genomic2cycle15F3=setEBV(Genomic2cycle15F3, GSModel); cor(Genomic2cycle15F3@gv, Genomic2cycle15F3@ebv)
Genomic2cycle15F3S=selectInd(Genomic2cycle15F3, nInd=6000, use='ebv'); mean(Genomic2cycle15F3@gv); mean(Genomic2cycle15F3S@gv)
Genomic2cycle15F4=self(Genomic2cycle15F3S, nProgeny=1); Genomic2cycle15F5=self(Genomic2cycle15F4, nProgeny=1)

#Trials 
Genomic2cycle15Pre=setPheno(Genomic2cycle15F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic2cycle15Pre, traits=1, use='pheno')
Genomic2cycle15Pre=setEBV(Genomic2cycle15Pre, GSModel); Genomic2cycle15PreS=selectInd(Genomic2cycle15Pre, nInd=300, use='ebv')
Genomic2cycle15Adv=setPheno(Genomic2cycle15PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic2cycle15Adv, traits=1, use='pheno')
Genomic2cycle15Adv=setEBV(Genomic2cycle15Adv, GSModel); Genomic2cycle15AdvS=selectInd(Genomic2cycle15Adv, nInd=50, use='ebv')
Genomic2cycle15Eli=setPheno(Genomic2cycle15AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic2cycle15Eli, traits=1, use='pheno')
Genomic2cycle15Eli=setEBV(Genomic2cycle15Eli, GSModel); Genomic2cycle15EliS=selectInd(Genomic2cycle15Eli, nInd=10, use='ebv')
Genomic2cycle15Eli1=selectInd(Genomic2cycle15Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic2cycle14PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic2cycle15CS=selectInd(F3, 150, use='ebv')
Genomic2cycle15F3Acc=cor(F3@gv, F3@ebv); Genomic2cycle15F3Acc
mean(Par@gv); mean(Genomic2cycle15CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic2cycle15F3, Genomic2cycle15Pre, Genomic2cycle15EliS)
names(datasets)=c('Genomic2cycle15F3', 'Genomic2cycle15Pre', 'Genomic2cycle15EliS')
Genomic2cycle15=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic2cycle15)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic2cycle15
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic2cycle15PreQtl=al(Genomic2cycle15Pre, SP); Genomic2cycle15EliSQtl= al(Genomic2cycle15EliS, SP)

################################################################################
#Genomic2 exporting estimates and Qtl

Genomic2cycle1to15PreQtl=cbind(Burnincycle1PreQtl, Burnincycle2PreQtl, Burnincycle3PreQtl, Burnincycle4PreQtl, Burnincycle5PreQtl, Genomic2cycle6PreQtl, Genomic2cycle7PreQtl, Genomic2cycle8PreQtl, Genomic2cycle9PreQtl, Genomic2cycle10PreQtl, Genomic2cycle11PreQtl, Genomic2cycle12PreQtl, Genomic2cycle13PreQtl, Genomic2cycle14PreQtl, Genomic2cycle15PreQtl); write.csv(Genomic2cycle1to15PreQtl, 'r1_Genomic2PreQtl.csv')

Genomic2cycle1to15EliSQtl=cbind(Burnincycle1EliSQtl, Burnincycle2EliSQtl, Burnincycle3EliSQtl, Burnincycle4EliSQtl, Burnincycle5EliSQtl, Genomic2cycle6EliSQtl, Genomic2cycle7EliSQtl, Genomic2cycle8EliSQtl, Genomic2cycle9EliSQtl, Genomic2cycle10EliSQtl, Genomic2cycle11EliSQtl, Genomic2cycle12EliSQtl, Genomic2cycle13EliSQtl, Genomic2cycle14EliSQtl, Genomic2cycle15EliSQtl); write.csv(Genomic2cycle1to15EliSQtl, 'r1_Genomic2EliSQtl.csv')

Genomic2cycle1to15=rbind(Burnincycle1, Burnincycle2, Burnincycle3, Burnincycle4, Burnincycle5, Genomic2cycle6, Genomic2cycle7, Genomic2cycle8, Genomic2cycle9, Genomic2cycle10, Genomic2cycle11, Genomic2cycle12, Genomic2cycle13, Genomic2cycle14, Genomic2cycle15); write.csv(Genomic2cycle1to15, 'r1_Genomic2.csv')

################################################################################
#Genomic3cycle6
#Training, parentals and F1s
Tr=list(Burnincycle3PreS, Burnincycle2PreS); Tr=mergePops(Tr)
Par1=Burnincycle3PreS; Par1=selectInd(Par1, nInd=150, use='pheno')
Par2=Burnincycle4CS; Par2=selectInd(Par2, nInd=150, use='ebv')
Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic3cycle6Par=setEBV(Par, GSModel); Genomic3cycle6F1=randCross(Genomic3cycle6Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic3cycle6F2=self(Genomic3cycle6F1, nProgeny=1); Genomic3cycle6F3=self(Genomic3cycle6F2, nProgeny=25)
Genomic3cycle6F3=setEBV(Genomic3cycle6F3, GSModel); cor(Genomic3cycle6F3@gv, Genomic3cycle6F3@ebv)
Genomic3cycle6F3S=selectInd(Genomic3cycle6F3, nInd=6000, use='ebv'); mean(Genomic3cycle6F3@gv); mean(Genomic3cycle6F3S@gv)
Genomic3cycle6F4=self(Genomic3cycle6F3S, nProgeny=1); Genomic3cycle6F5=self(Genomic3cycle6F4, nProgeny=1)

#Trials 
Genomic3cycle6Pre=setPheno(Genomic3cycle6F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic3cycle6Pre, traits=1, use='pheno')
Genomic3cycle6Pre=setEBV(Genomic3cycle6Pre, GSModel); Genomic3cycle6PreS=selectInd(Genomic3cycle6Pre, nInd=300, use='ebv')
Genomic3cycle6Adv=setPheno(Genomic3cycle6PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic3cycle6Adv, traits=1, use='pheno')
Genomic3cycle6Adv=setEBV(Genomic3cycle6Adv, GSModel); Genomic3cycle6AdvS=selectInd(Genomic3cycle6Adv, nInd=50, use='ebv')
Genomic3cycle6Eli=setPheno(Genomic3cycle6AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic3cycle6Eli, traits=1, use='pheno')
Genomic3cycle6Eli=setEBV(Genomic3cycle6Eli, GSModel); Genomic3cycle6EliS=selectInd(Genomic3cycle6Eli, nInd=10, use='ebv')
Genomic3cycle6Eli1=selectInd(Genomic3cycle6Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Burnincycle4PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle6CS=selectInd(F3, 150, use='ebv')
Genomic3cycle6F3Acc=cor(F3@gv, F3@ebv); Genomic3cycle6F3Acc
mean(Par@gv); mean(Genomic3cycle6CS@gv)

ParRound2=Genomic3cycle6CS
F1=randCross(ParRound2, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12); F2=setPheno(F2, h2=0.10, reps=2)
GSModel=fastRRBLUP(F2, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle6CSRound2=selectInd(F3, 150, use='ebv')
Genomic3cycle6F3AccRound2=cor(F3@gv, F3@ebv); Genomic3cycle6F3AccRound2
mean(ParRound2@gv); mean(Genomic3cycle6CSRound2@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic3cycle6F3, Genomic3cycle6Pre, Genomic3cycle6EliS)
names(datasets)=c('Genomic3cycle6F3', 'Genomic3cycle6Pre', 'Genomic3cycle6EliS')
Genomic3cycle6=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic3cycle6)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic3cycle6
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic3cycle6PreQtl=al(Genomic3cycle6Pre, SP); Genomic3cycle6EliSQtl= al(Genomic3cycle6EliS, SP)

#Genomic3cycle7
#Training, parentals and F1s
Tr=list(Burnincycle4PreS, Burnincycle3PreS, Burnincycle2PreS); Tr=mergePops(Tr)
Par1=Burnincycle4PreS; Par1=selectInd(Par1, nInd=100, use='pheno')
Par2=Burnincycle5CS; Par2=selectInd(Par2, nInd=100, use='ebv')
Par3=Burnincycle4CSRound2; Par3=selectInd(Par3, nInd=100, use='ebv')
Par=list(Par1, Par2, Par3); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic3cycle7Par=setEBV(Par, GSModel); Genomic3cycle7F1=randCross(Genomic3cycle7Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic3cycle7F2=self(Genomic3cycle7F1, nProgeny=1); Genomic3cycle7F3=self(Genomic3cycle7F2, nProgeny=25)
Genomic3cycle7F3=setEBV(Genomic3cycle7F3, GSModel); cor(Genomic3cycle7F3@gv, Genomic3cycle7F3@ebv)
Genomic3cycle7F3S=selectInd(Genomic3cycle7F3, nInd=6000, use='ebv'); mean(Genomic3cycle7F3@gv); mean(Genomic3cycle7F3S@gv)
Genomic3cycle7F4=self(Genomic3cycle7F3S, nProgeny=1); Genomic3cycle7F5=self(Genomic3cycle7F4, nProgeny=1)

#Trials 
Genomic3cycle7Pre=setPheno(Genomic3cycle7F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic3cycle7Pre, traits=1, use='pheno')
Genomic3cycle7Pre=setEBV(Genomic3cycle7Pre, GSModel); Genomic3cycle7PreS=selectInd(Genomic3cycle7Pre, nInd=300, use='ebv')
Genomic3cycle7Adv=setPheno(Genomic3cycle7PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic3cycle7Adv, traits=1, use='pheno')
Genomic3cycle7Adv=setEBV(Genomic3cycle7Adv, GSModel); Genomic3cycle7AdvS=selectInd(Genomic3cycle7Adv, nInd=50, use='ebv')
Genomic3cycle7Eli=setPheno(Genomic3cycle7AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic3cycle7Eli, traits=1, use='pheno')
Genomic3cycle7Eli=setEBV(Genomic3cycle7Eli, GSModel); Genomic3cycle7EliS=selectInd(Genomic3cycle7Eli, nInd=10, use='ebv')
Genomic3cycle7Eli1=selectInd(Genomic3cycle7Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic3cycle6PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle7CS=selectInd(F3, 150, use='ebv')
Genomic3cycle7F3Acc=cor(F3@gv, F3@ebv); Genomic3cycle7F3Acc
mean(Par@gv); mean(Genomic3cycle7CS@gv)

ParRound2=Genomic3cycle7CS
F1=randCross(ParRound2, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12); F2=setPheno(F2, h2=0.10, reps=2)
GSModel=fastRRBLUP(F2, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle7CSRound2=selectInd(F3, 150, use='ebv')
Genomic3cycle7F3AccRound2=cor(F3@gv, F3@ebv); Genomic3cycle7F3AccRound2
mean(ParRound2@gv); mean(Genomic3cycle7CSRound2@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic3cycle7F3, Genomic3cycle7Pre, Genomic3cycle7EliS)
names(datasets)=c('Genomic3cycle7F3', 'Genomic3cycle7Pre', 'Genomic3cycle7EliS')
Genomic3cycle7=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic3cycle7)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic3cycle7
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic3cycle7PreQtl=al(Genomic3cycle7Pre, SP); Genomic3cycle7EliSQtl= al(Genomic3cycle7EliS, SP)

#Genomic3cycle8
#Training, parentals and F1s
Tr=list(Genomic3cycle6PreS, Burnincycle4PreS, Burnincycle3PreS); Tr=mergePops(Tr)
Par1=Genomic3cycle6PreS; Par1=selectInd(Par1, nInd=100, use='pheno')
Par2=Genomic3cycle6CS; Par2=selectInd(Par2, nInd=100, use='ebv')
Par3=Burnincycle5CSRound2; Par3=selectInd(Par3, nInd=100, use='ebv')
Par=list(Par1, Par2, Par3); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic3cycle8Par=setEBV(Par, GSModel); Genomic3cycle8F1=randCross(Genomic3cycle8Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic3cycle8F2=self(Genomic3cycle8F1, nProgeny=1); Genomic3cycle8F3=self(Genomic3cycle8F2, nProgeny=25)
Genomic3cycle8F3=setEBV(Genomic3cycle8F3, GSModel); cor(Genomic3cycle8F3@gv, Genomic3cycle8F3@ebv)
Genomic3cycle8F3S=selectInd(Genomic3cycle8F3, nInd=6000, use='ebv'); mean(Genomic3cycle8F3@gv); mean(Genomic3cycle8F3S@gv)
Genomic3cycle8F4=self(Genomic3cycle8F3S, nProgeny=1); Genomic3cycle8F5=self(Genomic3cycle8F4, nProgeny=1)

#Trials 
Genomic3cycle8Pre=setPheno(Genomic3cycle8F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic3cycle8Pre, traits=1, use='pheno')
Genomic3cycle8Pre=setEBV(Genomic3cycle8Pre, GSModel); Genomic3cycle8PreS=selectInd(Genomic3cycle8Pre, nInd=300, use='ebv')
Genomic3cycle8Adv=setPheno(Genomic3cycle8PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic3cycle8Adv, traits=1, use='pheno')
Genomic3cycle8Adv=setEBV(Genomic3cycle8Adv, GSModel); Genomic3cycle8AdvS=selectInd(Genomic3cycle8Adv, nInd=50, use='ebv')
Genomic3cycle8Eli=setPheno(Genomic3cycle8AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic3cycle8Eli, traits=1, use='pheno')
Genomic3cycle8Eli=setEBV(Genomic3cycle8Eli, GSModel); Genomic3cycle8EliS=selectInd(Genomic3cycle8Eli, nInd=10, use='ebv')
Genomic3cycle8Eli1=selectInd(Genomic3cycle8Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic3cycle7PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle8CS=selectInd(F3, 150, use='ebv')
Genomic3cycle8F3Acc=cor(F3@gv, F3@ebv); Genomic3cycle8F3Acc
mean(Par@gv); mean(Genomic3cycle8CS@gv)

ParRound2=Genomic3cycle8CS
F1=randCross(ParRound2, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12); F2=setPheno(F2, h2=0.10, reps=2)
GSModel=fastRRBLUP(F2, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle8CSRound2=selectInd(F3, 150, use='ebv')
Genomic3cycle8F3AccRound2=cor(F3@gv, F3@ebv); Genomic3cycle8F3AccRound2
mean(ParRound2@gv); mean(Genomic3cycle8CSRound2@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic3cycle8F3, Genomic3cycle8Pre, Genomic3cycle8EliS)
names(datasets)=c('Genomic3cycle8F3', 'Genomic3cycle8Pre', 'Genomic3cycle8EliS')
Genomic3cycle8=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic3cycle8)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic3cycle8
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic3cycle8PreQtl=al(Genomic3cycle8Pre, SP); Genomic3cycle8EliSQtl= al(Genomic3cycle8EliS, SP)

#Genomic3cycle9
#Training, parentals and F1s
Tr=list(Genomic3cycle7PreS, Genomic3cycle6PreS, Burnincycle4PreS); Tr=mergePops(Tr)
Par1=Genomic3cycle7PreS; Par1=selectInd(Par1, nInd=100, use='pheno')
Par2=Genomic3cycle7CS; Par2=selectInd(Par2, nInd=100, use='ebv')
Par3= Genomic3cycle6CSRound2; Par3=selectInd(Par3, nInd=100, use='ebv')
Par=list(Par1, Par2, Par3); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic3cycle9Par=setEBV(Par, GSModel); Genomic3cycle9F1=randCross(Genomic3cycle9Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic3cycle9F2=self(Genomic3cycle9F1, nProgeny=1); Genomic3cycle9F3=self(Genomic3cycle9F2, nProgeny=25)
Genomic3cycle9F3=setEBV(Genomic3cycle9F3, GSModel); cor(Genomic3cycle9F3@gv, Genomic3cycle9F3@ebv)
Genomic3cycle9F3S=selectInd(Genomic3cycle9F3, nInd=6000, use='ebv'); mean(Genomic3cycle9F3@gv); mean(Genomic3cycle9F3S@gv)
Genomic3cycle9F4=self(Genomic3cycle9F3S, nProgeny=1); Genomic3cycle9F5=self(Genomic3cycle9F4, nProgeny=1)

#Trials 
Genomic3cycle9Pre=setPheno(Genomic3cycle9F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic3cycle9Pre, traits=1, use='pheno')
Genomic3cycle9Pre=setEBV(Genomic3cycle9Pre, GSModel); Genomic3cycle9PreS=selectInd(Genomic3cycle9Pre, nInd=300, use='ebv')
Genomic3cycle9Adv=setPheno(Genomic3cycle9PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic3cycle9Adv, traits=1, use='pheno')
Genomic3cycle9Adv=setEBV(Genomic3cycle9Adv, GSModel); Genomic3cycle9AdvS=selectInd(Genomic3cycle9Adv, nInd=50, use='ebv')
Genomic3cycle9Eli=setPheno(Genomic3cycle9AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic3cycle9Eli, traits=1, use='pheno')
Genomic3cycle9Eli=setEBV(Genomic3cycle9Eli, GSModel); Genomic3cycle9EliS=selectInd(Genomic3cycle9Eli, nInd=10, use='ebv')
Genomic3cycle9Eli1=selectInd(Genomic3cycle9Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic3cycle8PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle9CS=selectInd(F3, 150, use='ebv')
Genomic3cycle9F3Acc=cor(F3@gv, F3@ebv); Genomic3cycle9F3Acc
mean(Par@gv); mean(Genomic3cycle9CS@gv)

ParRound2=Genomic3cycle9CS
F1=randCross(ParRound2, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12); F2=setPheno(F2, h2=0.10, reps=2)
GSModel=fastRRBLUP(F2, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle9CSRound2=selectInd(F3, 150, use='ebv')
Genomic3cycle9F3AccRound2=cor(F3@gv, F3@ebv); Genomic3cycle9F3AccRound2
mean(ParRound2@gv); mean(Genomic3cycle9CSRound2@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic3cycle9F3, Genomic3cycle9Pre, Genomic3cycle9EliS)
names(datasets)=c('Genomic3cycle9F3', 'Genomic3cycle9Pre', 'Genomic3cycle9EliS')
Genomic3cycle9=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic3cycle9)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic3cycle9
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic3cycle9PreQtl=al(Genomic3cycle9Pre, SP); Genomic3cycle9EliSQtl= al(Genomic3cycle9EliS, SP)

#Genomic3cycle10
#Training, parentals and F1s
Tr=list(Genomic3cycle8PreS, Genomic3cycle7PreS, Genomic3cycle6PreS); Tr=mergePops(Tr)
Par1=Genomic3cycle8PreS; Par1=selectInd(Par1, nInd=100, use='pheno')
Par2=Genomic3cycle8CS; Par2=selectInd(Par2, nInd=100, use='ebv')
Par3= Genomic3cycle7CSRound2; Par3=selectInd(Par3, nInd=100, use='ebv')
Par=list(Par1, Par2, Par3); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic3cycle10Par=setEBV(Par, GSModel); Genomic3cycle10F1=randCross(Genomic3cycle10Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic3cycle10F2=self(Genomic3cycle10F1, nProgeny=1); Genomic3cycle10F3=self(Genomic3cycle10F2, nProgeny=25)
Genomic3cycle10F3=setEBV(Genomic3cycle10F3, GSModel); cor(Genomic3cycle10F3@gv, Genomic3cycle10F3@ebv)
Genomic3cycle10F3S=selectInd(Genomic3cycle10F3, nInd=6000, use='ebv'); mean(Genomic3cycle10F3@gv); mean(Genomic3cycle10F3S@gv)
Genomic3cycle10F4=self(Genomic3cycle10F3S, nProgeny=1); Genomic3cycle10F5=self(Genomic3cycle10F4, nProgeny=1)

#Trials 
Genomic3cycle10Pre=setPheno(Genomic3cycle10F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic3cycle10Pre, traits=1, use='pheno')
Genomic3cycle10Pre=setEBV(Genomic3cycle10Pre, GSModel); Genomic3cycle10PreS=selectInd(Genomic3cycle10Pre, nInd=300, use='ebv')
Genomic3cycle10Adv=setPheno(Genomic3cycle10PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic3cycle10Adv, traits=1, use='pheno')
Genomic3cycle10Adv=setEBV(Genomic3cycle10Adv, GSModel); Genomic3cycle10AdvS=selectInd(Genomic3cycle10Adv, nInd=50, use='ebv')
Genomic3cycle10Eli=setPheno(Genomic3cycle10AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic3cycle10Eli, traits=1, use='pheno')
Genomic3cycle10Eli=setEBV(Genomic3cycle10Eli, GSModel); Genomic3cycle10EliS=selectInd(Genomic3cycle10Eli, nInd=10, use='ebv')
Genomic3cycle10Eli1=selectInd(Genomic3cycle10Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic3cycle9PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle10CS=selectInd(F3, 150, use='ebv')
Genomic3cycle10F3Acc=cor(F3@gv, F3@ebv); Genomic3cycle10F3Acc
mean(Par@gv); mean(Genomic3cycle10CS@gv)

ParRound2=Genomic3cycle10CS
F1=randCross(ParRound2, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12); F2=setPheno(F2, h2=0.10, reps=2)
GSModel=fastRRBLUP(F2, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle10CSRound2=selectInd(F3, 150, use='ebv')
Genomic3cycle10F3AccRound2=cor(F3@gv, F3@ebv); Genomic3cycle10F3AccRound2
mean(ParRound2@gv); mean(Genomic3cycle10CSRound2@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic3cycle10F3, Genomic3cycle10Pre, Genomic3cycle10EliS)
names(datasets)=c('Genomic3cycle10F3', 'Genomic3cycle10Pre', 'Genomic3cycle10EliS')
Genomic3cycle10=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic3cycle10)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic3cycle10
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic3cycle10PreQtl=al(Genomic3cycle10Pre, SP); Genomic3cycle10EliSQtl= al(Genomic3cycle10EliS, SP)

#Genomic3cycle11
#Training, parentals and F1s
Tr=list(Genomic3cycle9PreS, Genomic3cycle8PreS, Genomic3cycle7PreS); Tr=mergePops(Tr)
Par1=Genomic3cycle9PreS; Par1=selectInd(Par1, nInd=100, use='pheno')
Par2=Genomic3cycle9CS; Par2=selectInd(Par2, nInd=100, use='ebv')
Par3= Genomic3cycle8CSRound2; Par3=selectInd(Par3, nInd=100, use='ebv')
Par=list(Par1, Par2, Par3); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic3cycle11Par=setEBV(Par, GSModel); Genomic3cycle11F1=randCross(Genomic3cycle11Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic3cycle11F2=self(Genomic3cycle11F1, nProgeny=1); Genomic3cycle11F3=self(Genomic3cycle11F2, nProgeny=25)
Genomic3cycle11F3=setEBV(Genomic3cycle11F3, GSModel); cor(Genomic3cycle11F3@gv, Genomic3cycle11F3@ebv)
Genomic3cycle11F3S=selectInd(Genomic3cycle11F3, nInd=6000, use='ebv'); mean(Genomic3cycle11F3@gv); mean(Genomic3cycle11F3S@gv)
Genomic3cycle11F4=self(Genomic3cycle11F3S, nProgeny=1); Genomic3cycle11F5=self(Genomic3cycle11F4, nProgeny=1)

#Trials 
Genomic3cycle11Pre=setPheno(Genomic3cycle11F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic3cycle11Pre, traits=1, use='pheno')
Genomic3cycle11Pre=setEBV(Genomic3cycle11Pre, GSModel); Genomic3cycle11PreS=selectInd(Genomic3cycle11Pre, nInd=300, use='ebv')
Genomic3cycle11Adv=setPheno(Genomic3cycle11PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic3cycle11Adv, traits=1, use='pheno')
Genomic3cycle11Adv=setEBV(Genomic3cycle11Adv, GSModel); Genomic3cycle11AdvS=selectInd(Genomic3cycle11Adv, nInd=50, use='ebv')
Genomic3cycle11Eli=setPheno(Genomic3cycle11AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic3cycle11Eli, traits=1, use='pheno')
Genomic3cycle11Eli=setEBV(Genomic3cycle11Eli, GSModel); Genomic3cycle11EliS=selectInd(Genomic3cycle11Eli, nInd=10, use='ebv')
Genomic3cycle11Eli1=selectInd(Genomic3cycle11Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic3cycle10PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle11CS=selectInd(F3, 150, use='ebv')
Genomic3cycle11F3Acc=cor(F3@gv, F3@ebv); Genomic3cycle11F3Acc
mean(Par@gv); mean(Genomic3cycle11CS@gv)

ParRound2=Genomic3cycle11CS
F1=randCross(ParRound2, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12); F2=setPheno(F2, h2=0.10, reps=2)
GSModel=fastRRBLUP(F2, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle11CSRound2=selectInd(F3, 150, use='ebv')
Genomic3cycle11F3AccRound2=cor(F3@gv, F3@ebv); Genomic3cycle11F3AccRound2
mean(ParRound2@gv); mean(Genomic3cycle11CSRound2@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic3cycle11F3, Genomic3cycle11Pre, Genomic3cycle11EliS)
names(datasets)=c('Genomic3cycle11F3', 'Genomic3cycle11Pre', 'Genomic3cycle11EliS')
Genomic3cycle11=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic3cycle11)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic3cycle11
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic3cycle11PreQtl=al(Genomic3cycle11Pre, SP); Genomic3cycle11EliSQtl= al(Genomic3cycle11EliS, SP)

#Genomic3cycle12
#Training, parentals and F1s
Tr=list(Genomic3cycle10PreS, Genomic3cycle9PreS, Genomic3cycle8PreS); Tr=mergePops(Tr)
Par1=Genomic3cycle10PreS; Par1=selectInd(Par1, nInd=100, use='pheno')
Par2=Genomic3cycle10CS; Par2=selectInd(Par2, nInd=100, use='ebv')
Par3= Genomic3cycle9CSRound2; Par3=selectInd(Par3, nInd=100, use='ebv')
Par=list(Par1, Par2, Par3); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic3cycle12Par=setEBV(Par, GSModel); Genomic3cycle12F1=randCross(Genomic3cycle12Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic3cycle12F2=self(Genomic3cycle12F1, nProgeny=1); Genomic3cycle12F3=self(Genomic3cycle12F2, nProgeny=25)
Genomic3cycle12F3=setEBV(Genomic3cycle12F3, GSModel); cor(Genomic3cycle12F3@gv, Genomic3cycle12F3@ebv)
Genomic3cycle12F3S=selectInd(Genomic3cycle12F3, nInd=6000, use='ebv'); mean(Genomic3cycle12F3@gv); mean(Genomic3cycle12F3S@gv)
Genomic3cycle12F4=self(Genomic3cycle12F3S, nProgeny=1); Genomic3cycle12F5=self(Genomic3cycle12F4, nProgeny=1)

#Trials 
Genomic3cycle12Pre=setPheno(Genomic3cycle12F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic3cycle12Pre, traits=1, use='pheno')
Genomic3cycle12Pre=setEBV(Genomic3cycle12Pre, GSModel); Genomic3cycle12PreS=selectInd(Genomic3cycle12Pre, nInd=300, use='ebv')
Genomic3cycle12Adv=setPheno(Genomic3cycle12PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic3cycle12Adv, traits=1, use='pheno')
Genomic3cycle12Adv=setEBV(Genomic3cycle12Adv, GSModel); Genomic3cycle12AdvS=selectInd(Genomic3cycle12Adv, nInd=50, use='ebv')
Genomic3cycle12Eli=setPheno(Genomic3cycle12AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic3cycle12Eli, traits=1, use='pheno')
Genomic3cycle12Eli=setEBV(Genomic3cycle12Eli, GSModel); Genomic3cycle12EliS=selectInd(Genomic3cycle12Eli, nInd=10, use='ebv')
Genomic3cycle12Eli1=selectInd(Genomic3cycle12Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic3cycle11PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle12CS=selectInd(F3, 150, use='ebv')
Genomic3cycle12F3Acc=cor(F3@gv, F3@ebv); Genomic3cycle12F3Acc
mean(Par@gv); mean(Genomic3cycle12CS@gv)

ParRound2=Genomic3cycle12CS
F1=randCross(ParRound2, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12); F2=setPheno(F2, h2=0.10, reps=2)
GSModel=fastRRBLUP(F2, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle12CSRound2=selectInd(F3, 150, use='ebv')
Genomic3cycle12F3AccRound2=cor(F3@gv, F3@ebv); Genomic3cycle12F3AccRound2
mean(ParRound2@gv); mean(Genomic3cycle12CSRound2@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic3cycle12F3, Genomic3cycle12Pre, Genomic3cycle12EliS)
names(datasets)=c('Genomic3cycle12F3', 'Genomic3cycle12Pre', 'Genomic3cycle12EliS')
Genomic3cycle12=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic3cycle12)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic3cycle12
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic3cycle12PreQtl=al(Genomic3cycle12Pre, SP); Genomic3cycle12EliSQtl= al(Genomic3cycle12EliS, SP)

#Genomic3cycle13
#Training, parentals and F1s
Tr=list(Genomic3cycle11PreS, Genomic3cycle10PreS, Genomic3cycle9PreS); Tr=mergePops(Tr)
Par1=Genomic3cycle11PreS; Par1=selectInd(Par1, nInd=100, use='pheno')
Par2=Genomic3cycle11CS; Par2=selectInd(Par2, nInd=100, use='ebv')
Par3= Genomic3cycle10CSRound2; Par3=selectInd(Par3, nInd=100, use='ebv')
Par=list(Par1, Par2, Par3); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic3cycle13Par=setEBV(Par, GSModel); Genomic3cycle13F1=randCross(Genomic3cycle13Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic3cycle13F2=self(Genomic3cycle13F1, nProgeny=1); Genomic3cycle13F3=self(Genomic3cycle13F2, nProgeny=25)
Genomic3cycle13F3=setEBV(Genomic3cycle13F3, GSModel); cor(Genomic3cycle13F3@gv, Genomic3cycle13F3@ebv)
Genomic3cycle13F3S=selectInd(Genomic3cycle13F3, nInd=6000, use='ebv'); mean(Genomic3cycle13F3@gv); mean(Genomic3cycle13F3S@gv)
Genomic3cycle13F4=self(Genomic3cycle13F3S, nProgeny=1); Genomic3cycle13F5=self(Genomic3cycle13F4, nProgeny=1)

#Trials 
Genomic3cycle13Pre=setPheno(Genomic3cycle13F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic3cycle13Pre, traits=1, use='pheno')
Genomic3cycle13Pre=setEBV(Genomic3cycle13Pre, GSModel); Genomic3cycle13PreS=selectInd(Genomic3cycle13Pre, nInd=300, use='ebv')
Genomic3cycle13Adv=setPheno(Genomic3cycle13PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic3cycle13Adv, traits=1, use='pheno')
Genomic3cycle13Adv=setEBV(Genomic3cycle13Adv, GSModel); Genomic3cycle13AdvS=selectInd(Genomic3cycle13Adv, nInd=50, use='ebv')
Genomic3cycle13Eli=setPheno(Genomic3cycle13AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic3cycle13Eli, traits=1, use='pheno')
Genomic3cycle13Eli=setEBV(Genomic3cycle13Eli, GSModel); Genomic3cycle13EliS=selectInd(Genomic3cycle13Eli, nInd=10, use='ebv')
Genomic3cycle13Eli1=selectInd(Genomic3cycle13Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic3cycle12PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle13CS=selectInd(F3, 150, use='ebv')
Genomic3cycle13F3Acc=cor(F3@gv, F3@ebv); Genomic3cycle13F3Acc
mean(Par@gv); mean(Genomic3cycle13CS@gv)

ParRound2=Genomic3cycle13CS
F1=randCross(ParRound2, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12); F2=setPheno(F2, h2=0.10, reps=2)
GSModel=fastRRBLUP(F2, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle13CSRound2=selectInd(F3, 150, use='ebv')
Genomic3cycle13F3AccRound2=cor(F3@gv, F3@ebv); Genomic3cycle13F3AccRound2
mean(ParRound2@gv); mean(Genomic3cycle13CSRound2@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic3cycle13F3, Genomic3cycle13Pre, Genomic3cycle13EliS)
names(datasets)=c('Genomic3cycle13F3', 'Genomic3cycle13Pre', 'Genomic3cycle13EliS')
Genomic3cycle13=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic3cycle13)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic3cycle13
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic3cycle13PreQtl=al(Genomic3cycle13Pre, SP); Genomic3cycle13EliSQtl= al(Genomic3cycle13EliS, SP)

#Genomic3cycle14
#Training, parentals and F1s
Tr=list(Genomic3cycle12PreS, Genomic3cycle11PreS, Genomic3cycle10PreS); Tr=mergePops(Tr)
Par1=Genomic3cycle12PreS; Par1=selectInd(Par1, nInd=100, use='pheno')
Par2=Genomic3cycle12CS; Par2=selectInd(Par2, nInd=100, use='ebv')
Par3= Genomic3cycle11CSRound2; Par3=selectInd(Par3, nInd=100, use='ebv')
Par=list(Par1, Par2, Par3); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic3cycle14Par=setEBV(Par, GSModel); Genomic3cycle14F1=randCross(Genomic3cycle14Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic3cycle14F2=self(Genomic3cycle14F1, nProgeny=1); Genomic3cycle14F3=self(Genomic3cycle14F2, nProgeny=25)
Genomic3cycle14F3=setEBV(Genomic3cycle14F3, GSModel); cor(Genomic3cycle14F3@gv, Genomic3cycle14F3@ebv)
Genomic3cycle14F3S=selectInd(Genomic3cycle14F3, nInd=6000, use='ebv'); mean(Genomic3cycle14F3@gv); mean(Genomic3cycle14F3S@gv)
Genomic3cycle14F4=self(Genomic3cycle14F3S, nProgeny=1); Genomic3cycle14F5=self(Genomic3cycle14F4, nProgeny=1)

#Trials 
Genomic3cycle14Pre=setPheno(Genomic3cycle14F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic3cycle14Pre, traits=1, use='pheno')
Genomic3cycle14Pre=setEBV(Genomic3cycle14Pre, GSModel); Genomic3cycle14PreS=selectInd(Genomic3cycle14Pre, nInd=300, use='ebv')
Genomic3cycle14Adv=setPheno(Genomic3cycle14PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic3cycle14Adv, traits=1, use='pheno')
Genomic3cycle14Adv=setEBV(Genomic3cycle14Adv, GSModel); Genomic3cycle14AdvS=selectInd(Genomic3cycle14Adv, nInd=50, use='ebv')
Genomic3cycle14Eli=setPheno(Genomic3cycle14AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic3cycle14Eli, traits=1, use='pheno')
Genomic3cycle14Eli=setEBV(Genomic3cycle14Eli, GSModel); Genomic3cycle14EliS=selectInd(Genomic3cycle14Eli, nInd=10, use='ebv')
Genomic3cycle14Eli1=selectInd(Genomic3cycle14Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic3cycle13PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle14CS=selectInd(F3, 150, use='ebv')
Genomic3cycle14F3Acc=cor(F3@gv, F3@ebv); Genomic3cycle14F3Acc
mean(Par@gv); mean(Genomic3cycle14CS@gv)

ParRound2=Genomic3cycle14CS
F1=randCross(ParRound2, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12); F2=setPheno(F2, h2=0.10, reps=2)
GSModel=fastRRBLUP(F2, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle14CSRound2=selectInd(F3, 150, use='ebv')
Genomic3cycle14F3AccRound2=cor(F3@gv, F3@ebv); Genomic3cycle14F3AccRound2
mean(ParRound2@gv); mean(Genomic3cycle14CSRound2@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic3cycle14F3, Genomic3cycle14Pre, Genomic3cycle14EliS)
names(datasets)=c('Genomic3cycle14F3', 'Genomic3cycle14Pre', 'Genomic3cycle14EliS')
Genomic3cycle14=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic3cycle14)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic3cycle14
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic3cycle14PreQtl=al(Genomic3cycle14Pre, SP); Genomic3cycle14EliSQtl= al(Genomic3cycle14EliS, SP)

#Genomic3cycle15
#Training, parentals and F1s
Tr=list(Genomic3cycle13PreS, Genomic3cycle12PreS, Genomic3cycle11PreS); Tr=mergePops(Tr)
Par1=Genomic3cycle13PreS; Par1=selectInd(Par1, nInd=100, use='pheno')
Par2=Genomic3cycle13CS; Par2=selectInd(Par2, nInd=100, use='ebv')
Par3= Genomic3cycle12CSRound2; Par3=selectInd(Par3, nInd=100, use='ebv')
Par=list(Par1, Par2, Par3); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic3cycle15Par=setEBV(Par, GSModel); Genomic3cycle15F1=randCross(Genomic3cycle15Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic3cycle15F2=self(Genomic3cycle15F1, nProgeny=1); Genomic3cycle15F3=self(Genomic3cycle15F2, nProgeny=25)
Genomic3cycle15F3=setEBV(Genomic3cycle15F3, GSModel); cor(Genomic3cycle15F3@gv, Genomic3cycle15F3@ebv)
Genomic3cycle15F3S=selectInd(Genomic3cycle15F3, nInd=6000, use='ebv'); mean(Genomic3cycle15F3@gv); mean(Genomic3cycle15F3S@gv)
Genomic3cycle15F4=self(Genomic3cycle15F3S, nProgeny=1); Genomic3cycle15F5=self(Genomic3cycle15F4, nProgeny=1)

#Trials 
Genomic3cycle15Pre=setPheno(Genomic3cycle15F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic3cycle15Pre, traits=1, use='pheno')
Genomic3cycle15Pre=setEBV(Genomic3cycle15Pre, GSModel); Genomic3cycle15PreS=selectInd(Genomic3cycle15Pre, nInd=300, use='ebv')
Genomic3cycle15Adv=setPheno(Genomic3cycle15PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic3cycle15Adv, traits=1, use='pheno')
Genomic3cycle15Adv=setEBV(Genomic3cycle15Adv, GSModel); Genomic3cycle15AdvS=selectInd(Genomic3cycle15Adv, nInd=50, use='ebv')
Genomic3cycle15Eli=setPheno(Genomic3cycle15AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic3cycle15Eli, traits=1, use='pheno')
Genomic3cycle15Eli=setEBV(Genomic3cycle15Eli, GSModel); Genomic3cycle15EliS=selectInd(Genomic3cycle15Eli, nInd=10, use='ebv')
Genomic3cycle15Eli1=selectInd(Genomic3cycle15Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic3cycle14PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle15CS=selectInd(F3, 150, use='ebv')
Genomic3cycle15F3Acc=cor(F3@gv, F3@ebv); Genomic3cycle15F3Acc
mean(Par@gv); mean(Genomic3cycle15CS@gv)

ParRound2=Genomic3cycle15CS
F1=randCross(ParRound2, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12); F2=setPheno(F2, h2=0.10, reps=2)
GSModel=fastRRBLUP(F2, use='pheno'); F3=setEBV(F3, GSModel); Genomic3cycle15CSRound2=selectInd(F3, 150, use='ebv')
Genomic3cycle15F3AccRound2=cor(F3@gv, F3@ebv); Genomic3cycle15F3AccRound2
mean(ParRound2@gv); mean(Genomic3cycle15CSRound2@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic3cycle15F3, Genomic3cycle15Pre, Genomic3cycle15EliS)
names(datasets)=c('Genomic3cycle15F3', 'Genomic3cycle15Pre', 'Genomic3cycle15EliS')
Genomic3cycle15=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic3cycle15)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic3cycle15
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic3cycle15PreQtl=al(Genomic3cycle15Pre, SP); Genomic3cycle15EliSQtl= al(Genomic3cycle15EliS, SP)

################################################################################
#Genomic3 exporting estimates and Qtl

Genomic3cycle1to15PreQtl=cbind(Burnincycle1PreQtl, Burnincycle2PreQtl, Burnincycle3PreQtl, Burnincycle4PreQtl, Burnincycle5PreQtl, Genomic3cycle6PreQtl, Genomic3cycle7PreQtl, Genomic3cycle8PreQtl, Genomic3cycle9PreQtl, Genomic3cycle10PreQtl, Genomic3cycle11PreQtl, Genomic3cycle12PreQtl, Genomic3cycle13PreQtl, Genomic3cycle14PreQtl, Genomic3cycle15PreQtl); write.csv(Genomic3cycle1to15PreQtl, 'r1_Genomic3PreQtl.csv')

Genomic3cycle1to15EliSQtl=cbind(Burnincycle1EliSQtl, Burnincycle2EliSQtl, Burnincycle3EliSQtl, Burnincycle4EliSQtl, Burnincycle5EliSQtl, Genomic3cycle6EliSQtl, Genomic3cycle7EliSQtl, Genomic3cycle8EliSQtl, Genomic3cycle9EliSQtl, Genomic3cycle10EliSQtl, Genomic3cycle11EliSQtl, Genomic3cycle12EliSQtl, Genomic3cycle13EliSQtl, Genomic3cycle14EliSQtl, Genomic3cycle15EliSQtl); write.csv(Genomic3cycle1to15EliSQtl, 'r1_Genomic3EliSQtl.csv')

Genomic3cycle1to15=rbind(Burnincycle1, Burnincycle2, Burnincycle3, Burnincycle4, Burnincycle5, Genomic3cycle6, Genomic3cycle7, Genomic3cycle8, Genomic3cycle9, Genomic3cycle10, Genomic3cycle11, Genomic3cycle12, Genomic3cycle13, Genomic3cycle14, Genomic3cycle15); write.csv(Genomic3cycle1to15, 'r1_Genomic3.csv')

################################################################################
#Genomic4cycle6
#Training, parentals and F1s
Tr=list(Burnincycle3PreS, Burnincycle2PreS); Tr=mergePops(Tr)
Par1=Burnincycle3PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Burnincycle4CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic4cycle6Par=setEBV(Par, GSModel); Genomic4cycle6F1=randCross(Genomic4cycle6Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic4cycle6F2=self(Genomic4cycle6F1, nProgeny=1); Genomic4cycle6F3=self(Genomic4cycle6F2, nProgeny=25)
Genomic4cycle6F3=setEBV(Genomic4cycle6F3, GSModel); cor(Genomic4cycle6F3@gv, Genomic4cycle6F3@ebv)
Genomic4cycle6F3S=selectInd(Genomic4cycle6F3, nInd=6000, use='ebv'); mean(Genomic4cycle6F3@gv); mean(Genomic4cycle6F3S@gv)
Genomic4cycle6F4=self(Genomic4cycle6F3S, nProgeny=1); Genomic4cycle6F5=self(Genomic4cycle6F4, nProgeny=1)

#Trials 
Genomic4cycle6Pre=setPheno(Genomic4cycle6F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic4cycle6Pre, traits=1, use='pheno')
Genomic4cycle6Pre=setEBV(Genomic4cycle6Pre, GSModel); Genomic4cycle6PreS=selectInd(Genomic4cycle6Pre, nInd=300, use='ebv')
Genomic4cycle6Adv=setPheno(Genomic4cycle6PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic4cycle6Adv, traits=1, use='pheno')
Genomic4cycle6Adv=setEBV(Genomic4cycle6Adv, GSModel); Genomic4cycle6AdvS=selectInd(Genomic4cycle6Adv, nInd=50, use='ebv')
Genomic4cycle6Eli=setPheno(Genomic4cycle6AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic4cycle6Eli, traits=1, use='pheno')
Genomic4cycle6Eli=setEBV(Genomic4cycle6Eli, GSModel); Genomic4cycle6EliS=selectInd(Genomic4cycle6Eli, nInd=10, use='ebv')
Genomic4cycle6Eli1=selectInd(Genomic4cycle6Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Burnincycle4PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic4cycle6CS=selectInd(F3, 150, use='ebv'); Genomic4cycle6CS=editGenomeTopQtl(Genomic4cycle6CS, 150, 10)
Genomic4cycle6F3Acc=cor(F3@gv, F3@ebv); Genomic4cycle6F3Acc
mean(Par@gv); mean(Genomic4cycle6CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic4cycle6F3, Genomic4cycle6Pre, Genomic4cycle6EliS)
names(datasets)=c('Genomic4cycle6F3', 'Genomic4cycle6Pre', 'Genomic4cycle6EliS')
Genomic4cycle6=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic4cycle6)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic4cycle6
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic4cycle6PreQtl=al(Genomic4cycle6Pre, SP); Genomic4cycle6EliSQtl= al(Genomic4cycle6EliS, SP)

#Genomic4cycle7
#Training, parentals and F1s
Tr=list(Burnincycle4PreS, Burnincycle3PreS); Tr=mergePops(Tr)
Par1=Burnincycle4PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Burnincycle5CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic4cycle7Par=setEBV(Par, GSModel); Genomic4cycle7F1=randCross(Genomic4cycle7Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic4cycle7F2=self(Genomic4cycle7F1, nProgeny=1); Genomic4cycle7F3=self(Genomic4cycle7F2, nProgeny=25)
Genomic4cycle7F3=setEBV(Genomic4cycle7F3, GSModel); cor(Genomic4cycle7F3@gv, Genomic4cycle7F3@ebv)
Genomic4cycle7F3S=selectInd(Genomic4cycle7F3, nInd=6000, use='ebv'); mean(Genomic4cycle7F3@gv); mean(Genomic4cycle7F3S@gv)
Genomic4cycle7F4=self(Genomic4cycle7F3S, nProgeny=1); Genomic4cycle7F5=self(Genomic4cycle7F4, nProgeny=1)

#Trials 
Genomic4cycle7Pre=setPheno(Genomic4cycle7F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic4cycle7Pre, traits=1, use='pheno')
Genomic4cycle7Pre=setEBV(Genomic4cycle7Pre, GSModel); Genomic4cycle7PreS=selectInd(Genomic4cycle7Pre, nInd=300, use='ebv')
Genomic4cycle7Adv=setPheno(Genomic4cycle7PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic4cycle7Adv, traits=1, use='pheno')
Genomic4cycle7Adv=setEBV(Genomic4cycle7Adv, GSModel); Genomic4cycle7AdvS=selectInd(Genomic4cycle7Adv, nInd=50, use='ebv')
Genomic4cycle7Eli=setPheno(Genomic4cycle7AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic4cycle7Eli, traits=1, use='pheno')
Genomic4cycle7Eli=setEBV(Genomic4cycle7Eli, GSModel); Genomic4cycle7EliS=selectInd(Genomic4cycle7Eli, nInd=10, use='ebv')
Genomic4cycle7Eli1=selectInd(Genomic4cycle7Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic4cycle6PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic4cycle7CS=selectInd(F3, 150, use='ebv'); Genomic4cycle7CS=editGenomeTopQtl(Genomic4cycle7CS, 150, 10)

Genomic4cycle7F3Acc=cor(F3@gv, F3@ebv); Genomic4cycle7F3Acc
mean(Par@gv); mean(Genomic4cycle7CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic4cycle7F3, Genomic4cycle7Pre, Genomic4cycle7EliS)
names(datasets)=c('Genomic4cycle7F3', 'Genomic4cycle7Pre', 'Genomic4cycle7EliS')
Genomic4cycle7=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic4cycle7)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic4cycle7
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic4cycle7PreQtl=al(Genomic4cycle7Pre, SP); Genomic4cycle7EliSQtl= al(Genomic4cycle7EliS, SP)

#Genomic4cycle8
#Training, parentals and F1s
Tr=list(Genomic4cycle6PreS, Burnincycle4PreS); Tr=mergePops(Tr)
Par1=Genomic4cycle6PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic4cycle6CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic4cycle8Par=setEBV(Par, GSModel); Genomic4cycle8F1=randCross(Genomic4cycle8Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic4cycle8F2=self(Genomic4cycle8F1, nProgeny=1); Genomic4cycle8F3=self(Genomic4cycle8F2, nProgeny=25)
Genomic4cycle8F3=setEBV(Genomic4cycle8F3, GSModel); cor(Genomic4cycle8F3@gv, Genomic4cycle8F3@ebv)
Genomic4cycle8F3S=selectInd(Genomic4cycle8F3, nInd=6000, use='ebv'); mean(Genomic4cycle8F3@gv); mean(Genomic4cycle8F3S@gv)
Genomic4cycle8F4=self(Genomic4cycle8F3S, nProgeny=1); Genomic4cycle8F5=self(Genomic4cycle8F4, nProgeny=1)

#Trials 
Genomic4cycle8Pre=setPheno(Genomic4cycle8F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic4cycle8Pre, traits=1, use='pheno')
Genomic4cycle8Pre=setEBV(Genomic4cycle8Pre, GSModel); Genomic4cycle8PreS=selectInd(Genomic4cycle8Pre, nInd=300, use='ebv')
Genomic4cycle8Adv=setPheno(Genomic4cycle8PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic4cycle8Adv, traits=1, use='pheno')
Genomic4cycle8Adv=setEBV(Genomic4cycle8Adv, GSModel); Genomic4cycle8AdvS=selectInd(Genomic4cycle8Adv, nInd=50, use='ebv')
Genomic4cycle8Eli=setPheno(Genomic4cycle8AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic4cycle8Eli, traits=1, use='pheno')
Genomic4cycle8Eli=setEBV(Genomic4cycle8Eli, GSModel); Genomic4cycle8EliS=selectInd(Genomic4cycle8Eli, nInd=10, use='ebv')
Genomic4cycle8Eli1=selectInd(Genomic4cycle8Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic4cycle7PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic4cycle8CS=selectInd(F3, 150, use='ebv'); Genomic4cycle8CS=editGenomeTopQtl(Genomic4cycle8CS, 150, 10)

Genomic4cycle8F3Acc=cor(F3@gv, F3@ebv); Genomic4cycle8F3Acc
mean(Par@gv); mean(Genomic4cycle8CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic4cycle8F3, Genomic4cycle8Pre, Genomic4cycle8EliS)
names(datasets)=c('Genomic4cycle8F3', 'Genomic4cycle8Pre', 'Genomic4cycle8EliS')
Genomic4cycle8=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic4cycle8)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic4cycle8
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic4cycle8PreQtl=al(Genomic4cycle8Pre, SP); Genomic4cycle8EliSQtl= al(Genomic4cycle8EliS, SP)

#Genomic4cycle9
#Training, parentals and F1s
Tr=list(Genomic4cycle7PreS, Genomic4cycle6PreS); Tr=mergePops(Tr)
Par1=Genomic4cycle7PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic4cycle7CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic4cycle9Par=setEBV(Par, GSModel); Genomic4cycle9F1=randCross(Genomic4cycle9Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic4cycle9F2=self(Genomic4cycle9F1, nProgeny=1); Genomic4cycle9F3=self(Genomic4cycle9F2, nProgeny=25)
Genomic4cycle9F3=setEBV(Genomic4cycle9F3, GSModel); cor(Genomic4cycle9F3@gv, Genomic4cycle9F3@ebv)
Genomic4cycle9F3S=selectInd(Genomic4cycle9F3, nInd=6000, use='ebv'); mean(Genomic4cycle9F3@gv); mean(Genomic4cycle9F3S@gv)
Genomic4cycle9F4=self(Genomic4cycle9F3S, nProgeny=1); Genomic4cycle9F5=self(Genomic4cycle9F4, nProgeny=1)

#Trials 
Genomic4cycle9Pre=setPheno(Genomic4cycle9F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic4cycle9Pre, traits=1, use='pheno')
Genomic4cycle9Pre=setEBV(Genomic4cycle9Pre, GSModel); Genomic4cycle9PreS=selectInd(Genomic4cycle9Pre, nInd=300, use='ebv')
Genomic4cycle9Adv=setPheno(Genomic4cycle9PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic4cycle9Adv, traits=1, use='pheno')
Genomic4cycle9Adv=setEBV(Genomic4cycle9Adv, GSModel); Genomic4cycle9AdvS=selectInd(Genomic4cycle9Adv, nInd=50, use='ebv')
Genomic4cycle9Eli=setPheno(Genomic4cycle9AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic4cycle9Eli, traits=1, use='pheno')
Genomic4cycle9Eli=setEBV(Genomic4cycle9Eli, GSModel); Genomic4cycle9EliS=selectInd(Genomic4cycle9Eli, nInd=10, use='ebv')
Genomic4cycle9Eli1=selectInd(Genomic4cycle9Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic4cycle8PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic4cycle9CS=selectInd(F3, 150, use='ebv'); Genomic4cycle9CS=editGenomeTopQtl(Genomic4cycle9CS, 150, 10)

Genomic4cycle9F3Acc=cor(F3@gv, F3@ebv); Genomic4cycle9F3Acc
mean(Par@gv); mean(Genomic4cycle9CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic4cycle9F3, Genomic4cycle9Pre, Genomic4cycle9EliS)
names(datasets)=c('Genomic4cycle9F3', 'Genomic4cycle9Pre', 'Genomic4cycle9EliS')
Genomic4cycle9=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic4cycle9)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic4cycle9
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic4cycle9PreQtl=al(Genomic4cycle9Pre, SP); Genomic4cycle9EliSQtl= al(Genomic4cycle9EliS, SP)

#Genomic4cycle10
#Training, parentals and F1s
Tr=list(Genomic4cycle8PreS, Genomic4cycle7PreS); Tr=mergePops(Tr)
Par1=Genomic4cycle8PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic4cycle8CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic4cycle10Par=setEBV(Par, GSModel); Genomic4cycle10F1=randCross(Genomic4cycle10Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic4cycle10F2=self(Genomic4cycle10F1, nProgeny=1); Genomic4cycle10F3=self(Genomic4cycle10F2, nProgeny=25)
Genomic4cycle10F3=setEBV(Genomic4cycle10F3, GSModel); cor(Genomic4cycle10F3@gv, Genomic4cycle10F3@ebv)
Genomic4cycle10F3S=selectInd(Genomic4cycle10F3, nInd=6000, use='ebv'); mean(Genomic4cycle10F3@gv); mean(Genomic4cycle10F3S@gv)
Genomic4cycle10F4=self(Genomic4cycle10F3S, nProgeny=1); Genomic4cycle10F5=self(Genomic4cycle10F4, nProgeny=1)

#Trials 
Genomic4cycle10Pre=setPheno(Genomic4cycle10F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic4cycle10Pre, traits=1, use='pheno')
Genomic4cycle10Pre=setEBV(Genomic4cycle10Pre, GSModel); Genomic4cycle10PreS=selectInd(Genomic4cycle10Pre, nInd=300, use='ebv')
Genomic4cycle10Adv=setPheno(Genomic4cycle10PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic4cycle10Adv, traits=1, use='pheno')
Genomic4cycle10Adv=setEBV(Genomic4cycle10Adv, GSModel); Genomic4cycle10AdvS=selectInd(Genomic4cycle10Adv, nInd=50, use='ebv')
Genomic4cycle10Eli=setPheno(Genomic4cycle10AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic4cycle10Eli, traits=1, use='pheno')
Genomic4cycle10Eli=setEBV(Genomic4cycle10Eli, GSModel); Genomic4cycle10EliS=selectInd(Genomic4cycle10Eli, nInd=10, use='ebv')
Genomic4cycle10Eli1=selectInd(Genomic4cycle10Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic4cycle9PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic4cycle10CS=selectInd(F3, 150, use='ebv'); Genomic4cycle10CS=editGenomeTopQtl(Genomic4cycle10CS, 150, 10)
Genomic4cycle10F3Acc=cor(F3@gv, F3@ebv); Genomic4cycle10F3Acc
mean(Par@gv); mean(Genomic4cycle10CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic4cycle10F3, Genomic4cycle10Pre, Genomic4cycle10EliS)
names(datasets)=c('Genomic4cycle10F3', 'Genomic4cycle10Pre', 'Genomic4cycle10EliS')
Genomic4cycle10=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic4cycle10)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic4cycle10
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic4cycle10PreQtl=al(Genomic4cycle10Pre, SP); Genomic4cycle10EliSQtl= al(Genomic4cycle10EliS, SP)

#Genomic4cycle11
#Training, parentals and F1s
Tr=list(Genomic4cycle9PreS, Genomic4cycle8PreS); Tr=mergePops(Tr)
Par1=Genomic4cycle9PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic4cycle9CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic4cycle11Par=setEBV(Par, GSModel); Genomic4cycle11F1=randCross(Genomic4cycle11Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic4cycle11F2=self(Genomic4cycle11F1, nProgeny=1); Genomic4cycle11F3=self(Genomic4cycle11F2, nProgeny=25)
Genomic4cycle11F3=setEBV(Genomic4cycle11F3, GSModel); cor(Genomic4cycle11F3@gv, Genomic4cycle11F3@ebv)
Genomic4cycle11F3S=selectInd(Genomic4cycle11F3, nInd=6000, use='ebv'); mean(Genomic4cycle11F3@gv); mean(Genomic4cycle11F3S@gv)
Genomic4cycle11F4=self(Genomic4cycle11F3S, nProgeny=1); Genomic4cycle11F5=self(Genomic4cycle11F4, nProgeny=1)

#Trials 
Genomic4cycle11Pre=setPheno(Genomic4cycle11F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic4cycle11Pre, traits=1, use='pheno')
Genomic4cycle11Pre=setEBV(Genomic4cycle11Pre, GSModel); Genomic4cycle11PreS=selectInd(Genomic4cycle11Pre, nInd=300, use='ebv')
Genomic4cycle11Adv=setPheno(Genomic4cycle11PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic4cycle11Adv, traits=1, use='pheno')
Genomic4cycle11Adv=setEBV(Genomic4cycle11Adv, GSModel); Genomic4cycle11AdvS=selectInd(Genomic4cycle11Adv, nInd=50, use='ebv')
Genomic4cycle11Eli=setPheno(Genomic4cycle11AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic4cycle11Eli, traits=1, use='pheno')
Genomic4cycle11Eli=setEBV(Genomic4cycle11Eli, GSModel); Genomic4cycle11EliS=selectInd(Genomic4cycle11Eli, nInd=10, use='ebv')
Genomic4cycle11Eli1=selectInd(Genomic4cycle11Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic4cycle10PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic4cycle11CS=selectInd(F3, 150, use='ebv'); Genomic4cycle11CS=editGenomeTopQtl(Genomic4cycle11CS, 150, 10)
Genomic4cycle11F3Acc=cor(F3@gv, F3@ebv); Genomic4cycle11F3Acc
mean(Par@gv); mean(Genomic4cycle11CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic4cycle11F3, Genomic4cycle11Pre, Genomic4cycle11EliS)
names(datasets)=c('Genomic4cycle11F3', 'Genomic4cycle11Pre', 'Genomic4cycle11EliS')
Genomic4cycle11=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic4cycle11)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic4cycle11
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic4cycle11PreQtl=al(Genomic4cycle11Pre, SP); Genomic4cycle11EliSQtl= al(Genomic4cycle11EliS, SP)

#Genomic4cycle12
#Training, parentals and F1s
Tr=list(Genomic4cycle10PreS, Genomic4cycle9PreS); Tr=mergePops(Tr)
Par1=Genomic4cycle10PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic4cycle10CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic4cycle12Par=setEBV(Par, GSModel); Genomic4cycle12F1=randCross(Genomic4cycle12Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic4cycle12F2=self(Genomic4cycle12F1, nProgeny=1); Genomic4cycle12F3=self(Genomic4cycle12F2, nProgeny=25)
Genomic4cycle12F3=setEBV(Genomic4cycle12F3, GSModel); cor(Genomic4cycle12F3@gv, Genomic4cycle12F3@ebv)
Genomic4cycle12F3S=selectInd(Genomic4cycle12F3, nInd=6000, use='ebv'); mean(Genomic4cycle12F3@gv); mean(Genomic4cycle12F3S@gv)
Genomic4cycle12F4=self(Genomic4cycle12F3S, nProgeny=1); Genomic4cycle12F5=self(Genomic4cycle12F4, nProgeny=1)

#Trials 
Genomic4cycle12Pre=setPheno(Genomic4cycle12F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic4cycle12Pre, traits=1, use='pheno')
Genomic4cycle12Pre=setEBV(Genomic4cycle12Pre, GSModel); Genomic4cycle12PreS=selectInd(Genomic4cycle12Pre, nInd=300, use='ebv')
Genomic4cycle12Adv=setPheno(Genomic4cycle12PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic4cycle12Adv, traits=1, use='pheno')
Genomic4cycle12Adv=setEBV(Genomic4cycle12Adv, GSModel); Genomic4cycle12AdvS=selectInd(Genomic4cycle12Adv, nInd=50, use='ebv')
Genomic4cycle12Eli=setPheno(Genomic4cycle12AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic4cycle12Eli, traits=1, use='pheno')
Genomic4cycle12Eli=setEBV(Genomic4cycle12Eli, GSModel); Genomic4cycle12EliS=selectInd(Genomic4cycle12Eli, nInd=10, use='ebv')
Genomic4cycle12Eli1=selectInd(Genomic4cycle12Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic4cycle11PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic4cycle12CS=selectInd(F3, 150, use='ebv'); Genomic4cycle12CS=editGenomeTopQtl(Genomic4cycle12CS, 150, 10)
Genomic4cycle12F3Acc=cor(F3@gv, F3@ebv); Genomic4cycle12F3Acc
mean(Par@gv); mean(Genomic4cycle12CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic4cycle12F3, Genomic4cycle12Pre, Genomic4cycle12EliS)
names(datasets)=c('Genomic4cycle12F3', 'Genomic4cycle12Pre', 'Genomic4cycle12EliS')
Genomic4cycle12=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic4cycle12)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic4cycle12
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic4cycle12PreQtl=al(Genomic4cycle12Pre, SP); Genomic4cycle12EliSQtl= al(Genomic4cycle12EliS, SP)

#Genomic4cycle13
#Training, parentals and F1s
Tr=list(Genomic4cycle11PreS, Genomic4cycle10PreS); Tr=mergePops(Tr)
Par1=Genomic4cycle11PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic4cycle11CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic4cycle13Par=setEBV(Par, GSModel); Genomic4cycle13F1=randCross(Genomic4cycle13Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic4cycle13F2=self(Genomic4cycle13F1, nProgeny=1); Genomic4cycle13F3=self(Genomic4cycle13F2, nProgeny=25)
Genomic4cycle13F3=setEBV(Genomic4cycle13F3, GSModel); cor(Genomic4cycle13F3@gv, Genomic4cycle13F3@ebv)
Genomic4cycle13F3S=selectInd(Genomic4cycle13F3, nInd=6000, use='ebv'); mean(Genomic4cycle13F3@gv); mean(Genomic4cycle13F3S@gv)
Genomic4cycle13F4=self(Genomic4cycle13F3S, nProgeny=1); Genomic4cycle13F5=self(Genomic4cycle13F4, nProgeny=1)

#Trials 
Genomic4cycle13Pre=setPheno(Genomic4cycle13F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic4cycle13Pre, traits=1, use='pheno')
Genomic4cycle13Pre=setEBV(Genomic4cycle13Pre, GSModel); Genomic4cycle13PreS=selectInd(Genomic4cycle13Pre, nInd=300, use='ebv')
Genomic4cycle13Adv=setPheno(Genomic4cycle13PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic4cycle13Adv, traits=1, use='pheno')
Genomic4cycle13Adv=setEBV(Genomic4cycle13Adv, GSModel); Genomic4cycle13AdvS=selectInd(Genomic4cycle13Adv, nInd=50, use='ebv')
Genomic4cycle13Eli=setPheno(Genomic4cycle13AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic4cycle13Eli, traits=1, use='pheno')
Genomic4cycle13Eli=setEBV(Genomic4cycle13Eli, GSModel); Genomic4cycle13EliS=selectInd(Genomic4cycle13Eli, nInd=10, use='ebv')
Genomic4cycle13Eli1=selectInd(Genomic4cycle13Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic4cycle12PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic4cycle13CS=selectInd(F3, 150, use='ebv'); Genomic4cycle13CS=editGenomeTopQtl(Genomic4cycle13CS, 150, 10)
Genomic4cycle13F3Acc=cor(F3@gv, F3@ebv); Genomic4cycle13F3Acc
mean(Par@gv); mean(Genomic4cycle13CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic4cycle13F3, Genomic4cycle13Pre, Genomic4cycle13EliS)
names(datasets)=c('Genomic4cycle13F3', 'Genomic4cycle13Pre', 'Genomic4cycle13EliS')
Genomic4cycle13=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic4cycle13)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic4cycle13
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic4cycle13PreQtl=al(Genomic4cycle13Pre, SP); Genomic4cycle13EliSQtl= al(Genomic4cycle13EliS, SP)

#Genomic4cycle14
#Training, parentals and F1s
Tr=list(Genomic4cycle12PreS, Genomic4cycle11PreS); Tr=mergePops(Tr)
Par1=Genomic4cycle12PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic4cycle12CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic4cycle14Par=setEBV(Par, GSModel); Genomic4cycle14F1=randCross(Genomic4cycle14Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic4cycle14F2=self(Genomic4cycle14F1, nProgeny=1); Genomic4cycle14F3=self(Genomic4cycle14F2, nProgeny=25)
Genomic4cycle14F3=setEBV(Genomic4cycle14F3, GSModel); cor(Genomic4cycle14F3@gv, Genomic4cycle14F3@ebv)
Genomic4cycle14F3S=selectInd(Genomic4cycle14F3, nInd=6000, use='ebv'); mean(Genomic4cycle14F3@gv); mean(Genomic4cycle14F3S@gv)
Genomic4cycle14F4=self(Genomic4cycle14F3S, nProgeny=1); Genomic4cycle14F5=self(Genomic4cycle14F4, nProgeny=1)

#Trials 
Genomic4cycle14Pre=setPheno(Genomic4cycle14F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic4cycle14Pre, traits=1, use='pheno')
Genomic4cycle14Pre=setEBV(Genomic4cycle14Pre, GSModel); Genomic4cycle14PreS=selectInd(Genomic4cycle14Pre, nInd=300, use='ebv')
Genomic4cycle14Adv=setPheno(Genomic4cycle14PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic4cycle14Adv, traits=1, use='pheno')
Genomic4cycle14Adv=setEBV(Genomic4cycle14Adv, GSModel); Genomic4cycle14AdvS=selectInd(Genomic4cycle14Adv, nInd=50, use='ebv')
Genomic4cycle14Eli=setPheno(Genomic4cycle14AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic4cycle14Eli, traits=1, use='pheno')
Genomic4cycle14Eli=setEBV(Genomic4cycle14Eli, GSModel); Genomic4cycle14EliS=selectInd(Genomic4cycle14Eli, nInd=10, use='ebv')
Genomic4cycle14Eli1=selectInd(Genomic4cycle14Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic4cycle13PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic4cycle14CS=selectInd(F3, 150, use='ebv'); Genomic4cycle14CS=editGenomeTopQtl(Genomic4cycle14CS, 150, 10)
Genomic4cycle14F3Acc=cor(F3@gv, F3@ebv); Genomic4cycle14F3Acc
mean(Par@gv); mean(Genomic4cycle14CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic4cycle14F3, Genomic4cycle14Pre, Genomic4cycle14EliS)
names(datasets)=c('Genomic4cycle14F3', 'Genomic4cycle14Pre', 'Genomic4cycle14EliS')
Genomic4cycle14=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic4cycle14)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic4cycle14
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic4cycle14PreQtl=al(Genomic4cycle14Pre, SP); Genomic4cycle14EliSQtl= al(Genomic4cycle14EliS, SP)

#Genomic4cycle15
#Training, parentals and F1s
Tr=list(Genomic4cycle13PreS, Genomic4cycle12PreS); Tr=mergePops(Tr)
Par1=Genomic4cycle13PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic4cycle13CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic4cycle15Par=setEBV(Par, GSModel); Genomic4cycle15F1=randCross(Genomic4cycle15Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic4cycle15F2=self(Genomic4cycle15F1, nProgeny=1); Genomic4cycle15F3=self(Genomic4cycle15F2, nProgeny=25)
Genomic4cycle15F3=setEBV(Genomic4cycle15F3, GSModel); cor(Genomic4cycle15F3@gv, Genomic4cycle15F3@ebv)
Genomic4cycle15F3S=selectInd(Genomic4cycle15F3, nInd=6000, use='ebv'); mean(Genomic4cycle15F3@gv); mean(Genomic4cycle15F3S@gv)
Genomic4cycle15F4=self(Genomic4cycle15F3S, nProgeny=1); Genomic4cycle15F5=self(Genomic4cycle15F4, nProgeny=1)

#Trials 
Genomic4cycle15Pre=setPheno(Genomic4cycle15F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic4cycle15Pre, traits=1, use='pheno')
Genomic4cycle15Pre=setEBV(Genomic4cycle15Pre, GSModel); Genomic4cycle15PreS=selectInd(Genomic4cycle15Pre, nInd=300, use='ebv')
Genomic4cycle15Adv=setPheno(Genomic4cycle15PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic4cycle15Adv, traits=1, use='pheno')
Genomic4cycle15Adv=setEBV(Genomic4cycle15Adv, GSModel); Genomic4cycle15AdvS=selectInd(Genomic4cycle15Adv, nInd=50, use='ebv')
Genomic4cycle15Eli=setPheno(Genomic4cycle15AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic4cycle15Eli, traits=1, use='pheno')
Genomic4cycle15Eli=setEBV(Genomic4cycle15Eli, GSModel); Genomic4cycle15EliS=selectInd(Genomic4cycle15Eli, nInd=10, use='ebv')
Genomic4cycle15Eli1=selectInd(Genomic4cycle15Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic4cycle14PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic4cycle15CS=selectInd(F3, 150, use='ebv'); Genomic4cycle15CS=editGenomeTopQtl(Genomic4cycle15CS, 150, 10)
Genomic4cycle15F3Acc=cor(F3@gv, F3@ebv); Genomic4cycle15F3Acc
mean(Par@gv); mean(Genomic4cycle15CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic4cycle15F3, Genomic4cycle15Pre, Genomic4cycle15EliS)
names(datasets)=c('Genomic4cycle15F3', 'Genomic4cycle15Pre', 'Genomic4cycle15EliS')
Genomic4cycle15=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic4cycle15)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic4cycle15
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic4cycle15PreQtl=al(Genomic4cycle15Pre, SP); Genomic4cycle15EliSQtl= al(Genomic4cycle15EliS, SP)

################################################################################
#Genomic4 exporting estimates and Qtl

Genomic4cycle1to15PreQtl=cbind(Burnincycle1PreQtl, Burnincycle2PreQtl, Burnincycle3PreQtl, Burnincycle4PreQtl, Burnincycle5PreQtl, Genomic4cycle6PreQtl, Genomic4cycle7PreQtl, Genomic4cycle8PreQtl, Genomic4cycle9PreQtl, Genomic4cycle10PreQtl, Genomic4cycle11PreQtl, Genomic4cycle12PreQtl, Genomic4cycle13PreQtl, Genomic4cycle14PreQtl, Genomic4cycle15PreQtl); write.csv(Genomic4cycle1to15PreQtl, 'r1_Genomic4PreQtl.csv')

Genomic4cycle1to15EliSQtl=cbind(Burnincycle1EliSQtl, Burnincycle2EliSQtl, Burnincycle3EliSQtl, Burnincycle4EliSQtl, Burnincycle5EliSQtl, Genomic4cycle6EliSQtl, Genomic4cycle7EliSQtl, Genomic4cycle8EliSQtl, Genomic4cycle9EliSQtl, Genomic4cycle10EliSQtl, Genomic4cycle11EliSQtl, Genomic4cycle12EliSQtl, Genomic4cycle13EliSQtl, Genomic4cycle14EliSQtl, Genomic4cycle15EliSQtl); write.csv(Genomic4cycle1to15EliSQtl, 'r1_Genomic4EliSQtl.csv')

Genomic4cycle1to15=rbind(Burnincycle1, Burnincycle2, Burnincycle3, Burnincycle4, Burnincycle5, Genomic4cycle6, Genomic4cycle7, Genomic4cycle8, Genomic4cycle9, Genomic4cycle10, Genomic4cycle11, Genomic4cycle12, Genomic4cycle13, Genomic4cycle14, Genomic4cycle15); write.csv(Genomic4cycle1to15, 'r1_Genomic4.csv')

################################################################################
#Genomic5cycle6
#Training, parentals and F1s
Tr=list(Burnincycle3PreS, Burnincycle2PreS); Tr=mergePops(Tr)
Par1=Burnincycle3PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Burnincycle4CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic5cycle6Par=setEBV(Par, GSModel); Genomic5cycle6F1=randCross(Genomic5cycle6Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic5cycle6F2=self(Genomic5cycle6F1, nProgeny=1); Genomic5cycle6F3=self(Genomic5cycle6F2, nProgeny=25)
Genomic5cycle6F3=setEBV(Genomic5cycle6F3, GSModel); cor(Genomic5cycle6F3@gv, Genomic5cycle6F3@ebv)
Genomic5cycle6F3S=selectInd(Genomic5cycle6F3, nInd=6000, use='ebv'); mean(Genomic5cycle6F3@gv); mean(Genomic5cycle6F3S@gv)
Genomic5cycle6F4=self(Genomic5cycle6F3S, nProgeny=1); Genomic5cycle6F5=self(Genomic5cycle6F4, nProgeny=1)

#Trials 
Genomic5cycle6Pre=setPheno(Genomic5cycle6F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic5cycle6Pre, traits=1, use='pheno')
Genomic5cycle6Pre=setEBV(Genomic5cycle6Pre, GSModel); Genomic5cycle6PreS=selectInd(Genomic5cycle6Pre, nInd=300, use='ebv')
Genomic5cycle6Adv=setPheno(Genomic5cycle6PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic5cycle6Adv, traits=1, use='pheno')
Genomic5cycle6Adv=setEBV(Genomic5cycle6Adv, GSModel); Genomic5cycle6AdvS=selectInd(Genomic5cycle6Adv, nInd=50, use='ebv')
Genomic5cycle6Eli=setPheno(Genomic5cycle6AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic5cycle6Eli, traits=1, use='pheno')
Genomic5cycle6Eli=setEBV(Genomic5cycle6Eli, GSModel); Genomic5cycle6EliS=selectInd(Genomic5cycle6Eli, nInd=10, use='ebv')
Genomic5cycle6Eli1=selectInd(Genomic5cycle6Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Burnincycle4PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic5cycle6CS=selectInd(F3, 150, use='ebv'); Genomic5cycle6CS=editGenomeTopQtl(Genomic5cycle6CS, 150, 50)
Genomic5cycle6F3Acc=cor(F3@gv, F3@ebv); Genomic5cycle6F3Acc
mean(Par@gv); mean(Genomic5cycle6CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic5cycle6F3, Genomic5cycle6Pre, Genomic5cycle6EliS)
names(datasets)=c('Genomic5cycle6F3', 'Genomic5cycle6Pre', 'Genomic5cycle6EliS')
Genomic5cycle6=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic5cycle6)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic5cycle6
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic5cycle6PreQtl=al(Genomic5cycle6Pre, SP); Genomic5cycle6EliSQtl= al(Genomic5cycle6EliS, SP)

#Genomic5cycle7
#Training, parentals and F1s
Tr=list(Burnincycle4PreS, Burnincycle3PreS); Tr=mergePops(Tr)
Par1=Burnincycle4PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Burnincycle5CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic5cycle7Par=setEBV(Par, GSModel); Genomic5cycle7F1=randCross(Genomic5cycle7Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic5cycle7F2=self(Genomic5cycle7F1, nProgeny=1); Genomic5cycle7F3=self(Genomic5cycle7F2, nProgeny=25)
Genomic5cycle7F3=setEBV(Genomic5cycle7F3, GSModel); cor(Genomic5cycle7F3@gv, Genomic5cycle7F3@ebv)
Genomic5cycle7F3S=selectInd(Genomic5cycle7F3, nInd=6000, use='ebv'); mean(Genomic5cycle7F3@gv); mean(Genomic5cycle7F3S@gv)
Genomic5cycle7F4=self(Genomic5cycle7F3S, nProgeny=1); Genomic5cycle7F5=self(Genomic5cycle7F4, nProgeny=1)

#Trials 
Genomic5cycle7Pre=setPheno(Genomic5cycle7F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic5cycle7Pre, traits=1, use='pheno')
Genomic5cycle7Pre=setEBV(Genomic5cycle7Pre, GSModel); Genomic5cycle7PreS=selectInd(Genomic5cycle7Pre, nInd=300, use='ebv')
Genomic5cycle7Adv=setPheno(Genomic5cycle7PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic5cycle7Adv, traits=1, use='pheno')
Genomic5cycle7Adv=setEBV(Genomic5cycle7Adv, GSModel); Genomic5cycle7AdvS=selectInd(Genomic5cycle7Adv, nInd=50, use='ebv')
Genomic5cycle7Eli=setPheno(Genomic5cycle7AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic5cycle7Eli, traits=1, use='pheno')
Genomic5cycle7Eli=setEBV(Genomic5cycle7Eli, GSModel); Genomic5cycle7EliS=selectInd(Genomic5cycle7Eli, nInd=10, use='ebv')
Genomic5cycle7Eli1=selectInd(Genomic5cycle7Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic5cycle6PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic5cycle7CS=selectInd(F3, 150, use='ebv'); Genomic5cycle7CS=editGenomeTopQtl(Genomic5cycle7CS, 150, 50)

Genomic5cycle7F3Acc=cor(F3@gv, F3@ebv); Genomic5cycle7F3Acc
mean(Par@gv); mean(Genomic5cycle7CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic5cycle7F3, Genomic5cycle7Pre, Genomic5cycle7EliS)
names(datasets)=c('Genomic5cycle7F3', 'Genomic5cycle7Pre', 'Genomic5cycle7EliS')
Genomic5cycle7=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic5cycle7)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic5cycle7
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic5cycle7PreQtl=al(Genomic5cycle7Pre, SP); Genomic5cycle7EliSQtl= al(Genomic5cycle7EliS, SP)

#Genomic5cycle8
#Training, parentals and F1s
Tr=list(Genomic5cycle6PreS, Burnincycle4PreS); Tr=mergePops(Tr)
Par1=Genomic5cycle6PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic5cycle6CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic5cycle8Par=setEBV(Par, GSModel); Genomic5cycle8F1=randCross(Genomic5cycle8Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic5cycle8F2=self(Genomic5cycle8F1, nProgeny=1); Genomic5cycle8F3=self(Genomic5cycle8F2, nProgeny=25)
Genomic5cycle8F3=setEBV(Genomic5cycle8F3, GSModel); cor(Genomic5cycle8F3@gv, Genomic5cycle8F3@ebv)
Genomic5cycle8F3S=selectInd(Genomic5cycle8F3, nInd=6000, use='ebv'); mean(Genomic5cycle8F3@gv); mean(Genomic5cycle8F3S@gv)
Genomic5cycle8F4=self(Genomic5cycle8F3S, nProgeny=1); Genomic5cycle8F5=self(Genomic5cycle8F4, nProgeny=1)

#Trials 
Genomic5cycle8Pre=setPheno(Genomic5cycle8F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic5cycle8Pre, traits=1, use='pheno')
Genomic5cycle8Pre=setEBV(Genomic5cycle8Pre, GSModel); Genomic5cycle8PreS=selectInd(Genomic5cycle8Pre, nInd=300, use='ebv')
Genomic5cycle8Adv=setPheno(Genomic5cycle8PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic5cycle8Adv, traits=1, use='pheno')
Genomic5cycle8Adv=setEBV(Genomic5cycle8Adv, GSModel); Genomic5cycle8AdvS=selectInd(Genomic5cycle8Adv, nInd=50, use='ebv')
Genomic5cycle8Eli=setPheno(Genomic5cycle8AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic5cycle8Eli, traits=1, use='pheno')
Genomic5cycle8Eli=setEBV(Genomic5cycle8Eli, GSModel); Genomic5cycle8EliS=selectInd(Genomic5cycle8Eli, nInd=10, use='ebv')
Genomic5cycle8Eli1=selectInd(Genomic5cycle8Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic5cycle7PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic5cycle8CS=selectInd(F3, 150, use='ebv'); Genomic5cycle8CS=editGenomeTopQtl(Genomic5cycle8CS, 150, 50)

Genomic5cycle8F3Acc=cor(F3@gv, F3@ebv); Genomic5cycle8F3Acc
mean(Par@gv); mean(Genomic5cycle8CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic5cycle8F3, Genomic5cycle8Pre, Genomic5cycle8EliS)
names(datasets)=c('Genomic5cycle8F3', 'Genomic5cycle8Pre', 'Genomic5cycle8EliS')
Genomic5cycle8=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic5cycle8)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic5cycle8
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic5cycle8PreQtl=al(Genomic5cycle8Pre, SP); Genomic5cycle8EliSQtl= al(Genomic5cycle8EliS, SP)

#Genomic5cycle9
#Training, parentals and F1s
Tr=list(Genomic5cycle7PreS, Genomic5cycle6PreS); Tr=mergePops(Tr)
Par1=Genomic5cycle7PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic5cycle7CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic5cycle9Par=setEBV(Par, GSModel); Genomic5cycle9F1=randCross(Genomic5cycle9Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic5cycle9F2=self(Genomic5cycle9F1, nProgeny=1); Genomic5cycle9F3=self(Genomic5cycle9F2, nProgeny=25)
Genomic5cycle9F3=setEBV(Genomic5cycle9F3, GSModel); cor(Genomic5cycle9F3@gv, Genomic5cycle9F3@ebv)
Genomic5cycle9F3S=selectInd(Genomic5cycle9F3, nInd=6000, use='ebv'); mean(Genomic5cycle9F3@gv); mean(Genomic5cycle9F3S@gv)
Genomic5cycle9F4=self(Genomic5cycle9F3S, nProgeny=1); Genomic5cycle9F5=self(Genomic5cycle9F4, nProgeny=1)

#Trials 
Genomic5cycle9Pre=setPheno(Genomic5cycle9F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic5cycle9Pre, traits=1, use='pheno')
Genomic5cycle9Pre=setEBV(Genomic5cycle9Pre, GSModel); Genomic5cycle9PreS=selectInd(Genomic5cycle9Pre, nInd=300, use='ebv')
Genomic5cycle9Adv=setPheno(Genomic5cycle9PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic5cycle9Adv, traits=1, use='pheno')
Genomic5cycle9Adv=setEBV(Genomic5cycle9Adv, GSModel); Genomic5cycle9AdvS=selectInd(Genomic5cycle9Adv, nInd=50, use='ebv')
Genomic5cycle9Eli=setPheno(Genomic5cycle9AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic5cycle9Eli, traits=1, use='pheno')
Genomic5cycle9Eli=setEBV(Genomic5cycle9Eli, GSModel); Genomic5cycle9EliS=selectInd(Genomic5cycle9Eli, nInd=10, use='ebv')
Genomic5cycle9Eli1=selectInd(Genomic5cycle9Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic5cycle8PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic5cycle9CS=selectInd(F3, 150, use='ebv'); Genomic5cycle9CS=editGenomeTopQtl(Genomic5cycle9CS, 150, 50)

Genomic5cycle9F3Acc=cor(F3@gv, F3@ebv); Genomic5cycle9F3Acc
mean(Par@gv); mean(Genomic5cycle9CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic5cycle9F3, Genomic5cycle9Pre, Genomic5cycle9EliS)
names(datasets)=c('Genomic5cycle9F3', 'Genomic5cycle9Pre', 'Genomic5cycle9EliS')
Genomic5cycle9=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic5cycle9)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic5cycle9
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic5cycle9PreQtl=al(Genomic5cycle9Pre, SP); Genomic5cycle9EliSQtl= al(Genomic5cycle9EliS, SP)

#Genomic5cycle10
#Training, parentals and F1s
Tr=list(Genomic5cycle8PreS, Genomic5cycle7PreS); Tr=mergePops(Tr)
Par1=Genomic5cycle8PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic5cycle8CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic5cycle10Par=setEBV(Par, GSModel); Genomic5cycle10F1=randCross(Genomic5cycle10Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic5cycle10F2=self(Genomic5cycle10F1, nProgeny=1); Genomic5cycle10F3=self(Genomic5cycle10F2, nProgeny=25)
Genomic5cycle10F3=setEBV(Genomic5cycle10F3, GSModel); cor(Genomic5cycle10F3@gv, Genomic5cycle10F3@ebv)
Genomic5cycle10F3S=selectInd(Genomic5cycle10F3, nInd=6000, use='ebv'); mean(Genomic5cycle10F3@gv); mean(Genomic5cycle10F3S@gv)
Genomic5cycle10F4=self(Genomic5cycle10F3S, nProgeny=1); Genomic5cycle10F5=self(Genomic5cycle10F4, nProgeny=1)

#Trials 
Genomic5cycle10Pre=setPheno(Genomic5cycle10F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic5cycle10Pre, traits=1, use='pheno')
Genomic5cycle10Pre=setEBV(Genomic5cycle10Pre, GSModel); Genomic5cycle10PreS=selectInd(Genomic5cycle10Pre, nInd=300, use='ebv')
Genomic5cycle10Adv=setPheno(Genomic5cycle10PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic5cycle10Adv, traits=1, use='pheno')
Genomic5cycle10Adv=setEBV(Genomic5cycle10Adv, GSModel); Genomic5cycle10AdvS=selectInd(Genomic5cycle10Adv, nInd=50, use='ebv')
Genomic5cycle10Eli=setPheno(Genomic5cycle10AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic5cycle10Eli, traits=1, use='pheno')
Genomic5cycle10Eli=setEBV(Genomic5cycle10Eli, GSModel); Genomic5cycle10EliS=selectInd(Genomic5cycle10Eli, nInd=10, use='ebv')
Genomic5cycle10Eli1=selectInd(Genomic5cycle10Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic5cycle9PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic5cycle10CS=selectInd(F3, 150, use='ebv'); Genomic5cycle10CS=editGenomeTopQtl(Genomic5cycle10CS, 150, 50)
Genomic5cycle10F3Acc=cor(F3@gv, F3@ebv); Genomic5cycle10F3Acc
mean(Par@gv); mean(Genomic5cycle10CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic5cycle10F3, Genomic5cycle10Pre, Genomic5cycle10EliS)
names(datasets)=c('Genomic5cycle10F3', 'Genomic5cycle10Pre', 'Genomic5cycle10EliS')
Genomic5cycle10=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic5cycle10)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic5cycle10
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic5cycle10PreQtl=al(Genomic5cycle10Pre, SP); Genomic5cycle10EliSQtl= al(Genomic5cycle10EliS, SP)

#Genomic5cycle11
#Training, parentals and F1s
Tr=list(Genomic5cycle9PreS, Genomic5cycle8PreS); Tr=mergePops(Tr)
Par1=Genomic5cycle9PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic5cycle9CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic5cycle11Par=setEBV(Par, GSModel); Genomic5cycle11F1=randCross(Genomic5cycle11Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic5cycle11F2=self(Genomic5cycle11F1, nProgeny=1); Genomic5cycle11F3=self(Genomic5cycle11F2, nProgeny=25)
Genomic5cycle11F3=setEBV(Genomic5cycle11F3, GSModel); cor(Genomic5cycle11F3@gv, Genomic5cycle11F3@ebv)
Genomic5cycle11F3S=selectInd(Genomic5cycle11F3, nInd=6000, use='ebv'); mean(Genomic5cycle11F3@gv); mean(Genomic5cycle11F3S@gv)
Genomic5cycle11F4=self(Genomic5cycle11F3S, nProgeny=1); Genomic5cycle11F5=self(Genomic5cycle11F4, nProgeny=1)

#Trials 
Genomic5cycle11Pre=setPheno(Genomic5cycle11F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic5cycle11Pre, traits=1, use='pheno')
Genomic5cycle11Pre=setEBV(Genomic5cycle11Pre, GSModel); Genomic5cycle11PreS=selectInd(Genomic5cycle11Pre, nInd=300, use='ebv')
Genomic5cycle11Adv=setPheno(Genomic5cycle11PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic5cycle11Adv, traits=1, use='pheno')
Genomic5cycle11Adv=setEBV(Genomic5cycle11Adv, GSModel); Genomic5cycle11AdvS=selectInd(Genomic5cycle11Adv, nInd=50, use='ebv')
Genomic5cycle11Eli=setPheno(Genomic5cycle11AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic5cycle11Eli, traits=1, use='pheno')
Genomic5cycle11Eli=setEBV(Genomic5cycle11Eli, GSModel); Genomic5cycle11EliS=selectInd(Genomic5cycle11Eli, nInd=10, use='ebv')
Genomic5cycle11Eli1=selectInd(Genomic5cycle11Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic5cycle10PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic5cycle11CS=selectInd(F3, 150, use='ebv'); Genomic5cycle11CS=editGenomeTopQtl(Genomic5cycle11CS, 150, 50)
Genomic5cycle11F3Acc=cor(F3@gv, F3@ebv); Genomic5cycle11F3Acc
mean(Par@gv); mean(Genomic5cycle11CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic5cycle11F3, Genomic5cycle11Pre, Genomic5cycle11EliS)
names(datasets)=c('Genomic5cycle11F3', 'Genomic5cycle11Pre', 'Genomic5cycle11EliS')
Genomic5cycle11=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic5cycle11)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic5cycle11
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic5cycle11PreQtl=al(Genomic5cycle11Pre, SP); Genomic5cycle11EliSQtl= al(Genomic5cycle11EliS, SP)

#Genomic5cycle12
#Training, parentals and F1s
Tr=list(Genomic5cycle10PreS, Genomic5cycle9PreS); Tr=mergePops(Tr)
Par1=Genomic5cycle10PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic5cycle10CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic5cycle12Par=setEBV(Par, GSModel); Genomic5cycle12F1=randCross(Genomic5cycle12Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic5cycle12F2=self(Genomic5cycle12F1, nProgeny=1); Genomic5cycle12F3=self(Genomic5cycle12F2, nProgeny=25)
Genomic5cycle12F3=setEBV(Genomic5cycle12F3, GSModel); cor(Genomic5cycle12F3@gv, Genomic5cycle12F3@ebv)
Genomic5cycle12F3S=selectInd(Genomic5cycle12F3, nInd=6000, use='ebv'); mean(Genomic5cycle12F3@gv); mean(Genomic5cycle12F3S@gv)
Genomic5cycle12F4=self(Genomic5cycle12F3S, nProgeny=1); Genomic5cycle12F5=self(Genomic5cycle12F4, nProgeny=1)

#Trials 
Genomic5cycle12Pre=setPheno(Genomic5cycle12F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic5cycle12Pre, traits=1, use='pheno')
Genomic5cycle12Pre=setEBV(Genomic5cycle12Pre, GSModel); Genomic5cycle12PreS=selectInd(Genomic5cycle12Pre, nInd=300, use='ebv')
Genomic5cycle12Adv=setPheno(Genomic5cycle12PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic5cycle12Adv, traits=1, use='pheno')
Genomic5cycle12Adv=setEBV(Genomic5cycle12Adv, GSModel); Genomic5cycle12AdvS=selectInd(Genomic5cycle12Adv, nInd=50, use='ebv')
Genomic5cycle12Eli=setPheno(Genomic5cycle12AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic5cycle12Eli, traits=1, use='pheno')
Genomic5cycle12Eli=setEBV(Genomic5cycle12Eli, GSModel); Genomic5cycle12EliS=selectInd(Genomic5cycle12Eli, nInd=10, use='ebv')
Genomic5cycle12Eli1=selectInd(Genomic5cycle12Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic5cycle11PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic5cycle12CS=selectInd(F3, 150, use='ebv'); Genomic5cycle12CS=editGenomeTopQtl(Genomic5cycle12CS, 150, 50)
Genomic5cycle12F3Acc=cor(F3@gv, F3@ebv); Genomic5cycle12F3Acc
mean(Par@gv); mean(Genomic5cycle12CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic5cycle12F3, Genomic5cycle12Pre, Genomic5cycle12EliS)
names(datasets)=c('Genomic5cycle12F3', 'Genomic5cycle12Pre', 'Genomic5cycle12EliS')
Genomic5cycle12=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic5cycle12)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic5cycle12
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic5cycle12PreQtl=al(Genomic5cycle12Pre, SP); Genomic5cycle12EliSQtl= al(Genomic5cycle12EliS, SP)

#Genomic5cycle13
#Training, parentals and F1s
Tr=list(Genomic5cycle11PreS, Genomic5cycle10PreS); Tr=mergePops(Tr)
Par1=Genomic5cycle11PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic5cycle11CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic5cycle13Par=setEBV(Par, GSModel); Genomic5cycle13F1=randCross(Genomic5cycle13Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic5cycle13F2=self(Genomic5cycle13F1, nProgeny=1); Genomic5cycle13F3=self(Genomic5cycle13F2, nProgeny=25)
Genomic5cycle13F3=setEBV(Genomic5cycle13F3, GSModel); cor(Genomic5cycle13F3@gv, Genomic5cycle13F3@ebv)
Genomic5cycle13F3S=selectInd(Genomic5cycle13F3, nInd=6000, use='ebv'); mean(Genomic5cycle13F3@gv); mean(Genomic5cycle13F3S@gv)
Genomic5cycle13F4=self(Genomic5cycle13F3S, nProgeny=1); Genomic5cycle13F5=self(Genomic5cycle13F4, nProgeny=1)

#Trials 
Genomic5cycle13Pre=setPheno(Genomic5cycle13F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic5cycle13Pre, traits=1, use='pheno')
Genomic5cycle13Pre=setEBV(Genomic5cycle13Pre, GSModel); Genomic5cycle13PreS=selectInd(Genomic5cycle13Pre, nInd=300, use='ebv')
Genomic5cycle13Adv=setPheno(Genomic5cycle13PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic5cycle13Adv, traits=1, use='pheno')
Genomic5cycle13Adv=setEBV(Genomic5cycle13Adv, GSModel); Genomic5cycle13AdvS=selectInd(Genomic5cycle13Adv, nInd=50, use='ebv')
Genomic5cycle13Eli=setPheno(Genomic5cycle13AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic5cycle13Eli, traits=1, use='pheno')
Genomic5cycle13Eli=setEBV(Genomic5cycle13Eli, GSModel); Genomic5cycle13EliS=selectInd(Genomic5cycle13Eli, nInd=10, use='ebv')
Genomic5cycle13Eli1=selectInd(Genomic5cycle13Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic5cycle12PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic5cycle13CS=selectInd(F3, 150, use='ebv'); Genomic5cycle13CS=editGenomeTopQtl(Genomic5cycle13CS, 150, 50)
Genomic5cycle13F3Acc=cor(F3@gv, F3@ebv); Genomic5cycle13F3Acc
mean(Par@gv); mean(Genomic5cycle13CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic5cycle13F3, Genomic5cycle13Pre, Genomic5cycle13EliS)
names(datasets)=c('Genomic5cycle13F3', 'Genomic5cycle13Pre', 'Genomic5cycle13EliS')
Genomic5cycle13=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic5cycle13)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic5cycle13
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic5cycle13PreQtl=al(Genomic5cycle13Pre, SP); Genomic5cycle13EliSQtl= al(Genomic5cycle13EliS, SP)

#Genomic5cycle14
#Training, parentals and F1s
Tr=list(Genomic5cycle12PreS, Genomic5cycle11PreS); Tr=mergePops(Tr)
Par1=Genomic5cycle12PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic5cycle12CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic5cycle14Par=setEBV(Par, GSModel); Genomic5cycle14F1=randCross(Genomic5cycle14Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic5cycle14F2=self(Genomic5cycle14F1, nProgeny=1); Genomic5cycle14F3=self(Genomic5cycle14F2, nProgeny=25)
Genomic5cycle14F3=setEBV(Genomic5cycle14F3, GSModel); cor(Genomic5cycle14F3@gv, Genomic5cycle14F3@ebv)
Genomic5cycle14F3S=selectInd(Genomic5cycle14F3, nInd=6000, use='ebv'); mean(Genomic5cycle14F3@gv); mean(Genomic5cycle14F3S@gv)
Genomic5cycle14F4=self(Genomic5cycle14F3S, nProgeny=1); Genomic5cycle14F5=self(Genomic5cycle14F4, nProgeny=1)

#Trials 
Genomic5cycle14Pre=setPheno(Genomic5cycle14F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic5cycle14Pre, traits=1, use='pheno')
Genomic5cycle14Pre=setEBV(Genomic5cycle14Pre, GSModel); Genomic5cycle14PreS=selectInd(Genomic5cycle14Pre, nInd=300, use='ebv')
Genomic5cycle14Adv=setPheno(Genomic5cycle14PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic5cycle14Adv, traits=1, use='pheno')
Genomic5cycle14Adv=setEBV(Genomic5cycle14Adv, GSModel); Genomic5cycle14AdvS=selectInd(Genomic5cycle14Adv, nInd=50, use='ebv')
Genomic5cycle14Eli=setPheno(Genomic5cycle14AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic5cycle14Eli, traits=1, use='pheno')
Genomic5cycle14Eli=setEBV(Genomic5cycle14Eli, GSModel); Genomic5cycle14EliS=selectInd(Genomic5cycle14Eli, nInd=10, use='ebv')
Genomic5cycle14Eli1=selectInd(Genomic5cycle14Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic5cycle13PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic5cycle14CS=selectInd(F3, 150, use='ebv'); Genomic5cycle14CS=editGenomeTopQtl(Genomic5cycle14CS, 150, 50)
Genomic5cycle14F3Acc=cor(F3@gv, F3@ebv); Genomic5cycle14F3Acc
mean(Par@gv); mean(Genomic5cycle14CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic5cycle14F3, Genomic5cycle14Pre, Genomic5cycle14EliS)
names(datasets)=c('Genomic5cycle14F3', 'Genomic5cycle14Pre', 'Genomic5cycle14EliS')
Genomic5cycle14=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic5cycle14)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic5cycle14
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic5cycle14PreQtl=al(Genomic5cycle14Pre, SP); Genomic5cycle14EliSQtl= al(Genomic5cycle14EliS, SP)

#Genomic5cycle15
#Training, parentals and F1s
Tr=list(Genomic5cycle13PreS, Genomic5cycle12PreS); Tr=mergePops(Tr)
Par1=Genomic5cycle13PreS; Par1=selectInd(Par1, nInd=150, use='pheno'); Par2=Genomic5cycle13CS; Par=list(Par1, Par2); Par=mergePops(Par)
GSModel=fastRRBLUP(Tr, use='pheno'); Genomic5cycle15Par=setEBV(Par, GSModel); Genomic5cycle15F1=randCross(Genomic5cycle15Par, 1225, balance=true)

#Nurseries and genomic selection F3
Genomic5cycle15F2=self(Genomic5cycle15F1, nProgeny=1); Genomic5cycle15F3=self(Genomic5cycle15F2, nProgeny=25)
Genomic5cycle15F3=setEBV(Genomic5cycle15F3, GSModel); cor(Genomic5cycle15F3@gv, Genomic5cycle15F3@ebv)
Genomic5cycle15F3S=selectInd(Genomic5cycle15F3, nInd=6000, use='ebv'); mean(Genomic5cycle15F3@gv); mean(Genomic5cycle15F3S@gv)
Genomic5cycle15F4=self(Genomic5cycle15F3S, nProgeny=1); Genomic5cycle15F5=self(Genomic5cycle15F4, nProgeny=1)

#Trials 
Genomic5cycle15Pre=setPheno(Genomic5cycle15F5, h2=0.30, reps=5); GSModel=fastRRBLUP(Genomic5cycle15Pre, traits=1, use='pheno')
Genomic5cycle15Pre=setEBV(Genomic5cycle15Pre, GSModel); Genomic5cycle15PreS=selectInd(Genomic5cycle15Pre, nInd=300, use='ebv')
Genomic5cycle15Adv=setPheno(Genomic5cycle15PreS, h2=0.40, reps=40); GSModel=fastRRBLUP(Genomic5cycle15Adv, traits=1, use='pheno')
Genomic5cycle15Adv=setEBV(Genomic5cycle15Adv, GSModel); Genomic5cycle15AdvS=selectInd(Genomic5cycle15Adv, nInd=50, use='ebv')
Genomic5cycle15Eli=setPheno(Genomic5cycle15AdvS, h2=0.40, reps=80); GSModel=fastRRBLUP(Genomic5cycle15Eli, traits=1, use='pheno')
Genomic5cycle15Eli=setEBV(Genomic5cycle15Eli, GSModel); Genomic5cycle15EliS=selectInd(Genomic5cycle15Eli, nInd=10, use='ebv')
Genomic5cycle15Eli1=selectInd(Genomic5cycle15Eli, nInd=1, use='ebv')

#Rapid cycling
Par=Genomic5cycle14PreS
F1=randCross(Par, 1225, balance=true)
F2=self(F1, nProgeny=1); F3=self(F2, nProgeny=12)
GSModel=fastRRBLUP(Par, use='pheno'); F3=setEBV(F3, GSModel); Genomic5cycle15CS=selectInd(F3, 150, use='ebv'); Genomic5cycle15CS=editGenomeTopQtl(Genomic5cycle15CS, 150, 50)
Genomic5cycle15F3Acc=cor(F3@gv, F3@ebv); Genomic5cycle15F3Acc
mean(Par@gv); mean(Genomic5cycle15CS@gv)

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), ebv(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Genomic5cycle15F3, Genomic5cycle15Pre, Genomic5cycle15EliS)
names(datasets)=c('Genomic5cycle15F3', 'Genomic5cycle15Pre', 'Genomic5cycle15EliS')
Genomic5cycle15=do.call(rbind, lapply(datasets, calc_stats)); colnames(Genomic5cycle15)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Genomic5cycle15
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Genomic5cycle15PreQtl=al(Genomic5cycle15Pre, SP); Genomic5cycle15EliSQtl= al(Genomic5cycle15EliS, SP)

################################################################################
#Genomic5 exporting estimates and Qtl

Genomic5cycle1to15PreQtl=cbind(Burnincycle1PreQtl, Burnincycle2PreQtl, Burnincycle3PreQtl, Burnincycle4PreQtl, Burnincycle5PreQtl, Genomic5cycle6PreQtl, Genomic5cycle7PreQtl, Genomic5cycle8PreQtl, Genomic5cycle9PreQtl, Genomic5cycle10PreQtl, Genomic5cycle11PreQtl, Genomic5cycle12PreQtl, Genomic5cycle13PreQtl, Genomic5cycle14PreQtl, Genomic5cycle15PreQtl); write.csv(Genomic5cycle1to15PreQtl, 'r1_Genomic5PreQtl.csv')

Genomic5cycle1to15EliSQtl=cbind(Burnincycle1EliSQtl, Burnincycle2EliSQtl, Burnincycle3EliSQtl, Burnincycle4EliSQtl, Burnincycle5EliSQtl, Genomic5cycle6EliSQtl, Genomic5cycle7EliSQtl, Genomic5cycle8EliSQtl, Genomic5cycle9EliSQtl, Genomic5cycle10EliSQtl, Genomic5cycle11EliSQtl, Genomic5cycle12EliSQtl, Genomic5cycle13EliSQtl, Genomic5cycle14EliSQtl, Genomic5cycle15EliSQtl); write.csv(Genomic5cycle1to15EliSQtl, 'r1_Genomic5EliSQtl.csv')

Genomic5cycle1to15=rbind(Burnincycle1, Burnincycle2, Burnincycle3, Burnincycle4, Burnincycle5, Genomic5cycle6, Genomic5cycle7, Genomic5cycle8, Genomic5cycle9, Genomic5cycle10, Genomic5cycle11, Genomic5cycle12, Genomic5cycle13, Genomic5cycle14, Genomic5cycle15); write.csv(Genomic5cycle1to15, 'r1_Genomic5.csv')

####################################################################################################
#Phenotypic6cycle6

#Parentals, F1s and nurseries
Phenotypic6cycle6Par= Burnincycle3PreS; Phenotypic6cycle6F1=randCross(Phenotypic6cycle6Par, 1225, balance=true)
Phenotypic6cycle6F2=self(Phenotypic6cycle6F1, nProgeny=1); Phenotypic6cycle6F3=self(Phenotypic6cycle6F2, nProgeny=25)
Phenotypic6cycle6F4=self(Phenotypic6cycle6F3, nProgeny=1); Phenotypic6cycle6F5=self(Phenotypic6cycle6F4, nProgeny=1)

#Trials 
Phenotypic6cycle6Pro=setPheno(Phenotypic6cycle6F5, h2=0.10, reps=1); Phenotypic6cycle6ProS=selectWithinFam(Phenotypic6cycle6Pro, nInd=5, use='pheno')
Phenotypic6cycle6Pre=setPheno(Phenotypic6cycle6ProS, h2=0.30, reps=10); Phenotypic6cycle6PreS=selectInd(Phenotypic6cycle6Pre, nInd=300, use='pheno')
Phenotypic6cycle6Adv=setPheno(Phenotypic6cycle6PreS, h2=0.40, reps=40); Phenotypic6cycle6AdvS=selectInd(Phenotypic6cycle6Adv, nInd=50, use='pheno')
Phenotypic6cycle6Eli=setPheno(Phenotypic6cycle6AdvS, h2=0.40, reps=80); Phenotypic6cycle6EliS=selectInd(Phenotypic6cycle6Eli, nInd=10, use='pheno')
Phenotypic6cycle6Eli1=selectInd(Phenotypic6cycle6Eli, nInd=1, use='pheno')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), pheno(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Phenotypic6cycle6Pro, Phenotypic6cycle6Pre, Phenotypic6cycle6EliS)
names(datasets)=c('Phenotypic6cycle6Pro', 'Phenotypic6cycle6Pre', 'Phenotypic6cycle6EliS')
Phenotypic6cycle6=do.call(rbind, lapply(datasets, calc_stats)); colnames(Phenotypic6cycle6)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Phenotypic6cycle6
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Phenotypic6cycle6PreQtl=al(Phenotypic6cycle6Pre, SP); Phenotypic6cycle6EliSQtl= al(Phenotypic6cycle6EliS, SP)

####################################################################################################
#Phenotypic6cycle7

#Parentals, F1s and nurseries
Phenotypic6cycle7Par= Burnincycle4PreS; Phenotypic6cycle7F1=randCross(Phenotypic6cycle7Par, 1225, balance=true)
Phenotypic6cycle7F2=self(Phenotypic6cycle7F1, nProgeny=1); Phenotypic6cycle7F3=self(Phenotypic6cycle7F2, nProgeny=25)
Phenotypic6cycle7F4=self(Phenotypic6cycle7F3, nProgeny=1); Phenotypic6cycle7F5=self(Phenotypic6cycle7F4, nProgeny=1)

#Trials 
Phenotypic6cycle7Pro=setPheno(Phenotypic6cycle7F5, h2=0.10, reps=1); Phenotypic6cycle7ProS=selectWithinFam(Phenotypic6cycle7Pro, nInd=5, use='pheno')
Phenotypic6cycle7Pre=setPheno(Phenotypic6cycle7ProS, h2=0.30, reps=10); Phenotypic6cycle7PreS=selectInd(Phenotypic6cycle7Pre, nInd=300, use='pheno')
Phenotypic6cycle7Adv=setPheno(Phenotypic6cycle7PreS, h2=0.40, reps=40); Phenotypic6cycle7AdvS=selectInd(Phenotypic6cycle7Adv, nInd=50, use='pheno')
Phenotypic6cycle7Eli=setPheno(Phenotypic6cycle7AdvS, h2=0.40, reps=80); Phenotypic6cycle7EliS=selectInd(Phenotypic6cycle7Eli, nInd=10, use='pheno')
Phenotypic6cycle7Eli1=selectInd(Phenotypic6cycle7Eli, nInd=1, use='pheno')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), pheno(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Phenotypic6cycle7Pro, Phenotypic6cycle7Pre, Phenotypic6cycle7EliS)
names(datasets)=c('Phenotypic6cycle7Pro', 'Phenotypic6cycle7Pre', 'Phenotypic6cycle7EliS')
Phenotypic6cycle7=do.call(rbind, lapply(datasets, calc_stats)); colnames(Phenotypic6cycle7)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Phenotypic6cycle7
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Phenotypic6cycle7PreQtl=al(Phenotypic6cycle7Pre, SP); Phenotypic6cycle7EliSQtl= al(Phenotypic6cycle7EliS, SP)

####################################################################################################
#Phenotypic6cycle8

#Parentals, F1s and nurseries
Phenotypic6cycle8Par=Burnincycle5PreS; Phenotypic6cycle8F1=randCross(Phenotypic6cycle8Par, 1225, balance=true)
Phenotypic6cycle8F2=self(Phenotypic6cycle8F1, nProgeny=1); Phenotypic6cycle8F3=self(Phenotypic6cycle8F2, nProgeny=25)
Phenotypic6cycle8F4=self(Phenotypic6cycle8F3, nProgeny=1); Phenotypic6cycle8F5=self(Phenotypic6cycle8F4, nProgeny=1)

#Trials 
Phenotypic6cycle8Pro=setPheno(Phenotypic6cycle8F5, h2=0.10, reps=1); Phenotypic6cycle8ProS=selectWithinFam(Phenotypic6cycle8Pro, nInd=5, use='pheno')
Phenotypic6cycle8Pre=setPheno(Phenotypic6cycle8ProS, h2=0.30, reps=10); Phenotypic6cycle8PreS=selectInd(Phenotypic6cycle8Pre, nInd=300, use='pheno')
Phenotypic6cycle8Adv=setPheno(Phenotypic6cycle8PreS, h2=0.40, reps=40); Phenotypic6cycle8AdvS=selectInd(Phenotypic6cycle8Adv, nInd=50, use='pheno')
Phenotypic6cycle8Eli=setPheno(Phenotypic6cycle8AdvS, h2=0.40, reps=80); Phenotypic6cycle8EliS=selectInd(Phenotypic6cycle8Eli, nInd=10, use='pheno')
Phenotypic6cycle8Eli1=selectInd(Phenotypic6cycle8Eli, nInd=1, use='pheno')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), pheno(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Phenotypic6cycle8Pro, Phenotypic6cycle8Pre, Phenotypic6cycle8EliS)
names(datasets)=c('Phenotypic6cycle8Pro', 'Phenotypic6cycle8Pre', 'Phenotypic6cycle8EliS')
Phenotypic6cycle8=do.call(rbind, lapply(datasets, calc_stats)); colnames(Phenotypic6cycle8)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Phenotypic6cycle8
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Phenotypic6cycle8PreQtl=al(Phenotypic6cycle8Pre, SP); Phenotypic6cycle8EliSQtl= al(Phenotypic6cycle8EliS, SP)

####################################################################################################
#Phenotypic6cycle9

#Parentals, F1s and nurseries
Phenotypic6cycle9Par=Phenotypic6cycle6PreS; Phenotypic6cycle9F1=randCross(Phenotypic6cycle9Par, 1225, balance=true)
Phenotypic6cycle9F2=self(Phenotypic6cycle9F1, nProgeny=1); Phenotypic6cycle9F3=self(Phenotypic6cycle9F2, nProgeny=25)
Phenotypic6cycle9F4=self(Phenotypic6cycle9F3, nProgeny=1); Phenotypic6cycle9F5=self(Phenotypic6cycle9F4, nProgeny=1)

#Trials 
Phenotypic6cycle9Pro=setPheno(Phenotypic6cycle9F5, h2=0.10, reps=1); Phenotypic6cycle9ProS=selectWithinFam(Phenotypic6cycle9Pro, nInd=5, use='pheno')
Phenotypic6cycle9Pre=setPheno(Phenotypic6cycle9ProS, h2=0.30, reps=10); Phenotypic6cycle9PreS=selectInd(Phenotypic6cycle9Pre, nInd=300, use='pheno')
Phenotypic6cycle9Adv=setPheno(Phenotypic6cycle9PreS, h2=0.40, reps=40); Phenotypic6cycle9AdvS=selectInd(Phenotypic6cycle9Adv, nInd=50, use='pheno')
Phenotypic6cycle9Eli=setPheno(Phenotypic6cycle9AdvS, h2=0.40, reps=80); Phenotypic6cycle9EliS=selectInd(Phenotypic6cycle9Eli, nInd=10, use='pheno')
Phenotypic6cycle9Eli1=selectInd(Phenotypic6cycle9Eli, nInd=1, use='pheno')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), pheno(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Phenotypic6cycle9Pro, Phenotypic6cycle9Pre, Phenotypic6cycle9EliS)
names(datasets)=c('Phenotypic6cycle9Pro', 'Phenotypic6cycle9Pre', 'Phenotypic6cycle9EliS')
Phenotypic6cycle9=do.call(rbind, lapply(datasets, calc_stats)); colnames(Phenotypic6cycle9)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Phenotypic6cycle9
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Phenotypic6cycle9PreQtl=al(Phenotypic6cycle9Pre, SP); Phenotypic6cycle9EliSQtl= al(Phenotypic6cycle9EliS, SP)

####################################################################################################
#Phenotypic6cycle10

#Parentals, F1s and nurseries
Phenotypic6cycle10Par=Phenotypic6cycle7PreS; Phenotypic6cycle10F1=randCross(Phenotypic6cycle10Par, 1225, balance=true)
Phenotypic6cycle10F2=self(Phenotypic6cycle10F1, nProgeny=1); Phenotypic6cycle10F3=self(Phenotypic6cycle10F2, nProgeny=25)
Phenotypic6cycle10F4=self(Phenotypic6cycle10F3, nProgeny=1); Phenotypic6cycle10F5=self(Phenotypic6cycle10F4, nProgeny=1)

#Trials 
Phenotypic6cycle10Pro=setPheno(Phenotypic6cycle10F5, h2=0.10, reps=1); Phenotypic6cycle10ProS=selectWithinFam(Phenotypic6cycle10Pro, nInd=5, use='pheno')
Phenotypic6cycle10Pre=setPheno(Phenotypic6cycle10ProS, h2=0.30, reps=10); Phenotypic6cycle10PreS=selectInd(Phenotypic6cycle10Pre, nInd=300, use='pheno')
Phenotypic6cycle10Adv=setPheno(Phenotypic6cycle10PreS, h2=0.40, reps=40); Phenotypic6cycle10AdvS=selectInd(Phenotypic6cycle10Adv, nInd=50, use='pheno')
Phenotypic6cycle10Eli=setPheno(Phenotypic6cycle10AdvS, h2=0.40, reps=80); Phenotypic6cycle10EliS=selectInd(Phenotypic6cycle10Eli, nInd=10, use='pheno')
Phenotypic6cycle10Eli1=selectInd(Phenotypic6cycle10Eli, nInd=1, use='pheno')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), pheno(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Phenotypic6cycle10Pro, Phenotypic6cycle10Pre, Phenotypic6cycle10EliS)
names(datasets)=c('Phenotypic6cycle10Pro', 'Phenotypic6cycle10Pre', 'Phenotypic6cycle10EliS')
Phenotypic6cycle10=do.call(rbind, lapply(datasets, calc_stats)); colnames(Phenotypic6cycle10)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Phenotypic6cycle10
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Phenotypic6cycle10PreQtl=al(Phenotypic6cycle10Pre, SP); Phenotypic6cycle10EliSQtl= al(Phenotypic6cycle10EliS, SP)

####################################################################################################
#Phenotypic6cycle11

#Parentals, F1s and nurseries
Phenotypic6cycle11Par=Phenotypic6cycle8PreS; Phenotypic6cycle11F1=randCross(Phenotypic6cycle11Par, 1225, balance=true)
Phenotypic6cycle11F2=self(Phenotypic6cycle11F1, nProgeny=1); Phenotypic6cycle11F3=self(Phenotypic6cycle11F2, nProgeny=25)
Phenotypic6cycle11F4=self(Phenotypic6cycle11F3, nProgeny=1); Phenotypic6cycle11F5=self(Phenotypic6cycle11F4, nProgeny=1)

#Trials 
Phenotypic6cycle11Pro=setPheno(Phenotypic6cycle11F5, h2=0.10, reps=1); Phenotypic6cycle11ProS=selectWithinFam(Phenotypic6cycle11Pro, nInd=5, use='pheno')
Phenotypic6cycle11Pre=setPheno(Phenotypic6cycle11ProS, h2=0.30, reps=10); Phenotypic6cycle11PreS=selectInd(Phenotypic6cycle11Pre, nInd=300, use='pheno')
Phenotypic6cycle11Adv=setPheno(Phenotypic6cycle11PreS, h2=0.40, reps=40); Phenotypic6cycle11AdvS=selectInd(Phenotypic6cycle11Adv, nInd=50, use='pheno')
Phenotypic6cycle11Eli=setPheno(Phenotypic6cycle11AdvS, h2=0.40, reps=80); Phenotypic6cycle11EliS=selectInd(Phenotypic6cycle11Eli, nInd=10, use='pheno')
Phenotypic6cycle11Eli1=selectInd(Phenotypic6cycle11Eli, nInd=1, use='pheno')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), pheno(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Phenotypic6cycle11Pro, Phenotypic6cycle11Pre, Phenotypic6cycle11EliS)
names(datasets)=c('Phenotypic6cycle11Pro', 'Phenotypic6cycle11Pre', 'Phenotypic6cycle11EliS')
Phenotypic6cycle11=do.call(rbind, lapply(datasets, calc_stats)); colnames(Phenotypic6cycle11)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Phenotypic6cycle11
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Phenotypic6cycle11PreQtl=al(Phenotypic6cycle11Pre, SP); Phenotypic6cycle11EliSQtl= al(Phenotypic6cycle11EliS, SP)

####################################################################################################
#Phenotypic6cycle12

#Parentals, F1s and nurseries
Phenotypic6cycle12Par=Phenotypic6cycle9PreS; Phenotypic6cycle12F1=randCross(Phenotypic6cycle12Par, 1225, balance=true)
Phenotypic6cycle12F2=self(Phenotypic6cycle12F1, nProgeny=1); Phenotypic6cycle12F3=self(Phenotypic6cycle12F2, nProgeny=25)
Phenotypic6cycle12F4=self(Phenotypic6cycle12F3, nProgeny=1); Phenotypic6cycle12F5=self(Phenotypic6cycle12F4, nProgeny=1)

#Trials 
Phenotypic6cycle12Pro=setPheno(Phenotypic6cycle12F5, h2=0.10, reps=1); Phenotypic6cycle12ProS=selectWithinFam(Phenotypic6cycle12Pro, nInd=5, use='pheno')
Phenotypic6cycle12Pre=setPheno(Phenotypic6cycle12ProS, h2=0.30, reps=10); Phenotypic6cycle12PreS=selectInd(Phenotypic6cycle12Pre, nInd=300, use='pheno')
Phenotypic6cycle12Adv=setPheno(Phenotypic6cycle12PreS, h2=0.40, reps=40); Phenotypic6cycle12AdvS=selectInd(Phenotypic6cycle12Adv, nInd=50, use='pheno')
Phenotypic6cycle12Eli=setPheno(Phenotypic6cycle12AdvS, h2=0.40, reps=80); Phenotypic6cycle12EliS=selectInd(Phenotypic6cycle12Eli, nInd=10, use='pheno')
Phenotypic6cycle12Eli1=selectInd(Phenotypic6cycle12Eli, nInd=1, use='pheno')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), pheno(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Phenotypic6cycle12Pro, Phenotypic6cycle12Pre, Phenotypic6cycle12EliS)
names(datasets)=c('Phenotypic6cycle12Pro', 'Phenotypic6cycle12Pre', 'Phenotypic6cycle12EliS')
Phenotypic6cycle12=do.call(rbind, lapply(datasets, calc_stats)); colnames(Phenotypic6cycle12)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Phenotypic6cycle12
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Phenotypic6cycle12PreQtl=al(Phenotypic6cycle12Pre, SP); Phenotypic6cycle12EliSQtl= al(Phenotypic6cycle12EliS, SP)

####################################################################################################
#Phenotypic6cycle13

#Parentals, F1s and nurseries
Phenotypic6cycle13Par=Phenotypic6cycle10PreS; Phenotypic6cycle13F1=randCross(Phenotypic6cycle13Par, 1225, balance=true)
Phenotypic6cycle13F2=self(Phenotypic6cycle13F1, nProgeny=1); Phenotypic6cycle13F3=self(Phenotypic6cycle13F2, nProgeny=25)
Phenotypic6cycle13F4=self(Phenotypic6cycle13F3, nProgeny=1); Phenotypic6cycle13F5=self(Phenotypic6cycle13F4, nProgeny=1)

#Trials 
Phenotypic6cycle13Pro=setPheno(Phenotypic6cycle13F5, h2=0.10, reps=1); Phenotypic6cycle13ProS=selectWithinFam(Phenotypic6cycle13Pro, nInd=5, use='pheno')
Phenotypic6cycle13Pre=setPheno(Phenotypic6cycle13ProS, h2=0.30, reps=10); Phenotypic6cycle13PreS=selectInd(Phenotypic6cycle13Pre, nInd=300, use='pheno')
Phenotypic6cycle13Adv=setPheno(Phenotypic6cycle13PreS, h2=0.40, reps=40); Phenotypic6cycle13AdvS=selectInd(Phenotypic6cycle13Adv, nInd=50, use='pheno')
Phenotypic6cycle13Eli=setPheno(Phenotypic6cycle13AdvS, h2=0.40, reps=80); Phenotypic6cycle13EliS=selectInd(Phenotypic6cycle13Eli, nInd=10, use='pheno')
Phenotypic6cycle13Eli1=selectInd(Phenotypic6cycle13Eli, nInd=1, use='pheno')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), pheno(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Phenotypic6cycle13Pro, Phenotypic6cycle13Pre, Phenotypic6cycle13EliS)
names(datasets)=c('Phenotypic6cycle13Pro', 'Phenotypic6cycle13Pre', 'Phenotypic6cycle13EliS')
Phenotypic6cycle13=do.call(rbind, lapply(datasets, calc_stats)); colnames(Phenotypic6cycle13)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Phenotypic6cycle13
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Phenotypic6cycle13PreQtl=al(Phenotypic6cycle13Pre, SP); Phenotypic6cycle13EliSQtl= al(Phenotypic6cycle13EliS, SP)

####################################################################################################
#Phenotypic6cycle14

#Parentals, F1s and nurseries
Phenotypic6cycle14Par=Phenotypic6cycle11PreS; Phenotypic6cycle14F1=randCross(Phenotypic6cycle14Par, 1225, balance=true)
Phenotypic6cycle14F2=self(Phenotypic6cycle14F1, nProgeny=1); Phenotypic6cycle14F3=self(Phenotypic6cycle14F2, nProgeny=25)
Phenotypic6cycle14F4=self(Phenotypic6cycle14F3, nProgeny=1); Phenotypic6cycle14F5=self(Phenotypic6cycle14F4, nProgeny=1)

#Trials 
Phenotypic6cycle14Pro=setPheno(Phenotypic6cycle14F5, h2=0.10, reps=1); Phenotypic6cycle14ProS=selectWithinFam(Phenotypic6cycle14Pro, nInd=5, use='pheno')
Phenotypic6cycle14Pre=setPheno(Phenotypic6cycle14ProS, h2=0.30, reps=10); Phenotypic6cycle14PreS=selectInd(Phenotypic6cycle14Pre, nInd=300, use='pheno')
Phenotypic6cycle14Adv=setPheno(Phenotypic6cycle14PreS, h2=0.40, reps=40); Phenotypic6cycle14AdvS=selectInd(Phenotypic6cycle14Adv, nInd=50, use='pheno')
Phenotypic6cycle14Eli=setPheno(Phenotypic6cycle14AdvS, h2=0.40, reps=80); Phenotypic6cycle14EliS=selectInd(Phenotypic6cycle14Eli, nInd=10, use='pheno')
Phenotypic6cycle14Eli1=selectInd(Phenotypic6cycle14Eli, nInd=1, use='pheno')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), pheno(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Phenotypic6cycle14Pro, Phenotypic6cycle14Pre, Phenotypic6cycle14EliS)
names(datasets)=c('Phenotypic6cycle14Pro', 'Phenotypic6cycle14Pre', 'Phenotypic6cycle14EliS')
Phenotypic6cycle14=do.call(rbind, lapply(datasets, calc_stats)); colnames(Phenotypic6cycle14)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Phenotypic6cycle14
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Phenotypic6cycle14PreQtl=al(Phenotypic6cycle14Pre, SP); Phenotypic6cycle14EliSQtl= al(Phenotypic6cycle14EliS, SP)

####################################################################################################
#Phenotypic6cycle15

#Parentals, F1s and nurseries
Phenotypic6cycle15Par=Phenotypic6cycle12PreS; Phenotypic6cycle15F1=randCross(Phenotypic6cycle15Par, 1225, balance=true)
Phenotypic6cycle15F2=self(Phenotypic6cycle15F1, nProgeny=1); Phenotypic6cycle15F3=self(Phenotypic6cycle15F2, nProgeny=25)
Phenotypic6cycle15F4=self(Phenotypic6cycle15F3, nProgeny=1); Phenotypic6cycle15F5=self(Phenotypic6cycle15F4, nProgeny=1)

#Trials 
Phenotypic6cycle15Pro=setPheno(Phenotypic6cycle15F5, h2=0.10, reps=1); Phenotypic6cycle15ProS=selectWithinFam(Phenotypic6cycle15Pro, nInd=5, use='pheno')
Phenotypic6cycle15Pre=setPheno(Phenotypic6cycle15ProS, h2=0.30, reps=10); Phenotypic6cycle15PreS=selectInd(Phenotypic6cycle15Pre, nInd=300, use='pheno')
Phenotypic6cycle15Adv=setPheno(Phenotypic6cycle15PreS, h2=0.40, reps=40); Phenotypic6cycle15AdvS=selectInd(Phenotypic6cycle15Adv, nInd=50, use='pheno')
Phenotypic6cycle15Eli=setPheno(Phenotypic6cycle15AdvS, h2=0.40, reps=80); Phenotypic6cycle15EliS=selectInd(Phenotypic6cycle15Eli, nInd=10, use='pheno')
Phenotypic6cycle15Eli1=selectInd(Phenotypic6cycle15Eli, nInd=1, use='pheno')

#Estimates and Qtl
calc_stats=function(data) {corr=cor(gv(data), pheno(data))
data.frame(n=data@nInd, gv=mean(data@gv), ebv=mean(data@ebv), pheno=mean(data@pheno), varG=varG(data), corr=corr)} 
datasets=list(Phenotypic6cycle15Pro, Phenotypic6cycle15Pre, Phenotypic6cycle15EliS)
names(datasets)=c('Phenotypic6cycle15Pro', 'Phenotypic6cycle15Pre', 'Phenotypic6cycle15EliS')
Phenotypic6cycle15=do.call(rbind, lapply(datasets, calc_stats)); colnames(Phenotypic6cycle15)=c('n', 'gv', 'ebv', 'pheno', 'varG', 'corr'); Phenotypic6cycle15
al=function(data, simParam) {Qtl=t(as.data.frame(pullQtlGeno(data, simParam=simParam))); data.frame(al0=rowSums(Qtl==0), al1=rowSums(Qtl==1), al2=rowSums(Qtl==2))}
Phenotypic6cycle15PreQtl=al(Phenotypic6cycle15Pre, SP); Phenotypic6cycle15EliSQtl= al(Phenotypic6cycle15EliS, SP)

####################################################################################################
#Phenotypic6 exporting estimates and Qtl

Phenotypic6cycle1to15PreQtl <- cbind(Burnincycle1PreQtl, Burnincycle2PreQtl, Burnincycle3PreQtl, Burnincycle4PreQtl, Burnincycle5PreQtl, Phenotypic6cycle6PreQtl, Phenotypic6cycle7PreQtl, Phenotypic6cycle8PreQtl, Phenotypic6cycle9PreQtl, Phenotypic6cycle10PreQtl, Phenotypic6cycle11PreQtl, Phenotypic6cycle12PreQtl, Phenotypic6cycle13PreQtl, Phenotypic6cycle14PreQtl, Phenotypic6cycle15PreQtl); write.csv(Phenotypic6cycle1to15PreQtl, 'r1_Phenotypic6PreQtl.csv')

Phenotypic6cycle1to15EliSQtl <- cbind(Burnincycle1EliSQtl, Burnincycle2EliSQtl, Burnincycle3EliSQtl, Burnincycle4EliSQtl, Burnincycle5EliSQtl, Phenotypic6cycle6EliSQtl, Phenotypic6cycle7EliSQtl, Phenotypic6cycle8EliSQtl, Phenotypic6cycle9EliSQtl, Phenotypic6cycle10EliSQtl, Phenotypic6cycle11EliSQtl, Phenotypic6cycle12EliSQtl, Phenotypic6cycle13EliSQtl, Phenotypic6cycle14EliSQtl, Phenotypic6cycle15EliSQtl); write.csv(Phenotypic6cycle1to15EliSQtl, 'r1_Phenotypic6EliSQtl.csv')

Phenotypic6cycle1to15 <- rbind(Burnincycle1, Burnincycle2, Burnincycle3, Burnincycle4, Burnincycle5, Phenotypic6cycle6, Phenotypic6cycle7, Phenotypic6cycle8, Phenotypic6cycle9, Phenotypic6cycle10, Phenotypic6cycle11, Phenotypic6cycle12, Phenotypic6cycle13, Phenotypic6cycle14, Phenotypic6cycle15); write.csv(Phenotypic6cycle1to15, 'r1_Phenotypic6.csv')

Genomic1cycle15
Genomic2cycle15
Genomic3cycle15
Genomic4cycle15
Genomic5cycle15
Phenotypic6cycle15

end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)