###### SINGLE FUNCTIONS ######

###### Parameters estimation ######
ParamEstSeq <- function(DataMatrix,omega)
{
  T=ncol(DataMatrix)
  NExp<-nrow(DataMatrix)
  
  sigmax<-c()
  mi<-c()
  
  for (i in 1:NExp)
  {
    mi[i]<-0
    sigmax[i]<-mad(DataMatrix[i,])^2
  }
  
  smu<-sqrt(omega*sigmax)
  sepsilon<-sqrt((1-omega)*sigmax)
  Results<-list()
  Results$mi<-mi
  Results$smu<-smu
  Results$sepsilon<-sepsilon
  
  Results
}



###### MUK estimation ######
MukEst <- function(DataMatrix,mw)
{
  NExp<-dim(DataMatrix)[1]
  
  if (NExp==1)
  {
    muk<-rbind(seq(-1,1,by=0.1))
  }
  
  if (NExp>1)
  {
    DataMatrix[which(DataMatrix>1)]<-1
    DataMatrix[which(DataMatrix< -1)]<- -1
    
    DataMatrixA<-c()
    for (i in 1:NExp)
    {
      DataMatrixA<-rbind(DataMatrixA,SMA(DataMatrix[i,], n=mw))
    }
    
    DataMatrixA<-DataMatrixA[,2:length(DataMatrixA[1,])]
    
    
    binsize=0.2
    binVec<-c(seq(-1,-0.2,by=binsize),0,seq(0.2,1,by=binsize))
    binmean<-c(seq(-0.9,-0.3,by=binsize),0,0,seq(0.3,0.9,by=binsize))
    
    
    DataQuant<-DataMatrixA
    
    for (i in 1:(length(binVec)-1))
    {
      DataQuant[which(DataMatrixA>binVec[i] & DataMatrixA<=binVec[i+1])]<-binmean[i]
    }
    
    muk<-unique(DataQuant,  MARGIN = 2)
    muk<-muk[,-1]
  }
  muk
}



###### Segmentation core (homogeneous) ######
JointSeg <- function(DataMatrix,eta,omega,muk,mi,smu,sepsilon)
{
  T=ncol(DataMatrix)
  K0<-ncol(muk)
  etav<-log(rep(1,K0)*(1/K0))
  NExp<-nrow(DataMatrix)
  
  P<-matrix(data=0,nrow=K0,ncol=K0)
  G<-matrix(data=0,nrow=K0,ncol=K0)
  emission<-matrix(data=0,nrow=K0,ncol=T)
  out<-.Fortran("transemis",as.vector(muk),as.vector(mi),as.double(eta),as.matrix(DataMatrix),as.integer(K0),as.integer(NExp),as.vector(smu),as.vector(sepsilon),as.integer(T),as.matrix(G),as.matrix(P),as.matrix(emission))
  
  P<-out[[11]]
  emission<-out[[12]]
  
  psi<-matrix(data=0,nrow=K0,ncol=T)
  path<-c(as.integer(rep(0,T)))
  out2<-.Fortran("bioviterbi",as.vector(etav),as.matrix(P),as.matrix(emission),as.integer(T),as.integer(K0),as.vector(path),as.matrix(psi))
  s<-out2[[6]]
  
  sortResult <- SortState(s)
  TotalPredBrek<-sortResult[[3]]
  TotalPredBrek
}



###### Segmentation core (inhomogeneous) ######
JointSegIn <- function(DataMatrix,eta,omega,muk,mi,smu,sepsilon,Pos,stepeta)
{
  CovPos<-diff(Pos)
  CovPosNorm<-CovPos/stepeta
  etavec<-eta+((1-eta)*exp(log(eta)/CovPosNorm))
  
  NCov<-length(etavec)
  K0<-ncol(muk)
  etav<-log(rep(1,K0)*(1/K0))
  T=ncol(DataMatrix)
  NExp<-nrow(DataMatrix)
  
  P<-matrix(data=0,nrow=K0,ncol=(K0*NCov))
  G<-matrix(data=0,nrow=K0,ncol=K0)
  emission<-matrix(data=0,nrow=K0,ncol=T)
  out<-.Fortran("transemisi",as.vector(muk),as.vector(mi),as.double(etavec),as.integer(NCov),as.matrix(DataMatrix),as.integer(K0),as.integer(NExp),as.vector(smu),as.vector(sepsilon),as.integer(T),as.matrix(G),as.matrix(P),as.matrix(emission))
  
  P<-out[[12]]
  emission<-out[[13]]
  
  psi<-matrix(data=0,nrow=K0,ncol=T)
  path<-c(as.integer(rep(0,T)))
  out2<-.Fortran("bioviterbii",as.vector(etav),as.matrix(P),as.matrix(emission),as.integer(T),as.integer(K0),as.vector(path),as.matrix(psi))
  s<-out2[[6]]
  
  sortResult <- SortState(s)
  TotalPredBreak<-sortResult[[3]]
  TotalPredBreak
}



###### State sorting ######
SortState <- function(s)
{
  l<-1
  seg<-c()
  brek<-c()
  t<-1
  for (k in 1:(length(s)-1))
  {
    if (s[k]!=s[k+1])
    {
      brek[t]<-k
      t<-t+1
      if (length(which(seg==s[k]))==0)
      {
        seg[l]<-s[k]
        l<-l+1
      }
    }
  }
  
  brek<-c(0,brek,length(s))
  
  if (length(which(seg==s[length(s)]))==0)
  {
    seg<-c(seg,s[length(s)])
  }
  
  s0<-c()
  
  for (k in 1:length(seg))
  {
    s0[which(s==seg[k])]<-k
  }
  
  SortResult<-list()
  SortResult[[1]]<-s0
  SortResult[[2]]<-seg
  SortResult[[3]]<-brek
  SortResult
}



###### Small shifts filtering ######
FilterSeg <- function(TotalPredBreak,FW)
{
  controllength<-diff(TotalPredBreak)
  
  indF<-which(controllength<=FW)
  
  if (length(indF)!=0)
  {
    if (indF[1]==1)
    {
      indF[1]<-2
      indF<-unique(indF)
      TotalPredBreak1<-TotalPredBreak[-(indF)]
    }
    if (indF[1]!=1)
    {
      TotalPredBreak1<-TotalPredBreak[-(indF)]
    }
  }
  
  if (length(indF)==0)
  {
    TotalPredBreak1<-TotalPredBreak
  }
  TotalPredBreak1
}



###### Calculating segment mean after segmentation ######
SegResults <- function(DataSeq,TotalPredBreak)
{
  TotalPred<-c()
  NExp<-nrow(DataSeq)
  
  for (j in 1:NExp)
  {
    s<-rep(0,ncol(DataSeq))
    for (i in 1:(length(TotalPredBreak)-1))
    {
      s[(TotalPredBreak[i]+1):TotalPredBreak[i+1]]<-median(DataSeq[j,(TotalPredBreak[i]+1):TotalPredBreak[i+1]])
    }
    TotalPred<-rbind(TotalPred,s)
  }
  
  Result<-TotalPred
  Result
}



###### COMPLETE PIPELINES ######

###### Homogeneous segmentation ######
SLM <- function(log2r_data, omega, eta, FW)
{
  mw<-1
  LogDataNorm<-rbind(log2r_data)
  ParamList<-ParamEstSeq(LogDataNorm,omega)
  mi<-ParamList$mi
  smu<-ParamList$smu
  sepsilon<-ParamList$sepsilon
  muk<-MukEst(LogDataNorm,mw)
  PredBreak<-JointSeg(LogDataNorm,eta,omega,muk,mi,smu,sepsilon)
  PredBreak1<-FilterSeg(PredBreak,FW)
  DataSeg<-SegResults(LogDataNorm,PredBreak1)
  
  DataSeg
}



###### Heterogeneous segmentation ######
HSLM <- function(log2r_data, pos_data, omega, eta, stepeta, FW)
{
  mw<-1
  LogDataNorm<-rbind(log2r_data)
  ParamList<-ParamEstSeq(LogDataNorm,omega)
  mi<-ParamList$mi
  smu<-ParamList$smu
  sepsilon<-ParamList$sepsilon
  muk<-MukEst(LogDataNorm,mw)
  PredBreak<-JointSegIn(LogDataNorm,eta,omega,muk,mi,smu,sepsilon,pos_data,stepeta)
  PredBreak1<-FilterSeg(PredBreak,FW)
  DataSeg<-SegResults(LogDataNorm,PredBreak1)
  
  DataSeg
}
