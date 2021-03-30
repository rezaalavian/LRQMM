lrqmm_m<-function(id,sire,dam,X,Y,cova=NULL,alpha=0,tau=0.5,Factor=FALSE,maxTries=3000,interval=30){
  start.time <- Sys.time()
  data<-data.frame(id,sire,dam,X,Y) #preparation data befor using Pedigree function
  Ped<-GeneticsPed::Pedigree(x=data,subject="id",ascendant=c("sire","dam"))
  Y=as.matrix(Ped[5]);n<-dim(Y)[1]#preparation response variable
  Ped$id<-factor(Ped$id)
  X<-as.matrix(X);XN<-as.matrix(as.numeric(sub(",", ".", row.names(table(X)), fixed = TRUE)))
  if(Factor==TRUE && dim(X)[2]>=2){
    X1<-Matrix::sparse.model.matrix(~(factor(X[,1]))-1)
    for(i in 2:dim(X)[2]){
      X2<-Matrix::sparse.model.matrix(~(factor(X[,i]))-1)
      X1<-cbind(X1,X2)}
    X<-X1}
  if(Factor==TRUE&&dim(X)[2]==1){X<-factor(X);X<-Matrix::sparse.model.matrix(~X-1)}
  if(Factor==FALSE){X<-Matrix::Matrix(X,sparse = TRUE)}#preparation fixed effect
  ped<-Ped;Ped=GeneticsPed::extend(Ped);m<-dim(Ped)[1];recordPed<-as.matrix(as.numeric(sub(",", ".", Ped[,1], fixed = TRUE)))
  Z<-cbind(Matrix::Matrix(0,n,m-n,sparse=T),Matrix::sparse.model.matrix(~factor(ped$id)-1))
  Ped<-Ped[order(kinship2::kindepth(Ped[, 1], Ped[, 2], Ped[, 3]),decreasing = FALSE),]
  Ped<-as.data.frame(Ped)
  Ainv<-MCMCglmm::inverseA(Ped[,1:3])$`Ainv`;Ainv<-Matrix::Matrix(Ainv,sparse=TRUE)
  Z.t.inv<-Z #inverse of t(Z) is equal to Z.
  random<-Z+((Z.t.inv%*%Ainv)*alpha);dimZ2<-dim(Z)[2]
  E<-cbind(X,random)
  if(!is.null(cova)){cova<-Matrix::Matrix(cova,sparse=TRUE);E<-cbind(E,cova)};E<-as.matrix(E)
  SVD<-SVDmat(E,maxTries,interval)
  model<-quantreg::rq.fit(SparseM::as.matrix.csr(SVD$u),Y,method="sfn",tau=tau)
  coef<-model$coef;coef<-SVD$v%*%(solve(diag(SVD$d))%*%as.matrix(coef))
  coef<-STDE(coef,Y=Y,E=E,SVD=SVD,tau = tau,n = n)
  resi<-matrix(model$res); colnames(resi)<-c("residuals' value")
  fix.effect<-coef[1:dim(X)[2]]# estimate fixed effect(s)
  ifelse(Factor==TRUE,fix.effect.d<-cbind(XN,as.matrix(coef[1:dim(X)[2],1:2])),fix.effect.d<-as.matrix(coef[1:dim(X)[2],1:2]))
  fix.effect.d<-as.data.frame(fix.effect.d);  colnames(fix.effect.d)[1:2]<-c("fix.effect'names","fix.effect'value");row.names(fix.effect.d)<-c()
  random.effect<-coef[(dim(X)[2]+1):(dim(X)[2]+dim(random)[2])]# estimate random effects
  random.effect.d<-as.data.frame(cbind(recordPed,as.matrix(coef[(dim(X)[2]+1):(dim(X)[2]+dim(random)[2]),1:2])))
  colnames(random.effect.d)[1:2]<-c("Record's Ped","random.effects");row.names(random.effect.d)<-c()
  if(!is.null(cova)){cova.effect<-coef[(dim(X)[2]+dim(random)[2]+1):(dim(E)[2]),1:2]}
  MAE<-mean(abs(resi))
  end.time <- Sys.time()
  if(is.null(cova)){ans<-list(fix.effect=fix.effect.d,random.effect=random.effect.d,residuals=resi,Time_between_start_to_end=end.time-start.time)}
  if(!is.null(cova)){ans<-list(fix.effect=fix.effect.d,cova.effect=cova.effect,random.effect=random.effect.d,residuals=resi,Time_between_start_to_end=end.time-start.time)}
  append(ans,ans$summary<-c("estimated effects in quantile:"=tau, "MAE:"=MAE, "Var(response):"=stats::var(Y), "Var(pedigree's random.effect):"=stats::var(random.effect), "Var(record's random.effect):"=stats::var(random.effect[(m-n+1):m]) ,"observaions:"=n,"pedigree's length:"=m, "fix.effect.lavel:"=dim(X)[2],"random.effect.lavel:"=dimZ2))
  return(sapply(ans,round,4))}

