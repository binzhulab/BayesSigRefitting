SigRefitting <- function(BSR.obj, sig.select=0.8) {

  check_BSR.object(BSR.obj) 
  check_sig(sig.select)
  ret <- bsr_refit(BSR.obj$objects$select$MC.gamma, BSR.obj$objects$V, 
                   BSR.obj$objects$Sig.W, list(sig.select=sig.select))
  ret 
}

bsr_refit <- function(v1, V, w.s, obj) {

  v.s.b<-c(V)
  p<-dim(w.s)[2]  # full model size
  n<-dim(V)[2] # sample size

  nc.v1  <- ncol(v1)
  nc.w.s <- ncol(w.s)
  if (nc.v1 == nc.w.s + 2) {
    v1.m <- as.matrix(v1[, -c(1, nc.v1), drop=FALSE])
  } else if (nc.v1 == nc.w.s + 1) {
    v1.m <- as.matrix(v1[, -nc.v1, drop=FALSE])
  } else {
    stop("ERROR with object dimensions")
  }

  colnames(v1.m)<-colnames(w.s) 
  v1.PP<-v1[,dim(v1)[2]]-min(v1[,dim(v1)[2]])
  v1.rf<-apply(v1.m*v1.PP, 2, sum)/sum(v1.PP)
  v1.rf<-apply(v1.m, 2, mean)
  Selected.Sig<-which(v1.rf >= obj$sig.select)
  if (!length(Selected.Sig)) {
    stop("No signatures selected, try adjusting sig.select option")
  }
  W         <- w.s[,Selected.Sig, drop=FALSE]

  #### Signature refitting
  MC.gamma.Selected<-rep(0,p)
  MC.gamma.Selected[Selected.Sig]<-1

  # signature refitting based on the selected signatures
  Refit.Results <- log.p.gamma.v1(MC.gamma.Selected, w.s, v.s.b) 

  ### compute the signature activity matrix
  H<-matrix(Refit.Results[[1]][-1], length(Selected.Sig), n)
  row.names(H)<-colnames(W)
  colnames(H)<-colnames(V)
  
  list(selection=Selected.Sig, W=W, H=H)
}


#### Signature Refitting Function-------------------------------------------
log.p.gamma.v1<-function(b, w.s, v.s.b){
 
  #n <- ncol(v.s.b)
  n <- length(v.s.b)/96

  g<-which(b!=0)
  # diagonal block matrix of w 
  MC.gamma.List<-list()
  for(i in 1:n){
    MC.gamma.List[[i]]<-w.s[,g]
  }
  w.b<-as.matrix(bdiag(MC.gamma.List))
  if(sum(g)>0){
    
    # lasso poisson regression
    cvfit<-glmnet(w.b,v.s.b,lower.limits=0,
                  family=poisson(link = "identity"),
                  nlambda=10,nfolds=5,standardize = FALSE)
    
    # computing the signature activities and mean of mutation
    beta.l<-as.vector(coef(cvfit,min(cvfit$lambda)))
    lambda.gamma<-w.b%*%beta.l[-1]+10^(-7)
    
    # computing the BIC
    BIC<-(-2)*sum(dpois(v.s.b,lambda = as.vector(lambda.gamma), log = TRUE))+n*length(g)*log(96*n)
    
  }else{
    beta.l<-0
    BIC<-0
  }
  return(list(beta.l,BIC)) 
}
