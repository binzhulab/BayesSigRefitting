BayesSigRef <- function(V, Sig.W="SBS", n.iter=3000, tol=0.05, n.core=1,
                        sig.select=0.8) {

  n.iter <- check_int(n.iter, "n.iter", min=2)
  n.core <- check_int(n.core, "n.core")
  check_tol(tol)
  check_sig(sig.select)
  obj    <- list(n.iter=n.iter, tol=tol, n.core=n.core, sig.select=sig.select)
  check_mat(V, "V")
  check_Sig.W(Sig.W)

  ret <- bsr_main(V, Sig.W, obj) 
  ret
}

bsr_main <- function(V, Sig.W, obj) {

  select       <- bsr_select(V, Sig.W, obj)
  Sig.W        <- select$Sig.W
  select$Sig.W <- NULL
  refit        <- try(bsr_refit(select$MC.gamma, V, Sig.W, obj))
  ret          <- bsr_getRetObject(select, refit, V, Sig.W, obj) 
}

bsr_getRetObject <- function(select, refit, V, Sig.W, obj) {

  # Get selection vector
  if (!("try-error" %in% class(refit))) {
    W         <- refit$W
    H         <- refit$H
    selection <- refit$selection
    if (length(selection) == ncol(W)) names(selection) <- colnames(W)
  } else {
    W         <- NULL
    H         <- NULL
    selection <- NULL
  }

  objects <- list(V=V, Sig.W=Sig.W, options=obj, select=select)
  ret <- list(selection=selection, refitting=list(W=W, H=H), 
              objects=objects)
}

bsr_select <- function(V, w.s, obj) {

  if (isString(w.s)) w.s <- getDefSig.W()
  v.s.b <- c(V)

  #### Shotgun Stochastic Search------------------------------------------------

  ncores <- obj$n.core 
  MC.rep <- obj$n.iter

  n<-dim(V)[2] ## simulated sample size 
  p<-dim(w.s)[2]  # full model size
  MC.Sd<-rep(0,MC.rep) # matrix to save the change of the probability
  MC.Pi.Max<-rep(0,MC.rep) # matrix to save the max marginal probability of each iteration

  # generating the initial model
  MC.gamma.Min0<-rep(0,p)
  s.index<-sample(1:p,5,replace = FALSE)
  MC.gamma.Min0[s.index]<-1

  # computing the marginal likelihood of the initial model
  Pi.gamma.Max<-log.p.gamma.v(MC.gamma.Min0, w.s, v.s.b)[[1]] 

  conv  <- 0
  core1 <- ncores == 1
  for (i in 1:MC.rep) {
    if (core1) cat(paste0("Iteration ", i, "\n")) 
    if(i==1){ 
    
      MC.gamma.Min<-MC.gamma.Min0 
      MC.a<-which(MC.gamma.Min!=0)
      MC.b<-which(MC.gamma.Min==0)
      MC.gamma.T<-matrix(rep(0,3*(p+1)),3,p+1)
      Pi.check<-rep(0,3)
    
      ## Models obtained by deleting signatures
      #  Generating the candidate model within the set
      K<-length(MC.a)
      MC.gamma.D<-matrix(rep(0,K*(p+1)),K,p+1)
      for (j in 1:K) {
        MC.gamma.MinJ<-MC.gamma.Min
        MC.gamma.MinJ[MC.a[j]]<-1-MC.gamma.MinJ[MC.a[j]]
        MC.gamma.D[j,-(p+1)]<-MC.gamma.MinJ
      }
    
      # Computing the model scores
      List.D<-as.list(as.data.frame(t(MC.gamma.D[,-(p+1)])))
      log.p.D<-mclapply(List.D, log.p.gamma.v, w.s=w.s, v.s.b=v.s.b, mc.cores=ncores) 
      for (j in 1:K) {
        MC.gamma.D[j,p+1]<-log.p.D[[j]][[1]]
      }
    
      # Sampling one model from the deleting model set
      if(length(MC.a)==1){if(MC.gamma.D[,p+1]=="-Inf"){MC.gamma.D[,p+1]<--10^10}}
      Log.Pi.D<-exp(MC.gamma.D[,p+1]-max(MC.gamma.D[,p+1]))
      Prob.D<-Log.Pi.D/sum(Log.Pi.D) # normalized model scores in the neighborhood
      z.D<-sample(1:K,1,prob = Prob.D)
      MC.gamma.T[1,]<-MC.gamma.D[z.D,]
    
      ## Models obtained by adding signatures
      # Generating the candidate models within the set
      M<-length(MC.b)
      MC.gamma.A<-matrix(rep(0,M*(p+1)),M,p+1)
      for (m in 1:M) {
        MC.gamma.MinJ<-MC.gamma.Min
        MC.gamma.MinJ[MC.b[m]]<-1-MC.gamma.MinJ[MC.b[m]]
        MC.gamma.A[m,-(p+1)]<-MC.gamma.MinJ   
      }
    
      # Compute model scores
      List.A<-as.list(as.data.frame(t(MC.gamma.A[,-(p+1)])))
      log.p.A<-mclapply(List.A,log.p.gamma.v, w.s=w.s, v.s.b=v.s.b, mc.cores=ncores) 
      for (m in 1:M) {
        MC.gamma.A[m,p+1]<-log.p.A[[m]][[1]]
      }
    
      # Sampling one model from the adding model set
      Log.Pi.A<-MC.gamma.A[,p+1]
      Log.Pi.A<-exp(MC.gamma.A[,p+1]-max(MC.gamma.A[,p+1]))
      Prob.A<-Log.Pi.A/sum(Log.Pi.A) # normalized model scores in the neighborhood
      z.A<-sample(1:M,1,prob = Prob.A)
      MC.gamma.T[2,]<-MC.gamma.A[z.A,]
    
      ## Models obtained by swapping signatures
      # Generating the candidate models within the set
      S<-length(MC.a)*length(MC.b)
      MC.gamma.S<-matrix(rep(0,S*(p+1)),S,p+1)
      for (r in 1:length(MC.a)) {
        check.h<-MC.gamma.Min
        check.h[MC.a[r]]<-1-check.h[MC.a[r]]
        for (m in 1:length(MC.b)){
          check.hh<-check.h
          check.hh[MC.b[m]]<-1-check.hh[MC.b[m]]
          MC.gamma.S[(r-1)*length(MC.b)+m,-(p+1)]<-check.hh
        }
      }
    
      # Computing model scores
      List.S<-as.list(as.data.frame(t(MC.gamma.S[,-(p+1)])))
      log.p.S<-mclapply(List.S,log.p.gamma.v, w.s=w.s, v.s.b=v.s.b, mc.cores=ncores) 
      for (rm in 1:S) {
        MC.gamma.S[rm,p+1]<-log.p.S[[rm]][[1]]
      }
    
      # Sampling one model from the swapping model set   
      Log.Pi.S<-MC.gamma.S[,p+1]
      Log.Pi.S<-exp(MC.gamma.S[,p+1]-max(MC.gamma.S[,p+1]))
      Prob.S<-Log.Pi.S/sum(Log.Pi.S) # normalized model scores in the neighborhood
      z.S<-sample(1:S,1,prob = Prob.S)
      MC.gamma.T[3,]<-MC.gamma.S[z.S,]
    
      ### Updata Gamma set
      MC.gamma.C<-rbind(MC.gamma.D,MC.gamma.A,MC.gamma.S) # combine the models
      MC.gamma.C<-MC.gamma.C[order(MC.gamma.C[,p+1],decreasing = TRUE),]
      #MC.gamma<-MC.gamma.C[1:100,] 
      MC.gamma<-MC.gamma.C[1:min(100, nrow(MC.gamma.C)), , drop=FALSE] 
    }else{ 
    
      MC.gamma.Min<-check 
      MC.a<-which(MC.gamma.Min!=0)
      MC.b<-which(MC.gamma.Min==0)
      MC.gamma.T<-matrix(rep(0,3*(p+1)),3,p+1)
      Pi.check<-rep(0,3)
    
      ## Models obtained by deleting signatures
      K<-length(MC.a)
      MC.gamma.D<-matrix(rep(0,K*(p+1)),K,p+1)
      for (j in 1:K) {
        MC.gamma.MinJ<-MC.gamma.Min
        MC.gamma.MinJ[MC.a[j]]<-1-MC.gamma.MinJ[MC.a[j]]
        MC.gamma.D[j,-(p+1)]<-MC.gamma.MinJ
      }
    
      ## Models obtained by adding signatures
      M<-length(MC.b)
      MC.gamma.A<-matrix(rep(0,M*(p+1)),M,p+1)
      for (m in 1:M) {
        MC.gamma.MinJ<-MC.gamma.Min
        MC.gamma.MinJ[MC.b[m]]<-1-MC.gamma.MinJ[MC.b[m]]
        MC.gamma.A[m,-(p+1)]<-MC.gamma.MinJ
      }
    
      ## Models obtained by swapping signatures
      S<-length(MC.a)*length(MC.b)
      MC.gamma.S<-matrix(rep(0,S*(p+1)),S,p+1)
      for (r in 1:length(MC.a)) {
        check.h<-MC.gamma.Min
        check.h[MC.a[r]]<-1-check.h[MC.a[r]]
        for (m in 1:length(MC.b)){
          check.hh<-check.h
          check.hh[MC.b[m]]<-1-check.hh[MC.b[m]]
          MC.gamma.S[(r-1)*length(MC.b)+m,-(p+1)]<-check.hh
        }
      }
    
      ## mathcing with the all the models in the total model set   
      MC.gamma.B<-rbind(MC.gamma.D,MC.gamma.A,MC.gamma.S) # combine the 3 sets of models
      Index.C<-row.match(as.data.frame(MC.gamma.B[,-(p+1)]),as.data.frame(MC.gamma.C[,-(p+1)]))
      In.C<-Index.C[which(Index.C!="NA")]
      In.B<-which(Index.C!="NA")
      MC.gamma.B[In.B,p+1]<-MC.gamma.C[In.C,p+1]
    
      ## Computing model score for all the models with no match in the total model set
      MC.B<-which(MC.gamma.B[,p+1]==0)
      if(length(MC.B)!=0){
        List.B<-as.list(as.data.frame(t(MC.gamma.B[MC.B,-(p+1)])))
        log.p.B<-mclapply(List.B,log.p.gamma.v, w.s=w.s, v.s.b=v.s.b, mc.cores=ncores) 
        for (mb in 1:length(List.B)){
          MC.gamma.B[MC.B[mb],p+1]<-log.p.B[[mb]][[1]]
        }
      }
    
      ## Sampling one model from the deleting models
      MC.gamma.D[,p+1]<-MC.gamma.B[1:length(MC.a),p+1]
      if(length(MC.a)==1){if(MC.gamma.D[,p+1]=="-Inf"){MC.gamma.D[,p+1]<--10^10}}
      # normalizing model scores in the neighborhood
      Log.Pi.D<-exp(MC.gamma.D[,p+1]-max(MC.gamma.D[,p+1]))
      Prob.D<-Log.Pi.D/sum(Log.Pi.D) 
      # sampling one model according to the normalized model scores
      z.D<-sample(1:K,1,prob = Prob.D)
      MC.gamma.T[1,]<-MC.gamma.D[z.D,]
    
      ## Sampling one model from the adding model set
      MC.gamma.A[,p+1]<-MC.gamma.B[(K+1):p,p+1]
      # normalized model scores in the neighborhood
      Log.Pi.A<-exp(MC.gamma.A[,p+1]-max(MC.gamma.A[,p+1]))
      Prob.A<-Log.Pi.A/sum(Log.Pi.A)
      # sampling one model according to the normalized model scores
      z.A<-sample(1:M,1,prob = Prob.A)
      MC.gamma.T[2,]<-MC.gamma.A[z.A,]
    
      ## Sampling one model from the swapping model set
      dim.B<-dim(MC.gamma.B)[1]
      MC.gamma.S[,p+1]<-MC.gamma.B[(p+1):dim.B,p+1]
      # normalized model scores in the neighborhood
      Log.Pi.S<-exp(MC.gamma.S[,p+1]-max(MC.gamma.S[,p+1]))
      Prob.S<-Log.Pi.S/sum(Log.Pi.S) 
      # sampling one model according to the normalized model scores
      z.S<-sample(1:S,1,prob = Prob.S)
      MC.gamma.T[3,]<-MC.gamma.S[z.S,]
    
      ## Updating Gamma model set & total model set   
      MC.gamma.BB<-rbind(MC.gamma.B,MC.gamma)
      MC.gamma.BB<-MC.gamma.BB[order(MC.gamma.BB[,p+1],decreasing = TRUE),]

      #MC.gamma<-MC.gamma.BB[!duplicated(MC.gamma.BB),][1:100,]
      tmp      <- MC.gamma.BB[!duplicated(MC.gamma.BB), , drop=FALSE]
      MC.gamma <- tmp[1:min(100, nrow(tmp)), , drop=FALSE] 
      tmp      <- NULL
 
      MC.gamma.C<-rbind(MC.gamma.C,MC.gamma.B)
      MC.gamma.C<-MC.gamma.C[!duplicated(MC.gamma.C),]
    
    }
  
    ## Sampling one model from the final model set 
    Log.Pi.T<-exp(MC.gamma.T[,p+1]-max(MC.gamma.T[,p+1]))
    Prob.T<-Log.Pi.T/sum(Log.Pi.T)
    z.T<-sample(1:3,1,prob = Prob.T)
    check<-MC.gamma.T[z.T,-(p+1)]
  
    ## Monitoring convergence and output Gamma set
    if(i==1){
      tol<-1 # the change of the total of the marginal likelihood 
      Sum.Pi.gamma<-sum(MC.gamma[,p+1])
    }else{
      tol<-abs((sum(MC.gamma[,p+1])-Sum.Pi.gamma)/Sum.Pi.gamma) 
      Sum.Pi.gamma<-sum(MC.gamma[,p+1])
    }
  
    Pi.gamma.Max<-max(MC.gamma[,p+1]) # updating the max marginal likelihood
    MC.Sd[i]<-tol # saving the change of the total marginal likelihood
    MC.Pi.Max[i]<-Pi.gamma.Max # saving the max marginal likelihood
    # stop the algorithm when the convergence is achieved
    if(tol<=0.01){ 
      conv <- 1
      break 
      #return(MC.gamma,MC.Sd,MC.Pi.Max)
    }else{ 
      Pi.gamma.Max<-max(MC.gamma[,p+1]) 
    }
  } 
  msg <- ""
  if (!conv) {
    msg <- "Algorithm did not converge, ran out of iterations"
    warning(msg)
  }

  MC.Sd     <- MC.Sd[1:i]
  MC.Pi.Max <- MC.Pi.Max[1:i]
  gc()

  list(MC.gamma=MC.gamma, MC.Sd=MC.Sd, MC.Pi.Max=MC.Pi.Max, Sig.W=w.s, msg=msg)
}

#### Posterior Probability Function-------------------------------------------
log.p.gamma.v<-function(b, w.s, v.s.b){

  n <- length(v.s.b)/96
  g <- which(b!=0)
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
    # computing the value of lambda and hyperparameter s
    beta.l<-as.vector(coef(cvfit,min(cvfit$lambda)))
    lambda.gamma<-w.b%*%beta.l[-1]+10^(-7)
    s<-1/(cvfit$lambda.min)
    # computing the value of hessian matrix
    ev<-eigen(t(w.s[,g])%*%w.s[,g])
    max.en<-median(as.vector((v.s.b+0.0000001)/(lambda.gamma*lambda.gamma)))
    det.gamma.1<-prod((ev$values)*max.en)
    det.gamma<-n*log(det.gamma.1)
    # computing the marginal likelihood
    out<-sum(dpois(v.s.b,lambda = lambda.gamma, log = TRUE))-
         sum(n*length(g)*log(2*s))-sum((1/s)*sum(abs(beta.l)))-
         0.5*det.gamma+((n*length(g))/2)*log(2*pi)
  }else{
    out<-sum(v.s.b*log(0.0000001)-factorial(v.s.b))
    beta.l<-0
  }
  return(list(out,beta.l,BIC)) 
}

getDefSig.W <- function() {

    Sig.W <- NULL
    dir <- system.file("data", package="BayesSigRefitting", mustWork=TRUE)
    f   <- file.path(dir, "Sig.W.rda")
    load(f)
    Sig.W
}
