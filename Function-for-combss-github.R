noise.compute.from.SNR <- function(SNR,beta,sigma) {
  sigma.2 <- t(matrix(beta,ncol=1))%*%sigma%*%matrix(beta,ncol=1)/SNR
  return(sigma.2)
}

cov.X <- function(p,c){
  mat.cov <- diag(p)
  for (i in 1:p){
    for (j in i:p){
      mat.cov[i,j] <- c**(abs(i-j))
      mat.cov[j,i] <- mat.cov[i,j]
    }
  }
  return(mat.cov)
}


map.w.to.t <- function(w) 1-exp(-w**2)
map.t.to.w <- function(t) (-log(1-t))**(0.5)


map_to_combss_object <- function(x,tau=0.5,n,X,y){
  t <- x[[1]]
  s <- as.vector((t>tau))
  beta.hat.s <- rep(0,dim(X)[1])
  
  if(sum(s)==0|sum(s)>n){beta.hat.s <- 0}else{
    beta.hat.s <- solve((t(X[,s])%*%X[,s]))%*%t(X[,s])%*%matrix(y,ncol=1)}
  
  index.cov <- which(s==TRUE)
  result <- list(index.cov=index.cov,s=s,t=t,Niter=x[[4]],beta.hat.s=beta.hat.s,convergence=x[[3]])
  class(result) <- "COMBSS"
  return(result)
}



predict.COMBSS <- function(object,newdata){
  if(sum(object$s)==0|sum(object$beta.hat.s)==0){y.pred <- rep(0,dim(newdata)[1])}else{
    Xpred <- newdata[,object$s] 
    y.pred <- Xpred%*%object$beta.hat.s}
  return(y.pred)
}



ADAM.COMBSS <- function(X,y,delta=n,lambda,tau=0.5,Niter=200,alpha=0.001,psy=c(0.9,0.999),epoch=10,tol=0.0001,CG=TRUE,trunc=0.001){
  n <- dim(X)[1]
  p <- dim(X)[2]
  indice <- 1:p
  Xnew <- X[,indice]
  XtXnew  <- (1/n)*crossprod(Xnew)
  Xtynew <- (1/n)*crossprod(Xnew,matrix(y,ncol=1))
  
  t0 <- rep(0.5,p)
  w <- map.t.to.w(t0)
  c <- 10e-8
  #w <- w - alpha*u/sqrt(v+c)
  u <- rep(0,p)
  v <- rep(0,p)
  u.tilde <- rep(0,p)
  v.tilde <- rep(0,p)
  t.new <- t0
  i <- 0
  counter <- 0
  pnew <- p
  wnew <- w
  p.new.iter <- p
  while (i<Niter&counter< epoch){
  
    if(pnew<p.new.iter){
      Xnew <- X[,indice]
      pnew <- length(indice)
      XtXnew  <- (1/n)*crossprod(Xnew)
      Xtynew <- (1/n)*crossprod(Xnew,matrix(y,ncol=1))
    }
    wnew <- w[indice]
    p.new.iter <- pnew
    i <- i+1
    t.old <- t.new 
    gradw <- gradient.Combss.Im(w=wnew,p=pnew,X=Xnew,XtXn=XtXnew,Xtyn=Xtynew,delta=n,lambda=lambda,CG=CG,n=n)
  
    u[indice] <-psy[1]*u[indice]-(1-psy[1])*gradw
    v[indice] <-psy[2]*v[indice]+(1-psy[2])*(gradw*gradw)
    
    u.tilde[indice] <- u[indice]/(1-psy[1]**i)
    v.tilde[indice] <- v[indice]/(1-psy[2]**i)
    w[indice] <- w[indice] + alpha*u.tilde[indice]/sqrt(v.tilde[indice]+c)
    
    t.new <- map.w.to.t(w)
    if(!is.null(trunc)){
    if(any(t.new[indice]<trunc)){
      indice <- which(t.new>trunc)
      t.new[-c(indice)] <- 0
      pnew <- length(indice)
    }
    if(is.null(pnew)) break
    }
    
    maxnorm <- max(abs(t.old-t.new))
    if(maxnorm < tol){counter <- counter + 1} else counter <-0
  }
  
  
  if(i<Niter) convergence <- TRUE else convergence <- FALSE
  t <- map.w.to.t(w)
  s <- as.vector((t>tau))
  
  if(sum(s)==0|sum(s)>n) beta.hat.s <- 0 else{
    beta.hat.s <- solve((t(X[,s])%*%X[,s]))%*%t(X[,s])%*%matrix(y,ncol=1)}
  
  
  
  result <- list(index.cov=which(s==TRUE),s=s,t=t,w=w,Niter=i,beta.hat.s=beta.hat.s,convergence=convergence)
  class(result) <- "COMBSS"
  return(result)
}




gradient.Combss.Im <- function(w,p,X,XtXn=XtXn,Xtyn=Xtyn,delta,lambda=2,CG,n){
  t <- as.vector(map.w.to.t(w))
  B <- t((t%d*%XtXn))
  #Lt <- t%d*%B

  A <- t%d*%B
  if(p==1) Lt<-A +(1/n)*(delta*(1-t**2)) else Lt<-A +(1/n)*(delta*(diag(1-t**2)))
  
  
  if(CG==TRUE){
    u <- t%d*%Xtyn # u 
    if(p<n){
      beta.tilde <- cgsolve(Lt,u,maxIter = 100) 
    }else{
      rt<-(n/delta)*(1/(1-t**2)) 
      temp0 <- rt*(t**2)
      temp0 <- rt*t
      B <- temp0%d*%t(X)
      
      Lttilde <- (1/n)*X%*%(t%d*%B) ### Ltilde_t-I
      diag(Lttilde ) <-diag(Lttilde )+1  ### B is  Ltilde_t
      
      temp <- X%*%matrix(t*rt*u,ncol=1)
  
      beta.tilde <- rt*u-(1/n)*B%*%(cgsolve(Lttilde ,temp,maxIter = 100))
      
    }
  }else{
    Lt.inv <- solve(Lt)
    beta.tilde <- Lt.inv%*%(t%d*%Xtyn)
  }
  t.beta.tilde <- matrix(t*beta.tilde,ncol=1)
  at <-  XtXn%*%t.beta.tilde-Xtyn
  bt <- at -(delta/n)*(t.beta.tilde)
  
  
  if(CG==TRUE){
    if(p<n){
      ct <- cgsolve(Lt,matrix(t*at,ncol=1),maxIter = 100)}else{
        #time <- proc.time()
        u <- t*at
        rt<-(n/delta)*(1/(1-t**2)) 
        temp0 <- rt*(t**2)
        temp0 <- rt*t
        B <- temp0%d*%t(X)
        Lttilde <- (1/n)*X%*%(t%d*%B) ### Ltilde_t-I
        diag(Lttilde ) <-diag(Lttilde )+1  ### B is  Ltilde_t
        
        temp <- X%*%matrix(t*rt*u,ncol=1)
        ct <- rt*u-(1/n)*B%*%(cgsolve(Lttilde ,temp,maxIter = 100))
       
      }
  }else{
    ct <- Lt.inv%*%(matrix(t*at,ncol=1))
  }
  dt <- (XtXn-delta/n*diag(p))%*%(matrix(t*ct,ncol=1))
  zetat <- 2*(beta.tilde*(at-dt))-2*(bt*ct)
  gradient.g <- (zetat+lambda)*2*w*exp(-w**2)
  return(gradient.g)
}






