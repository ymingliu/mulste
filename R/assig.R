assig <- function(n_args){
  cargs <- vector("list", length(n_args))
  for(i in 1:length(n_args)) cargs[[i]] <- 1:n_args[i]
  t(expand.grid(cargs))
}

##------------------------------------------------------------------------##
##--------------------produce the B-spline functions----------------------##
bsbasefun <- function(X,K,degr){
  n = dim(X)[1]
  p = dim(X)[2]
  nk = K - degr
  u.k = seq(0, 1, length=nk+2)[-c(1,nk+2)]
  BS = NULL
  for(j in 1:p){
    Knots = as.numeric(quantile(X[,j], u.k))  
    BS0 = bs(X[,j], knots=Knots, intercept=TRUE, degree=degr)
    BS = cbind(BS,BS0[,-1])
  }
  BS = scale(BS,center = T, scale = F)
  id = seq(1,p*K,K)
  Z = NULL
  for(j in 1:K){
    Z = cbind(Z,BS[,id+j-1])
  }
  return(Z)
}

##------------------------------------------------------------------------##
##-----------------generate scenario I data in MARM----------------------##
# e.g. mydata=marm3.sim.fbs(200,5,100,10,rep(1:4,each=25))
marm3.sim.fbs <- function(n,q,p,s,group=NULL,r10=2,r20=2,r30=2,isfixedR=0,D3=NULL,
                           K=6,degr=3,sigma2=NULL,seed_id=NULL,
                           r1_index=NULL,r2_index=NULL,r3_index=NULL,D0=NULL){
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s<1) stop("s must be not smaller than 1")
  if(is.null(group)) stop("the parameter group must be entered")
  gunique <- unique(group)
  ng = length(gunique)
  pg = p/ng
  if(is.null(sigma2)) sigma2 = 0.1
  if(is.null(seed_id)) seed_id=1000
  
  if(is.null(D3)) {
    set.seed(2)
    S3 <- matrix(runif(r10*r20*r30,5,10),nrow = r30)
    T1 <- matrix(rnorm(s*r10),nrow = s)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r20),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(q*r30),nrow = q)
    U3 <- qr.Q(qr(T1))
    D3 <- U3%*%S3%*%t(kronecker(U2,U1))
  }

  set.seed(seed_id)
  X <- matrix(runif(n*p), nrow = n)
  
  f0 = matrix(0,n,q)
  for(j in 1:ng){
    Xj <- X[,group==gunique[j]]
    X1 <- Xj[,1:s]
    f0 = f0 + bsbasefun(X1,K,degr)%*%t(D3)
  }
  eps <- matrix(rnorm(n*q),n,q)
  Y <- f0 + eps*sigma2
  
  # initialize
  if(is.null(r1_index)) r1_index = 1:min(4,s)
  if(is.null(r2_index)) r2_index = 1:min(4,K)
  if(is.null(r3_index)) r3_index = 1:min(4,q)
  if(is.null(D0)) {
    set.seed(10)
    if(isfixedR){
      T1 = matrix(rnorm(pg*r10),pg,r10)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r20),K,r20)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r30),q,r30)
      C <- qr.Q(qr(T1))
      S = matrix(runif(r10*r20*r30),r30,r10*r20)
      SABC = list(S=S,A=A,B=B,C=C)
      D0 = list()
      for(j in 1:ng) D0[[j]] = SABC
    }else{
      r1_max = max(r1_index) 
      r2_max = max(r2_index) 
      r3_max = max(r3_index) 
      K_max = K
      T1 = matrix(rnorm(pg*r1_max),pg,r1_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_max),K_max,r2_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r3_max),q,r3_max)
      C <- qr.Q(qr(T1))
      S = matrix(runif(r1_max*r2_max*r3_max),r3_max,r1_max*r2_max)
      SABC = list(S=S,A=A,B=B,C=C)
      D0 = list()
      for(j in 1:ng) D0[[j]] = SABC
    }
  }
  
  return(list(Y=Y,X=X,f0=f0,group=group,D3=D3,n=n,q=q,p=p,s=s,r10=r10,r20=r20,r30=r30,K=K,degr=degr,sigma2=sigma2,seed_id=seed_id,
              r1_index=r1_index,r2_index=r2_index,r3_index=r3_index,D0=D0))
}

##------------------------------------------------------------------------##
##----------------generate scenario II data in MARM----------------------##
# e.g. mydata=marm3.sim.fsin(200,5,100,10,rep(1:4,each=25))
marm3.sim.fsin <- function(n,q,p,s,group=NULL,r10=2,r20=2,r30=2,isfixedR=0,D2=NULL,
                            K=6,degr=3,sigma2=NULL,seed_id=NULL,
                            r1_index=NULL,r2_index=NULL,r3_index=NULL,D0=NULL){
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s<1) stop("s must be not smaller than 1")
  if(is.null(group)) stop("the parameter group must be entered")
  gunique <- unique(group)
  ng = length(gunique)
  pg = p/ng
  if(is.null(sigma2)) sigma2 = 0.1
  if(is.null(seed_id)) seed_id=1000
  
  if(is.null(D2)) {
    set.seed(2)
    S3 <- matrix(runif(r10*r20*r30,5,10),nrow = r30)
    T1 <- matrix(rnorm(s*r10),nrow = s)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r20),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(q*r30),nrow = q)
    U3 <- qr.Q(qr(T1))
    D3 <- U3%*%S3%*%t(kronecker(U2,U1))
    D2 <- TransferModalUnfoldingsT(D3,3,2,c(s,K,q))
  }
  
  set.seed(seed_id)
  X <- matrix(runif(n*p), nrow = n)
  
  f0 = matrix(0,n,q)
  for(j in 1:ng){
    Xj <- X[,group==gunique[j]]
    X1 <- Xj[,1:s]
    basefuns1 <- sin(2*pi*X1)
    basefuns2 <- cos(pi*X1)
    f0 <- f0 + basefuns1%*%matrix(D2[1,],nrow=s)+basefuns2%*%matrix(D2[2,],nrow=s)
  }
  eps <- matrix(rnorm(n*q),n,q)
  Y <- f0 + sigma2*eps
  
  # initialize
  if(is.null(r1_index)) r1_index = 1:min(4,s)
  if(is.null(r2_index)) r2_index = 1:min(4,K)
  if(is.null(r3_index)) r3_index = 1:min(4,q)
  if(is.null(D0)) {
    set.seed(10)
    if(isfixedR){
      T1 = matrix(rnorm(pg*r10),pg,r10)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r20),K,r20)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r30),q,r30)
      C <- qr.Q(qr(T1))
      S = matrix(runif(r10*r20*r30),r30,r10*r20)
      SABC = list(S=S,A=A,B=B,C=C)
      D0 = list()
      for(j in 1:ng) D0[[j]] = SABC
    }else{
      r1_max = max(r1_index) 
      r2_max = max(r2_index) 
      r3_max = max(r3_index) 
      K_max = K
      T1 = matrix(rnorm(pg*r1_max),pg,r1_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_max),K_max,r2_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r3_max),q,r3_max)
      C <- qr.Q(qr(T1))
      S = matrix(runif(r1_max*r2_max*r3_max),r3_max,r1_max*r2_max)
      SABC = list(S=S,A=A,B=B,C=C)
      D0 = list()
      for(j in 1:ng) D0[[j]] = SABC
    }
  }
  
  return(list(Y=Y,X=X,f0=f0,group=group,D3=D3,D2=D2,n=n,q=q,p=p,s=s,r10=r10,r20=r20,r30=r30,K=K,degr=degr,sigma2=sigma2,seed_id=seed_id,
              r1_index=r1_index,r2_index=r2_index,r3_index=r3_index,D0=D0))
}

##------------------------------------------------------------------------##
##--------------generate scenario I data in structural MARM--------------##
# e.g. mydata=marm4.sim.fbs(200,5,100,10,rep(1:4,each=25))
marm4.sim.fbs <- function(n,q,p,s,group=NULL,r10=2,r20=2,r30=2,r40=2,isfixedR=0,D44=NULL,
                           K=6,degr=3,sigma2=NULL,seed_id=NULL,
                           r1_index=NULL,r2_index=NULL,r3_index=NULL,r4_index=NULL,D0=NULL){
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s<1) stop("s must be not smaller than 1")
  if(is.null(group)) stop("the parameter group must be entered")
  gunique <- unique(group)
  ng = length(gunique)
  pg = p/ng
  if(ng<r30) stop("ng must be not smaller than r30")
  if(is.null(sigma2)) sigma2 = 0.1
  if(is.null(seed_id)) seed_id=1000
  
  if(is.null(D44)) {
    set.seed(2)
    T1 <- matrix(rnorm(s*r10),nrow = s)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r20),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(ng*r30),nrow = ng)
    U3 <- as.matrix(qr.Q(qr(T1)))
    T1 <- matrix(rnorm(q*r40),nrow = q)
    U4 <- qr.Q(qr(T1)) 
    S4 <- matrix(runif(r10*r20*r30*r40,5,10),nrow = r30)
    D44 <- U4%*%S4%*%t(kronecker(U3,kronecker(U2,U1)))
  }

  set.seed(seed_id)
  X <- matrix(runif(n*p), nrow = n)
  
  f0 = matrix(0,n,q)
  for(j in 1:ng){
    D3 = D44[,((j-1)*s*K+1):(j*s*K)]
    Xj <- X[,group==gunique[j]]
    X1 <- Xj[,1:s]
    f0 = f0 + bsbasefun(X1,K,degr)%*%t(D3)
  }
  eps <- matrix(rnorm(n*q),n,q)
  Y <- f0 + eps*sigma2
  
  # initialize
  if(is.null(r1_index)) r1_index = 1:min(4,s)
  if(is.null(r2_index)) r2_index = 1:min(4,K)
  if(is.null(r3_index)) r3_index = 1:min(4,ng)
  if(is.null(r4_index)) r4_index = 1:min(4,q)
  if(is.null(D0)) {
    set.seed(10)
    if(isfixedR){
      T1 = matrix(rnorm(pg*r10),pg,r10)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r20),K,r20)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(ng*r30),ng,r30)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r40),q,r40)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r10*r20*r30*r40),r40,r10*r20*r30)
      D0 = list(S=S,A=A,B=B,C=C,D=D)
    }else{
      r1_max = max(r1_index) 
      r2_max = max(r2_index) 
      r3_max = max(r3_index) 
      r4_max = max(r4_index) 
      K_max = K
      T1 = matrix(rnorm(pg*r1_max),pg,r1_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_max),K_max,r2_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(ng*r3_max),ng,r3_max)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r4_max),q,r4_max)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r1_max*r2_max*r3_max*r4_max),r4_max,r1_max*r2_max*r3_max)
      D0 = list(S=S,A=A,B=B,C=C,D=D)
    }
  }
  
  return(list(Y=Y,X=X,f0=f0,group=group,D44=D44,n=n,q=q,p=p,s=s,r10=r10,r20=r20,r30=r30,r40=r40,K=K,degr=degr,sigma2=sigma2,seed_id=seed_id,
              r1_index=r1_index,r2_index=r2_index,r3_index=r3_index,r4_index=r4_index,D0=D0))
}

##------------------------------------------------------------------------##
##-------------generate scenario II data in structural MARM--------------##
# e.g. mydata=marm4.sim.fsin(200,5,100,10,rep(1:4,each=25))
marm4.sim.fsin <- function(n,q,p,s,group=NULL,r10=2,r20=2,r30=2,r40=2,isfixedR=0,D42=NULL,
                            K=6,degr=3,sigma2=NULL,seed_id=NULL,
                            r1_index=NULL,r2_index=NULL,r3_index=NULL,r4_index=NULL,D0=NULL){
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s<1) stop("s must be not smaller than 1")
  if(is.null(group)) stop("the parameter group must be entered")
  gunique <- unique(group)
  ng = length(gunique)
  pg = p/ng
  if(ng<r30) stop("ng must be not smaller than r30")
  if(is.null(sigma2)) sigma2 = 0.1
  if(is.null(seed_id)) seed_id=1000
  
  if(is.null(D42)) {
    set.seed(2)
    T1 <- matrix(rnorm(s*r10),nrow = s)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r20),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(ng*r30),nrow = ng)
    U3 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(q*r40),nrow = q)
    U4 <- qr.Q(qr(T1)) 
    S4 <- matrix(runif(r10*r20*r30*r40,5,10),nrow = r30)
    D44 <- U4%*%S4%*%t(kronecker(U3,kronecker(U2,U1)))
    D42 = TransferModalUnfoldingsT(D44,4,2,c(s,K,ng,q))
  }
  
  set.seed(seed_id)
  X <- matrix(runif(n*p), nrow = n)
  
  f0 = matrix(0,n,q)
  id = NULL
  for(j in 1:q) id = c(id, c(1:s)+(j-1)*s*ng)
  for(j in 1:ng){
    D2 = D42[,id+(j-1)*s]
    Xj <- X[,group==gunique[j]]
    X1 <- Xj[,1:s]
    basefuns1 <- sin(2*pi*X1)
    basefuns2 <- cos(pi*X1)
    f0 <- f0 + basefuns1%*%matrix(D2[1,],nrow=s)+basefuns2%*%matrix(D2[2,],nrow=s)
  }
  eps <- matrix(rnorm(n*q),n,q)
  Y <- f0 + sigma2*eps
  
  # initialize
  if(is.null(r1_index)) r1_index = 1:min(4,s)
  if(is.null(r2_index)) r2_index = 1:min(4,K)
  if(is.null(r3_index)) r3_index = 1:min(4,ng)
  if(is.null(r4_index)) r4_index = 1:min(4,q)
  if(is.null(D0)) {
    set.seed(10)
    if(isfixedR){
      T1 = matrix(rnorm(pg*r10),pg,r10)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r20),K,r20)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(ng*r30),ng,r30)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r40),q,r40)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r10*r20*r30*r40),r40,r10*r20*r30)
      D0 = list(S=S,A=A,B=B,C=C,D=D)
    }else{
      r1_max = max(r1_index) 
      r2_max = max(r2_index) 
      r3_max = max(r3_index) 
      r4_max = max(r4_index) 
      K_max = K
      T1 = matrix(rnorm(pg*r1_max),pg,r1_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_max),K_max,r2_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(ng*r3_max),ng,r3_max)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r4_max),q,r4_max)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r1_max*r2_max*r3_max*r4_max),r4_max,r1_max*r2_max*r3_max)
      D0 = list(S=S,A=A,B=B,C=C,D=D)
    }
  }
  
  return(list(Y=Y,X=X,f0=f0,group=group,D44=D44,D42=D42,n=n,q=q,p=p,s=s,r10=r10,r20=r20,r30=r30,r40=r40,K=K,degr=degr,sigma2=sigma2,seed_id=seed_id,
              r1_index=r1_index,r2_index=r2_index,r3_index=r3_index,r4_index=r4_index,D0=D0))
}



##------------------------------------------------------------------------##
##--------------generate scenario I data in composed model----------------##
# e.g. mydata=marmComposed.sim.fbs(200,5,100,15,10,1,rep(1:4,each=25))
marmComposed.sim.fbs <- function(n,q,p,s1,s2,G1=NULL,group=NULL,r310=2,r320=2,r330=2,r410=2,r420=2,r430=2,r440=2,isfixedR=0,D3=NULL,D44=NULL,
                                 K=6,degr=3,sigma2=NULL,seed_id=NULL,
                                 r1_t3_index=NULL,r2_t3_index=NULL,r3_t3_index=NULL,r1_t4_index=NULL,r2_t4_index=NULL,r3_t4_index=NULL,r4_t4_index=NULL,D0_t3=NULL,D0_t4=NULL){
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s1<1) stop("s1 must be not smaller than 1")
  if(s2<1) stop("s2 must be not smaller than 1")
  if(is.null(G1)) stop("the parameter G1 must be entered")
  if(is.null(group)) stop("the parameter group must be entered")
  gunique <- unique(group)
  G = length(gunique)
  G2 = G-G1
  pg1 = sum(group==gunique[G1])
  pg2 = sum(group==gunique[G1+1])
  p1 = pg1*G1
  p2 = pg2*G2
  group1 = group[which(group<=G1)]
  group2 = group[which(group>G1)]-G1
  gunique1 <- unique(group1)
  gunique2 <- unique(group2)
  if(is.null(sigma2)) sigma2 = 0.1
  if(is.null(seed_id)) seed_id=1000
  
  if(is.null(D3) | is.null(D44)) {
    set.seed(2)
    S3 <- matrix(runif(r310*r320*r330,5,10),nrow = r330)
    T1 <- matrix(rnorm(s1*r310),nrow = s1)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r320),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(q*r330),nrow = q)
    U3 <- qr.Q(qr(T1))
    D3 <- U3%*%S3%*%t(kronecker(U2,U1))

    T1 <- matrix(rnorm(s2*r410),nrow = s2)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r420),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(G2*r430),nrow = G2)
    U3 <- as.matrix(qr.Q(qr(T1)))
    T1 <- matrix(rnorm(q*r440),nrow = q)
    U4 <- qr.Q(qr(T1)) 
    S4 <- matrix(runif(r410*r420*r430*r440,5,10),nrow = r430)
    D44 <- U4%*%S4%*%t(kronecker(U3,kronecker(U2,U1)))
  }
  
  
  set.seed(seed_id)
  X1 <- matrix(runif(n*p1), nrow = n)
  X2 <- matrix(runif(n*p2), nrow = n)

  
  f01 = matrix(0,n,q)
  for(j in 1:G1){
    X1j <- X1[,group1==gunique1[j]]
    X11 <- X1j[,1:s1]
    f01 = f01 + bsbasefun(X11,K,degr)%*%t(D3)
  }
  len = s2
  f02 = matrix(0,n,q)
  for(j in 1:G2){
    D23 = D44[,((j-1)*len*K+1):(j*len*K)]
    X2j <- X2[,group2==gunique2[j]]
    X21 <- X2j[,1:s2]
    f02 = f02 + bsbasefun(X21,K,degr)%*%t(D23)
  }
  eps <- matrix(rnorm(n*q),n,q)
  Y <- f01 + f02 + eps*sigma2
  
  
  # initialize
  K_index = K
  if(is.null(r1_t3_index)) r1_t3_index = 1:min(4,s1)
  if(is.null(r2_t3_index)) r2_t3_index = 1:min(4,K)
  if(is.null(r3_t3_index)) r3_t3_index = 1:min(4,q)
  if(is.null(r1_t4_index)) r1_t4_index = 1:min(4,s2)
  if(is.null(r2_t4_index)) r2_t4_index = 1:min(4,K)
  if(is.null(r3_t4_index)) r3_t4_index = 1:min(4,G2) 
  if(is.null(r4_t4_index)) r4_t4_index = 1:min(4,q)
  r_index = list(r1_t3_index = r1_t3_index, r2_t3_index = r2_t3_index, r3_t3_index = r3_t3_index, 
                 r1_t4_index = r1_t4_index, r2_t4_index = r2_t4_index, r3_t4_index = r3_t4_index, r4_t4_index = r4_t4_index)
  if(is.null(D0_t3) | is.null(D0_t4)) {
    set.seed(10)
    if(isfixedR){
      T1 = matrix(rnorm(pg1*r310),pg1,r310)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r320),K,r320)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r330),q,r330) 
      C <- qr.Q(qr(T1))
      S = matrix(runif(r310*r320*r330),r330,r310*r320)
      SABC = list(S=S,A=A,B=B,C=C)
      D0_t3 = list()
      for(j in 1:G1) D0_t3[[j]] = SABC
      
      T1 = matrix(rnorm(pg2*r410),pg2,r410)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r420),K,r420)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(G2*r430),G2,r430)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r440),q,r440)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r410*r420*r430*r440),r440,r410*r420*r430)
      D0_t4 = list(S=S,A=A,B=B,C=C,D=D)
    }else{
      r1_t3_max = max(r1_t3_index) 
      r2_t3_max = max(r2_t3_index) 
      r3_t3_max = max(r3_t3_index) 
      K_max = max(K_index)
      T1 = matrix(rnorm(pg1*r1_t3_max),pg1,r1_t3_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_t3_max),K_max,r2_t3_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r3_t3_max),q,r3_t3_max) #note! r3_t3_max = r4_t4_max
      C <- qr.Q(qr(T1))
      S = matrix(runif(r1_t3_max*r2_t3_max*r3_t3_max),r3_t3_max,r1_t3_max*r2_t3_max)
      SABC = list(S=S,A=A,B=B,C=C)
      D0_t3 = list()
      for(j in 1:G1) D0_t3[[j]] = SABC
      
      r1_t4_max = max(r1_t4_index) 
      r2_t4_max = max(r2_t4_index) 
      r3_t4_max = max(r3_t4_index) 
      r4_t4_max = max(r4_t4_index) 
      K_max = max(K_index)
      T1 = matrix(rnorm(pg2*r1_t4_max),pg2,r1_t4_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_t4_max),K_max,r2_t4_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(G2*r3_t4_max),G2,r3_t4_max)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r4_t4_max),q,r4_t4_max)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r1_t4_max*r2_t4_max*r3_t4_max*r4_t4_max),r4_t4_max,r1_t4_max*r2_t4_max*r3_t4_max)
      D0_t4 = list(S=S,A=A,B=B,C=C,D=D)
    }
  }
  
  return(list(Y=Y,X1=X1,X2=X2,f01=f01,f02=f02,group=group,D3=D3,D44=D44,n=n,q=q,p=p,s1=s1,s2=s2,r310=r310,r320=r320,r330=r330,r410=r410,r420=r420,r430=r430,r440=r440,
              K=K,degr=degr,sigma2=sigma2,seed_id=seed_id,is.fabs=1,
              r1_t3_index = r1_t3_index, r2_t3_index = r2_t3_index, r3_t3_index = r3_t3_index, 
              r1_t4_index = r1_t4_index, r2_t4_index = r2_t4_index, r3_t4_index = r3_t4_index, r4_t4_index = r4_t4_index,
              r_index=r_index,D0_t3=D0_t3,D0_t4=D0_t4))
}




##------------------------------------------------------------------------##
##--------------generate scenario I data in composed model----------------##
# e.g. mydata=marmComposed.sim.fsin(200,5,100,15,10,1,rep(1:4,each=25))
marmComposed.sim.fsin <- function(n,q,p,s1,s2,G1=NULL,group=NULL,r310=2,r320=2,r330=2,r410=2,r420=2,r430=2,r440=2,isfixedR=0,D2=NULL,D42=NULL,
                                  K=6,degr=3,sigma2=NULL,seed_id=NULL,
                                  r1_t3_index=NULL,r2_t3_index=NULL,r3_t3_index=NULL,r1_t4_index=NULL,r2_t4_index=NULL,r3_t4_index=NULL,r4_t4_index=NULL,D0_t3=NULL,D0_t4=NULL){
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s1<1) stop("s1 must be not smaller than 1")
  if(s2<1) stop("s2 must be not smaller than 1")
  if(is.null(G1)) stop("the parameter G1 must be entered")
  if(is.null(group)) stop("the parameter group must be entered")
  gunique <- unique(group)
  G = length(gunique)
  G2 = G-G1
  pg1 = sum(group==gunique[G1])
  pg2 = sum(group==gunique[G1+1])
  p1 = pg1*G1
  p2 = pg2*G2
  group1 = group[which(group<=G1)]
  group2 = group[which(group>G1)]-G1
  gunique1 <- unique(group1)
  gunique2 <- unique(group2)
  if(is.null(sigma2)) sigma2 = 0.1
  if(is.null(seed_id)) seed_id=1000
  
  if(is.null(D2) | is.null(D42)) {
    set.seed(2)
    S3 <- matrix(runif(r310*r320*r330,5,10),nrow = r330)
    T1 <- matrix(rnorm(s1*r310),nrow = s1)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r320),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(q*r330),nrow = q)
    U3 <- qr.Q(qr(T1))
    D3 <- U3%*%S3%*%t(kronecker(U2,U1))
    D2 <- TransferModalUnfoldingsT(D3,3,2,c(s1,K,q))
    
    T1 <- matrix(rnorm(s2*r410),nrow = s2)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r420),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(G2*r430),nrow = G2)
    U3 <- as.matrix(qr.Q(qr(T1)))
    T1 <- matrix(rnorm(q*r440),nrow = q)
    U4 <- qr.Q(qr(T1)) 
    S4 <- matrix(runif(r410*r420*r430*r440,5,10),nrow = r430)
    D44 <- U4%*%S4%*%t(kronecker(U3,kronecker(U2,U1)))
    D42 <- TransferModalUnfoldingsT(D44,4,2,c(s2,K,G2,q))
  }
  
  set.seed(seed_id)
  X1 <- matrix(runif(n*p1), nrow = n)
  X2 <- matrix(runif(n*p2), nrow = n)
  
  f01 = matrix(0,n,q)
  for(j in 1:G1){
    X1j <- X1[,group1==gunique1[j]]
    X11 <- X1j[,1:s1]
    basefuns11 <- sin(2*pi*X11)
    basefuns12 <- cos(pi*X11)
    f01 <- f01 + basefuns11%*%matrix(D2[1,],nrow=s1)+basefuns12%*%matrix(D2[2,],nrow=s1)
  }
  len = s2
  f02 = matrix(0,n,q)
  id = NULL
  for(j in 1:q) id = c(id, c(1:len)+(j-1)*len*G2)
  for(j in 1:G2){
    D22 = D42[,id+(j-1)*len]
    X2j <- X2[,group2==gunique2[j]]
    X21 <- X2j[,1:s2]
    basefuns21 <- sin(2*pi*X21)
    basefuns22 <- cos(pi*X21)
    f02 <- f02 + basefuns21%*%matrix(D22[1,],nrow=s2)+basefuns22%*%matrix(D22[2,],nrow=s2)
  }
  eps <- matrix(rnorm(n*q),n,q)
  Y = f01 + f02 + eps*sigma2
  
  # initialize
  K_index = K
  if(is.null(r1_t3_index)) r1_t3_index = 1:min(4,s1)
  if(is.null(r2_t3_index)) r2_t3_index = 1:min(4,K)
  if(is.null(r3_t3_index)) r3_t3_index = 1:min(4,q)
  if(is.null(r1_t4_index)) r1_t4_index = 1:min(4,s2)
  if(is.null(r2_t4_index)) r2_t4_index = 1:min(4,K)
  if(is.null(r3_t4_index)) r3_t4_index = 1:min(4,G2) 
  if(is.null(r4_t4_index)) r4_t4_index = 1:min(4,q)
  r_index = list(r1_t3_index = r1_t3_index, r2_t3_index = r2_t3_index, r3_t3_index = r3_t3_index, 
                 r1_t4_index = r1_t4_index, r2_t4_index = r2_t4_index, r3_t4_index = r3_t4_index, r4_t4_index = r4_t4_index)
  if(is.null(D0_t3) | is.null(D0_t4)) {
    set.seed(10)
    if(isfixedR){
      T1 = matrix(rnorm(pg1*r310),pg1,r310)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r320),K,r320)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r330),q,r330) 
      C <- qr.Q(qr(T1))
      S = matrix(runif(r310*r320*r330),r330,r310*r320)
      SABC = list(S=S,A=A,B=B,C=C)
      D0_t3 = list()
      for(j in 1:G1) D0_t3[[j]] = SABC
      
      T1 = matrix(rnorm(pg2*r410),pg2,r410)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r420),K,r420)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(G2*r430),G2,r430)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r440),q,r440)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r410*r420*r430*r440),r440,r410*r420*r430)
      D0_t4 = list(S=S,A=A,B=B,C=C,D=D)
    }else{
      r1_t3_max = max(r1_t3_index) 
      r2_t3_max = max(r2_t3_index) 
      r3_t3_max = max(r3_t3_index) 
      K_max = max(K_index)
      T1 = matrix(rnorm(pg1*r1_t3_max),pg1,r1_t3_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_t3_max),K_max,r2_t3_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r3_t3_max),q,r3_t3_max) #note! r3_t3_max = r4_t4_max
      C <- qr.Q(qr(T1))
      S = matrix(runif(r1_t3_max*r2_t3_max*r3_t3_max),r3_t3_max,r1_t3_max*r2_t3_max)
      SABC = list(S=S,A=A,B=B,C=C)
      D0_t3 = list()
      for(j in 1:G1) D0_t3[[j]] = SABC
      
      r1_t4_max = max(r1_t4_index) 
      r2_t4_max = max(r2_t4_index) 
      r3_t4_max = max(r3_t4_index) 
      r4_t4_max = max(r4_t4_index) 
      K_max = max(K_index)
      T1 = matrix(rnorm(pg2*r1_t4_max),pg2,r1_t4_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_t4_max),K_max,r2_t4_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(G2*r3_t4_max),G2,r3_t4_max)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r4_t4_max),q,r4_t4_max)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r1_t4_max*r2_t4_max*r3_t4_max*r4_t4_max),r4_t4_max,r1_t4_max*r2_t4_max*r3_t4_max)
      D0_t4 = list(S=S,A=A,B=B,C=C,D=D)
    }
  }
  
  return(list(Y=Y,X1=X1,X2=X2,f01=f01,f02=f02,group=group,D2=D2,D42=D42,n=n,q=q,p=p,s1=s1,s2=s2,r310=r310,r320=r320,r330=r330,r410=r410,r420=r420,r430=r430,r440=r440,
              K=K,degr=degr,sigma2=sigma2,seed_id=seed_id,is.fabs=0,
              r1_t3_index = r1_t3_index, r2_t3_index = r2_t3_index, r3_t3_index = r3_t3_index, 
              r1_t4_index = r1_t4_index, r2_t4_index = r2_t4_index, r3_t4_index = r3_t4_index, r4_t4_index = r4_t4_index,
              r_index=r_index,D0_t3=D0_t3,D0_t4=D0_t4))
}


















