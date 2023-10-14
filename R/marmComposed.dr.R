
##--------------Estimation with Penalty by CV----------------------##
marmComposed.dr <- 
  function(Y,X1,X2,G1=NULL,group=NULL,is.fabs=1,K_index=NULL,r_index=NULL,method="BIC",ncv=10,penalty="LASSO",lambda=NULL,
           D0_t3=NULL,D0_t4=NULL,intercept=TRUE,nlam=50,degr=3,lam_min=0.01,eps=1e-4,max_step=20,eps1=1e-4,max_step1_t3=20,max_step1_t4=20,gamma=2,dfmax1=NULL,dfmax2=NULL,alpha=1,vnorm_ratio=1){
    n <- nrow(Y)
    q <- ncol(Y)
    nx1 <- ncol(X1)
    nx2 <- ncol(X2)
    nx = nx1+nx2
    if(is.null(G1)) stop("the parameter G1 must be entered.")
    if(is.null(group)) stop("the parameter group must be entered.")
    gunique <- unique(group)
    G = length(gunique)
    G2 = G-G1
    group1 = group[which(group<=G1)]
    group2 = group[which(group>G1)]-G1
    gunique1 <- unique(group1)
    gunique2 <- unique(group2)
    p = rep(0,G)
    for(g in 1:G) p[g] = sum(group==gunique[g])
    p1 = p[1:G1]
    p2 = p[G1+1] #note!
    K1 <- 6
    if(degr>min(6,K1-1)-1) stop("K must be larger than degree+1 !")
    if(is.null(K_index)) K_index = min(6,K1-1):max(8,K1+1)
    r1_t3_index = r_index$r1_t3_index
    r2_t3_index = r_index$r2_t3_index
    r3_t3_index = r_index$r3_t3_index
    if(is.null(r1_t3_index)) r1_t3_index = 1:min(floor(log(n)),min(p1))
    if(is.null(r2_t3_index)) r2_t3_index = 1:min(K_index)
    if(is.null(r3_t3_index)) r3_t3_index = 1:min(floor(log(n)),q)
    r1_t4_index = r_index$r1_t4_index
    r2_t4_index = r_index$r2_t4_index
    r3_t4_index = r_index$r3_t4_index
    r4_t4_index = r_index$r4_t4_index
    if(is.null(r1_t4_index)) r1_t4_index = 1:min(floor(log(n)),p2)
    if(is.null(r2_t4_index)) r2_t4_index = 1:min(K_index)
    if(is.null(r3_t4_index)) r3_t4_index = 1:min(floor(log(n)),G2)
    if(is.null(r4_t4_index)) r4_t4_index = 1:min(floor(log(n)),q)
    if (penalty == "LASSO") pen <- 1
    if (penalty == "MCP")   pen <- 2 
    if (penalty=="SCAD"){    
      gamma <- 3
      pen <- 3
    }  
    if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MCP penalty")
    if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
    
    if (is.null(dfmax1)) dfmax1 = nx1 + 1
    if (is.null(dfmax2)) dfmax2 = p2 + 1
    
    # initial A,B,C,D,S
    if(is.null(D0_t3)){
      set.seed(1)
      r1_t3_max = max(r1_t3_index) 
      r2_t3_max = max(r2_t3_index) 
      r3_t3_max = max(r3_t3_index) 
      K_max = max(K_index)
      B = rbind(diag(r2_t3_max), matrix(0,K_max-r2_t3_max,r2_t3_max))
      C = rbind(diag(r3_t3_max), matrix(0,q-r3_t3_max,r3_t3_max))
      S = matrix(rnorm(r1_t3_max*r2_t3_max*r3_t3_max),r3_t3_max,r1_t3_max*r2_t3_max)
      D0_t3 = list(); 
      for(j in 1:G1){
        A = rbind(diag(r1_t3_max), matrix(0,p1[j]-r1_t3_max,r1_t3_max))
        SABC = list(S=S,A=A,B=B,C=C)
        D0_t3[[j]] = SABC
      }
    }
    if(is.null(D0_t4)){
      set.seed(1)
      r1_t4_max = max(r1_t4_index) 
      r2_t4_max = max(r2_t4_index) 
      r3_t4_max = max(r3_t4_index) 
      r4_t4_max = max(r4_t4_index) 
      K_max = max(K_index)
      A = rbind(diag(r1_t4_max), matrix(0,p2-r1_t4_max,r1_t4_max))
      B = rbind(diag(r2_t4_max), matrix(0,K_max-r2_t4_max,r2_t4_max))
      C = rbind(diag(r3_t4_max), matrix(0,G2-r3_t4_max,r3_t4_max))
      D = rbind(diag(r4_t4_max), matrix(0,q-r4_t4_max,r4_t4_max))
      S = matrix(rnorm(r1_t4_max*r2_t4_max*r3_t4_max*r4_t4_max),r4_t4_max,r1_t4_max*r2_t4_max*r3_t4_max)
      D0_t4 = list(S=S,A=A,B=B,C=C,D=D)
    }
    else{
      A = D0_t4$A #note!
      B = D0_t4$B
      C = D0_t4$C
      D = D0_t4$D
      S = D0_t4$S
    }
    
    opts = list(eps=eps,eps1=eps1,max_step=max_step,max_step1=-1,max_step1_t3=max_step1_t3,max_step1_t4=max_step1_t4,n=n,r1=2,r2=2,r3=2,r4=2,p1=p1,p2=p2,q=q,G=G,G1=G1,G2=G2,degr=degr,K=max(K_index),nx=-1,nx1=nx1,nx2=nx2,is.fabs=is.fabs)
    
    if (is.null(lambda)) {
      is_setlam = 1
      if (nlam < 1||is.null(nlam)) stop("nlambda must be at least 1")
      setlam = c(1,lam_min,alpha,nlam)
      Sinit = list()
      Ainit = list()
      Binit = list()
      Cinit = list()
      Z_t3 = list()
      for(i in 1:G1){
        Sinit[[i]] = D0_t3[[i]]$S
        Ainit[[i]] = D0_t3[[i]]$A
        Binit[[i]] = D0_t3[[i]]$B
        Cinit[[i]] = D0_t3[[i]]$C
        Z_t3[[i]] = bsbasefun(X1[,group1==gunique1[i]],max(K_index),degr) #note!
        Zbar = colMeans(Z_t3[[i]])
        Z_t3[[i]] = Z_t3[[i]] - matrix(rep(Zbar,each=n),n)
      }
      Z_t4 = NULL;    for(i in 1:G2)  Z_t4 = cbind(Z_t4, bsbasefun(X2[,group2==gunique2[i]],max(K_index),degr))
      Zbar = colMeans(Z_t4)
      Z_t4 = Z_t4 - matrix(rep(Zbar,each=n),n)
      
      Ybar = colMeans(Y)
      Y1 = Y - matrix(rep(Ybar,each=n),n)
      fit_lambda = setuplambdaV(Y1,Z_t3,Sinit,Ainit,Binit,Cinit,Z_t4,S,A,B,C,D,nx1,G1,nlam,setlam)
      lambda = fit_lambda$lambda
      vnorm_ratio0 = fit_lambda$vnorm_ratio
    }
    else {
      is_setlam = 0
      nlam = length(lambda)
      setlam = c(1,lam_min,alpha,nlam)
    }
    
    opts_pen = list(gamma=gamma,dfmax=-1,dfmax1=dfmax1,dfmax2=dfmax2,pen=pen,nlam=nlam,lam_max=1,lam_min=lam_min,alpha=alpha,isPenColumn=1,vnorm_ratio=vnorm_ratio,vnorm_ratio0=vnorm_ratio0) 
    
    if((max(r1_t4_index)>dim(A)[2])|(max(r2_t4_index)>dim(B)[2])|(max(r3_t4_index)>dim(C)[2])|(max(r4_t4_index)>dim(D)[2]))
      stop("maximum number of index sequence of r1, r2, r3 and r4 must not be larger than A, B, C and D, respectively !")
    #---------------- The selection by CV  ---------------------#  
    if(method=="BIC") fit_dr = marmComposed.bic(Y,X1,X2,G1,group,K_index,r_index,lambda,D0_t3,D0_t4,intercept,opts,opts_pen)
    #if(method=="CV") fit_dr = gmam_sparse_composed_cv(Y,X,group,ncv,K_index,r1_index,r2_index,r3_index,lambda,D0,intercept,opts,opts_pen)
    
    return(fit_dr)
  }