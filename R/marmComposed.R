marmComposed <- 
  function(Y,X1,X2,G1=NULL,group=NULL,is.fabs=1,K=6,r31=NULL,r32=NULL,r33=NULL,r41=NULL,r42=NULL,r43=NULL,r44=NULL,method="BIC",ncv=10,penalty="LASSO",lambda=NULL,
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
    pv1 = p[1:G1]
    p1 = p[G1]
    p2 = p[G1+1] #note!
    if(degr>K-1) stop("K must be larger than degree+1 !")
    if(is.null(r31)) r31 <- 2
    if(is.null(r32)) r32 <- 2
    if(is.null(r33)) r33 <- 2
    if(is.null(r41)) r41 <- 2
    if(is.null(r42)) r42 <- 2
    if(is.null(r43)) r43 <- 2
    if(is.null(r44)) r44 <- 2
    r_index = c(r31=r31,r32=r32,r33=r33,r41=r41,r42=r42,r43=r43,r44=r44)
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
      B = rbind(diag(r32), matrix(0,K-r32,r32))
      C = rbind(diag(r33), matrix(0,q-r33,r33))
      S = matrix(rnorm(r31*r32*r33),r33,r31*r32)
      D0_t3 = list(); 
      for(j in 1:G1){
        A = rbind(diag(r31), matrix(0,pv1[j]-r31,r31))
        SABC = list(S=S,A=A,B=B,C=C)
        D0_t3[[j]] = SABC
      }
    }
    if(is.null(D0_t4)){
      set.seed(1)
      A = rbind(diag(r41), matrix(0,p2-r41,r41))
      B = rbind(diag(r42), matrix(0,K-r42,r42))
      C = rbind(diag(r43), matrix(0,G2-r43,r43))
      D = rbind(diag(r44), matrix(0,q-r44,r44))
      S = matrix(rnorm(r41*r42*r43*r44),r44,r41*r42*r43)
      D0_t4 = list(S=S,A=A,B=B,C=C,D=D)
    }
    else{
      A = D0_t4$A #note!
      B = D0_t4$B
      C = D0_t4$C
      D = D0_t4$D
      S = D0_t4$S
    }
    #note!
    opts = list(eps=eps,eps1=eps1,max_step=max_step,max_step1=-1,max_step1_t3=max_step1_t3,max_step1_t4=max_step1_t4,n=n,r1=r41,r2=r42,r3=r43,r4=r44,pv1=pv1,p1=p1,p2=p2,q=q,G=G,G1=G1,G2=G2,degr=degr,K=K,nx=-1,nx1=nx1,nx2=nx2,is.fabs=is.fabs)
    Sinit = list()
    Ainit = list()
    Binit = list()
    Cinit = list()
    Z_t3 = list()
    Zbar_t3 = NULL
    for(i in 1:G1){
      Sinit[[i]] = D0_t3[[i]]$S
      Ainit[[i]] = D0_t3[[i]]$A
      Binit[[i]] = D0_t3[[i]]$B
      Cinit[[i]] = D0_t3[[i]]$C
      Z_t3[[i]] = bsbasefun(X1[,group1==gunique1[i]],K,degr) #note!
      Zbar = colMeans(Z_t3[[i]])
      Z_t3[[i]] = Z_t3[[i]] - matrix(rep(Zbar,each=n),n)
      Zbar_t3 = cbind(Zbar_t3, Zbar) 
    }
    Z_t4 = NULL;    for(i in 1:G2)  Z_t4 = cbind(Z_t4, bsbasefun(X2[,group2==gunique2[i]],K,degr))
    Zbar_t4 = colMeans(Z_t4) 
    Z_t4 = Z_t4 - matrix(rep(Zbar_t4,each=n),n)
    
    Ybar = colMeans(Y)
    Y1 = Y - matrix(rep(Ybar,each=n),n)
    if (is.null(lambda)) {
      is_setlam = 1
      if (nlam < 1||is.null(nlam)) stop("nlambda must be at least 1")
      setlam = c(1,lam_min,alpha,nlam)
      
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
    
    if((r41>dim(A)[2])|(r42>dim(B)[2])|(r43>dim(C)[2])|(r44>dim(D)[2]))
      stop("r41, r42, r43 and r44 of the fourth-order tensor must not be larger than its components A, B, C and D, respectively !")
    #---------------- The selection by BIC and CV  ---------------------#  
    #if(method=="BIC") fit_dr = marmComposed_bic(Y,X1,X2,G1,group,K,r_index,lambda,D0_t3,D0_t4,intercept,opts,opts_pen)
    #if(method=="CV") fit_dr = gmam_sparse_composed_cv(Y,X,group,ncv,K_index,r1_index,r2_index,r3_index,lambda,D0,intercept,opts,opts_pen)
    
    if(method=="BIC"){
      if(is.fabs){
        fit = EstPenColumnComposed1(Y1,Z_t3,Sinit,Ainit,Binit,Cinit,Z_t4,S,A,B,C,D,lambda,opts,opts_pen) 
      }else{
        fit = EstPenColumnComposed2(Y1,Z_t3,Sinit,Ainit,Binit,Cinit,Z_t4,S,A,B,C,D,lambda,opts,opts_pen)
      }
      
      df = NULL  #note!
      for (ll in 1:nlam) {
        df1 = 0
        for (g in 1:G1) {
          if(fit$df_t3[g,ll]<r31) {
            df1 = df1 + r31*r32*r33 + fit$df_t3[g,ll]*r31 + K*r32 + q*r33 - r32^2 - r33^2
          }else {
            df1 = df1 + r31*r32*r33 + fit$df_t3[g,ll]*r31 + K*r32 + q*r33 - r31^2 - r32^2 - r33^2
          }
        }
        if(fit$df_t4[ll]<r41) {
          df1 = df1 + r41*r42*r43*r44 + fit$df_t4[ll]*r41 + K*r42 + G2*r43 + q*r44 - r42^2 - r43^2 - r44^2
        }else {
          df1 = df1 + r41*r42*r43*r44 + fit$df_t4[ll]*r41 + K*r42 + G2*r43 + q*r44 - r41^2 - r42^2 - r43^2 - r44^2
        }
        df = c(df, df1)
      }
      
      bic = log(fit$likhd_t3/(n*q)) + log(n*q)*df/(n*q)
      selected = which.min(bic)
      lambda_opt = lambda[selected]
      opts_pen$nlam = length(lambda[1:selected])
      if(is.fabs){
        fit = EstPenColumnComposed1(Y1,Z_t3,Sinit,Ainit,Binit,Cinit,Z_t4,S,A,B,C,D,lambda[1:selected],opts,opts_pen) 
      }else{
        fit = EstPenColumnComposed2(Y1,Z_t3,Sinit,Ainit,Binit,Cinit,Z_t4,S,A,B,C,D,lambda[1:selected],opts,opts_pen)
      }
      
      activeX_t3 = fit$betapath_t3[,selected]
      Snew_t3 = fit$Snew_t3[[selected]]
      Anew_t3 = fit$Anew_t3[[selected]]
      Bnew_t3 = fit$Bnew_t3[[selected]]
      Cnew_t3 = fit$Cnew_t3[[selected]]
      Dn_t3 = NULL
      for (g in 1:G1) Dn_t3 = cbind(Dn_t3, Cnew_t3[[g]] %*% Snew_t3[[g]] %*%t(kronecker(Bnew_t3[[g]], Anew_t3[[g]])))
      activeX_t4 = fit$betapath_t4[,selected]
      Anew_t4=matrix(fit$Apath_t4[,selected],nrow=p2)
      Bnew_t4=matrix(fit$Bpath_t4[,selected],nrow=K)
      Cnew_t4=matrix(fit$Cpath_t4[,selected],nrow=G2)
      Dnew_t4=matrix(fit$Dpath_t4[,selected],nrow=q)
      Snew_t4=matrix(fit$Spath_t4[,selected],nrow=r44)
      Dn_t4 = Dnew_t4 %*% Snew_t4 %*%t(kronecker(Cnew_t4, kronecker(Bnew_t4,Anew_t4)))
      
      if(intercept)  mu = Ybar-Dn_t3%*%as.vector(Zbar_t3)-Dn_t4%*%Zbar_t4
      else mu = rep(0,q)
    }
    if(method=="CV"&&nlam>1){
      len_cv = ceiling(n/ncv)
      RSS = rep(0,nlam)
      for(jj in 1:ncv){
        cv.id = ((jj-1)*len_cv+1):(jj*len_cv)
        if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
        Ytrain = Y1[-cv.id,]
        X1train = X1[-cv.id,]
        X2train = X2[-cv.id,]
        Ytest = Y1[cv.id,]
        X1test = X1[cv.id,]
        X2test = X2[cv.id,]
        Ztrain_t3 = list()
        Ztest_t3 = list()
        for(i in 1:G1) {
          Ztrain_t3[[i]] = as.matrix(bsbasefun(X1train[,group1==gunique1[i]],K,degr))
          Ztest_t3[[i]] = as.matrix(bsbasefun(X1test[,group1==gunique1[i]],K,degr))
        }
        Ztrain_t4 = NULL;  for(i in 1:G2)  Ztrain_t4 = cbind(Ztrain_t4, bsbasefun(X2train[,group2==gunique2[i]],K,degr))
        Ztest_t4 = NULL;  for(i in 1:G2)  Ztest_t4 = cbind(Ztest_t4, bsbasefun(X2test[,group2==gunique2[i]],K,degr))
        
        if(is.fabs){
          fit = EstPenColumnComposed1CV(Ytrain,Ztrain_t3,Ztrain_t4,Ytest,Ztest_t3,Ztest_t4,Sinit,Ainit,Binit,Cinit,S,A,B,C,D,lambda,opts,opts_pen) 
        }else{
          fit = EstPenColumnComposed2CV(Ytrain,Ztrain_t3,Ztrain_t4,Ytest,Ztest_t3,Ztest_t4,Sinit,Ainit,Binit,Cinit,S,A,B,C,D,lambda,opts,opts_pen)
        }
        RSS = RSS + fit$likhd
      } 
      selected = which.min(RSS)
      lambda_opt = lambda[selected]
      
      opts_pen$nlam = length(lambda[1:selected])
      if(is.fabs){
        fit = EstPenColumnComposed1(Y1,Z_t3,Sinit,Ainit,Binit,Cinit,Z_t4,S,A,B,C,D,lambda[1:selected],opts,opts_pen) 
      }else{
        fit = EstPenColumnComposed2(Y1,Z_t3,Sinit,Ainit,Binit,Cinit,Z_t4,S,A,B,C,D,lambda[1:selected],opts,opts_pen) 
      }
      
      activeX_t3 = fit$betapath_t3[,selected]
      Snew_t3 = fit$Snew_t3[[selected]]
      Anew_t3 = fit$Anew_t3[[selected]]
      Bnew_t3 = fit$Bnew_t3[[selected]]
      Cnew_t3 = fit$Cnew_t3[[selected]]
      Dn_t3 = NULL
      for (g in 1:G1) Dn_t3 = cbind(Dn_t3, Cnew_t3[[g]] %*% Snew_t3[[g]] %*%t(kronecker(Bnew_t3[[g]], Anew_t3[[g]])))
      activeX_t4 = fit$betapath_t4[,selected]
      Anew_t4=matrix(fit$Apath_t4[,selected],nrow=p2)
      Bnew_t4=matrix(fit$Bpath_t4[,selected],nrow=K)
      Cnew_t4=matrix(fit$Cpath_t4[,selected],nrow=G2)
      Dnew_t4=matrix(fit$Dpath_t4[,selected],nrow=q)
      Snew_t4=matrix(fit$Spath_t4[,selected],nrow=r44)
      Dn_t4 = Dnew_t4 %*% Snew_t4 %*%t(kronecker(Cnew_t4, kronecker(Bnew_t4,Anew_t4)))
      
      if(intercept)  mu = Ybar-Dn_t3%*%as.vector(Zbar_t3)-Dn_t4%*%Zbar_t4
      else mu = rep(0,q)
    }
    if(method=="CV"&&nlam==1){
      if(is.fabs){
        fit = EstPenColumnComposed1(Y1,Z_t3,Sinit,Ainit,Binit,Cinit,Z_t4,S,A,B,C,D,lambda,opts,opts_pen) 
      }else{
        fit = EstPenColumnComposed2(Y1,Z_t3,Sinit,Ainit,Binit,Cinit,Z_t4,S,A,B,C,D,lambda,opts,opts_pen) 
      }
      selected = 1
      lambda_opt = lambda
      
      activeX_t3 = fit$betapath_t3[,selected]
      Snew_t3 = fit$Snew_t3[[selected]]
      Anew_t3 = fit$Anew_t3[[selected]]
      Bnew_t3 = fit$Bnew_t3[[selected]]
      Cnew_t3 = fit$Cnew_t3[[selected]]
      Dn_t3 = NULL
      for (g in 1:G1) Dn_t3 = cbind(Dn_t3, Cnew_t3[[g]] %*% Snew_t3[[g]] %*%t(kronecker(Bnew_t3[[g]], Anew_t3[[g]])))
      activeX_t4 = fit$betapath_t4[,selected]
      Anew_t4=matrix(fit$Apath_t4[,selected],nrow=p2)
      Bnew_t4=matrix(fit$Bpath_t4[,selected],nrow=K)
      Cnew_t4=matrix(fit$Cpath_t4[,selected],nrow=G2)
      Dnew_t4=matrix(fit$Dpath_t4[,selected],nrow=q)
      Snew_t4=matrix(fit$Spath_t4[,selected],nrow=r44)
      Dn_t4 = Dnew_t4 %*% Snew_t4 %*%t(kronecker(Cnew_t4, kronecker(Bnew_t4,Anew_t4)))
      
      if(intercept)  mu = Ybar-Dn_t3%*%as.vector(Zbar_t3)-Dn_t4%*%Zbar_t4
      else mu = rep(0,q)
    }
    
    return(list(D_t3 = Dn_t3,
                D_t4 = Dn_t4,
                mu = mu,
                S_t3.opt = Snew_t3,
                A_t3.opt = Anew_t3,
                B_t3.opt = Bnew_t3,
                C_t3.opt = Cnew_t3,
                S_t4.opt = Snew_t4,
                A_t4.opt = Anew_t4,
                B_t4.opt = Bnew_t4,
                C_t4.opt = Cnew_t4,
                D_t4.opt = Dnew_t4,
                lambda.seq = lambda,
                lambda.opt = lambda_opt,
                rss = fit$likhd_t3[selected],
                df_t3 = fit$df_t3,
                df_t4 = fit$df_t4,
                activeX_t3 = activeX_t3,
                activeX_t4 = activeX_t4,
                opts = opts,
                opts_pen = opts_pen))
  }























