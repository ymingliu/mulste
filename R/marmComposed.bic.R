
##--------------Estimation with Penalty by BIC----------------------##
marmComposed.bic <- 
  function(Y,X1,X2,G1,group,K_index,r_index,lambda,D0_t3,D0_t4,intercept,opts,opts_pen){
    is.fabs = opts$is.fabs
    n = opts$n
    p1 = opts$p1
    p2 = opts$p2
    q = opts$q
    nlam = length(lambda)
    degr = opts$degr
    
    gunique <- unique(group)
    G = length(gunique)
    G2 = G-G1
    group1 = group[which(group<=G1)]
    group2 = group[which(group>G1)]-G1
    gunique1 <- unique(group1)
    gunique2 <- unique(group2)
    
    r1_index = r_index$r1_t4_index #note!
    r2_index = r_index$r2_t4_index
    r3_index = r_index$r3_t4_index
    r4_index = r_index$r4_t4_index
    
    Ybar = colMeans(Y)
    Y1 = Y - matrix(rep(Ybar,each=n),n)
    Sinit = list()
    Ainit = list()
    Binit = list()
    Cinit = list()
    Z_t3 = list()
    RSS = NULL
    
    for(K in K_index){
      for(i in 1:G1) {
        Z_t3[[i]] = as.matrix(bsbasefun(X1[,group1==gunique1[i]],K,degr))
        Zbar = colMeans(Z_t3[[i]])
        Z_t3[[i]] = Z_t3[[i]] - matrix(rep(Zbar,each=n),n)
      }
      Z_t4 = NULL;  for(i in 1:G2)  Z_t4 = cbind(Z_t4, bsbasefun(X2[,group2==gunique2[i]],K,degr))
      Zbar = colMeans(Z_t4)
      Z_t4 = Z_t4 - matrix(rep(Zbar,each=n),n)
      
      opts$K = K
      for(r4 in r4_index) {
        opts$r4 = r4
        for(r3 in r3_index){
          opts$r3 = r3
          for(r2 in r2_index){
            opts$r2 = r2
            for(r1 in r1_index){
              opts$r1 = r1
              for(g in 1:G1){
                Sinit[[g]] = as.matrix(D0_t3[[g]]$S[1:r4,1:(r1*r2)]) #note!
                Ainit[[g]] = as.matrix(D0_t3[[g]]$A[,1:r1])
                Binit[[g]] = as.matrix(D0_t3[[g]]$B[,1:r2])
                Cinit[[g]] = as.matrix(D0_t3[[g]]$C[,1:r4]) #note!
                if(dim(Sinit[[g]])[1]!=dim(Cinit[[g]])[2]) Sinit[[g]] = t(Sinit[[g]])
              }
              S = as.matrix(D0_t4$S[1:r4,1:(r1*r2*r3)])
              A = as.matrix(D0_t4$A[,1:r1])
              B = as.matrix(D0_t4$B[,1:r2])
              C = as.matrix(D0_t4$C[,1:r3])
              D = as.matrix(D0_t4$D[,1:r4])
              if(dim(S)[1]!=dim(D)[2]) S=t(S)
              
              if(is.fabs){
                fit = EstPenColumnComposed1(Y1,Z_t3,Sinit,Ainit,Binit,Cinit,Z_t4,S,A,B,C,D,lambda,opts,opts_pen) 
              }else{
                fit = EstPenColumnComposed2(Y1,Z_t3,Sinit,Ainit,Binit,Cinit,Z_t4,S,A,B,C,D,lambda,opts,opts_pen)
              }
              
              df = NULL  #note!
              for (ll in 1:nlam) {
                df1 = 0
                for (g in 1:G1) {
                  if(fit$df_t3[g,ll]<r1) {
                    df1 = df1 + r1*r2*r4 + fit$df_t3[g,ll]*r1 + K*r2 + q*r4 - r2^2 - r4^2
                  }else {
                    df1 = df1 + r1*r2*r4 + fit$df_t3[g,ll]*r1 + K*r2 + q*r4 - r1^2 - r2^2 - r4^2
                  }
                }
                if(fit$df_t4[ll]<r1) {
                  df1 = df1 + r1*r2*r3*r4 + fit$df_t4[ll]*r1 + K*r2 + G2*r3 + q*r4 - r2^2 - r3^2 - r4^2
                }else {
                  df1 = df1 + r1*r2*r3*r4 + fit$df_t4[ll]*r1 + K*r2 + G2*r3 + q*r4 - r1^2 - r2^2 - r3^2 - r4^2
                }
                df = c(df, df1)
              }

              RSS = cbind(RSS,log(fit$likhd_t3/(n*q)) + log(n*q)*df/(n*q)) 
            }
          }
        }
      }
    }
    selected = which.min(RSS)
    qj = ceiling(selected/nlam)
    qj1 = selected%%nlam
    if(qj1==0) qj1=nlam
    
    lambda_opt = lambda[qj1]
    opt = assig(c(length(r1_index),length(r2_index),length(r3_index),length(r4_index),length(K_index)))[,qj]
    r1_opt = r1_index[opt[1]]
    r2_opt = r2_index[opt[2]]
    r3_opt = r3_index[opt[3]]
    r4_opt = r4_index[opt[4]]
    K_opt = K_index[opt[5]]
    
    #---------------- The estimation after selection ---------------------#
    opts$r1 = r1_opt
    opts$r2 = r2_opt
    opts$r3 = r3_opt
    opts$r4 = r4_opt
    opts$K = K_opt
    
    Z_t3 = list()
    Zbar_t3 = NULL
    for(i in 1:G1){
      Sinit[[i]] = as.matrix(D0_t3[[i]]$S[1:r4_opt,1:(r1_opt*r2_opt)])
      Ainit[[i]] = as.matrix(D0_t3[[i]]$A[,1:r1_opt])
      Binit[[i]] = as.matrix(D0_t3[[i]]$B[,1:r2_opt])
      Cinit[[i]] = as.matrix(D0_t3[[i]]$C[,1:r4_opt])
      if(dim(Sinit[[i]])[1]!=dim(Cinit[[i]])[2]) Sinit[[i]] = t(Sinit[[i]])
      Z_t3[[i]] = as.matrix(bsbasefun(X1[,group1==gunique1[i]],K_opt,degr))
      Zbar1 = colMeans(Z_t3[[i]])
      Z_t3[[i]] = Z_t3[[i]] - matrix(rep(Zbar1,each=n),n)
      Zbar_t3 = cbind(Zbar_t3, Zbar1) 
    }
    Z_t4 = NULL;  for(i in 1:G2)  Z_t4 = cbind(Z_t4, bsbasefun(X2[,group2==gunique2[i]],K_opt,degr))
    Zbar_t4 = colMeans(Z_t4) 
    Z_t4 = Z_t4 - matrix(rep(Zbar_t4,each=n),n)
    
    S = as.matrix(D0_t4$S[1:r4_opt,1:(r1_opt*r2_opt*r3_opt)])
    A = as.matrix(D0_t4$A[,1:r1_opt])
    B = as.matrix(D0_t4$B[,1:r2_opt])
    C = as.matrix(D0_t4$C[,1:r3_opt])
    D = as.matrix(D0_t4$D[,1:r4_opt])
    if(dim(S)[1]!=dim(D)[2]) S=t(S)
    
    opts_pen$nlam = length(lambda[1:qj1])
    opts_pen$lam_max = max(lambda[1:qj1])
    opts_pen$lam_min = min(lambda[1:qj1])
    if(is.fabs){
      fit_opt = EstPenColumnComposed1(Y1,Z_t3,Sinit,Ainit,Binit,Cinit,Z_t4,S,A,B,C,D,lambda[1:qj1],opts,opts_pen) 
    }else{
      fit_opt = EstPenColumnComposed2(Y1,Z_t3,Sinit,Ainit,Binit,Cinit,Z_t4,S,A,B,C,D,lambda[1:qj1],opts,opts_pen) 
    }
    
    activeX_t3 = fit_opt$betapath_t3[,qj1]
    Snew_t3 = fit_opt$Snew_t3[[qj1]]
    Anew_t3 = fit_opt$Anew_t3[[qj1]]
    Bnew_t3 = fit_opt$Bnew_t3[[qj1]]
    Cnew_t3 = fit_opt$Cnew_t3[[qj1]]
    Dn_t3 = NULL
    for (g in 1:G1) Dn_t3 = cbind(Dn_t3, Cnew_t3[[g]] %*% Snew_t3[[g]] %*%t(kronecker(Bnew_t3[[g]], Anew_t3[[g]])))
    activeX_t4 = fit_opt$betapath_t4[,qj1]
    Anew_t4=matrix(fit_opt$Apath_t4[,qj1],nrow=p2)
    Bnew_t4=matrix(fit_opt$Bpath_t4[,qj1],nrow=K_opt)
    Cnew_t4=matrix(fit_opt$Cpath_t4[,qj1],nrow=G2)
    Dnew_t4=matrix(fit_opt$Dpath_t4[,qj1],nrow=q)
    Snew_t4=matrix(fit_opt$Spath_t4[,qj1],nrow=r4_opt)
    Dn_t4 = Dnew_t4 %*% Snew_t4 %*%t(kronecker(Cnew_t4, kronecker(Bnew_t4,Anew_t4)))
    
    if(intercept)  mu = Ybar-Dn_t3%*%as.vector(Zbar_t3)-Dn_t4%*%Zbar_t4
    else mu = rep(0,q)
    
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
                rk_t3_opt = c(r1_opt,r2_opt,r4_opt,K_opt),
                rk_t4_opt = c(r1_opt,r2_opt,r3_opt,r4_opt,K_opt),
                lambda.seq = lambda,
                lambda.opt = lambda_opt,
                rss = fit_opt$likhd_t3[qj1],
                df_t3 = fit_opt$df_t3,
                df_t4 = fit_opt$df_t4,
                activeX_t3 = activeX_t3,
                activeX_t4 = activeX_t4,
                opts = opts,
                opts_pen = opts_pen))
}
