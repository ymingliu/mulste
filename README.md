# mulste
Multi-omics data integration with multi-view learning via composed tensors.
 
  In the multivariate additive model (MARM) for a multi-view data, the B-splines are applied to approximate the component functions. We treat the coefficients as multiple third-order tensors (MARM) or a fourth-order tensor. The composed model (COMARM) can be used when the number of covariates in each view is not equal. With the tensor low-rankness, the Tucker decomposition and group sparse penalty (lasso, mcp or scad) reduce the number of parameters. An alternative updating algorithm based on the coordinate descent strategy is used to estimate the core tensors and factor matrices, and further additive functions.
# Installation

    #install Rtools 3.5 (http://cran.r-project.org/bin/windows/Rtools)
    #install.packages("devtools")
    #install.packages("Rcpp")
    library(devtools)
    install_github("xliusufe/mulste")

# Usage

   - [x] [mulste-manual](https://github.com/xliusufe/mulste/blob/master/inst/mulste-manual.pdf) ------------ Details of the usage of the package.
# Examples
    
    library(mulste)
    # Example 1
    # The usage of function "marm3.dr"
    n <- 200; q <- 5; p <- 100; s <- 3; ng = 4
    group <- rep(1:ng,each=p/ng)
    mydata <- marm3.sim.fbs(n,q,p,s,group)
    fit <- with(mydata, marm3.dr(Y,X,group,K,r1_index,r2_index,r3_index,D0=D0,nlam=5))
    
    # Example 2
    # The usage of function "marm4.dr"
    n <- 200; q <- 5; p <- 100; s <- 3; ng = 4
    group <- rep(1:ng,each=p/ng)
    mydata <- marm4.sim.fbs(n,q,p,s,group)
    fit <- with(mydata, marm4.dr(Y,X,group,K,r1_index,r2_index,r3_index,r4_index,D0=D0,nlam=5))
    
    # Example 3
    # The usage of function "marmComposed.dr"
    n <- 200; q <- 5; p <- 100; s1 <- 5; s2 <- 3; G1 <- 1; ng = 4
    group <- rep(1:ng,each=p/ng)
    mydata <- marmComposed.sim.fbs(n,q,p,s1,s2,G1,group)
    fit <- with(mydata, marmComposed.dr(Y,X1,X2,G1,group,is.fabs,K,r_index,D0_t3=D0_t3,D0_t4=D0_t4,nlam=5))
    
    # Example 4
    # The usage of function "marm3"
    n <- 200; q <- 5; p <- 100; s <- 3; ng <- 4
    group <- rep(1:ng,each=p/ng)
    mydata <- marm3.sim.fbs(n,q,p,s,group,isfixedR=1)
    fit <- with(mydata, marm3(Y,X,group,K,r10,r20,r30,D0=D0,nlam=5))
    
    # Example 5
    # The usage of function "marm4"
    n <- 200; q <- 5; p <- 100; s <- 3; ng <- 4
    group <- rep(1:ng,each=p/ng)
    mydata <- marm4.sim.fbs(n,q,p,s,group,isfixedR=1)
    fit <- with(mydata, marm4(Y,X,group,K,r10,r20,r30,r40,D0=D0,nlam=5))
    
    # Example 6
    # The usage of function "marmComposed"
    n <- 200; q <- 5; p <- 100; s1 <- 5; s2 <- 3; G1 <- 1; ng <- 4
    group <- rep(1:ng,each=p/ng)
    mydata <- marmComposed.sim.fbs(n,q,p,s1,s2,G1,group,isfixedR=1)
    fit <- with(mydata, marmComposed(Y,X1,X2,G1,group,is.fabs,K,r310,r320,r330,r410,r420,r430,r440,D0_t3=D0_t3,D0_t4=D0_t4,nlam=5))
    

 
 # References
Multi-omics data integration with multi-view learning via composed tensors. Manuscript.

# Development
The R-package is developed by Xu Liu (liu.xu@sufe.edu.cn) and Yiming Liu.
