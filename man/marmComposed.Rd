\name{marmComposed}
\alias{marmComposed-function}
\alias{marmComposed}
\docType{package}
\title{
  Fit composed model (COMARM) with sparsity assumption and fixed ranks.
}
\description{
  Fit a composed model for multi-view data using B-splines with given ranks (\eqn{r_{1g}, r_{2g}, r_{3g}} and \eqn{r_{1}, r_{2}, r_{3}, r_{4}}). Multiple third-order coefficient tensors and a fourth-order coefficient tensor can be estimated by this function. The group sparse penalty such as LASSO, MCP or SCAD and the coordinate descent algorithm are used to yield a sparsity estimator. The BIC or cross-validation method are used to search the optimal regularization parameter.
}

\usage{
marmComposed <-
  function(Y,X1,X2,G1=NULL,group=NULL,is.fabs=1,K=6,r31=NULL,r32=NULL,r33=NULL,r41=NULL,r42=NULL,r43=NULL,r44=NULL,method="BIC",ncv=10,penalty="LASSO",lambda=NULL,D0_t3=NULL,D0_t4=NULL,intercept=TRUE,nlam=50,degr=3,lam_min=0.01,eps=1e-4,max_step=20,eps1=1e-4,max_step1_t3=20,max_step1_t4=20,gamma=2,dfmax1=NULL,dfmax2=NULL,alpha=1,vnorm_ratio=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}

  \item{X1}{A \eqn{n\times p1} numeric design matrix for the composed model, where \eqn{p1=\sum_{1\leq g\leq G1}p_g}}

  \item{X2}{A \eqn{n\times p2} numeric design matrix for the composed model, where \eqn{p2=\sum_{G1+1 \leq g \leq ng}p_g} and ng is the number of views.}

  \item{G1}{The number of views that we do not consider their intergroup correlation.}

  \item{group}{A \eqn{p} vector of the grouping index of predictors, e.g., \eqn{group=c(1,1,1,2,2,2)} means there are \eqn{6} perdictors in the model, and the first three perdictors are in the same group and the last three perdictors are in another one. By default, we set \eqn{group=rep(1,p)}.}

  \item{K}{The number of B-spline basis functions, that is the plus of both degrees of basis functioins and the number of knots. Default is \code{6}, which means cubic spline.}

  \item{r31}{The first dimension of single value matrix of the third-order tensor. Default is \code{2}.}

  \item{r32}{The second dimension of single value matrix of the third-order tensor. Default is \code{2}.}

  \item{r33}{The third dimension of single value matrix of the third-order tensor. Default is \code{2}.}

  \item{r41}{The first dimension of single value matrix of the fourth-order tensor. Default is \code{2}.}

  \item{r42}{The second dimension of single value matrix of the fourth-order tensor. Default is \code{2}.}

  \item{r43}{The third dimension of single value matrix of the fourth-order tensor. Default is \code{2}.}

  \item{r44}{The fourth dimension of single value matrix of the fourth-order tensor. Default is \code{2}.}

  \item{method}{The method to be applied to select regularization parameters.  Either \code{BIC} (default), or \code{CV}.}

  \item{ncv}{The number of cross-validation folds. Default is \code{10}. If \code{method} is not \code{CV}, \code{ncv} is useless.}

  \item{penalty}{The penalty to be applied to the model. Either \code{LASSO} (the default), \code{MCP} or \code{SCAD}.}

  \item{lambda}{A user-specified sequence of lambda values.  By default, a sequence of values of length \code{nlam} is computed, equally spaced on the log scale.}

  \item{D0_t3}{A user-specified list of initialized values for the third-order tensor, including \code{ng} sub-lists where ng is the number of groups. For each sub-list, it has four initialized matrix \code{S_{(3)}} (called \code{S}), \code{A}, \code{B} and \code{C}. By default, a list of initialization satisfying fixed ranks is computed by random.}

  \item{D0_t4}{A user-specified list of initialized values for the fourth-order tensor, including five initialized matrix \eqn{S_{(4)}}  (called \code{S}), \code{A}, \code{B}, \code{C} and \code{D}. By default, a list of initialization satisfying fixed ranks is computed by random.}

  \item{intercept}{A logical value indicating whether the intercept is fitted. Default is \code{TRUE} or set to zero by \code{FALSE}.}

  \item{degr}{The number of knots of B-spline base function. Default is \code{3}.}

  \item{nlam}{The number of lambda values. Default is \code{20}.}

  \item{lam_min}{The smallest value for lambda, as a fraction of lambda.max.  Default is \code{0.01}.}

  \item{eps}{Convergence threshhold. The algorithm iterates until the relative change in any coefficient is less than \code{eps}. Default is \code{1e-4}.}

  \item{max_step}{Maximum number of iterations. Default is \code{20}.}

  \item{eps1}{Convergence threshhold. The Coordinate descent method algorithm iterates until the relative change in any coefficient is less than \code{eps1}. Default is \code{1e-4}.}

  \item{max_step1_t3}{The maximum iterates number of coordinate descent method during the estimation of the third-order tensor. Default is \code{20}.}

  \item{max_step1_t4}{The maximum iterates number of coordinate descent method during the estimation of the fourth-order tensor. Default is \code{20}.}

  \item{gamma}{The tuning parameter of the MCP/SCAD penalty.}

  \item{dfmax1}{Upper bound for the number of nonzero coefficients during the estimation of the third-order tensor. Default is no upper bound.  However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.}

  \item{dfmax2}{Upper bound for the number of nonzero coefficients during the estimation of the fourth-order tensor. Default is no upper bound.  However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.}

  \item{alpha}{Tuning parameter for the Mnet estimator which controls the relative contributions from the LASSO, MCP or SCAD penalty and the ridge, or L2 penalty.  \code{alpha=1} is equivalent to LASSO, MCP or SCAD penalty, while \code{alpha=0} would be equivalent to ridge regression. However, \code{alpha=0} is not supported; \code{alpha} may be arbitrarily small, but not exactly 0.}

  \item{vnorm_ratio}{The ratio between \eqn{lambda_1} and \eqn{lambda_2}, that is, \eqn{lambda_2}=vnorm_ratio\eqn{\times lambda_1}. Default is 1.}

}

\details{
  This function gives \code{pq} functional coefficients' estimators of composed model. Multiple third-order tensors with multiple ranks (\eqn{r_{1g}, r_{2g}, r_{3g}}) and a fourth-order tensor with multiple ranks (\eqn{r_1, r_2, r_3, r_4}) need to be estimated. We fix these ranks and use an alternative updating algorithm to update its core tensor and factor matrices based on Tucker decomposition. Group LASSO, SCAD or MCP penalty is applied on the row of each factor matrix \eqn{A^g} and \eqn{A} to achieve variable selection.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{D_t3}{Estimator of coefficients corresponding to the third-order tensor \eqn{D_{(3)}=(D^1_{(3)},...,D^{G1}_{(3)})}.}

  \item{D_t4}{Estimator of coefficients corresponding to the fourth-order tensor \eqn{D_{(4)}}.}

  \item{mu}{Estimator of intercept \eqn{\mu}.}

  \item{S_t3.opt}{A length-\eqn{G1} list including estimator of the core tensor \eqn{S_{(3)}} of each third-order coefficient tensor.}

  \item{A_t3.opt}{A length-\eqn{G1} list including estimator of the factor matrix \eqn{A} of each third-order coefficient tensor.}

  \item{B_t3.opt}{A length-\eqn{G1} list including estimator of the factor matrix \eqn{B} of each third-order coefficient tensor.}

  \item{C_t3.opt}{A length-\eqn{G1} list including estimator of the factor matrix \eqn{C} of each third-order coefficient tensor.}

  \item{S_t4.opt}{A length-\eqn{ng-G1} list including estimator of the core tensor \eqn{S_{(3)}} of each fourth-order coefficient tensor.}

  \item{A_t4.opt}{A length-\eqn{ng-G1} list including estimator of the factor matrix \eqn{A} of each fourth-order coefficient tensor.}

  \item{B_t4.opt}{A length-\eqn{ng-G1} list including estimator of the factor matrix \eqn{B} of each fourth-order coefficient tensor.}

  \item{C_t4.opt}{A length-\eqn{ng-G1} list including estimator of the factor matrix \eqn{C} of each fourth-order coefficient tensor.}

  \item{D_t4.opt}{A length-\eqn{ng-G1} list including estimator of the factor matrix \eqn{D} of each fourth-order coefficient tensor.}

  \item{lambda.seq}{The sequence of regularization parameter values in the path.}

  \item{lambda_opt}{The value of \code{lambda} with the minimum \code{BIC} or \code{CV} value.}

  \item{rss}{Residual sum of squares (RSS).}

  \item{df_t3}{Degrees of freedom for all third-order tensors.}

  \item{df_t4}{Degrees of freedom for the fourth-order tensor.}

  \item{activeX_t3}{The active set of \eqn{X1}. A length-\eqn{p1} vector.}

  \item{activeX_t4}{The active set of \eqn{X2}. A length-\eqn{p2} vector.}

  \item{opts}{Other related parameters used in algorithm. Some of them are set by default.}

  \item{opts_pen}{Other related parameters used in algorithm (especially parameters in peanlty). Some of them are set by default.}

  %\item{...}{ Other options for CompositeQuantile.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}
\references{
  Multi-omics data integration with multi-view learning via composed tensors.
}
\keyword{ Group sparsity; Tensor low-rankness; Tucker decomposition; Multivariate additive regression; Multi-view data. }
\seealso{
  marmComposed.dr
}

\examples{
library(mulste)
n <- 200; q <- 5; p <- 100; s1 <- 5; s2 <- 3; G1 <- 1; ng <- 4
group <- rep(1:ng,each=p/ng)
mydata <- marmComposed.sim.fbs(n,q,p,s1,s2,G1,group,isfixedR=1)
fit <- with(mydata, marmComposed(Y,X1,X2,G1,group,is.fabs,K,r310,r320,r330,r410,r420,r430,r440,D0_t3=D0_t3,D0_t4=D0_t4,nlam=5))
}
