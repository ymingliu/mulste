\name{mulste-package}
\alias{mulste-package}
\alias{mulste}
\docType{package}
\title{
  Multi-omics data integration with multi-view learning via composed tensors
}
\description{
  For an integrative multi-view multivariate additive model (MARM), the B-splines are applied to approximate the component functions. We treat the coefficients as multiple third-order tensors (MARM) or even a fourth-order tensor (structural MARM). The composed model (COMARM) can be used when the number of covariates in each view is not equal. With the tensor low-rankness, the Tucker decomposition and group sparse penalty (lasso, mcp or scad) reduce the number of parameters. An alternative updating algorithm based on the coordinate descent strategy is used to estimate the core tensors and factor matrices, and further additive functions.

}
\details{
  This pacakge includes six main functions and six generating functions. \code{marm3} and \code{marm3.dr} yield the estimator of MARM. The difference is the former requires the fixed ranks and the latter can search the optimal ranks and regularization parameter simultaneously by BIC or CV method. \code{marm3.sim.fbs} and \code{marm3.sim.fsin} generate data of scenario I and II respectively. Scenario I assumes that the true functions are exactly residing in the space of B-spline basis functions. Scenario II assumes that the true functions are some linear combination of \eqn{sin(2\pi x)} and \eqn{cos(\pi x)}. Similarly, \code{marm4} and \code{marm4.dr} yield the estimator of structural MARM. \code{marm4.sim.fbs} and \code{marm4.sim.fsin} are two generating functions of structural MARM with scenario I and II settings. \code{marmComposed} and \code{marmComposed.dr} yield the estimator of composed model. \code{marmComposed.sim.fbs} and \code{marmComposed.sim.fsin} are two generating functions of composed model with scenario I and II settings. They all have the same assumptions as MARM.
}
\author{
Xu Liu

Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
}
\references{
  Multi-omics data integration with multi-view learning via composed tensors.
}
\keyword{ Group sparsity; Tensor low-rankness; Tucker decomposition; Multivariate additive regression; Multi-view data. }
