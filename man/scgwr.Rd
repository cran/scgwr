\name{scgwr}
\alias{scgwr}
\title{Scalable Geographically Weighted Regression}
\usage{
scgwr( coords, y, x, knn = 100, kernel = "gau", p = 4 )
}
\description{
This function estimates a scalable geographically weighted regression (GWR) model.
}

\arguments{
\item{coords}{Matrix of spatial point coordinates (N x 2)}
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables (N x K). Default is NULL}
\item{knn}{Number of nearest-neighbors being geographically weighted. Default is 100. Larger knn is better for larger samples (see Murakami er al., 2019)}
\item{kernel}{Kernel to model spatial heterogeneity. Gaussian kernel ("gau") and exponential kernel ("exp") are available}
\item{p}{Degree of the polynomial to approximate the kernel function. Default is 4}
}
\value{
\item{b}{Matrix of estimated coefficients (N x K)}
\item{bse}{Matrix of the standard errors for the coefficients (N x k)}
\item{t}{Matrix of the t-values for the coefficients (N x K)}
\item{p}{Matrix of the p-values for the coefficients (N x K)}
\item{pa}{Matrix of the p-values adjusted to address the multiple testing problem using the approach of de Silva and Fotheringham (2016) (N x K)}
\item{par}{Estimated model parameters includeing a scale parameter and a shrinkage parameter if penalty = TRUE (see Murakami et al., 2018)}
\item{e}{Error statistics. It includes sum of squared errors (SSE), residual standard error (resid_SE), R-squared (R2), adjusted R2 (adjR2), log-likelihood (logLik), corrected Akaike information criterion (AICc), and the cross-validation (CV) score measured by root mean squared error (RMSE) (CV_score(RMSE))}
\item{pred}{Vector of predicted values (N x 1)}
\item{resid}{Vector of residuals (N x 1)}
}
\references{

Murakami, D., Tsutsumida, N., Yoshida, T., Nakaya, T., and Lu, B. (2019) Scalable GWR: A linear-time algorithm for large-scale geographically weighted regression with polynomial kernels. <arXiv:1905.00266>.

da Silva, A.R., and Fotheringham, A.S. (2016) The multiple testing issue in geographically weighted regression. Geographical Analysis, 48(3), 233-247.
}

\examples{
require(spData)
data(boston)
coords <- boston.c[ , c("LON", "LAT") ]
y      <- log(boston.c[,"MEDV"])
x      <- boston.c[, c("CRIM", "ZN", "INDUS", "CHAS", "AGE")]
res    <- scgwr( coords = coords, y = y, x = x )
res$b
res$bse
res$p
res$e

}
