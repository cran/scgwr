\name{scgwr_p}
\alias{scgwr_p}
\title{Parallel implementation of scalable geographically weighted regression}
\usage{
scgwr_p( coords, y, x = NULL, knn = 100, kernel = "gau",
       p = 4, approach = "CV", nsamp = NULL, cl = NULL)
}
\description{
Parallel implementation of scalable geographically weighted regression for large samples
}

\arguments{
\item{coords}{Matrix of spatial point coordinates (N x 2)}
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables (N x K). Default is NULL}
\item{knn}{Number of nearest-neighbors being geographically weighted. Default is 100. Larger knn is better for larger samples (see Murakami er al., 2019)}
\item{kernel}{Kernel to model spatial heterogeneity. Gaussian kernel ("gau") and exponential kernel ("exp") are available}
\item{p}{Degree of the polynomial to approximate the kernel function. Default is 4}
\item{approach}{If "CV", leave-one-out cross-validation is used for the model calibration. If "AICc", the corrected Akaike Information Criterion is minimized for the calibation. Default is "CV"}
\item{nsamp}{Number of samples used to approximate the cross-validation. The samples are randomly selected. If the value is large enough (e.g., 10,000), error due to the sampling is quite small owing to the central limit theorem. The value must be smaller than the sample size. Default is NULL}
\item{cl}{Number of cores used for the parallel computation. If cl = NULL, which is the default, the number of available cores is detected and used}
}
\value{
\item{b}{Matrix of estimated coefficients (N x K)}
\item{bse}{Matrix of the standard errors for the coefficients (N x k)}
\item{t}{Matrix of the t-values for the coefficients (N x K)}
\item{p}{Matrix of the p-values for the coefficients (N x K)}
\item{par}{Estimated model parameters includeing a scale parameter and a shrinkage parameter if penalty = TRUE (see Murakami et al., 2018)}
\item{e}{Error statistics. It includes sum of squared errors (SSE), residual standard error (resid_SE), R-squared (R2), adjusted R2 (adjR2), log-likelihood (logLik), corrected Akaike information criterion (AICc), and the cross-validation (CV) score measured by root mean squared error (RMSE) (CV_score(RMSE))}
\item{pred}{Vector of predicted values (N x 1)}
\item{resid}{Vector of residuals (N x 1)}
\item{other}{Other objects internally used}
}
\references{
Murakami, D., Tsutsumida, N., Yoshida, T., Nakaya, T., and Lu, B. (2019) Scalable GWR: A linear-time algorithm for large-scale geographically weighted regression with polynomial kernels. <arXiv:1905.00266>.
}
\seealso{
\code{\link{scgwr}}, \code{\link{predict0}}
}
\examples{
# require(spData);require(sp)
# data(house)
# dat   <- data.frame(coordinates(house), house@data[,c("price","age","rooms","beds","syear")])
# coords<- dat[ ,c("long","lat")]
# y	    <- log(dat[,"price"])
# x     <- dat[,c("age","rooms","beds","syear")]

# Parallel estimation
# res1  <- scgwr_p( coords = coords, y = y, x = x )
# res1

# Parallel estimation + Approximate cross-validation using 10000 samples
# res2  <- scgwr_p( coords = coords, y = y, x = x, nsamp = 10000 )
# res2
}



