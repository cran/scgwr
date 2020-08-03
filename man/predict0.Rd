\name{predict0}
\alias{predict0}
\title{Spatial prediction using the scalable GWR model}
\usage{
predict0( mod, coords0, x0 = NULL )
}
\description{
This function predicts explained variables and spatially varying coefficients at unobserved sites using the scalable GWR model.
}

\arguments{
\item{mod}{Output from the scgwr function}
\item{coords0}{Matrix of spatial point coordinates at predicted sites (N0 x 2)}
\item{x0}{Matrix of explanatory variables at predicted sites (N0 x K). If NULL, explained variables are not predicted (only spatially varying coefficients are predicted). Default is NULL}
}
\value{
\item{pred}{Vector of predicted values (N0 x 1)}
\item{b}{Matrix of estimated coefficients (N0 x K)}
\item{bse}{Matrix of the standard errors for the coefficients (N0 x k)}
\item{t}{Matrix of the t-values for the coefficients (N0 x K)}
\item{p}{Matrix of the p-values for the coefficients (N0 x K)}
}

\examples{
require(spData)
data(boston)

id_obs  <-sample(dim(boston.c)[1],400)

######################### data at observed sites
y       <- log(boston.c[id_obs,"MEDV"])
x       <- boston.c[id_obs, c("CRIM", "INDUS","ZN","NOX","AGE")]
coords  <- boston.c[id_obs , c("LON", "LAT") ]

######################### data at predicted sites
x0      <- boston.c[-id_obs, c("CRIM", "INDUS","ZN","NOX", "AGE")]
coords0 <- boston.c[-id_obs , c("LON", "LAT") ]

mod     <- scgwr( coords = coords, y = y, x = x )
pred0   <- predict0( mod=mod, coords0=coords0, x0=x0)

pred    <- pred0$pred # predicted value
b       <- pred0$b    # spatially varying coefficients
b[1:5,]

bse     <- pred0$bse  # standard error of the coefficients
bt      <- pred0$t    # t-values
bp      <- pred0$p    # p-values

}
