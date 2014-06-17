\name{spcr}
\alias{spcr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit a sparse principal component regression (SPCR)}
\description{This function computes a principal component regression model via sparse regularization.}
\usage{spcr(x, y, k, lambda.B, lambda.gamma, w=0.1, xi=0.01, 
	adaptive=FALSE, center=TRUE, scale=FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data matrix.
}
  \item{y}{
A response vector.
}
  \item{k}{
The number of principal components.
}
  \item{lambda.B}{
The regularization parameter for the parameter \eqn{B}.
}
  \item{lambda.gamma}{
The regularization parameter for the coefficient vector \eqn{\gamma}.
}
  \item{w}{
Weight parameter with \eqn{0\le w \le 1}. The default is 0.1.
}
  \item{xi}{
The elastic net mixing parameter with \eqn{0\le \alpha \le 1}. The default is 0.01.
}
 \item{adaptive}{If \code{"TRUE"}, the adaptive SPCR is used.}
 \item{center}{If \code{"TRUE"}, the data matrix is centered.}
 \item{scale}{If \code{"TRUE"}, the data matrix is scaled.}
}
\value{
\item{loadings.B}{B}
\item{gamma}{the coefficient}
\item{gamma0}{intercept}
\item{loadings.A}{A}
}
\references{
%% ~put references to the literature/web site here ~
%Friedman, J. (2008).  Fast sparse regression and classification. \eqn{Technical report}, Standford University.\cr
Kawano, S., Fujisawa, H., Takada, T. and Shiroishi, T. (2014). 
\emph{Sparse principal component regression with adaptive loading}.
  \url{http://arxiv.org/abs/1402.6455}
}
\author{
Shuichi Kawano\cr
\email{skawano@ms.osakafu-u.ac.jp}
}
%\note{
%  ~~further notes~~
%}
% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{\code{cv.spcr}}
\examples{
#data
n <- 100
np <- 5
set.seed(4)
nu0 <- c(-1, 1)
x <- matrix( rnorm(np*n), n, np )
e <- rnorm(n)
y <- nu0[1]*x[ ,1] + nu0[2]*x[ ,2] + e

#fit
spcr.fit <- spcr(x=x, y=y, k=2, lambda.B=6, lambda.gamma=2)
spcr.fit

#fit (adaptive SPCR)
adaspcr.fit <- spcr(x=x, y=y, k=2, lambda.B=6, lambda.gamma=2, adaptive=TRUE)
adaspcr.fit
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
