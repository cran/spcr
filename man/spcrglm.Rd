\name{spcrglm}
\alias{spcrglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit a sparse principal component regression for generalized linear models (SPCR-glm)}
\description{This function computes a principal component regression for generalized linear models via sparse regularization.}
\usage{spcrglm(x, y, k, family=c("binomial","poisson","multinomial"), lambda.B, 
	lambda.gamma, w=0.1, xi=0.01, adaptive=FALSE, q=1, center=TRUE, scale=FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data matrix.
}
  \item{y}{
A response data.
}
  \item{k}{
The number of principal components.
}
  \item{family}{
Response type.
}
  \item{lambda.B}{
The regularization parameter for the parameter \eqn{B}.
}
  \item{lambda.gamma}{
The regularization parameter for the coefficient vector \eqn{\gamma}.
}
  \item{w}{
Weight parameter with \eqn{w \ge 0}. The default is 0.1.
}
  \item{xi}{
The elastic net mixing parameter with \eqn{0\le \alpha \le 1}. The default is 0.01.
}
 \item{adaptive}{If \code{"TRUE"}, the adaptive SPCR-glm (aSPCR-glm) is used.}
  \item{q}{
The tuning parameter that controls weights in aSPCR-glm. The default is 1.
}
 \item{center}{If \code{"TRUE"}, the data matrix is centered.}
 \item{scale}{If \code{"TRUE"}, the data matrix is scaled.}
}
\value{
\item{loadings.B}{the loading matrix B}
\item{gamma}{the coefficient}
\item{gamma0}{intercept}
\item{loadings.A}{the loading matrix A}
}
\references{
%% ~put references to the literature/web site here ~
%Friedman, J. (2008).  Fast sparse regression and classification. \eqn{Technical report}, Standford University.\cr
Kawano, S., Fujisawa, H., Takada, T. and Shiroishi, T. (2016). 
\emph{Sparse principal component regression for generalized linear models}.
  \url{http://arxiv.org/abs/1609.08886}
}
\author{
Shuichi Kawano\cr
\email{skawano@uec.ac.jp}
}
%\note{
%  ~~further notes~~
%}
% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{\code{cv.spcrglm}}
\examples{
# binomial
n <- 100
np <- 5
nu0 <- c(-1, 1)
set.seed(4)
x <- matrix( rnorm(np*n), n, np )
y <- rbinom(n,1,1-1/(1+exp(  (nu0[1]*x[ ,1] + nu0[2]*x[ ,2]  ))))
spcrglm.fit <- spcrglm(x=x, y=y, k=2, family="binomial", lambda.B=2, lambda.gamma=1)
spcrglm.fit

# Poisson
set.seed(4)
y <- rpois(n, exp(  (nu0[1]*x[ ,1] + nu0[2]*x[ ,2] ) ))spcrglm.fit <- spcrglm(x=x, y=y, k=2, family="poisson", lambda.B=2, lambda.gamma=1)
spcrglm.fit

# multinomial
set.seed(4)
y <- sample(1:4, n, replace=TRUE)
spcrglm.fit <- spcrglm(x=x, y=y, k=2, family="multinomial", lambda.B=2, lambda.gamma=2)
spcrglm.fit
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
