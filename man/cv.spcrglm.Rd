\name{cv.spcrglm}
\alias{cv.spcrglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for spcr-glm}
\description{This function performs cross-validation for SPCR-glm. \code{cv.spcrglm} enables us to determine two regularization parameters \eqn{\lambda_\beta} and \eqn{\lambda_\gamma} objectively.}
\usage{cv.spcrglm(x, y, k, family=c("binomial","poisson","multinomial"), 
	w=0.1, xi=0.01, nfolds=5, adaptive=FALSE, q=1, center=TRUE, 
	scale=FALSE, lambda.B.length=10, lambda.gamma.length=10,
	lambda.B=NULL, lambda.gamma=NULL)}
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
  \item{family}{
Response type.
}
  \item{w}{
Weight parameter with \eqn{w \ge 0}. The default is 0.1.
}
  \item{xi}{
The elastic net mixing parameter with \eqn{0\le \alpha \le 1}. The default is 0.01.
}
  \item{nfolds}{
The number of folds. The default is 5.
}
 \item{adaptive}{If \code{"TRUE"}, the adaptive SPCR-glm (aSPCR-glm) is used.}
 \item{q}{The tuning parameter that controls weights in aSPCR-glm. The default is 1.}
  \item{center}{If \code{"TRUE"}, the data matrix is centered.}
 \item{scale}{If \code{"TRUE"}, the data matrix is scaled.}
 \item{lambda.B.length}{
The number of candidates for the parameter \eqn{\lambda_\beta}. The default is 10.
}
 \item{lambda.gamma.length}{
The number of candidates for the parameter \eqn{\lambda_\gamma}. The default is 10.
}
 \item{lambda.B}{
Optional user-supplied candidates for the parameter \eqn{\lambda_\beta}.  The default is NULL. 
}
 \item{lambda.gamma}{
Optional user-supplied candidates for the parameter \eqn{\lambda_\gamma}.  The default is NULL. 
}
}
\value{
\item{lambda.gamma.seq}{The values of \code{lambda.gamma} in the fit.}
\item{lambda.B.seq}{The values of \code{lambda.B} in the fit.}
\item{CV.mat}{Matrix of the mean values of cross-validation. The row shows a sequence of \code{lambda.gamma}. The column shows a sequence of \code{lambda.B}.}
\item{lambda.gamma.cv}{The value of \code{lambda.gamma} selected by cross-validation.}
\item{lambda.B.cv}{The value of \code{lambda.B} selected by cross-validation.}
\item{cvm}{The minimum of the mean cross-validated error.}
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
\seealso{\code{spcrglm}}
\examples{
# binomial
n <- 100
np <- 3
nu0 <- c(-1, 1)
set.seed(4)
x <- matrix( rnorm(np*n), n, np )
y <- rbinom(n,1,1-1/(1+exp(  (nu0[1]*x[ ,1] + nu0[2]*x[ ,2]  ))))
cv.spcrglm.fit <- cv.spcrglm(x=x, y=y, k=1, family="binomial")
cv.spcrglm.fit

# Poisson
set.seed(5)
y <- rpois(n, 1)
cv.spcrglm.fit <- cv.spcrglm(x=x, y=y, k=1, family="poisson")
cv.spcrglm.fit

# multinomial
set.seed(4)
y <- sample(1:4, n, replace=TRUE)
cv.spcrglm.fit <- cv.spcrglm(x=x, y=y, k=1, family="multinomial")
cv.spcrglm.fit
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
