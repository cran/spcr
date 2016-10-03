cv.spcrglm <- function(x, y, k, family=c("binomial","poisson","multinomial"), w=0.1, xi=0.01, nfolds=5, adaptive=FALSE, q=1, center=TRUE, scale=FALSE, lambda.B.length=10, lambda.gamma.length=10, lambda.B=NULL, lambda.gamma=NULL){
	if( !is.matrix(x) ) stop("x must be a matrix.")
	if( mode(x)!="numeric" ) stop("x must be numeric.")
	if ( !is.vector(y) ) stop("y must be a vector.")
	if( mode(y)!="numeric" ) stop("y must be numeric.")
	if( k < 1 ) stop("k is an integer and larger than one.")
	
    n <- nrow(x)
    
    ini.lambda.beta <- ini.lambda.gamma <- ini.lambda( x=x, y=y, k=k, w=w, xi=xi )
    lambda.beta.candidate <- rev( seq( n*0.005, ini.lambda.beta, length=lambda.B.length ) )
    lambda.gamma.candidate <- rev( seq(n* 0.005, ini.lambda.gamma, length=lambda.gamma.length ) )
    if( is.null(lambda.B) != TRUE ) lambda.beta.candidate <- sort(lambda.B, decreasing=TRUE)
    if( is.null(lambda.gamma) != TRUE ) lambda.gamma.candidate <- sort(lambda.gamma, decreasing=TRUE)
    
    ### binomial
    if( family=="binomial" ){
        if( adaptive==FALSE ){
            spcr.object <- CV.SPCRLoG( x, y, k, xi, w, nfolds=nfolds, lambda.beta.candidate=lambda.beta.candidate, lambda.gamma.candidate=lambda.gamma.candidate,  center=center, scale=scale, adaptive=adaptive, q=q )
        }
        if( adaptive==TRUE ){
            spcr.object <- CV.SPCRLoG( x, y, k, xi, w, nfolds=nfolds, lambda.beta.candidate=lambda.beta.candidate, lambda.gamma.candidate=lambda.gamma.candidate, center=center, scale=scale, adaptive=adaptive, q=q )
        }
    }
	
    ### poisson
    if( family=="poisson" ){
        if( adaptive==FALSE ){
            spcr.object <- CV.SPCRPoi( x, y, k, xi, w, nfolds=nfolds, lambda.beta.candidate=lambda.beta.candidate, lambda.gamma.candidate=lambda.gamma.candidate, center=center, scale=scale, adaptive=adaptive, q=q )
        }
        if( adaptive==TRUE ){
            spcr.object <- CV.SPCRPoi( x, y, k, xi, w, nfolds=nfolds, lambda.beta.candidate=lambda.beta.candidate, lambda.gamma.candidate=lambda.gamma.candidate, center=center, scale=scale, adaptive=adaptive, q=q )
        }
    }
    
    ### multinomial
    if( family=="multinomial" ){

        Y = y
        unique.Y = unique(Y)
        y = matrix(0, nrow(x), length(unique.Y))
        for(i in 1:nrow(x)) for(j in 1:length(unique.Y)) if( Y[i] == unique.Y[j] ) y[i,j] <- 1

        if( adaptive==FALSE ){
            spcr.object <- CV.SPCRMultiLoG( x, y, k, xi, w, nfolds=nfolds, lambda.beta.candidate=lambda.beta.candidate, lambda.gamma.candidate=lambda.gamma.candidate, center=center, scale=scale, adaptive=adaptive, q=q )
        }
        if( adaptive==TRUE ){
            spcr.object <- CV.SPCRMultiLoG( x, y, k, xi, w, nfolds=nfolds, lambda.beta.candidate=lambda.beta.candidate, lambda.gamma.candidate=lambda.gamma.candidate, center=center, scale=scale, adaptive=adaptive, q=q )
        }
    }
  
    ans <- list( lambda.gamma.seq = spcr.object$lambda.gamma.candidate, lambda.B.seq = spcr.object$lambda.beta.candidate, CV.mat = - spcr.object$CV.mat, lambda.gamma.cv = spcr.object$lambda.gamma.cv, lambda.B.cv = spcr.object$lambda.beta.cv, cvm = - spcr.object$cvm, call=match.call() )
    class(ans) <- "cv.spcrglm"
    return( ans )
    
}
