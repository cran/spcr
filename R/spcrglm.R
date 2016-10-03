spcrglm <- function(x, y, k, family=c("binomial","poisson","multinomial"), lambda.B, lambda.gamma, w=0.1, xi=0.01, adaptive=FALSE, q=1, center=TRUE, scale=FALSE){
	if( !is.matrix(x) ) stop("x must be a matrix.")
	if( mode(x)!="numeric" ) stop("x must be numeric.")
	if ( !is.vector(y) ) stop("y must be a vector.")
	if( mode(y)!="numeric" ) stop("y must be numeric.")
	
	if( center==TRUE ) x <- sweep(x, 2, apply(x,2,mean))
	if( scale==TRUE ) x <- scale(x)
	
    ### binomial
    if( family=="binomial" ){
        if( adaptive==FALSE ){
            A <- as.matrix(eigen(var(x))$vectors[ ,1:k])
            gamma0 <- mean(y)
            gamma <- rep(0, k)
            Beta <- matrix( 0, nrow(A), k )
            spcr.object <- SPCRLoG( x, y, k, xi, w, A, gamma0, gamma, Beta, lambda.gamma, lambda.B )
            ans <- list( loadings.B=spcr.object$Beta, gamma=spcr.object$gamma, gamma0=spcr.object$gamma0, loadings.A=spcr.object$A, call=match.call() )
            class(ans) <- "spcrglm"
        }
        if( adaptive==TRUE ){
            A <- as.matrix(eigen(var(x))$vectors[ ,1:k])
            gamma0 <- mean(y)
            gamma <- rep(0, k)
            Beta <- matrix( 0, nrow(A), k )
            spcr.object <- SPCRLoG( x, y, k, xi, w, A, gamma0, gamma, Beta, lambda.gamma, lambda.B )
            Beta <- spcr.object$Beta
            gamma <- spcr.object$gamma
            gamma0 <- spcr.object$gamma0
            A <- spcr.object$A
            BetaWeight <- Beta/sum(abs(Beta))
            adaspcr.object <- adaSPCRLoG( x, y, k, q, xi, w, A, gamma0, gamma, Beta, lambda.gamma, lambda.B, BetaWeight )
            ans <- list( loadings.B=adaspcr.object$Beta, gamma=adaspcr.object$gamma, gamma0=adaspcr.object$gamma0, loadings.A=adaspcr.object$A, call=match.call()  )
            class(ans) <- "spcrglm"
        }
    }
    
    ### poisson
    if( family=="poisson" ){
        if( adaptive==FALSE ){
            A <- as.matrix(eigen(var(x))$vectors[ ,1:k])
            gamma0 <- mean(y)
            gamma <- rep(0, k)
            Beta <- matrix( 0, nrow(A), k )
            spcr.object <- SPCRPoi( x, y, k, xi, w, A, gamma0, gamma, Beta, lambda.gamma, lambda.B )
            ans <- list( loadings.B=spcr.object$Beta, gamma=spcr.object$gamma, gamma0=spcr.object$gamma0, loadings.A=spcr.object$A, call=match.call() )
            class(ans) <- "spcrglm"
        }
        if( adaptive==TRUE ){
            A <- as.matrix(eigen(var(x))$vectors[ ,1:k])
            gamma0 <- mean(y)
            gamma <- rep(0, k)
            Beta <- matrix( 0, nrow(A), k )
            spcr.object <- SPCRPoi( x, y, k, xi, w, A, gamma0, gamma, Beta, lambda.gamma, lambda.B )
            Beta <- spcr.object$Beta
            gamma <- spcr.object$gamma
            gamma0 <- spcr.object$gamma0
            A <- spcr.object$A
            BetaWeight <- Beta/sum(abs(Beta))
            adaspcr.object <- adaSPCRPoi( x, y, k, q, xi, w, A, gamma0, gamma, Beta, lambda.gamma, lambda.B, BetaWeight )
            ans <- list( loadings.B=adaspcr.object$Beta, gamma=adaspcr.object$gamma, gamma0=adaspcr.object$gamma0, loadings.A=adaspcr.object$A, call=match.call() )
            class(ans) <- "spcrglm"
        }
    }
    
    ### multinomial
    if( family=="multinomial" ){
        
        Y = y
        unique.Y = unique(Y)
        y = matrix(0, nrow(x), length(unique.Y))
        for(i in 1:nrow(x)) for(j in 1:length(unique.Y)) if( Y[i] == unique.Y[j] ) y[i,j] <- 1
        
        if( adaptive==FALSE ){
            A = as.matrix(eigen(var(x))$vectors[ ,1:k])
            gamma0 = apply(y, 2, mean)
            gamma = matrix(0, k, ncol(y))
            Beta = matrix(0, ncol(x), k)
            spcr.object <- SPCRMultiLoG( x, y, k, xi, w, A, gamma0, gamma, Beta, lambda.gamma, lambda.B )
            ans <- list( loadings.B=spcr.object$Beta, gamma=spcr.object$gamma, gamma0=spcr.object$gamma0, loadings.A=spcr.object$A, call=match.call() )
            class(ans) <- "spcrglm"
        }
        if( adaptive==TRUE ){
            A = as.matrix(eigen(var(x))$vectors[ ,1:k])
            gamma0 = apply(y, 2, mean)
            gamma = matrix(0, k, ncol(y))
            Beta = matrix(0, ncol(x), k)
            spcr.object <- SPCRMultiLoG( x, y, k, xi, w, A, gamma0, gamma, Beta, lambda.gamma, lambda.B )
            Beta <- spcr.object$Beta
            gamma <- spcr.object$gamma
            gamma0 <- spcr.object$gamma0
            A <- spcr.object$A
            BetaWeight <- Beta/sum(abs(Beta))
            adaspcr.object <- adaSPCRMultiLoG( x, y, k, q, xi, w, A, gamma0, gamma, Beta, lambda.gamma, lambda.B, BetaWeight )
            ans <- list( loadings.B=adaspcr.object$Beta, gamma=adaspcr.object$gamma, gamma0=adaspcr.object$gamma0, loadings.A=adaspcr.object$A, call=match.call() )
            class(ans) <- "spcrglm"
        }
    }
    
    return( ans )
    
}
