cv.spcr <- function(x, y, k, w=0.1, xi=0.01, nfolds=5, adaptive=FALSE, lambda.B.length=10, lambda.gamma.length=10){
	if( !is.matrix(x) ) stop("x must be a matrix.")
	if( mode(x)!="numeric" ) stop("x must be numeric.")
	if ( !is.vector(y) ) stop("y must be a vector.")
	if( mode(y)!="numeric" ) stop("y must be numeric.")
	if( k < 1 ) stop("k is an integer and larger than one.")
	
	n <- length(y)
	
	ini.lambda.beta <- ini.lambda.gamma <- ini.lambda( x=x, y=y, k=k, w=w, xi=xi )
	lambda.beta.candidate <- rev( seq( n*0.005, ini.lambda.beta, length=lambda.B.length ) )
	lambda.gamma.candidate <- rev( seq(n* 0.005, ini.lambda.gamma, length=lambda.gamma.length ) )
	
	A.ini <- as.matrix(eigen(var(x))$vectors[ ,1:k])
	gamma0.ini <- mean(y)
	gamma.ini <- rep(0, k)
	Beta.ini <- matrix( 0, nrow(A.ini), k )
	
	### CV_mat : estimated CV errors
	CV.mat <- matrix( 0, lambda.gamma.length, lambda.B.length )

	foldid <- sample(rep(seq(nfolds),length=n))
	x.all <- x
	y.all <- y
	
	for(i in seq(nfolds))
	{
		x <- x.all[ -(foldid==i), ]
		y <- y.all[ -(foldid==i) ]
		x.test.cv <- x.all[ (foldid==i), ]
		y.test.cv <- y.all[ (foldid==i) ]
		####### START Estimate parameters (gamma_0, gamma, A, Beta)
		for( itr.lambda.gamma in 1:lambda.gamma.length )
		{
			lambda.gamma <- lambda.gamma.candidate[itr.lambda.gamma]
			A <- A.ini
			gamma0 <- gamma0.ini
			gamma <- gamma.ini
			Beta <- Beta.ini
			
			for( itr.lambda.beta in 1:lambda.B.length )
			{
				lambda.beta <- lambda.beta.candidate[itr.lambda.beta]
				
				para_old <- c(gamma0, gamma, matrix(Beta, 1, nrow(Beta)*ncol(Beta)))
				para_new <- para_old + 10
				
				if( adaptive==FALSE ){
#					spcr.object = myfunc(x, y, A, Beta, gamma, gamma0, lambda.beta, lambda.gamma, xi, w)
					spcr.object <- .Call( "spcr", x, y, A, Beta, gamma, gamma0, lambda.beta, lambda.gamma, xi, w )
					Beta <- spcr.object[[1]]
					gamma <- spcr.object[[2]]
					gamma0 <- spcr.object[[3]]
					A <- spcr.object[[4]]
				} else {
					spcr.object <- .Call( "spcr", x, y, A, Beta, gamma, gamma0, lambda.beta, lambda.gamma, xi, w )
					Beta <- spcr.object[[1]]
					gamma <- spcr.object[[2]]
					gamma0 <- spcr.object[[3]]
					A <- spcr.object[[4]]
					BetaWeight <- Beta/sum(abs(Beta))		
					adaspcr.object <- .Call( "adaspcr", x, y, A, Beta, gamma, gamma0, lambda.beta, lambda.gamma, xi, w , BetaWeight)
					Beta <- adaspcr.object[[1]]
					gamma <- adaspcr.object[[2]]
					gamma0 <- adaspcr.object[[3]]
					A <- adaspcr.object[[4]]
				}
				
				para_old <- para_new
				para_new <- c(gamma0, gamma, c(Beta))
				if( mean(abs(para_new-para_old)) == 0 ) break
				
				### CV-error
				s_cv <- mean( ( y.test.cv - gamma0 - t(gamma) %*% t(Beta) %*% t(x.test.cv) )^2 )
				
				### Strock of CV-error
				CV.mat[ itr.lambda.gamma, itr.lambda.beta ] <- CV.mat[ itr.lambda.gamma, itr.lambda.beta ] + s_cv/nrow(x.test.cv)				
			}
		}
		CV.mat <- CV.mat/nfolds
	}
	
	### START Search of min CV
	minCandi.col <- whichiminCandi.col <- rep(0, nrow(CV.mat))
	for(i in 1:nrow(CV.mat))
	{
		whichiminCandi.col[i] <- which.min(CV.mat[i, ])
		minCandi.col[i] <- min(CV.mat[i, ])
	}
	
	minCandi.row <- whichiminCandi.row <- rep(0, ncol(CV.mat))
	for(i in 1:ncol(CV.mat))
	{
		whichiminCandi.row[i] <- which.min(CV.mat[ ,i])
		minCandi.row[i] <- min(CV.mat[ ,i])
	}
	### END Search of min CV
	
	### Selected tuning parameters by CV
	lambda.gamma.cv <- lambda.gamma.candidate[ whichiminCandi.row[ which.min(minCandi.row) ] ]
	lambda.beta.cv <- lambda.beta.candidate[ whichiminCandi.col[ which.min(minCandi.col) ] ]
	
	ans <- list( lambda.gamma.seq=lambda.gamma.candidate, lambda.B.seq=lambda.beta.candidate, CV.mat=CV.mat, lambda.gamma.cv=lambda.gamma.cv, lambda.B.cv=lambda.beta.cv, cvm=min(minCandi.row), call=match.call() )
	class(ans) <- "cv.spcr"
	ans		
}
