####### SPCR for Poisson regression (non-adaptive)
SPCRPoi <- function(x, y, k, xi, w, A, gamma0, gamma, Beta, lambda_gamma, lambda_beta){
	
	PARA_old <- c(gamma0, gamma, c(Beta))
	PARA_new <- PARA_old + 10

	while( max(abs(PARA_new-PARA_old)[PARA_new-PARA_old != 0]) > 1e-4 )
	{
		para_old <- c(gamma0, gamma, c(Beta))
		para_new <- para_old + 10
		
		eta <- exp( gamma0 + x%*%Beta%*%gamma )
		z <- gamma0 + x%*%Beta%*%gamma + ( y - eta )/eta
		
		while( max(abs(para_new-para_old)[para_new-para_old != 0]) > 1e-4 )
		{
			
			y_star <- x %*% A	
			
			### Estimate Beta
			for( l in 1:ncol(x) )
			{
				for( j in 1:k )
				{			
					Z <- z - gamma0 - x[,-l] %*% Beta[-l, ] %*% gamma  - x[, l] * sum( gamma[-j]*Beta[l,-j] )	# by fujisawa ; largely modified			
					Y_star_j <- y_star[, j] - x[, -l] %*% Beta[-l, j]	# by fujisawa			
					s <- sum( x[, l] * ( eta*Z*gamma[j] + 2*w*Y_star_j ) )	# by fujisawa			
					Beta[ l, j ] <- softsh( s, lambda_beta*(1-xi) )/( gamma[j]^2*sum( eta*x[ ,l]^2 ) + 2*w*sum( x[ ,l]^2 ) + 2*lambda_beta*xi )
				}
			}
			
			### Estimate gamma
	#		print("estimates of gamma")
		#	W <- sqrt( diag(var( x %*% Beta )) )  ### weights for adaptive lasso
			x_star <- x %*% Beta	# by fujisawa		
			for( l in 1:k )
			{
		#		x_star = x %*% Beta	# by fujisawa		
				z_star_2 <- z - gamma0 - as.matrix(x_star[, -l]) %*% as.matrix(gamma[-l])	# by fujisawa		
				s = sum( eta * x_star[, l] * z_star_2 )	# by fujisawa
				gamma[ l ] <- softsh( s, lambda_gamma)/(sum(eta*x_star[ ,l]^2))   ### ordinal lasso
		#		gamma[ l ] <- softsh( s, 0.5*lambda_gamma/W[l])/(sum(x_star[ ,l]^2))   ### adaptive lasso
				if( gamma[ l ] == "NaN" ) gamma[ l ] <- 0		
			}
			
			### Estimate gamma0
            ws = eta * (z - x  %*% Beta %*% gamma)		# by fujisawa
            gamma0 = sum(ws)/sum(eta)		# by fujisawa
			
			### Estimate A
			SVD <- svd( ( t(x) %*% x ) %*% Beta )
			A <- SVD$u %*% t(SVD$v)
			
			para_old <- para_new
			para_new <- c(gamma0, gamma, c(Beta))
			if( mean(abs(para_new-para_old)) == 0 ) break
			
		}
		
		PARA_old <- PARA_new
		PARA_new <- c(gamma0, gamma, c(Beta))
		if( mean(abs(PARA_new-PARA_old)) == 0 ) break
				
	}
	list( gamma0=gamma0, gamma=gamma, Beta=Beta, A=A )
}

####### SPCR for Logistic regression (adaptive)
adaSPCRPoi <- function(x, y, k, q=1, xi, w, A, gamma0, gamma, Beta, lambda_gamma, lambda_beta, BetaWeight){
	
	PARA_old <- c(gamma0, gamma, c(Beta))
	PARA_new <- PARA_old + 10
	
	while( max(abs(PARA_new-PARA_old)[PARA_new-PARA_old != 0]) > 1e-4 )
	{
		para_old <- c(gamma0, gamma, c(Beta))
		para_new <- para_old + 10
		
        eta <- exp( gamma0 + x%*%Beta%*%gamma )
        z <- gamma0 + x%*%Beta%*%gamma + ( y - eta )/eta
        
		while( max(abs(para_new-para_old)[para_new-para_old != 0]) > 1e-4 )
		{
			
			y_star <- x %*% A	
			
			### Estimate Beta
			for( l in 1:ncol(x) )
			{
				for( j in 1:k )
				{			
					Z <- z - gamma0 - x[,-l] %*% Beta[-l, ] %*% gamma  - x[, l] * sum( gamma[-j]*Beta[l,-j] )	# by fujisawa ; largely modified			
					Y_star_j <- y_star[, j] - x[, -l] %*% Beta[-l, j]	# by fujisawa			
					s <- sum( x[, l] * ( eta*Z*gamma[j] + 2*w*Y_star_j ) )	# by fujisawa			
					Beta[ l, j ] <- softsh( s, lambda_beta*(1-xi)/( abs(BetaWeight[ l, j ])^q + 1e-7 ) )/( gamma[j]^2*sum( eta*x[ ,l]^2 ) + 2*w*sum( x[ ,l]^2 ) + 2*lambda_beta*xi )
				}
			}
			
			### Estimate gamma
	#		print("estimates of gamma")
		#	W <- sqrt( diag(var( x %*% Beta )) )  ### weights for adaptive lasso
			x_star <- x %*% Beta	# by fujisawa		
			for( l in 1:k )
			{
		#		x_star = x %*% Beta	# by fujisawa		
				z_star_2 <- z - gamma0 - as.matrix(x_star[, -l]) %*% as.matrix(gamma[-l])	# by fujisawa		
				s = sum( eta * x_star[, l] * z_star_2 )	# by fujisawa
				gamma[ l ] <- softsh( s, lambda_gamma)/(sum(eta*x_star[ ,l]^2))   ### ordinal lasso
		#		gamma[ l ] <- softsh( s, 0.5*lambda_gamma/W[l])/(sum(x_star[ ,l]^2))   ### adaptive lasso
				if( gamma[ l ] == "NaN" ) gamma[ l ] <- 0		
			}
			
			### Estimate gamma0
            ws = eta * (z - x  %*% Beta %*% gamma)		# by fujisawa
            gamma0 = sum(ws)/sum(eta)		# by fujisawa
			
			### Estimate A
			SVD <- svd( ( t(x) %*% x ) %*% Beta )
			A <- SVD$u %*% t(SVD$v)
			
			para_old <- para_new
			para_new <- c(gamma0, gamma, c(Beta))
			if( mean(abs(para_new-para_old)) == 0 ) break
			
		}
		
		PARA_old <- PARA_new
		PARA_new <- c(gamma0, gamma, c(Beta))
		if( mean(abs(PARA_new-PARA_old)) == 0 ) break
		
	}
	list( gamma0=gamma0, gamma=gamma, Beta=Beta, A=A )
}

####### Cross-Validation for SPCRLoG
CV.SPCRPoi <- function(x, y, k, xi, w, nfolds=5, lambda.beta.candidate, lambda.gamma.candidate, center=TRUE, scale=FALSE, adaptive=FALSE, q=1, lambda.beta.seq=NULL, lambda.gamma.seq=NULL){

    n <- nrow(x)
    
    ####### Initialization of parameters (A, gamma0, gamma, Beta)
    A.ini <- as.matrix(eigen(var(x))$vectors[ ,1:k])
    gamma0.ini <- mean(y)
    gamma.ini <- rep(0, k)
    Beta.ini <- matrix( 0, nrow(A.ini), k )

    ### CV_mat : estimated CV errors
    CV.mat <- matrix( 0, length(lambda.gamma.candidate), length(lambda.beta.candidate) )
	
	foldid <- sample(rep(seq(nfolds),length=n))
	x.all <- x
	y.all <- y
	
	for(i in seq(nfolds))
	{
		num.foldid <- which(foldid==i)
		x <- x.all[ -num.foldid, ]
		y <- y.all[ -num.foldid ]
		x.test.cv <- x.all[ num.foldid, ]
		y.test.cv <- y.all[ num.foldid ]
		
	if( center==TRUE ){
		x_ori  <- x
		x <- sweep(x_ori, 2, apply(x_ori,2,mean))
		x.test.cv <- sweep(x.test.cv, 2, apply(x_ori,2,mean))
	}
	if( scale==TRUE ){
		x_ori  <- x
		x <- scale(x_ori)
		x.test.cv <- sweep(sweep(x.test.cv, 2, apply(x_ori, 2, mean)), 2, apply(x_ori, 2, sd), FUN="/")
	}
		
		####### START Estimate parameters (gamma_0, gamma, A, Beta)
		for( itr.lambda.gamma in 1:length(lambda.gamma.candidate) )
		{
			lambda.gamma <- lambda.gamma.candidate[itr.lambda.gamma]
			A <- A.ini
			gamma0 <- gamma0.ini
			gamma <- gamma.ini
			Beta <- Beta.ini
			
			for( itr.lambda.beta in 1:length(lambda.beta.candidate) )
			{
				lambda.beta <- lambda.beta.candidate[itr.lambda.beta]
				
				if( adaptive==FALSE ){
#					spcr.object = myfunc(x, y, A, Beta, gamma, gamma0, lambda.beta, lambda.gamma, xi, w)
#					spcr.object <- .Call( "spcr", x, y, A, Beta, gamma, gamma0, lambda.beta, lambda.gamma, xi, w )
					spcr.object <- SPCRPoi( x=x, y=y, A=A, k=k, gamma0=gamma0, gamma=gamma, Beta=Beta, lambda_beta=lambda.beta, lambda_gamma=lambda.gamma, xi=xi, w=w )
					Beta <- spcr.object$Beta
					gamma <- spcr.object$gamma
					gamma0 <- spcr.object$gamma0
					A <- spcr.object$A
				} else {
					spcr.object <- SPCRPoi( x=x, y=y, A=A, k=k, gamma0=gamma0, gamma=gamma, Beta=Beta, lambda_beta=lambda.beta, lambda_gamma=lambda.gamma, xi=xi, w=w )
					Beta <- spcr.object$Beta
					gamma <- spcr.object$gamma
					gamma0 <- spcr.object$gamma0
					A <- spcr.object$A
#					BetaWeight <- Beta
					if( sum(abs(Beta))==0 ) BetaWeight <- Beta
					if( sum(abs(Beta))!=0 ) BetaWeight <- Beta/sum(abs(Beta))
					adaspcr.object <- adaSPCRPoi( x=x, y=y, A=A, k=k, q=q, gamma0=gamma0, gamma=gamma, Beta=Beta, lambda_beta=lambda.beta, lambda_gamma=lambda.gamma, xi=xi, w=w, BetaWeight=BetaWeight )
					Beta <- adaspcr.object$Beta
					gamma <- adaspcr.object$gamma
					gamma0 <- adaspcr.object$gamma0
					A <- adaspcr.object$A
				}
				
				### CV-error
				s_cv <- ( gamma0*sum(y.test.cv) + t(y.test.cv)%*%x.test.cv%*%Beta%*%gamma - sum( exp(gamma0 + x.test.cv%*%Beta%*%gamma)) - sum( log( factorial(y.test.cv) ) ) )/length(y.test.cv)
				
				### Strock of CV-error
				CV.mat[ itr.lambda.gamma, itr.lambda.beta ] <- CV.mat[ itr.lambda.gamma, itr.lambda.beta ] + s_cv
			}
		}
	}
	CV.mat <- CV.mat/nfolds
	
	
	### START Search of max CV
	maxCandi.col <- whichimaxCandi.col <- rep(0, nrow(CV.mat))
	for(i in 1:nrow(CV.mat))
	{
		whichimaxCandi.col[i] <- which.max(CV.mat[i, ])
		maxCandi.col[i] <- max(CV.mat[i, ])
	}
	whichimaxCandi.row <- which.max( CV.mat[ , whichimaxCandi.col[ which.max(maxCandi.col) ]] )
	maxCandi.row <- max( CV.mat[ , whichimaxCandi.col[ which.max(maxCandi.col) ]] )		
	### END Search of max CV
	
	list(CV.mat = CV.mat, lambda.beta.candidate = lambda.beta.candidate, lambda.gamma.candidate = lambda.gamma.candidate, lambda.gamma.cv = lambda.gamma.candidate[ whichimaxCandi.row ], lambda.beta.cv = lambda.beta.candidate[ whichimaxCandi.col[ which.max(maxCandi.col) ] ], cvm=max(maxCandi.row))
}
