####### SPCR for Logistic regression (non-adaptive)
SPCRMultiLoG <- function(x, y, k, xi, w, A, gamma0, gamma, Beta, lambda_gamma, lambda_beta){
	
	PARA_old <- c(gamma0, c(gamma), c(Beta))
	PARA_new <- PARA_old + 10

##	ITR <- 0
	
	while( max(abs(PARA_new-PARA_old)[PARA_new-PARA_old != 0]) > 1e-3 )
	{
		
		### Recentering
		gammaMed = apply(gamma, 1, median)
		gamma = gamma - gammaMed
		gamma0 = gamma0 - mean(gamma0)
		
		### Update Lq
		wa = rep(1, nrow(x))%*%t(gamma0) + x%*%Beta%*%gamma
		ww = apply( wa, 1, max )
		wb = wa - ww
		ptilde <- exp(wb) / apply(exp(wb), 1, sum)
		eta <- ptilde*( 1 - ptilde )
		etaz = eta * ( rep(1, nrow(x))%*%t(gamma0) + x%*%Beta%*%gamma ) + (y - ptilde)

        ### Estimate Beta
        y_star <- x %*% A	
        x_star <- x %*% Beta
        for( l in 1:ncol(x) )
        {
            for( j in 1:k )
            {			
                etaZ <- etaz - eta * ( rep(1, nrow(x))%*%t(gamma0) + x[,-l] %*% Beta[-l, ] %*% gamma + as.matrix(x[, l]) %*% as.matrix(t(Beta[l,-j])) %*% gamma[-j, ] )	# by fujisawa ; largely modified
                Y_star_j <- y_star[, j] - x[, -l] %*% Beta[-l, j]	# by fujisawa			
                s <- sum( x[, l] * ( (etaZ)%*%gamma[j, ] + 2*w*Y_star_j ) )	# by fujisawa			
                Beta[ l, j ] <- softsh( s, lambda_beta*(1-xi) )/( (gamma[j, ]^2)%*%t(eta)%*%(x[ ,l]^2) + 2*w*sum( x[ ,l]^2 ) + 2*lambda_beta*xi )
            }
        }
        
        ### Estimate A
        SVD <- svd( ( t(x) %*% x ) %*% Beta )
        A <- SVD$u %*% t(SVD$v)
        
        ### Estimate gamma & gamma0
        x_star <- x %*% Beta

        for(g in 1:ncol(y))
        {	
            ### Estimate gamma
            for( l in 1:k )
            {					
                etaz_star_2 <- etaz[ ,g] - eta[,g] * ( gamma0[g] + as.matrix(x_star[, -l]) %*% as.matrix(gamma[-l,g]) )	# by fujisawa
                s = sum( x_star[, l] * etaz_star_2 )	# by fujisawa
                ww = sum(eta[ ,g]*(x_star[ ,l]^2))
                if( ww  == 0 ){ 
                    gamma[ l, g ] = 0
                }else{
                    gamma[ l, g ] <- softsh( s, lambda_gamma)/(1e-7+sum(eta[ ,g]*(x_star[ ,l]^2)))   ### ordinal lasso
                }
            }
            
            ### Estimate gamma0
            ws =  etaz[ ,g] - eta[ ,g] * ( x %*% Beta %*% gamma[ ,g] )		# by fujisawa
            gamma0[g] = sum(ws)/sum(eta[ ,g])		# by fujisawa
        }

        PARA_old <- PARA_new
		PARA_new <- c(gamma0, c(gamma), c(Beta))
		if( mean(abs(PARA_new-PARA_old)) == 0 ) break
				
	}
	list( gamma0=gamma0, gamma=gamma, Beta=Beta, A=A )
}

####### SPCR for Logistic regression (adaptive)
adaSPCRMultiLoG <- function(x, y, k, q=1, xi, w, A, gamma0, gamma, Beta, lambda_gamma, lambda_beta, BetaWeight){
    
    PARA_old <- c(gamma0, c(gamma), c(Beta))
    PARA_new <- PARA_old + 10
    
    ##	ITR <- 0
    
    while( max(abs(PARA_new-PARA_old)[PARA_new-PARA_old != 0]) > 1e-3 )
    {
        
        ### Recentering
        gammaMed = apply(gamma, 1, median)
        gamma = gamma - gammaMed
        gamma0 = gamma0 - mean(gamma0)
        
        ### Update Lq
        wa = rep(1, nrow(x))%*%t(gamma0) + x%*%Beta%*%gamma
        ww = apply( wa, 1, max )
        wb = wa - ww
        ptilde <- exp(wb) / apply(exp(wb), 1, sum)
        eta <- ptilde*( 1 - ptilde )
        etaz = eta * ( rep(1, nrow(x))%*%t(gamma0) + x%*%Beta%*%gamma ) + (y - ptilde)
        
        ### Estimate Beta
        y_star <- x %*% A
        x_star <- x %*% Beta
        for( l in 1:ncol(x) )
        {
            for( j in 1:k )
            {
                etaZ <- etaz - eta * ( rep(1, nrow(x))%*%t(gamma0) + x[,-l] %*% Beta[-l, ] %*% gamma + as.matrix(x[, l]) %*% as.matrix(t(Beta[l,-j])) %*% gamma[-j, ] )	# by fujisawa ; largely modified
                Y_star_j <- y_star[, j] - x[, -l] %*% Beta[-l, j]	# by fujisawa
                s <- sum( x[, l] * ( (etaZ)%*%gamma[j, ] + 2*w*Y_star_j ) )	# by fujisawa
                Beta[ l, j ] <- softsh( s, lambda_beta*(1-xi)/( abs(BetaWeight[ l, j ])^q + 1e-7 ) )/( (gamma[j, ]^2)%*%t(eta)%*%(x[ ,l]^2) + 2*w*sum( x[ ,l]^2 ) + 2*lambda_beta*xi )
            }
        }
        
        ### Estimate A
        SVD <- svd( ( t(x) %*% x ) %*% Beta )
        A <- SVD$u %*% t(SVD$v)
        
        ### Estimate gamma & gamma0
        x_star <- x %*% Beta
        
        for(g in 1:ncol(y))
        {
            ### Estimate gamma
            for( l in 1:k )
            {
                etaz_star_2 <- etaz[ ,g] - eta[,g] * ( gamma0[g] + as.matrix(x_star[, -l]) %*% as.matrix(gamma[-l,g]) )	# by fujisawa
                s = sum( x_star[, l] * etaz_star_2 )	# by fujisawa
                ww = sum(eta[ ,g]*(x_star[ ,l]^2))
                if( ww  == 0 ){
                    gamma[ l, g ] = 0
                }else{
                    gamma[ l, g ] <- softsh( s, lambda_gamma)/(1e-7+sum(eta[ ,g]*(x_star[ ,l]^2)))   ### ordinal lasso
                }
            }
            
            ### Estimate gamma0
            ws =  etaz[ ,g] - eta[ ,g] * ( x %*% Beta %*% gamma[ ,g] )		# by fujisawa
            gamma0[g] = sum(ws)/sum(eta[ ,g])		# by fujisawa
        }
        
        PARA_old <- PARA_new
        PARA_new <- c(gamma0, c(gamma), c(Beta))
        if( mean(abs(PARA_new-PARA_old)) == 0 ) break
        
    }
    list( gamma0=gamma0, gamma=gamma, Beta=Beta, A=A )
}

####### Cross-Validation for SPCRMultiLoG
CV.SPCRMultiLoG <- function(x, y, k, xi, w, nfolds=5, lambda.beta.candidate, lambda.gamma.candidate, center=TRUE, scale=FALSE, adaptive=FALSE, q=1){
    
    n <- nrow(x)
    
    ####### Initialization of parameters (A, gamma0, gamma, Beta)
    A.ini <- as.matrix(eigen(var(x))$vectors[ ,1:k])
    gamma0.ini <- apply(y, 2, mean)
    gamma.ini <- matrix(0, k, ncol(y))
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
        y <- y.all[ -num.foldid,  ]
        x.test.cv <- x.all[ num.foldid, ]
        y.test.cv <- y.all[ num.foldid,  ]
        
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
                    spcr.object <- SPCRMultiLoG( x=x, y=y, k=k, xi=xi, w=w, A=A, gamma0=gamma0, gamma=gamma, Beta=Beta, lambda_gamma=lambda.gamma, lambda_beta=lambda.beta )
                    Beta <- spcr.object$Beta
                    gamma <- spcr.object$gamma
                    gamma0 <- spcr.object$gamma0
                    A <- spcr.object$A
                } else {
                    spcr.object <- SPCRMultiLoG( x=x, y=y, k=k, xi=xi, w=w, A=A, gamma0=gamma0, gamma=gamma, Beta=Beta, lambda_gamma=lambda.gamma, lambda_beta=lambda.beta )
                    Beta <- spcr.object$Beta
                    gamma <- spcr.object$gamma
                    gamma0 <- spcr.object$gamma0
                    A <- spcr.object$A
                    if( sum(abs(Beta))==0 ) BetaWeight <- Beta
                    if( sum(abs(Beta))!=0 ) BetaWeight <- Beta/sum(abs(Beta))
                    adaspcr.object <- adaSPCRMultiLoG( x=x, y=y, A=A, k=k, q=q, gamma0=gamma0, gamma=gamma, Beta=Beta, lambda_beta=lambda.beta, lambda_gamma=lambda.gamma, xi=xi, w=w, BetaWeight=BetaWeight )
                    Beta <- adaspcr.object$Beta
                    gamma <- adaspcr.object$gamma
                    gamma0 <- adaspcr.object$gamma0
                    A <- adaspcr.object$A
                }
                
                ### CV-error
                s_cv <- ( sum( y.test.cv %*% gamma0 ) + sum( y.test.cv * ( x.test.cv %*% Beta %*% gamma ) ) - sum( log( 1 + exp( gamma0 + x.test.cv %*% Beta %*% gamma ) ) ) ) / nrow(y.test.cv)
                
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
