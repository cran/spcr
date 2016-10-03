//Ｒ用 Ｃ
//#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Parse.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>


//////////////////////////////////////////////////
//BLAS: y := alpha*A*x + beta*y
long int gl_incx = 1;
long int gl_incy = 1;
double alphaOne = 1.0;
double alphaminusOne = -1.0;
double betaOne = 1.0;
double betaZero = 0.0;
double Zero = 0.0;
char TRANSN = 'N';
char TRANST = 'T';
char UPPER = 'U';
char ATIL = 'A';


double softsh(double z, double eta)
{
	double sss;
	if( eta < fabs(z) ){
		if( z > 0.0 ) {
			sss = z - eta;
			return sss;
		} else {
			sss = z + eta;
			return sss;
		}
	} else {
		sss = 0.0;
		return sss;
	}
}


// spcr
SEXP spcr(SEXP ex_x, SEXP ex_y, SEXP ex_A, SEXP ex_Beta, SEXP ex_gamma, SEXP ex_gamma0, SEXP ex_lambda_beta, SEXP ex_lambda_gamma, SEXP ex_xi, SEXP ex_w)
{
	int n = INTEGER(GET_DIM(ex_x))[0];
	int p = INTEGER(GET_DIM(ex_x))[1];
	int k = INTEGER(GET_DIM(ex_A))[1];
	int pk = p*k;
	int pkkone = pk+k+1;
	int one=1;
	int i;
	int j;
	int l;
	int lwork;
	int info;
	
	double s;
	double s_sum;
	double wkopt;
	double* work;
	double diff_max;
	double residual_A;
	double residual_B;
	double Beta_old;
	double gamma_old;
	
    
	
	SEXP x_star;
	SEXP y_star;
	SEXP A;
	SEXP Beta;
	SEXP gamma;
	SEXP gamma0;
	SEXP x_Beta_ForGamma0;
	SEXP x_Beta_gamma_ForGamma0;
	SEXP ans;
	SEXP A_singular;
	SEXP u;	
	SEXP vt;
	SEXP iwork;
	SEXP tXX;
	SEXP tXX_Beta;
	SEXP para_old;
	SEXP para_new;
	SEXP para_diff;
	SEXP n_itr;
	SEXP X_inn;
	SEXP XinY;
	SEXP Xingamma0;
	SEXP XinY_star;
	SEXP unit_vec;
	SEXP BetainGamma;
	SEXP X_innBetainGamma;
	SEXP X_innBeta;
	SEXP X_star_inn;
	SEXP X_staringamma0;
	SEXP X_starinY;
	SEXP X_starinngamma;
	
	
	PROTECT(x_star = allocMatrix(REALSXP, n, k));
	PROTECT(y_star = allocMatrix(REALSXP, n, k));
	PROTECT(A = allocMatrix(REALSXP, p, k));
	PROTECT(Beta = allocMatrix(REALSXP, p, k));
	PROTECT(gamma = allocVector(REALSXP, k));
	PROTECT(gamma0 = allocVector(REALSXP, one));
	PROTECT(x_Beta_ForGamma0 = allocMatrix(REALSXP, n, k));
	PROTECT(x_Beta_gamma_ForGamma0 = allocVector(REALSXP, n));
	PROTECT(ans = allocVector(VECSXP, 5));
	PROTECT(A_singular = allocVector(REALSXP, p));
	PROTECT(u = allocMatrix(REALSXP, p, p));
	PROTECT(vt = allocMatrix(REALSXP, k, k));
	PROTECT(iwork = allocVector(INTSXP, 8*k));
	PROTECT(tXX = allocMatrix(REALSXP, p, p));
	PROTECT(tXX_Beta = allocMatrix(REALSXP, p, k));
	PROTECT(para_old = allocVector(REALSXP, pkkone));
	PROTECT(para_new = allocVector(REALSXP, pkkone));
	PROTECT(para_diff = allocVector(REALSXP, pkkone));
	PROTECT(n_itr = allocVector(INTSXP, one));
	PROTECT(X_inn = allocMatrix(REALSXP, p, p));
	PROTECT(XinY = allocVector(REALSXP, p));
	PROTECT(Xingamma0 = allocVector(REALSXP, p));
	PROTECT(XinY_star = allocMatrix(REALSXP, p, k));
	PROTECT(unit_vec = allocVector(REALSXP, n));
	PROTECT(BetainGamma = allocVector(REALSXP, p));
	PROTECT(X_innBetainGamma = allocVector(REALSXP, p));
	PROTECT(X_innBeta = allocMatrix(REALSXP, p, k));
	PROTECT(X_star_inn = allocMatrix(REALSXP, k, k));
	PROTECT(X_staringamma0 = allocVector(REALSXP, k));
	PROTECT(X_starinY = allocVector(REALSXP, k));
	PROTECT(X_starinngamma = allocVector(REALSXP, k));
	

	F77_CALL(dcopy)(&pk, REAL(ex_Beta), &one, REAL(Beta), &one);
	F77_CALL(dcopy)(&k, REAL(ex_gamma), &one, REAL(gamma), &one);
	F77_CALL(dcopy)(&one, REAL(ex_gamma0), &one, REAL(gamma0), &one);
	F77_CALL(dcopy)(&pk, REAL(ex_A), &one, REAL(A), &one);
	

	
	for (i=0; i<pk; i++) {
		REAL(para_new)[i] = REAL(Beta)[i];
	}
	for (i=pk; i<(pk+k); i++) {
		REAL(para_new)[i] = REAL(gamma)[(i-pk)];
	}
	REAL(para_new)[(pk+k)] = REAL(gamma0)[0];
	
	diff_max=1.0;	
	INTEGER(n_itr)[0] = 0;
	
	/// X_inn の作成
	F77_CALL(dgemm)(&TRANST, &TRANSN, &p, &p, &n, &alphaOne, REAL(ex_x), &n, REAL(ex_x), &n, &betaZero, REAL(X_inn), &p);
	
	/// XinY の作成
	F77_CALL(dgemv)(&TRANST, &n, &p, &alphaOne, REAL(ex_x), &n, REAL(ex_y), &one, &betaZero, REAL(XinY), &one);

	// パラメータの更新スタート
	while ( diff_max > 0.001 ) {
		F77_CALL(dcopy)(&pkkone, REAL(para_new), &one, REAL(para_old), &one);				
		
		/////////////////////////
		///// Estimate Beta /////
		/////////////////////////
		
		/// y_star の作成
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &n, &k, &p, &alphaOne, REAL(ex_x), &n, REAL(A), &p, &betaZero, REAL(y_star), &n);
		
		/// Xingamma0 の作成
		for (i=0; i<n; i++) {
			REAL(unit_vec)[i] = REAL(gamma0)[0];
		}
		F77_CALL(dgemv)(&TRANST, &n, &p, &alphaOne, REAL(ex_x), &n, REAL(unit_vec), &one, &betaZero, REAL(Xingamma0), &one);	
		
		/// XinY_star の作成
		F77_CALL(dgemm)(&TRANST, &TRANSN, &p, &k, &n, &alphaOne, REAL(ex_x), &n, REAL(y_star), &n, &betaZero, REAL(XinY_star), &p);
		
		/// X_innBetainGamma の作成
		F77_CALL(dgemv)(&TRANSN, &p, &k, &alphaOne, REAL(Beta), &p, REAL(gamma), &one, &betaZero, REAL(BetainGamma), &one);
		F77_CALL(dgemv)(&TRANSN, &p, &p, &alphaOne, REAL(X_inn), &p, REAL(BetainGamma), &one, &betaZero, REAL(X_innBetainGamma), &one);	
		
		/// X_innBeta の作成
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &k, &p, &alphaOne, REAL(X_inn), &p, REAL(Beta), &p, &betaZero, REAL(X_innBeta), &p);
		
		/// パラメータの更新 via Covariance Update
		for (j=0; j<k; j++) {	
			for (l=0; l<p; l++) {
				
				residual_A = REAL(XinY)[l] - REAL(Xingamma0)[l] - REAL(X_innBetainGamma)[l];
				residual_B = REAL(XinY_star)[l+p*j] - REAL(X_innBeta)[l+p*j];
				s = (1.0-REAL(ex_w)[0])*REAL(gamma)[j]*residual_A + REAL(ex_w)[0]*residual_B + REAL(Beta)[l+p*j]*REAL(X_inn)[l+p*l]*(REAL(ex_w)[0] + (1.0-REAL(ex_w)[0])*(pow(REAL(gamma)[j], 2.0)));
				Beta_old = REAL(Beta)[l+p*j];
				REAL(Beta)[l+p*j] = softsh(s, 0.5*REAL(ex_lambda_beta)[0]*(1.0-REAL(ex_xi)[0]))/(((1.0-REAL(ex_w)[0])*pow(REAL(gamma)[j], 2.0)+REAL(ex_w)[0])*REAL(X_inn)[l+p*l] + REAL(ex_lambda_beta)[0]*REAL(ex_xi)[0]);
				
				if( fabs(REAL(Beta)[l+p*j]) != 0.0 ){
					for (i=0; i<p; i++) {
						REAL(X_innBetainGamma)[i] = REAL(X_innBetainGamma)[i] + REAL(gamma)[j]*REAL(Beta)[l+p*j]*REAL(X_inn)[i+p*l] - REAL(gamma)[j]*Beta_old*REAL(X_inn)[i+p*l];
						REAL(X_innBeta)[i+p*j] = REAL(X_innBeta)[i+p*j] + REAL(Beta)[l+p*j]*REAL(X_inn)[i+p*l] - Beta_old*REAL(X_inn)[i+p*l];
					}
					
				}
				
				if( fabs(Beta_old) != 0.0 && fabs(REAL(Beta)[l+p*j]) == 0.0 ){
					for (i=0; i<p; i++) {
						REAL(X_innBetainGamma)[i] = REAL(X_innBetainGamma)[i] - REAL(gamma)[j]*Beta_old*REAL(X_inn)[i+p*l];
						REAL(X_innBeta)[i+p*j] = REAL(X_innBeta)[i+p*j] - Beta_old*REAL(X_inn)[i+p*l];
					}				
				}
				
			}
		}
			
		///////////////////////////
		///// Estimate gamma //////
		///////////////////////////

		/// x_star の作成
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &n, &k, &p, &alphaOne, REAL(ex_x), &n, REAL(Beta), &p, &betaZero, REAL(x_star), &n);
		
		/// X_star_inn の作成
		F77_CALL(dgemm)(&TRANST, &TRANSN, &k, &k, &n, &alphaOne, REAL(x_star), &n, REAL(x_star), &n, &betaZero, REAL(X_star_inn), &k);
		
		/// X_starinY の作成
		F77_CALL(dgemv)(&TRANST, &n, &k, &alphaOne, REAL(x_star), &n, REAL(ex_y), &one, &betaZero, REAL(X_starinY), &one);	
		
		/// X_staringamma0 の作成
		F77_CALL(dgemv)(&TRANST, &n, &k, &alphaOne, REAL(x_star), &n, REAL(unit_vec), &one, &betaZero, REAL(X_staringamma0), &one);	
		
		/// X_starinngamma の作成
		F77_CALL(dgemv)(&TRANSN, &k, &k, &alphaOne, REAL(X_star_inn), &k, REAL(gamma), &one, &betaZero, REAL(X_starinngamma), &one);	
		
		/// パラメータの更新 via Covariance Update
		for (l=0; l<k; l++) {
			residual_A = REAL(X_starinY)[l] - REAL(X_staringamma0)[l] - REAL(X_starinngamma)[l];
			s = residual_A + REAL(gamma)[l]*REAL(X_star_inn)[l+k*l];
			gamma_old = REAL(gamma)[l];
			REAL(gamma)[l] = softsh((1.0-REAL(ex_w)[0])*s, 0.5*REAL(ex_lambda_gamma)[0])/((1.0-REAL(ex_w)[0])*REAL(X_star_inn)[l+k*l]);
			
			if( isnan(REAL(gamma)[l]) ) REAL(gamma)[l] = 0.0;
			
			if( fabs(REAL(gamma)[l]) != 0.0 ){			
				for (i=0; i<k; i++) {
					REAL(X_starinngamma)[i] = REAL(X_starinngamma)[i] + REAL(gamma)[l]*REAL(X_star_inn)[i+k*l] - gamma_old*REAL(X_star_inn)[i+k*l];
				}
				
			}
			
			if( fabs(gamma_old) != 0.0 && fabs(REAL(gamma)[l]) == 0.0 ){
				for (i=0; i<k; i++) {
					REAL(X_starinngamma)[i] = REAL(X_starinngamma)[i] - gamma_old*REAL(X_star_inn)[i+k*l];
				}				
			}		
		}
		
	
		// Estimate gamma0
		//  x %*% Betaを作る過程 -> x_Beta_ForGamma0
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &n, &k, &p, &alphaOne, REAL(ex_x), &n, REAL(Beta), &p, &betaZero, REAL(x_Beta_ForGamma0), &n);	
		//  x %*% Beta %*% gamma を作る過程 -> x_Beta_gamma_ForGamma0
		F77_CALL(dgemv)(&TRANSN, &n, &k, &alphaOne, REAL(x_Beta_ForGamma0), &n, REAL(gamma), &one, &betaZero, REAL(x_Beta_gamma_ForGamma0), &one);
		
		s=0.0;
		for (i=0; i<n; i++) {
			s = s + REAL(ex_y)[i] - REAL(x_Beta_gamma_ForGamma0)[i];
		}
		REAL(gamma0)[0] = s/((double)n);
		
		

		// Estimate A
		//  t(x) %*% x を作る過程 -> tXX
		F77_CALL(dgemm)(&TRANST, &TRANSN, &p, &p, &n, &alphaOne, REAL(ex_x), &n, REAL(ex_x), &n, &betaZero, REAL(tXX), &p);	
		//  t(x) %*% x  %*% Beta を作る過程 -> tXX_Beta
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &k, &p, &alphaOne, REAL(tXX), &p, REAL(Beta), &p, &betaZero, REAL(tXX_Beta), &p);	
		// implement SVD
		lwork=-1;
		F77_CALL(dgesdd)(&ATIL, &p, &k, REAL(tXX_Beta), &p, REAL(A_singular), REAL(u), &p, REAL(vt), &k, &wkopt, &lwork, INTEGER(iwork), &info);
		lwork = (int) wkopt;
		work = (double*) malloc( lwork*sizeof(double) );
		/* Compute SVD */
		F77_CALL(dgesdd)(&ATIL, &p, &k, REAL(tXX_Beta), &p, REAL(A_singular), REAL(u), &p, REAL(vt), &k, work, &lwork, INTEGER(iwork), &info);
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &k, &k, &alphaOne, REAL(u), &p, REAL(vt), &k, &betaZero, REAL(A), &p);	
	
		
		
		// para_new を更新
		for (i=0; i<pk; i++) {
			REAL(para_new)[i] = REAL(Beta)[i];
		}
		for (i=pk; i<(pk+k); i++) {
			REAL(para_new)[i] = REAL(gamma)[(i-pk)];
		}
		REAL(para_new)[(pk+k)] = REAL(gamma0)[0];
		
		
		// para_diff を更新
		for (i=0; i<pkkone; i++) {
			REAL(para_diff)[i] = fabs( REAL(para_new)[i] - REAL(para_old)[i] );
		}
		
		// diff_max を更新
		diff_max=REAL(para_diff)[0];
		for (i=1; i<pkkone; i++) {
			if( diff_max < REAL(para_diff)[i] ) diff_max = REAL(para_diff)[i];
		}
		
		// 直前と平均的になにもなかったら break する
		s_sum = 0.0;
		for (i=0; i<pkkone; i++) {
			s_sum = s_sum + REAL(para_diff)[i];
		}
		s_sum = s_sum/((double) pkkone);
		if( s_sum == 0.0 ) break;
		
		// n_itr の更新
		INTEGER(n_itr)[0] = INTEGER(n_itr)[0] + 1;
		
		// メモリの解放
		free(work);	
		
	}		
	
	
	SET_VECTOR_ELT(ans, 0, Beta);
	SET_VECTOR_ELT(ans, 1, gamma);
	SET_VECTOR_ELT(ans, 2, gamma0);
	SET_VECTOR_ELT(ans, 3, A);
	SET_VECTOR_ELT(ans, 4, n_itr);
	UNPROTECT(31);
	return(ans);
		
}


// adaspcr
SEXP adaspcr(SEXP ex_x, SEXP ex_y, SEXP ex_A, SEXP ex_Beta, SEXP ex_gamma, SEXP ex_gamma0, SEXP ex_lambda_beta, SEXP ex_lambda_gamma, SEXP ex_xi, SEXP ex_w, SEXP ex_BetaWeight)
{
	int n = INTEGER(GET_DIM(ex_x))[0];
	int p = INTEGER(GET_DIM(ex_x))[1];
	int k = INTEGER(GET_DIM(ex_A))[1];
	int pk = p*k;
	int pkkone = pk+k+1;
	int one=1;
	int i;
	int j;
	int l;
	int lwork;
	int info;
	
	double s;
	double s_sum;
	double wkopt;
	double* work;
	double diff_max;
	double residual_A;
	double residual_B;
	double Beta_old;
	double gamma_old;
	
	
	SEXP x_star;
	SEXP y_star;
	SEXP A;
	SEXP Beta;
	SEXP gamma;
	SEXP gamma0;
	SEXP x_Beta_ForGamma0;
	SEXP x_Beta_gamma_ForGamma0;
	SEXP ans;
	SEXP A_singular;
	SEXP u;	
	SEXP vt;
	SEXP iwork;
	SEXP tXX;
	SEXP tXX_Beta;
	SEXP para_old;
	SEXP para_new;
	SEXP para_diff;
	SEXP n_itr;
	SEXP X_inn;
	SEXP XinY;
	SEXP Xingamma0;
	SEXP XinY_star;
	SEXP unit_vec;
	SEXP BetainGamma;
	SEXP X_innBetainGamma;
	SEXP X_innBeta;
	SEXP X_star_inn;
	SEXP X_staringamma0;
	SEXP X_starinY;
	SEXP X_starinngamma;
	SEXP BetaWeight;
	
	
	
	PROTECT(x_star = allocMatrix(REALSXP, n, k));
	PROTECT(y_star = allocMatrix(REALSXP, n, k));
	PROTECT(A = allocMatrix(REALSXP, p, k));
	PROTECT(Beta = allocMatrix(REALSXP, p, k));
	PROTECT(gamma = allocVector(REALSXP, k));
	PROTECT(gamma0 = allocVector(REALSXP, one));
	PROTECT(x_Beta_ForGamma0 = allocMatrix(REALSXP, n, k));
	PROTECT(x_Beta_gamma_ForGamma0 = allocVector(REALSXP, n));
	PROTECT(ans = allocVector(VECSXP, 5));
	PROTECT(A_singular = allocVector(REALSXP, p));
	PROTECT(u = allocMatrix(REALSXP, p, p));
	PROTECT(vt = allocMatrix(REALSXP, k, k));
	PROTECT(iwork = allocVector(INTSXP, 8*k));
	PROTECT(tXX = allocMatrix(REALSXP, p, p));
	PROTECT(tXX_Beta = allocMatrix(REALSXP, p, k));
	PROTECT(para_old = allocVector(REALSXP, pkkone));
	PROTECT(para_new = allocVector(REALSXP, pkkone));
	PROTECT(para_diff = allocVector(REALSXP, pkkone));
	PROTECT(n_itr = allocVector(INTSXP, one));
	PROTECT(X_inn = allocMatrix(REALSXP, p, p));
	PROTECT(XinY = allocVector(REALSXP, p));
	PROTECT(Xingamma0 = allocVector(REALSXP, p));
	PROTECT(XinY_star = allocMatrix(REALSXP, p, k));
	PROTECT(unit_vec = allocVector(REALSXP, n));
	PROTECT(BetainGamma = allocVector(REALSXP, p));
	PROTECT(X_innBetainGamma = allocVector(REALSXP, p));
	PROTECT(X_innBeta = allocMatrix(REALSXP, p, k));
	PROTECT(X_star_inn = allocMatrix(REALSXP, k, k));
	PROTECT(X_staringamma0 = allocVector(REALSXP, k));
	PROTECT(X_starinY = allocVector(REALSXP, k));
	PROTECT(X_starinngamma = allocVector(REALSXP, k));
	PROTECT(BetaWeight = allocMatrix(REALSXP, p, k));
	
	
	F77_CALL(dcopy)(&pk, REAL(ex_Beta), &one, REAL(Beta), &one);
	F77_CALL(dcopy)(&k, REAL(ex_gamma), &one, REAL(gamma), &one);
	F77_CALL(dcopy)(&one, REAL(ex_gamma0), &one, REAL(gamma0), &one);
	F77_CALL(dcopy)(&pk, REAL(ex_A), &one, REAL(A), &one);
	
	F77_CALL(dcopy)(&pk, REAL(ex_BetaWeight), &one, REAL(BetaWeight), &one);
	
	
	
	for (i=0; i<pk; i++) {
		REAL(para_new)[i] = REAL(Beta)[i];
	}
	for (i=pk; i<(pk+k); i++) {
		REAL(para_new)[i] = REAL(gamma)[(i-pk)];
	}
	REAL(para_new)[(pk+k)] = REAL(gamma0)[0];
	
	diff_max=1.0;	
	INTEGER(n_itr)[0] = 0;
	
	/// X_inn の作成
	F77_CALL(dgemm)(&TRANST, &TRANSN, &p, &p, &n, &alphaOne, REAL(ex_x), &n, REAL(ex_x), &n, &betaZero, REAL(X_inn), &p);
	
	/// XinY の作成
	F77_CALL(dgemv)(&TRANST, &n, &p, &alphaOne, REAL(ex_x), &n, REAL(ex_y), &one, &betaZero, REAL(XinY), &one);
	
	// パラメータの更新スタート
	while ( diff_max > 0.001 ) {
		F77_CALL(dcopy)(&pkkone, REAL(para_new), &one, REAL(para_old), &one);				
		
		/////////////////////////
		///// Estimate Beta /////
		/////////////////////////
		
		/// y_star の作成
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &n, &k, &p, &alphaOne, REAL(ex_x), &n, REAL(A), &p, &betaZero, REAL(y_star), &n);
		
		/// Xingamma0 の作成
		for (i=0; i<n; i++) {
			REAL(unit_vec)[i] = REAL(gamma0)[0];
		}
		F77_CALL(dgemv)(&TRANST, &n, &p, &alphaOne, REAL(ex_x), &n, REAL(unit_vec), &one, &betaZero, REAL(Xingamma0), &one);	
		
		/// XinY_star の作成
		F77_CALL(dgemm)(&TRANST, &TRANSN, &p, &k, &n, &alphaOne, REAL(ex_x), &n, REAL(y_star), &n, &betaZero, REAL(XinY_star), &p);
		
		/// X_innBetainGamma の作成
		F77_CALL(dgemv)(&TRANSN, &p, &k, &alphaOne, REAL(Beta), &p, REAL(gamma), &one, &betaZero, REAL(BetainGamma), &one);	
		F77_CALL(dgemv)(&TRANSN, &p, &p, &alphaOne, REAL(X_inn), &p, REAL(BetainGamma), &one, &betaZero, REAL(X_innBetainGamma), &one);	
		
		/// X_innBeta の作成
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &k, &p, &alphaOne, REAL(X_inn), &p, REAL(Beta), &p, &betaZero, REAL(X_innBeta), &p);
		
		/// パラメータの更新 via Covariance Update
		for (j=0; j<k; j++) {	
			for (l=0; l<p; l++) {
				
				residual_A = REAL(XinY)[l] - REAL(Xingamma0)[l] - REAL(X_innBetainGamma)[l];
				residual_B = REAL(XinY_star)[l+p*j] - REAL(X_innBeta)[l+p*j];
				s = (1.0-REAL(ex_w)[0])*REAL(gamma)[j]*residual_A + REAL(ex_w)[0]*residual_B + REAL(Beta)[l+p*j]*REAL(X_inn)[l+p*l]*(REAL(ex_w)[0] + (1.0-REAL(ex_w)[0])*(pow(REAL(gamma)[j], 2.0)));
				Beta_old = REAL(Beta)[l+p*j];
				REAL(Beta)[l+p*j] = softsh(s, 0.5*REAL(ex_lambda_beta)[0]*(1.0-REAL(ex_xi)[0])/(fabs(REAL(BetaWeight)[l+p*j])+0.0000001))/(((1.0-REAL(ex_w)[0])*pow(REAL(gamma)[j], 2.0)+REAL(ex_w)[0])*REAL(X_inn)[l+p*l] + REAL(ex_lambda_beta)[0]*REAL(ex_xi)[0]);
				
				if( fabs(REAL(Beta)[l+p*j]) != 0.0 ){
					for (i=0; i<p; i++) {
						REAL(X_innBetainGamma)[i] = REAL(X_innBetainGamma)[i] + REAL(gamma)[j]*REAL(Beta)[l+p*j]*REAL(X_inn)[i+p*l] - REAL(gamma)[j]*Beta_old*REAL(X_inn)[i+p*l];
						REAL(X_innBeta)[i+p*j] = REAL(X_innBeta)[i+p*j] + REAL(Beta)[l+p*j]*REAL(X_inn)[i+p*l] - Beta_old*REAL(X_inn)[i+p*l];
					}
					
				}
				
				if( fabs(Beta_old) != 0.0 && fabs(REAL(Beta)[l+p*j]) == 0.0 ){
					for (i=0; i<p; i++) {
						REAL(X_innBetainGamma)[i] = REAL(X_innBetainGamma)[i] - REAL(gamma)[j]*Beta_old*REAL(X_inn)[i+p*l];
						REAL(X_innBeta)[i+p*j] = REAL(X_innBeta)[i+p*j] - Beta_old*REAL(X_inn)[i+p*l];
					}				
				}
				
			}
		}
		
		///////////////////////////
		///// Estimate gamma //////
		///////////////////////////
		
		/// x_star の作成
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &n, &k, &p, &alphaOne, REAL(ex_x), &n, REAL(Beta), &p, &betaZero, REAL(x_star), &n);
		
		/// X_star_inn の作成
		F77_CALL(dgemm)(&TRANST, &TRANSN, &k, &k, &n, &alphaOne, REAL(x_star), &n, REAL(x_star), &n, &betaZero, REAL(X_star_inn), &k);
		
		/// X_starinY の作成
		F77_CALL(dgemv)(&TRANST, &n, &k, &alphaOne, REAL(x_star), &n, REAL(ex_y), &one, &betaZero, REAL(X_starinY), &one);	
		
		/// X_staringamma0 の作成
		F77_CALL(dgemv)(&TRANST, &n, &k, &alphaOne, REAL(x_star), &n, REAL(unit_vec), &one, &betaZero, REAL(X_staringamma0), &one);	
		
		/// X_starinngamma の作成
		F77_CALL(dgemv)(&TRANSN, &k, &k, &alphaOne, REAL(X_star_inn), &k, REAL(gamma), &one, &betaZero, REAL(X_starinngamma), &one);	
		
		/// パラメータの更新 via Covariance Update
		for (l=0; l<k; l++) {
			residual_A = REAL(X_starinY)[l] - REAL(X_staringamma0)[l] - REAL(X_starinngamma)[l];
			s = residual_A + REAL(gamma)[l]*REAL(X_star_inn)[l+k*l];
			gamma_old = REAL(gamma)[l];
			REAL(gamma)[l] = softsh((1.0-REAL(ex_w)[0])*s, 0.5*REAL(ex_lambda_gamma)[0])/((1.0-REAL(ex_w)[0])*REAL(X_star_inn)[l+k*l]);
			
			if( isnan(REAL(gamma)[l]) ) REAL(gamma)[l] = 0.0;
			
			if( fabs(REAL(gamma)[l]) != 0.0 ){			
				for (i=0; i<k; i++) {
					REAL(X_starinngamma)[i] = REAL(X_starinngamma)[i] + REAL(gamma)[l]*REAL(X_star_inn)[i+k*l] - gamma_old*REAL(X_star_inn)[i+k*l];
				}
				
			}
			
			if( fabs(gamma_old) != 0.0 && fabs(REAL(gamma)[l]) == 0.0 ){
				for (i=0; i<k; i++) {
					REAL(X_starinngamma)[i] = REAL(X_starinngamma)[i] - gamma_old*REAL(X_star_inn)[i+k*l];
				}				
			}		
		}
		
		
		// Estimate gamma0
		//  x %*% Betaを作る過程 -> x_Beta_ForGamma0
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &n, &k, &p, &alphaOne, REAL(ex_x), &n, REAL(Beta), &p, &betaZero, REAL(x_Beta_ForGamma0), &n);	
		//  x %*% Beta %*% gamma を作る過程 -> x_Beta_gamma_ForGamma0
		F77_CALL(dgemv)(&TRANSN, &n, &k, &alphaOne, REAL(x_Beta_ForGamma0), &n, REAL(gamma), &one, &betaZero, REAL(x_Beta_gamma_ForGamma0), &one);
		
		s=0.0;
		for (i=0; i<n; i++) {
			s = s + REAL(ex_y)[i] - REAL(x_Beta_gamma_ForGamma0)[i];
		}
		REAL(gamma0)[0] = s/((double)n);
		
		
		
		// Estimate A
		//  t(x) %*% x を作る過程 -> tXX
		F77_CALL(dgemm)(&TRANST, &TRANSN, &p, &p, &n, &alphaOne, REAL(ex_x), &n, REAL(ex_x), &n, &betaZero, REAL(tXX), &p);	
		//  t(x) %*% x  %*% Beta を作る過程 -> tXX_Beta
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &k, &p, &alphaOne, REAL(tXX), &p, REAL(Beta), &p, &betaZero, REAL(tXX_Beta), &p);	
		// implement SVD
		lwork=-1;
		F77_CALL(dgesdd)(&ATIL, &p, &k, REAL(tXX_Beta), &p, REAL(A_singular), REAL(u), &p, REAL(vt), &k, &wkopt, &lwork, INTEGER(iwork), &info);
		lwork = (int) wkopt;
		work = (double*) malloc( lwork*sizeof(double) );
		/* Compute SVD */
		F77_CALL(dgesdd)(&ATIL, &p, &k, REAL(tXX_Beta), &p, REAL(A_singular), REAL(u), &p, REAL(vt), &k, work, &lwork, INTEGER(iwork), &info);
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &k, &k, &alphaOne, REAL(u), &p, REAL(vt), &k, &betaZero, REAL(A), &p);	
		
		
		
		// para_new を更新
		for (i=0; i<pk; i++) {
			REAL(para_new)[i] = REAL(Beta)[i];
		}
		for (i=pk; i<(pk+k); i++) {
			REAL(para_new)[i] = REAL(gamma)[(i-pk)];
		}
		REAL(para_new)[(pk+k)] = REAL(gamma0)[0];
		
		
		
		
		
		// para_diff を更新
		for (i=0; i<pkkone; i++) {
			REAL(para_diff)[i] = fabs( REAL(para_new)[i] - REAL(para_old)[i] );
		}
		
		// diff_max を更新
		diff_max=REAL(para_diff)[0];
		for (i=1; i<pkkone; i++) {
			if( diff_max < REAL(para_diff)[i] ) diff_max = REAL(para_diff)[i];
		}
		
		// 直前と平均的になにもなかったら break する
		s_sum = 0.0;
		for (i=0; i<pkkone; i++) {
			s_sum = s_sum + REAL(para_diff)[i];
		}
		s_sum = s_sum/((double) pkkone);
		if( s_sum == 0.0 ) break;
		
		// n_itr の更新
		INTEGER(n_itr)[0] = INTEGER(n_itr)[0] + 1;
		
		// メモリの解放
		free(work);	
		
	}		
	
	
	SET_VECTOR_ELT(ans, 0, Beta);
	SET_VECTOR_ELT(ans, 1, gamma);
	SET_VECTOR_ELT(ans, 2, gamma0);
	SET_VECTOR_ELT(ans, 3, A);
	SET_VECTOR_ELT(ans, 4, n_itr);
	UNPROTECT(32);
	return(ans);
	
}




// inispcr
SEXP inispcr(SEXP ex_x, SEXP ex_y, SEXP ex_A, SEXP ex_Beta, SEXP ex_gamma, SEXP ex_gamma0, SEXP ex_lambda_beta, SEXP ex_lambda_gamma, SEXP ex_xi, SEXP ex_w)
{
	int n = INTEGER(GET_DIM(ex_x))[0];
	int p = INTEGER(GET_DIM(ex_x))[1];
	int k = INTEGER(GET_DIM(ex_A))[1];
	int pk = p*k;
	int pkkone = pk+k+1;
	int one=1;
	int i;
	int j;
	int l;
	int lwork;
	int info;
	
	double s;
	double wkopt;
	double* work;
	double residual_A;
	double residual_B;
	double Beta_old;
	double gamma_old;
	
    
	
	SEXP x_star;
	SEXP y_star;
	SEXP A;
	SEXP Beta;
	SEXP gamma;
	SEXP gamma0;
	SEXP x_Beta_ForGamma0;
	SEXP x_Beta_gamma_ForGamma0;
	SEXP ans;
	SEXP A_singular;
	SEXP u;	
	SEXP vt;
	SEXP iwork;
	SEXP tXX;
	SEXP tXX_Beta;
	SEXP para_new;
	SEXP para_diff;
	SEXP n_itr;
	SEXP X_inn;
	SEXP XinY;
	SEXP Xingamma0;
	SEXP XinY_star;
	SEXP unit_vec;
	SEXP BetainGamma;
	SEXP X_innBetainGamma;
	SEXP X_innBeta;
	SEXP X_star_inn;
	SEXP X_staringamma0;
	SEXP X_starinY;
	SEXP X_starinngamma;
	
	
	PROTECT(x_star = allocMatrix(REALSXP, n, k));
	PROTECT(y_star = allocMatrix(REALSXP, n, k));
	PROTECT(A = allocMatrix(REALSXP, p, k));
	PROTECT(Beta = allocMatrix(REALSXP, p, k));
	PROTECT(gamma = allocVector(REALSXP, k));
	PROTECT(gamma0 = allocVector(REALSXP, one));
	PROTECT(x_Beta_ForGamma0 = allocMatrix(REALSXP, n, k));
	PROTECT(x_Beta_gamma_ForGamma0 = allocVector(REALSXP, n));
	PROTECT(ans = allocVector(VECSXP, 5));
	PROTECT(A_singular = allocVector(REALSXP, p));
	PROTECT(u = allocMatrix(REALSXP, p, p));
	PROTECT(vt = allocMatrix(REALSXP, k, k));
	PROTECT(iwork = allocVector(INTSXP, 8*k));
	PROTECT(tXX = allocMatrix(REALSXP, p, p));
	PROTECT(tXX_Beta = allocMatrix(REALSXP, p, k));
	PROTECT(para_new = allocVector(REALSXP, pkkone));
	PROTECT(para_diff = allocVector(REALSXP, pkkone));
	PROTECT(n_itr = allocVector(INTSXP, one));
	PROTECT(X_inn = allocMatrix(REALSXP, p, p));
	PROTECT(XinY = allocVector(REALSXP, p));
	PROTECT(Xingamma0 = allocVector(REALSXP, p));
	PROTECT(XinY_star = allocMatrix(REALSXP, p, k));
	PROTECT(unit_vec = allocVector(REALSXP, n));
	PROTECT(BetainGamma = allocVector(REALSXP, p));
	PROTECT(X_innBetainGamma = allocVector(REALSXP, p));
	PROTECT(X_innBeta = allocMatrix(REALSXP, p, k));
	PROTECT(X_star_inn = allocMatrix(REALSXP, k, k));
	PROTECT(X_staringamma0 = allocVector(REALSXP, k));
	PROTECT(X_starinY = allocVector(REALSXP, k));
	PROTECT(X_starinngamma = allocVector(REALSXP, k));
	
	
	F77_CALL(dcopy)(&pk, REAL(ex_Beta), &one, REAL(Beta), &one);
	F77_CALL(dcopy)(&k, REAL(ex_gamma), &one, REAL(gamma), &one);
	F77_CALL(dcopy)(&one, REAL(ex_gamma0), &one, REAL(gamma0), &one);
	F77_CALL(dcopy)(&pk, REAL(ex_A), &one, REAL(A), &one);
	
	
	
	for (i=0; i<pk; i++) {
		REAL(para_new)[i] = REAL(Beta)[i];
	}
	for (i=pk; i<(pk+k); i++) {
		REAL(para_new)[i] = REAL(gamma)[(i-pk)];
	}
	REAL(para_new)[(pk+k)] = REAL(gamma0)[0];
	
	INTEGER(n_itr)[0] = 0;
	
	/// X_inn の作成
	F77_CALL(dgemm)(&TRANST, &TRANSN, &p, &p, &n, &alphaOne, REAL(ex_x), &n, REAL(ex_x), &n, &betaZero, REAL(X_inn), &p);
	
	/// XinY の作成
	F77_CALL(dgemv)(&TRANST, &n, &p, &alphaOne, REAL(ex_x), &n, REAL(ex_y), &one, &betaZero, REAL(XinY), &one);
	
	
	
	// パラメータの更新スタート			
	
	/////////////////////////
	///// Estimate Beta /////
	/////////////////////////
	
	/// y_star の作成
	F77_CALL(dgemm)(&TRANSN, &TRANSN, &n, &k, &p, &alphaOne, REAL(ex_x), &n, REAL(A), &p, &betaZero, REAL(y_star), &n);
	
	/// Xingamma0 の作成
	for (i=0; i<n; i++) {
		REAL(unit_vec)[i] = REAL(gamma0)[0];
	}
	F77_CALL(dgemv)(&TRANST, &n, &p, &alphaOne, REAL(ex_x), &n, REAL(unit_vec), &one, &betaZero, REAL(Xingamma0), &one);	
	
	/// XinY_star の作成
	F77_CALL(dgemm)(&TRANST, &TRANSN, &p, &k, &n, &alphaOne, REAL(ex_x), &n, REAL(y_star), &n, &betaZero, REAL(XinY_star), &p);
	
	/// X_innBetainGamma の作成
	F77_CALL(dgemv)(&TRANSN, &p, &k, &alphaOne, REAL(Beta), &p, REAL(gamma), &one, &betaZero, REAL(BetainGamma), &one);	
	F77_CALL(dgemv)(&TRANSN, &p, &p, &alphaOne, REAL(X_inn), &p, REAL(BetainGamma), &one, &betaZero, REAL(X_innBetainGamma), &one);	
	
	/// X_innBeta の作成
	F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &k, &p, &alphaOne, REAL(X_inn), &p, REAL(Beta), &p, &betaZero, REAL(X_innBeta), &p);
	
	/// パラメータの更新 via Covariance Update
	for (j=0; j<k; j++) {	
		for (l=0; l<p; l++) {
			
			residual_A = REAL(XinY)[l] - REAL(Xingamma0)[l] - REAL(X_innBetainGamma)[l];
			residual_B = REAL(XinY_star)[l+p*j] - REAL(X_innBeta)[l+p*j];
			s = (1.0-REAL(ex_w)[0])*REAL(gamma)[j]*residual_A + REAL(ex_w)[0]*residual_B + REAL(Beta)[l+p*j]*REAL(X_inn)[l+p*l]*(REAL(ex_w)[0] + (1.0-REAL(ex_w)[0])*(pow(REAL(gamma)[j], 2.0)));
			Beta_old = REAL(Beta)[l+p*j];
			REAL(Beta)[l+p*j] = softsh(s, 0.5*REAL(ex_lambda_beta)[0]*(1.0-REAL(ex_xi)[0]))/(((1.0-REAL(ex_w)[0])*pow(REAL(gamma)[j], 2.0)+REAL(ex_w)[0])*REAL(X_inn)[l+p*l] + REAL(ex_lambda_beta)[0]*REAL(ex_xi)[0]);
			
			if( fabs(REAL(Beta)[l+p*j]) != 0.0 ){
				for (i=0; i<p; i++) {
					REAL(X_innBetainGamma)[i] = REAL(X_innBetainGamma)[i] + REAL(gamma)[j]*REAL(Beta)[l+p*j]*REAL(X_inn)[i+p*l] - REAL(gamma)[j]*Beta_old*REAL(X_inn)[i+p*l];
					REAL(X_innBeta)[i+p*j] = REAL(X_innBeta)[i+p*j] + REAL(Beta)[l+p*j]*REAL(X_inn)[i+p*l] - Beta_old*REAL(X_inn)[i+p*l];
				}
				
			}
			
			if( fabs(Beta_old) != 0.0 && fabs(REAL(Beta)[l+p*j]) == 0.0 ){
				for (i=0; i<p; i++) {
					REAL(X_innBetainGamma)[i] = REAL(X_innBetainGamma)[i] - REAL(gamma)[j]*Beta_old*REAL(X_inn)[i+p*l];
					REAL(X_innBeta)[i+p*j] = REAL(X_innBeta)[i+p*j] - Beta_old*REAL(X_inn)[i+p*l];
				}				
			}
			
		}
	}
	
	///////////////////////////
	///// Estimate gamma //////
	///////////////////////////
	
	/// x_star の作成
	F77_CALL(dgemm)(&TRANSN, &TRANSN, &n, &k, &p, &alphaOne, REAL(ex_x), &n, REAL(Beta), &p, &betaZero, REAL(x_star), &n);
	
	/// X_star_inn の作成
	F77_CALL(dgemm)(&TRANST, &TRANSN, &k, &k, &n, &alphaOne, REAL(x_star), &n, REAL(x_star), &n, &betaZero, REAL(X_star_inn), &k);
	
	/// X_starinY の作成
	F77_CALL(dgemv)(&TRANST, &n, &k, &alphaOne, REAL(x_star), &n, REAL(ex_y), &one, &betaZero, REAL(X_starinY), &one);	
	
	/// X_staringamma0 の作成
	F77_CALL(dgemv)(&TRANST, &n, &k, &alphaOne, REAL(x_star), &n, REAL(unit_vec), &one, &betaZero, REAL(X_staringamma0), &one);	
	
	/// X_starinngamma の作成
	F77_CALL(dgemv)(&TRANSN, &k, &k, &alphaOne, REAL(X_star_inn), &k, REAL(gamma), &one, &betaZero, REAL(X_starinngamma), &one);	
	
	/// パラメータの更新 via Covariance Update
	for (l=0; l<k; l++) {
		residual_A = REAL(X_starinY)[l] - REAL(X_staringamma0)[l] - REAL(X_starinngamma)[l];
		s = residual_A + REAL(gamma)[l]*REAL(X_star_inn)[l+k*l];
		gamma_old = REAL(gamma)[l];
		REAL(gamma)[l] = softsh((1.0-REAL(ex_w)[0])*s, 0.5*REAL(ex_lambda_gamma)[0])/((1.0-REAL(ex_w)[0])*REAL(X_star_inn)[l+k*l]);
		
		if( isnan(REAL(gamma)[l]) ) REAL(gamma)[l] = 0.0;
		
		if( fabs(REAL(gamma)[l]) != 0.0 ){			
			for (i=0; i<k; i++) {
				REAL(X_starinngamma)[i] = REAL(X_starinngamma)[i] + REAL(gamma)[l]*REAL(X_star_inn)[i+k*l] - gamma_old*REAL(X_star_inn)[i+k*l];
			}
			
		}
		
		if( fabs(gamma_old) != 0.0 && fabs(REAL(gamma)[l]) == 0.0 ){
			for (i=0; i<k; i++) {
				REAL(X_starinngamma)[i] = REAL(X_starinngamma)[i] - gamma_old*REAL(X_star_inn)[i+k*l];
			}				
		}		
	}
	
	
	// Estimate gamma0
	//  x %*% Betaを作る過程 -> x_Beta_ForGamma0
	F77_CALL(dgemm)(&TRANSN, &TRANSN, &n, &k, &p, &alphaOne, REAL(ex_x), &n, REAL(Beta), &p, &betaZero, REAL(x_Beta_ForGamma0), &n);	
	//  x %*% Beta %*% gamma を作る過程 -> x_Beta_gamma_ForGamma0
	F77_CALL(dgemv)(&TRANSN, &n, &k, &alphaOne, REAL(x_Beta_ForGamma0), &n, REAL(gamma), &one, &betaZero, REAL(x_Beta_gamma_ForGamma0), &one);
	
	s=0.0;
	for (i=0; i<n; i++) {
		s = s + REAL(ex_y)[i] - REAL(x_Beta_gamma_ForGamma0)[i];
	}
	REAL(gamma0)[0] = s/((double)n);
	
	
	
	// Estimate A
	//  t(x) %*% x を作る過程 -> tXX
	F77_CALL(dgemm)(&TRANST, &TRANSN, &p, &p, &n, &alphaOne, REAL(ex_x), &n, REAL(ex_x), &n, &betaZero, REAL(tXX), &p);	
	//  t(x) %*% x  %*% Beta を作る過程 -> tXX_Beta
	F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &k, &p, &alphaOne, REAL(tXX), &p, REAL(Beta), &p, &betaZero, REAL(tXX_Beta), &p);	
	// implement SVD
	lwork=-1;
	F77_CALL(dgesdd)(&ATIL, &p, &k, REAL(tXX_Beta), &p, REAL(A_singular), REAL(u), &p, REAL(vt), &k, &wkopt, &lwork, INTEGER(iwork), &info);
	lwork = (int) wkopt;
	work = (double*) malloc( lwork*sizeof(double) );
	/* Compute SVD */
	F77_CALL(dgesdd)(&ATIL, &p, &k, REAL(tXX_Beta), &p, REAL(A_singular), REAL(u), &p, REAL(vt), &k, work, &lwork, INTEGER(iwork), &info);
	F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &k, &k, &alphaOne, REAL(u), &p, REAL(vt), &k, &betaZero, REAL(A), &p);	
	
	
	
	// para_new を更新
	for (i=0; i<pk; i++) {
		REAL(para_new)[i] = REAL(Beta)[i];
	}
	for (i=pk; i<(pk+k); i++) {
		REAL(para_new)[i] = REAL(gamma)[(i-pk)];
	}
	REAL(para_new)[(pk+k)] = REAL(gamma0)[0];
	
	
	// n_itr の更新
	INTEGER(n_itr)[0] = INTEGER(n_itr)[0] + 1;
	
	// メモリの解放
	free(work);	
	
	SET_VECTOR_ELT(ans, 0, Beta);
	SET_VECTOR_ELT(ans, 1, gamma);
	SET_VECTOR_ELT(ans, 2, gamma0);
	SET_VECTOR_ELT(ans, 3, A);
	SET_VECTOR_ELT(ans, 4, n_itr);
	UNPROTECT(30);
	return(ans);
	
}
