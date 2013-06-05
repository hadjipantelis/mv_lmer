#Code is based on a presentation by D.Bates in  University of Lausanne on July 2, 2009
#Depending on the provided nlev* variable we get the corresponding mixed effects regression
#This code is experimental and has known "shortcomings"; criticism will very appreciated.

PrfDvnce<- function(theta, Zt, X, y, XtX, nlevU = NA, nlevL = NA, nlev2 = NA, SimpleCase = NA, ML=0, 
  GiveBeta=0, GiveSigma=0, GiveGamma=0, GiveL= 0, GiveRX= 0, GiveUt=0, GiveRZX= 0, GiveLambda = 0, GiveFitted =0, GiveCorr=0){

  #theta 	: vector of parameters to optimize for. Usually ratios like sqrt( sigma_e^2 / sigma_g^2) and correlations (doubles)
	#Zt		: traspose of Z matrix. Usually an indicator matrix (matrix)
	#X			: Design matrix for the mean (matrix)
	#XtX		: crossproduct of X, provided to avoid unnecessary computation (matrix)
	#nlevU	: number of upper grouping (betwee group) in the case of multilevel regression (int)
	#nlevL 	: number of lower grouping (within group) in the case of multilevel regression (int)
	#nlev2 	: number of second grouping in the case of having multiple regressors (int)
	#SimpleCase: Do we actually have a vanilla case of single random effect? (1/0)
	#ML		: produce ML instead of REML profiled estimate (1/0)
	#GiveBeta : Provide Fixed Effects covariates (1/0)
	#GiveSigma: Provide Standard Deviation of the residuals (1/0)
	#GiveGamma: Provide Random Effects estimates (1/0)
	#GiveL     	: Provide  L Cholesky factor
	#GiveRX	: Provide RX submatrix
	#GiveRZX  : Provide RZX submatrix
	#GiveUt	: Provide Ut matrix
	#GiveLambda: Provide Lambda matrix
	#GiveFitted: Provide Fitted values		

	#Calculate TOTAL number of random effects levels
	nlev <-  dim(Zt)[1] 
	#Abs() variables to make sure all ratios are positive
	theta <- abs(theta)

	#Construct Lambda	

	if( !is.na(nlev2) ){
	#if you have a two crossed ranodm effects LMER
		#Set Lambda the diagonal matrix holding the ratio sqrt( sigma_e^2 / sigma_g^2)
		Lambda <- Diagonal(x=c(rep.int(theta[1], nlev-nlev2), rep.int(theta[2], nlev2 ) ))

	} else if ( !is.na(nlevU) &&  !is.na(nlevL) ){
	#if you have a two nested random efeects LMER	
		#Set diagonal matrix with the relative ratios of each grouping
		D <- Diagonal(theta[1:nlevU],n=nlevU) 
		#Set lower triangular view of the correlation matrix "cos-ing" to make sure everying is -1
		Pl <- diag(x= 1, nlevU )
		Pl[ lower.tri(Pl)] <-  cos(theta[ (nlevU +1): length(theta ) ]) * .2 #######THIS .3 IS JUST A MANUAL VALUE! NEEDS TO BE COMPUTED
		#Get lower triangular view of relative compacted precission matrix
		if (GiveCorr == 1) {return (Pl)}
		Kl = D %*% Pl ##  %*% D delte
		#Expamd Lambda to each appropriate dimensions
		Lambda <-  kronecker( Kl , Diagonal(x=1, nlevL)) 
	
	} else if ( !is.na(SimpleCase) ){ 
	#if you have a single randome effect LMER	
		#Set Lambda the diagonal matrix holding the ratio sqrt( sigma_e^2 / sigma_g^2)
		Lambda <- Diagonal(x=c(rep.int(theta[1], nlev)))

	} else {
		print('Something has not gone according to plan. Check your nlev* arguments')
		return ( 1/0 ) 
	}

	#Set Lambda the diagonal matrix holding the ratio sqrt( sigma_e^2 / sigma_g^2)
	if( GiveLambda ==1) { return( Lambda) }
	#Calculate the product U(theta) = Z* Lambda(theta) holding the spherical random-effects vector Ut
	Ut <- crossprod(Lambda, Zt)
	if( GiveUt==1) {return (Ut) }
	#Calculate L as the sparse LowerTriangular matrix satisfying t(U)*U+I = L*t(L)
	L <- Cholesky(tcrossprod(Ut), LDL = FALSE, Imult = 1) 
	#Left Cholesky factor of a positive-definite symmetric matrix.
	L <- update(L, Ut, mult = 1)
	if( GiveL==1) {return (L) }
	#Set cu as the solution to Lcu = PU'y or practically the 
	cu <- solve(L, solve(L, Ut %*% y, sys = "P"), sys = "L")
	RZX <- solve(L, solve(L, Ut %*% X, sys = "P"), sys = "L")
	if(GiveRZX==1) {return (RZX)}
	RX <- chol(XtX - crossprod(RZX))
	if(GiveRX==1) {return (RX)}
	cb <- solve(t(RX),crossprod(X,y)- crossprod(RZX, cu))
	beta <- solve(RX, cb)
	if (GiveBeta==1) { return(beta) }
	u <- solve(L,solve(L,cu - RZX %*% beta, sys="Lt"), sys="Pt")
	if (GiveGamma==1) {
		if (!is.na(SimpleCase) ){ return ( u * theta ) }
		if (!is.na(nlev2) ){return( u * diag(Lambda) ) }
		if (!is.na(nlevL) ){return( matrix( Lambda %*%  u %*% Lambda , nrow= nlevL)) } #This is almost surely wrong. Placeholder code.
	}
	#Calculated fitted values as sum of random and fixed effect influence
	yHat <- as.vector(crossprod(Ut, u) + X %*% beta)
	if (GiveFitted==1) { return( yHat ) }
	#Set Penalized Residual Sum of Square differences
	prss <- sum(c(y - yHat, as.vector(u))^2)  
	#Set number of samples and number of parameters
	n <- length(y); p <- ncol(RX)
	#Conditional Estimate for the variance:
	if( GiveSigma == 1) {  return (sqrt(prss/(n-p))) }
	if( ML ==1) { return( 2*determinant(L)$mod  +(n)*(1+log(2*pi*prss/(n)))) }
	return( 2*determinant(L)$mod+ 2* determinant(RX)$mod +(n-p)*(1+log(2*pi*prss/(n-p))))
}
