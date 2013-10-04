#Depending on the provided mask and nlev* variables we get the corresponding mixed effects regression
#Code is based on a presentation by D.Bates in  University of Lausanne on July 2, 2009
#This version accommodates one or two random effects only. It can work in a univariate and a multivariate setting.
#It allows specification of the correlation patterns between the mixed effects covariances.

PrfDvnce_Mask_GaCo <- function(theta, Zt, X, y, XtX, nlevU = 1, nlevL1 = NA, nlevL2 = NA, ML=0, GiveBeta=0, GiveSigma=0, GiveGamma=0, GiveL= 0, GiveRX= 0, GiveUt=0, GiveRZX= 0, GiveLambda = 0, GiveFitted =0, GiveCorrelations =0,GiveCovariances = 0, GiveResiduals=0, Mask = NA, GarCol=0){

  #theta 		: vector of parameters to optimize for. ratios like sqrt( sigma_e^2 / sigma_g^2)
	#Zt		: traspose of Z matrix. Usually an indicator matrix (matrix)
	#X		: Design matrix for the mean (matrix)
	#XtX		: crossproduct of X, provided to avoid unnecessary computation (matrix)
	#nlevU		: number of upper grouping (betwee group) in the case of multilevel regression (int)
	#nlevL 		: number of lower grouping (within group) in the case of multilevel regression (int)
	#ML		: produce ML instead of REML profiled estimate (1/0)
	#GiveBeta 	: Provide Fixed Effects covariates (1/0)
	#GiveSigma	: Provide Standard Deviation of the residuals (1/0)
	#GiveGamma	: Provide Random Effects estimates (1/0)
	#GiveL  	: Provide L Cholesky factor (1/0)
	#GiveRX		: Provide RX submatrix (1/0)
	#GiveRZX	: Provide RZX submatrix (1/0)
	#GiveUt		: Provide Ut matrix (1/0)
	#GiveLambda	: Provide Lambda matrix (1/0)
	#GiveFitted	: Provide Fitted values (1/0)
	#Mask		: Masking Matrix, that enforces the independence assumptions we want. (square bounded matrix)
	#GarCol		: Garbage collection on
	
	#Check if Mask exists
	if (typeof(Mask) != "double"){
		Mask = matrix( rep(1, nlevU^2), nrow = nlevU)} 
	#and if it is marginally reasonable
	if(  (max(Mask) > 1) | (min(Mask) <0) | (dim(Mask)[1] != nlevU) |  (dim(Mask)[2] != nlevU)  ) {
		print('Something has not gone according to plan. Check your Mask\'s properties.')
		return ( 1/0 ) 	}

	Checker=0; #Dummy variable to check if a valid Lambda was computed.

	#Constructing Lambda For A single Random Effect
	if ( !is.na(nlevU) &&  !is.na(nlevL1) &&  !is.na(nlevL2) && ( nlevL2 == 0 ) ){
	#if you have a nested random effects structure
 		middle <- sum(1:nlevU)
		
		#Built up the "wrong/full structure" covariance matrix using the Unconstrainted L1
		preKl1 <- diag(nlevU)
		preKl1[lower.tri(preKl1, diag = TRUE)] <- theta[1:middle]		
		fullK1 <- tcrossprod( preKl1 )	
		
		#Enforce the mask on it
		sparseK1 <- Mask * fullK1;
		#Check that it is PSD

		Kl1 <- try(   t(chol(sparseK1)), silent=T);
		if (!is.matrix(Kl1) ){
			eigenvalues1 <- eigen(only.values=T, sparseK1)$values
			#If you don't have a PSD make it at least PSD and then Tikhonov it.
			sparseK1 <- sparseK1 +  ((diag(nlevU)) * (abs(min(eigenvalues1)) + 10^(-6)) )}
			#Compute the actual lower triangular matrix we want without worrying about it's structure.

		Kl1 <- t(chol(sparseK1));

		FullSizeOfL<-nlevU * (nlevL1)
		Lambda <- Matrix(sparse=T, diag(FullSizeOfL));
		
		#Expamd Lambda to its appropriate dimensions
		Lambda[1:(nlevU *nlevL1),1:(nlevU *nlevL1)]  <- kronecker( Kl1 , Diagonal(x=1, nlevL1)) 
		Checker=1}	
			
	#Constructing Lambda	
	if ( !is.na(nlevU) &&  !is.na(nlevL1) &&  !is.na(nlevL2) && ( nlevL2 != 0 ) ){
	#if you have a nested random effects structure
 		middle <- sum(1:nlevU)
		
		#Built up the "wrong/full structure" covariance matrix using the Unconstrainted L1
		preKl1 <- diag(nlevU)
		preKl1[lower.tri(preKl1, diag = TRUE)] <- theta[1:middle]		
		fullK1 <- tcrossprod( preKl1 )	
		#Built up the "wrong/full structure" covariance matrix using the Unconstrainted L2
		preKl2 <- diag(nlevU)
		preKl2[lower.tri(preKl2, diag = TRUE)] <- theta[(1+middle): length(theta)]		
		fullK2 <- tcrossprod( preKl2 )	

		#Enforce the mask on it
		sparseK1 <- Mask * fullK1;
			
		#Check that it is PSD
		Kl1 <- try(   t(chol(sparseK1)), silent=T);
		if (!is.matrix(Kl1) ){
			eigenvalues1 <- eigen(only.values=T, sparseK1)$values
			#If you don't have a PSD make it at least PSD and then Tikhonov it. 
			sparseK1 <- sparseK1 +  ((diag(nlevU)) * (abs(min(eigenvalues1)) + 10^(-6)) )}
			#Compute the actual lower triangular matrix we want without worrying about it's structure.

		Kl1 <- t(chol(sparseK1));

		#Enforce the mask on it
		sparseK2 <- Mask * fullK2;
		
		#Check that it is PSD
		Kl2 <- try(   t(chol(sparseK2)), silent=T);
		if (!is.matrix(Kl2) ){
			eigenvalues2 <- eigen(only.values=T, sparseK2)$values
			#If you don't have a PSD make it at least PSD and then Tikhonov it. #This is bonkers.
			sparseK2 <- sparseK2 +  ((diag(nlevU)) * (abs(min(eigenvalues2)) + 10^(-6)) )}
			#Compute the actual lower triangular matrix we want without worrying about it's structure.

		Kl2 <- t(chol(sparseK2));


		FullSizeOfL<-nlevU * (nlevL2 + nlevL1)
		Lambda <- Matrix(sparse=T, diag(FullSizeOfL));
		
		#Expamd Lambda to its appropriate dimensions
		Lambda[1:(nlevU *nlevL1),1:(nlevU *nlevL1)]  <- kronecker( Kl1 , Diagonal(x=1, nlevL1)) 
		Lambda[(1+(nlevU *nlevL1)):FullSizeOfL, (1+(nlevU *nlevL1)): FullSizeOfL]  <- kronecker( Kl2 , Diagonal(x=1, nlevL2)) 
		Checker =1;} 

	if (Checker==0){
		print('Something has not gone according to plan. Check your nlev* arguments')
		return ( 1/0 ) } 
	#Set Lambda the diagonal matrix holding the ratio sqrt( sigma_e^2 / sigma_g^2)
	if( GiveLambda ==1) { 
		return( Lambda) }

	#Calculate the product U(theta) = Z* Lambda(theta) : spherical random-effects vector Ut
	Ut <- crossprod(Lambda, Zt)
	if( GiveUt==1) {
		return (Ut) }
	#Calculate L as the sparse LowerTriangular matrix satisfying t(U)*U+I = L*t(L)
	L <- Cholesky(tcrossprod(Ut), LDL = FALSE, Imult = 1) 
	#Left Cholesky factor of a positive-definite symmetric matrix.
	L <- update(L, Ut, mult = 1)
	if( GiveL==1) {
		return (L) }

	#Set cu as the solution to Lcu = PU'y 
	cu <- solve(L, solve(L, Ut %*% y, sys = "P"), sys = "L")
	#From this point on we are practically solving by parts the QR decomposition 
	RZX <- solve(L, solve(L, Ut %*% X, sys = "P"), sys = "L")
	if(GiveRZX==1) {
		return (RZX)} 
	dimRZX <- dim(RZX) 
	#Recasting RZX as a full matrix is important as the crossproduct of sparse matrices is quite slow
	RX <- chol(XtX - crossprod( matrix( RZX, nrow= dimRZX[1]) ))
	if(GiveRX==1) {
		return (RX)}
	cb <- solve(t(RX),crossprod(X,y)- crossprod(RZX, cu))
	beta <- solve(RX, cb)

	# the whitened realization of the mixed effects
	u <- solve(L,solve(L,cu - RZX %*% beta, sys="Lt"), sys="Pt") 
	if (GiveGamma==1) {
		UnstructuredGamma <- matrix( t(Lambda %*% u), nrow= 1)
		MV_1 <-  (matrix( UnstructuredGamma[1: (nlevU* nlevL1 )], ncol= nlevU))
		if (nlevL2!=0) {
		MV_2 <-  (matrix( UnstructuredGamma[ (1 + (nlevU* nlevL1 )) : length(UnstructuredGamma)], ncol= nlevU))		
		return( list(MV_1,MV_2)  ) }  #Double check it does what you think it does.
		else{
		return( list(MV_1)) } 	}
	rm(Lambda);  rm(RZX); if(GarCol==1) {	gc(); }	
	#Calculated fitted values as sum of random and fixed effect influence
	yHat <- as.vector(crossprod(Ut, u) + X %*% beta)
	if (GiveFitted==1) { 
		return( matrix(yHat,ncol=nlevU) ) }
	if (GiveResiduals==1) { 
		return( matrix(y - yHat,ncol=nlevU) ) } 
	#Set Penalized Residual Sum of Square differences
	prss <- sum(c(y - yHat, as.vector(u))^2)  
	#Set number of samples and number of parameters
	n <- length(y); p <- ncol(RX)
	#Conditional Estimate for the variances:
	if( GiveSigma == 1) {
		if(nlevU==1) { return (sqrt(prss/(n-p)) ) }
		else
		{		
		SigmaVector <- rep(0, nlevU +1)
		for (i in 1:nlevU){
			x1 <- 1 + (i-1) * (n/nlevU);
			x2 <- i * (n/nlevU)
			x3 <- 1 + (i-1) * nlevL1
			x4 <- i * nlevL1
			J <- sum( (c( y[x1:x2] - yHat[x1:x2], u[x3:x3]))^2)
			SigmaVector[i] = (  sqrt((nlevU*J)/(n-p)) ) # if you don't multiply by nlevU you need to change n & p
		}
		SigmaVector[i+1] = sqrt(prss/(n-p)) #this is the "fake" cross sample variance
		return (SigmaVector) }}
 
	if (GiveBeta==1) { 
		Betas <- matrix(beta, ncol=nlevU)
		if(nlevU==1) {  
			SE_of_Betas <- matrix( sqrt(diag(solve(crossprod(RX))))*sqrt(prss/(n-p)),ncol=nlevU)
			return(list(Betas,SE_of_Betas))
		}
		else 
		{
		SigmaVector <- rep(0, nlevU )
		for (i in 1:nlevU){
			x1 <- 1 + (i-1) * (n/nlevU);
			x2 <- i * (n/nlevU)
			x3 <- 1 + (i-1) * nlevL1
			x4 <- i * nlevL1
			J <- sum( (c( y[x1:x2] - yHat[x1:x2], u[x3:x3]))^2)
			SigmaVector[i] = (  sqrt((nlevU*J)/(n-p)) ) # if you don't multiply by nlevU you need to change n & p
		}
		
		
		SEB<- matrix( sqrt(diag(solve(crossprod(RX))))* rep( SigmaVector, each= (p/nlevU)) ,ncol=nlevU)
		return( list(Betas, SEB))
		}	 
	}
 
	if( (GiveCorrelations == 1) || ( GiveCovariances == 1) ) {  
		SigmaE <- sqrt(prss/(n-p))			#"fake" cross sample variance
		VCV1 <- tcrossprod( SigmaE *  Kl1)		#Variance covariance matrix 
		if( nlevL2!=0) {
			VCV2 <- tcrossprod( SigmaE *  Kl2)	#Variance covariance matrix 
			if(  GiveCovariances == 1) { return (list(VCV1, VCV2)) } }
		else{
			if(  GiveCovariances == 1) { return(list(VCV1))} }		
		CorMat1  <- cov2cor(VCV1)			#Correlation matrix
		if( nlevL2!=0) {
			CorMat2  <- cov2cor(VCV2)		#Correlation matrix
			return( list(CorMat1, CorMat2)) }
		else{
			return( list(CorMat1))} }

	if( ML ==1) { 
		ML_value = 2*determinant(L)$mod  +(n)*(1+log(2*pi*prss/(n)));  if(GarCol==1) {	gc(); }
		return( ML_value ) }

	REML_value = 2*determinant(L)$mod+ 2* determinant(RX)$mod +(n-p)*(1+log(2*pi*prss/(n-p)));  if(GarCol==1) {	gc(); }
	return( REML_value )
}
