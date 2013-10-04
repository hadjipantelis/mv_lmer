rm(list=ls()); 
library(lme4)
#==================
#Univariate Example
#================== 
	print("                                                                                                              ")
	print("                                                                                                              ")
	print("****Univariate Example****")
	print("(we effectively replicate all of lme4's lmer() results)")
#Initialization of basic constants:
	N = 4000; seed_=1; nlev1 = 50; nlev2 = 25; 
	set.seed(seed_); 

#Error
	E1 = rnorm(N); 
	E1 = scale(E1) * 1;

#Fixed Effects
	X1 = rnorm(N, mean=10, sd=2);  
	X2 = rt(n=N, df=3)+3; 
	X3 = rnorm(N, mean=-2, sd=8);
	X4 = rt(n=N, df=2)+2; 
	X5 = runif(n=N, max= 10, min=-10) ; 
	X = matrix( c(X1,X2, X3, X4, X5), nrow=N);  B = matrix(c(70,10,12,03,05), nrow=1);

#Random Effects
	Class1 =  sample(1:nlev1, N, replace=T) 	#Not perfectly balanced design
	Zsingle1 <- model.matrix( ~-1 + as.factor(  Class1 ) ) 
	G1 =  scale(rnorm(nlev1))*4

	Class2 =  sample(1:nlev2, N, replace=T) 	#Not perfectly balanced design
	Zsingle2 <- model.matrix( ~-1 + as.factor(  Class2 ) ) 
	G2 =  scale(rnorm(nlev2))
	
#Generative model
	A = X %*% t(B) + Zsingle1 %*% (G1) + Zsingle2 %*% (G2) + E1
	 
#Regression using LMER
	MV_Class1  <- as.factor(rep(Class1,1))
	MV_Class2  <- as.factor(rep(Class2,1))
	print( "Optimizing LMER....")
	MV_lmer_A=  lmer(as.vector( matrix(A, ncol=1) ) ~ -1 + (X) + (1|MV_Class1)   + (1|MV_Class2)  );
	print( "Optimization finished....")	

#LMER results	
	Q <-summary(MV_lmer_A)

#Load implementation
	source("PrfDvnce_Mask_GaCo.R")

#Construct R.E. Design Matrix for custom evaluation
	Z12 = matrix(0, nrow = N, ncol = (nlev2+ nlev1))
	Z12[ ,1:nlev1] = Zsingle1
	Z12[ , (1+nlev1) : (nlev2+ nlev1)] = Zsingle2 
	Z12 = Matrix(Z12, sparse=T)

#Compare Results if one had the same theta:
	#REML final value
	print( "==============================================================================================================")
	print( "REML value if one used the same theta as the ones from lme4")
	print( paste( "lme4 result  :",  as.numeric( Q$devcomp$cmp["REML"] ) , sep= " ")) 
	print( paste( "MVLMER result:", as.numeric( PrfDvnce_Mask_GaCo(theta=MV_lmer_A@theta, X= X, XtX= crossprod(X), Zt= t(Z12), y= matrix(A, ncol=1), nlevU=1,  nlevL1 = nlev1, nlevL2=nlev2) )))
	#Betas
	print( "==============================================================================================================")
	print( "Betas if one used the same theta as the ones from lme4")
	print( Q$coefficients )
	HatBeta= matrix(rep(0, 2*5), ncol=2)
	HatBeta[,1] = PrfDvnce_Mask_GaCo(theta=MV_lmer_A@theta, X= X, XtX= crossprod(X), Zt= t(Z12), y= matrix(A, ncol=1), nlevU=1,  nlevL1 = nlev1, nlevL2=nlev2, GiveBeta=1) [[1]]
	HatBeta[,2] = PrfDvnce_Mask_GaCo(theta=MV_lmer_A@theta, X= X, XtX= crossprod(X), Zt= t(Z12), y= matrix(A, ncol=1), nlevU=1,  nlevL1 = nlev1, nlevL2=nlev2, GiveBeta=1) [[2]]
	print(HatBeta)
	#Sigma_E
	print( "==============================================================================================================")
	print( "Sigma if one used the same theta as the ones from lme4")
	print( paste( "lme4 result  :",  as.numeric(Q$sigma), sep= " ")) 
	print( paste( "MVLMER result:", as.numeric( PrfDvnce_Mask_GaCo(theta=MV_lmer_A@theta, X= X, XtX= crossprod(X), Zt= t(Z12), y= matrix(A, ncol=1), nlevU=1,  nlevL1 = nlev1, nlevL2=nlev2, GiveSigma=1) )))
	#Sigma_G
	print( "==============================================================================================================")
	print( "Std.Dev of the random effects if one used the same theta as the ones from lme4")
	Vars=   c( as.numeric(attr( Q$varcor$MV_Class1, "stddev")), as.numeric(attr( Q$varcor$MV_Class2, "stddev")))
 	print( do.call (paste, as.list( unlist( list( "lme4 result  :", Vars)))))
 	Vars=   sqrt( unlist(PrfDvnce_Mask_GaCo(theta=MV_lmer_A@theta, X= X, XtX= crossprod(X), Zt= t(Z12), y= matrix(A, ncol=1), nlevU=1,  nlevL1 = nlev1, nlevL2=nlev2, GiveCovariances=1)))
	print(do.call (paste, as.list( unlist( list( "MVLMER result:", Vars)))) )

#Optimize for current implementation using off-the-shelf Simplex Solver
	print( "Optimizing MVLMER....")
	OPT <- optim( PrfDvnce_Mask_GaCo, par=c(1,1) , X= X, XtX= crossprod(X), Zt= t(Z12), y= matrix(A, ncol=1), nlevU=1,  nlevL1 = nlev1, nlevL2=nlev2)
	print( "Optimization finished....")	

	print( "==============================================================================================================")
	print( "REML value from MVLMER N-M optimization")
	print( paste( "MVLMER result:", as.numeric( OPT$value)))
	print( "==============================================================================================================")
	print( "Betas value from MVLMER N-M optimization")
	HatBeta= matrix(rep(0, 2*5), ncol=2)
	HatBeta[,1] = PrfDvnce_Mask_GaCo(theta= OPT$par, X= X, XtX= crossprod(X), Zt= t(Z12), y= matrix(A, ncol=1), nlevU=1,  nlevL1 = nlev1, nlevL2=nlev2, GiveBeta=1) [[1]]
	HatBeta[,2] = PrfDvnce_Mask_GaCo(theta= OPT$par, X= X, XtX= crossprod(X), Zt= t(Z12), y= matrix(A, ncol=1), nlevU=1,  nlevL1 = nlev1, nlevL2=nlev2, GiveBeta=1) [[2]]
	print(HatBeta)
	print( "==============================================================================================================")
	print( "Sigma value from MVLMER N-M optimization")
	print( paste( "MVLMER result:", as.numeric( PrfDvnce_Mask_GaCo(theta= OPT$par, X= X, XtX= crossprod(X), Zt= t(Z12), y= matrix(A, ncol=1), nlevU=1,  nlevL1 = nlev1, nlevL2=nlev2, GiveSigma=1) )))
	print( "==============================================================================================================")
	print( "Std.Dev of the random effects values from MVLMER N-M optimization")
	Vars=   sqrt( unlist(PrfDvnce_Mask_GaCo(theta= OPT$par, X= X, XtX= crossprod(X), Zt= t(Z12), y= matrix(A, ncol=1), nlevU=1,  nlevL1 = nlev1, nlevL2=nlev2, GiveCovariances=1)))
	print( do.call (paste, as.list( unlist( list( "MVLMER result:", Vars)))) )


#==================
#Multivariate Example
#==================
rm(list=ls()); 
library(lme4)
	print("                                                                                                              ")
	print("                                                                                                              ")
	print("****Multivariate Example****")
#Initialization of basic constants:
	N = 4000; seed_=1; nlev1 = 50; nlev2 = 25; 
	set.seed(seed_);  
	p = 4; #Dimensionality
	
	correlation1_ =0.00; 
	correlation2_ =0.20; correlation3_ =0.15; 
	correlation4_ =0.20; correlation5_ =0.25; correlation6_ =0; 
	
	customcorr = 1; # Use custom correlation structure

#Errors
	E1 = rnorm(N); E1 = scale(E1) * 1;
	E2 = rnorm(N); E2 = scale(E2) * 2;
	E3 = rnorm(N); E3 = scale(E3) * 3; 
	E4 = rnorm(N); E4 = scale(E4) * 1; 
	E = (matrix( c( E1, E2 , E3, E4), nrow= N))
	P <- diag(p) ; 				
	#Correlated Error Structure (not strong)
	if ( customcorr == 1 ) {
		P <- diag(p)* .99 + .01
		P[1,2] = .05; P[2,1] = P[1,2]; 
		P[1,3] = .10; P[3,1] = P[1,3]; 
		P[1,4] = .15; P[4,1] = P[1,4];
		P[2,4] = .10; P[4,2] = P[2,4];
		P[3,4] = .05; P[4,3] = P[3,4]; 
	}
	E = t( t(chol(P)) %*% t(E)) 
	
#Fixed Effects
	X1 = rnorm(N, mean=10, sd=2);  
	X2 = rt(n=N, df=3)+3; 
	X3 = rnorm(N, mean=-2, sd=8);
	X4 = rt(n=N, df=2)+2; 
	X5 = runif(n=N, max= 10, min=-10) ; 
	X = matrix( c(X1,X2, X3, X4, X5), nrow=N);  B = matrix(c(70,10,2,.5,6,4,7,9.5,19, 6,6,60, 8,16,24,1,2,3,4,5), nrow=5)
 
#Random Effects
	Class1 = sample(1:nlev1, N, replace=T)  	#Random Nvel Classes
	Zsingle1 <- model.matrix( ~-1 + as.factor(  Class1 ) ) 

	Class2 = sample(1:nlev2, N, replace=T ) 	#Random Nvel Classes
	Zsingle2 <- model.matrix( ~-1 + as.factor(  Class2 ) ) 

	P <- diag(p) ; 				#R.E. Correlation structure
	P[1,2] =correlation1_; P[2,1] = P[1,2];
	P[1,3] =correlation2_; P[3,1] = P[1,3];
	P[1,4] =correlation3_; P[4,1] = P[1,4];
	P[2,3] =correlation4_; P[3,2] = P[2,3];
	P[2,4] =correlation5_; P[4,2] = P[2,4];
	P[3,4] =correlation6_; P[4,3] = P[3,4];

	#Make Mask Matrix: 
	Mask <- matrix( nrow=4, as.double( !(P == 0)))

	D1 <- diag(c(4,6,8,10))			#Deviation/Variance magnitude
	K1 <- D1 %*% (Mask * P) %*% D1		#Covariance structure by scaling the masked correlation matrix
	R1 <- matrix( rnorm(nlev1* p), nrow=p)	#Random Normal Matrix
	#Make R perfect
	R1 <- t(scale(t(R1)))
	#Make Gamma
	Gamma1 = t(chol(K1)) %*% R1;  
	
	D2 <- diag(c(5,3,1,4))			#Deviation/Variance magnitude
	K2 <- D2 %*% (Mask * P) %*% D2		#Covariance structure by scaling the masked correlation matrix
	R2 <- matrix( rnorm(nlev2* p), nrow=p)	#Random Normal Matrix
	#Make R perfect
	R2 <- t(scale(t(R2)))
	#Make Gamma
	Gamma2 = t(chol(K2)) %*% R2;  

#Generative model
	A = X %*% B + Zsingle1 %*% t(Gamma1) + Zsingle2 %*% t(Gamma2) + (E)

#Load implementation
	source("PrfDvnce_Mask_GaCo.R")

#Construct Design Matrices for custom evaluation
	mv_x =  kronecker(Diagonal(p), X)
	mv_xxt = crossprod( mv_x ) 
	mv_z1 <- kronecker(Diagonal(p),  Zsingle1)
	mv_z2 <- kronecker(Diagonal(p),  Zsingle2)
	mv_z <-  kronecker(Diagonal(p), matrix( 0, nrow= N, ncol= nlev1+nlev2))
	mv_z[, 1:(p*nlev1)] <- mv_z1
	mv_z[, (1+(p*nlev1) ): (p*(nlev2+nlev1))] <- mv_z2


if(1==3){
#Invalid Regression using LMER 
	Indicator <- as.factor(rep(c(1,2,3,4), each=N) );  
	MV_Class1  <- as.factor(rep(Class1,p))
	MV_Class2  <- as.factor(rep(Class2,p))
	print( "Optimizing LMER.... (this may take a moment)")
	MV_lmer_B =  lmer(as.vector( matrix(A, ncol=1) ) ~ -1 + kronecker(diag(p),X) + ( 0+Indicator|MV_Class2)   + ( 0+Indicator|MV_Class1)  );
	print( "Optimization finished....")	

#LMER results
	Q <-summary(MV_lmer_B)

#Compare Results if one had the same theta:
	#REML final value
	print( "==============================================================================================================")
	print( "Print REML value if one used the same theta as the ones from lme4")
	print( paste( "lme4 result  :",  as.numeric( Q$devcomp$cmp["REML"] ) , sep= " ")) 
	print( paste( "MVLMER result:", as.numeric( PrfDvnce_Mask_GaCo(theta=MV_lmer_B@theta, X= mv_x, XtX=mv_xxt, Zt= t(mv_z), y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2) )))
	
	#Betas / #The standard errors SHOULD NOT be the same
	print( "==============================================================================================================")
	print( "Print Betas if one used the same theta as the ones from lme4; the Std.Error are not the same.")
	print( Q$coefficients )
	HatBeta= matrix(rep(0, p*2*5), ncol=2)
	HatBeta[,1] = matrix( ncol=1, PrfDvnce_Mask_GaCo(theta=MV_lmer_B@theta, X= mv_x, XtX=mv_xxt, Zt= t(mv_z), y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, GiveBeta=1) [[1]])
	HatBeta[,2] = matrix( ncol=1,  PrfDvnce_Mask_GaCo(theta=MV_lmer_B@theta, X= mv_x, XtX=mv_xxt, Zt= t(mv_z), y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, GiveBeta=1) [[2]])
	print(HatBeta)
	
	#Sigma_E /#Only the "sample-wide" sigma should be the same (last entry)
	print( "==============================================================================================================")
	print( "Print Sigma(s) if one used the same theta as the ones from lme4; MVLMER returns a sigma for each dimension, last entry is the \"sample\"-wide sigma ")
	print( paste( "lme4 result  :",  as.numeric(Q$sigma), sep= " "))
	Vars = as.numeric( PrfDvnce_Mask_GaCo(theta=MV_lmer_B@theta,X= mv_x, XtX=mv_xxt, Zt= t(mv_z) , y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, GiveSigma=1))
	print( do.call (paste, as.list( unlist( list( "MVLMER result:", Vars)))) ) 
	
	#Sigma_G /#Correlations and StdDevs
	print( "==============================================================================================================")
	print( "Print VCV-structure for each Rand.Effect if one used the same theta as the ones from lme4")
	print( "Correlations first and then Standard Deviations")

	print( "====================================================")
	print( "Correlation for Random Effect 1:")
	print( "LMER:")
	print( round( digits=3, attr( Q$varcor$MV_Class1, "correlation")))
	print( "MVLMER:")
	print( round( digits=3, PrfDvnce_Mask_GaCo(theta=MV_lmer_B@theta, X= mv_x, XtX=mv_xxt, Zt= t(mv_z), y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, GiveCorrelations=1)[[1]]))

	print( "Correlation for Random Effect 2:")
	print( "LMER:")
	print( round( digits=3, attr( Q$varcor$MV_Class2, "correlation")))
	print( "MVLMER:")
	print( round( digits=3, PrfDvnce_Mask_GaCo(theta=MV_lmer_B@theta, X= mv_x, XtX=mv_xxt, Zt= t(mv_z), y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, GiveCorrelations=1)[[2]]))
	
	print( "====================================================")
	print( "Std.Devs for Random Effect 1:")
	print( "LMER:")
	print( as.numeric( attr( Q$varcor$MV_Class1, "stddev")))
	print( "MVLMER:")
	print( sqrt(diag(PrfDvnce_Mask_GaCo(theta=MV_lmer_B@theta, X= mv_x, XtX=mv_xxt, Zt= t(mv_z), y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, GiveCovariances=1)[[1]])) )
	print( "Std.Devs for Random Effect 2:")	
	print( "LMER:")
	print( as.numeric(attr( Q$varcor$MV_Class2, "stddev")) )
	print( "MVLMER:")
	print( sqrt(diag(PrfDvnce_Mask_GaCo(theta=MV_lmer_B@theta, X= mv_x, XtX=mv_xxt, Zt= t(mv_z), y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, GiveCovariances=1)[[2]])) )
}
#Optimize for current implementation using off-the-shelf Simplex Solver	
	#We now "know" the approximate covariance between the random effects.
	print( "==============================================================================================================")
	print( "Optimizing MVLMER....")
	OPT <- optim( PrfDvnce_Mask_GaCo, par=rep(1,20), X= mv_x, XtX=mv_xxt, Zt= t(mv_z), y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, Mask=Mask)
	OPT2 <- optim( PrfDvnce_Mask_GaCo, par=OPT$par , X= mv_x, XtX=mv_xxt, Zt= t(mv_z), y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, Mask=Mask, method="BFGS")
	OPT <- optim( PrfDvnce_Mask_GaCo, par=OPT2$par , X= mv_x, XtX=mv_xxt, Zt= t(mv_z), y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, Mask=Mask )
	print( "Optimization finished....")
	print( "REML value from MVLMER N-M/BFGS optimization, one can use the Mask information here")
	print( paste( "MVLMER result:", as.numeric( OPT$value)))
	print( "==============================================================================================================")
	print( "Betas value from MVLMER N-M/BFGS optimization")
	HatBeta= matrix(rep(0, p*2*5), ncol=2)
	HatBeta[,1] = PrfDvnce_Mask_GaCo(theta= OPT$par, X= mv_x, XtX=mv_xxt, Zt= t(mv_z),  y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, GiveBeta=1, Mask=Mask) [[1]]
	HatBeta[,2] = PrfDvnce_Mask_GaCo(theta= OPT$par, X= mv_x, XtX=mv_xxt, Zt= t(mv_z),  y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, GiveBeta=1, Mask=Mask) [[2]]
	print(HatBeta)
	print( "==============================================================================================================")
	print( "Sigma value from MVLMER N-M/BFGS optimization")
	Vars = as.numeric( PrfDvnce_Mask_GaCo(theta=OPT$par,X= mv_x, XtX=mv_xxt, Zt= t(mv_z) , y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, GiveSigma=1, Mask=Mask))
	print( do.call (paste, as.list( unlist( list( "MVLMER result:", Vars)))) ) 
	print( "==============================================================================================================")
	print( "Print VCV-structure for each Rand.Effect from MVLMER N-M/BFGS optimization")
	print( "Correlations first and then Standard Deviations")

	print( "====================================================")
	print( "Correlation for Random Effect 1:")
	print( round( digits=3, PrfDvnce_Mask_GaCo(theta= OPT$par, X= mv_x, XtX=mv_xxt, Zt= t(mv_z), y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, GiveCorrelations=1, Mask=Mask)[[1]]))
	print( "Correlation for Random Effect 2:")
	print( round( digits=3, PrfDvnce_Mask_GaCo(theta= OPT$par, X= mv_x, XtX=mv_xxt, Zt= t(mv_z), y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, GiveCorrelations=1, Mask=Mask)[[2]]))
	
	print( "====================================================")
	print( "Std.Devs for Random Effect 1:")
	print( sqrt(diag(PrfDvnce_Mask_GaCo(theta= OPT$par, X= mv_x, XtX=mv_xxt, Zt= t(mv_z), y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, GiveCovariances=1, Mask=Mask)[[1]])) )
	print( "Std.Devs for Random Effect 2:")	
	print( sqrt(diag(PrfDvnce_Mask_GaCo(theta= OPT$par, X= mv_x, XtX=mv_xxt, Zt= t(mv_z), y= matrix(A, ncol=1), nlevU=p,  nlevL1 = nlev1, nlevL2=nlev2, GiveCovariances=1, Mask=Mask)[[2]])) )
