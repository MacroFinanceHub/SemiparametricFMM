##########################################################################
####  Get Spline Z matrix
####  First N rows are Z matrix for data including N observations
####	Next ngrid rows are for evaluating spline on grid of size ngrid of x
####  with endpoints a and b
##########################################################################

GetSpline<-function(x,numIntKnots=5,ngrid=71,a=20,b=90){

	xg <- seq(a,b,length=ngrid)

	library(splines)
	intKnots <- quantile(unique(x),seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))])

	names(intKnots) <- NULL
	B <- bs(x,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE)
	BTB <- crossprod(B,B) ; 
	cat('intKnots=',intKnots,'\n')

	#### Derivative of B-spline basis
	formBprime <- function(x,a,b,intKnots)
	{
   		allKnots <- c(rep(a,4),intKnots,rep(b,4)) 
   		K <- length(intKnots) ; L <- 3*(K+8)
   		Bd <- spline.des(allKnots,x,derivs=rep(1,length(x)),outer.ok=TRUE)$design     
   		return(Bd)
	}

	formOmega <- function(a,b,intKnots)
	{
   		allKnots <- c(rep(a,4),intKnots,rep(b,4)) 
   		K <- length(intKnots) ; L <- 3*(K+8)
   		xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+rep(allKnots,each=3)[-c(1,2,L)])/2
   		wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
   		Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),outer.ok=TRUE)$design  
   		Omega     <- t(Bdd*wts)%*%Bdd     
   		return(Omega)
	}

	Omega <- formOmega(a,b,intKnots)
	eigOmega <- eigen(Omega)

	# Obtain the matrix for linear transformation of $\\bB$ to $\\bZ$:
   
	indsZ <- 1:(numIntKnots+2)
	UZ <- eigOmega$vectors[,indsZ]
	LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))

	indsX <- (numIntKnots+3):(numIntKnots+4)
	UX <- eigOmega$vectors[,indsX]   
	L <- cbind( UX, LZ )
	stabCheck <- t(crossprod(L,t(crossprod(L,Omega))))          
	if (sum(stabCheck^2) > 1.0001*(numIntKnots+2))
    		print("WARNING: NUMERICAL INSTABILITY ARISING FROM SPECTRAL DECOMPOSITION")

	X <- cbind(rep(1,length(x)),x)
	Z <- B%*%LZ

	Bg <- bs(xg,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE)
	Zg <- Bg%*%LZ

	Bdg <- formBprime(xg,a,b,intKnots)
	Zdg <- Bdg%*%LZ

	Zres=rbind(Z,Zg,Zdg)

	return(list(Zres=Zres,B=B,Omega=Omega))
}
