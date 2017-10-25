KMTest.surv = function(Z=NULL, U, Delta, X = NULL, gamma = NULL, kernel, npert = 1000){
  #=======================================================================================================================#
  #=======================================================================================================================#
  # Version 3: 10.10.2015                                                                                                 #
  # Changes in Version 3                                                                                                  #
  # Minor changes to increase efficiency, shd give identical results as V2                                                #
  #                                                                                                                       #
  # Version 2: 09.30.2013                                                                                                 #
  # Changes in Version 2:                                                                                                 #
  # doesn't check if Z, U, Delta, X etc are matrices etc                                                                  #
  # kernel must be an n x n kernel matrix that is supplied                                                                #
  # Z is not needed anymore, and removed the weights argument                                                             #
  #                                                                                                                       #
  # Notes: This function tests for association between a set of SNPS and a censored survival outcome                      #
  #        using the survival kernel machine.  The IBS and weighted kernels assume that the SNPs are in additive mode,    # 
  #        but can easily be adapted for SNPs in dominant mode by either defining new kernel functions or entering the    #
  #        nXn kernel matrix directly.                                                                                    #
  #        The function also assumes that survival times are right-censored.                                              # 
  #        The function also assumes that there are no missing data, individuals with missing data should be removed      #
  #        prior to analysis.                                                                                             #
  #                                                                                                                       #
  #                                                                                                                       #
  # :::INPUT:::                                                                                                           #
  #                                                                                                                       #
  # X:      is a nXR  matrix of relevant COVARIATES with each row as a different individual and each column               #
  #         as a separate covariate measurement. If no additional covariates are present, X can be left unspecified or    #
  #         left as NULL. Note that each column of X has to be a numerical variable, non-numerical variables have to be   #
  #         recoded appropriately before analysis.                                                                        #                         
  # Z:      is a nXS  matrix of the relevant SNPS with each row as a different individual and each column                 #
  #         as a separate SNP.  We asume the genotypes are in additive mode: 0,1,2 for non-carriers, heterozygotes        #
  #         and homozygous carriers, respectively.                                                                        #
  # U:      is a nX1 vector containing the observed times. Note: U=min(C,T) where C = censoring time, T = survival time   #
  # Delta:  is a nX1 vector containing the event indicator.                                                               #
  # gamma: Unless X = NULL, gamma has to be supplied.     gamma <- coxph(Surv(U,Delta)~X)$coef (?????)                    #                         
  # kernel: defines the kernel matrix.  This may be defined in several ways:                                              #
  #         kernel can be the nXn Kernel matrix OR                                                                        #
  #         kernel can be:                                                                                                #
  #                       "linear":  Linear Kernel                                                                        #
  #                       "IBS": IBS Kernel                                                                               #
  #                       "IBS.weighted": weighted IBS (weights must be specified a priori)                               #
  # weights: is a vector of length S*n of prespecified weights for the weighted IBS Kernel                                #
  #          See Kern.FUN.weightedIBS() to for information on how weights are defined                                     #
  # npert:   is the number of perturbations used to calculate p-value (default =1000)                                     #
  #                                                                                                                       #
  #                                                                                                                       #                                                                                                                                                                                                                                    
  # :::OUTPUT:::                                                                                                          #
  #                                                                                                                       # 
  # n : no. of samples                                                                                                    #
  # nsnps: no. of SNPs = S                                                                                                #   
  # p.value.ptb: single p-value for testing the null of no association based on resampling                                #
  # pvalue.chisq: single p-value for testing the null of no association based on Satterthwaite approximation              #
  #               Note that the Satterthwaite approximation uses resampling to estimate the scale parameter and df,       #
  #               and changes with each run. Use set.seed() to obtain identical p-values.                                 #
  # df.chisq: the degrees of freedom of the chi-square statistic                                                          #
  # score : the unscaled score statistic                                                                                  #
  # wptb: uncentered realizations of score statistic under null                                                           #                                                                                                                      
  # Citation:                                                                                                             #
  #                                                                                                                       #
  #                                                                                                                       #
  # If you use the survival kernel machine, please cite the following two papers:                                         #
  #                                                                                                                       #
  # Cai T, Tonini G and Lin X. 2011. Kernel machine approach to testing the significance of genomic pathways.             #
  # Biometrics 67: 975-86.                                                                                                #
  #                                                                                                                       #
  # Lin X, Cai T, Wu M, Zhou Q, Liu G, Christiani D, Lin X. 2011. Kernel Machine SNP-set Analysis for Censored            #
  # Survival Outcomes in Genome-wide Association Studies. Genetic Epidemiology 35: 620-31.                                #
  #                                                                                                                       #
  #                                                                                                                       #
  #                                                                                                                       #
  #=======================================================================================================================#
  #=======================================================================================================================#
	
      #if (!is.null(X)) {
   	  #   if (class(X)!= "matrix") stop("X is not a matrix")
	  #   if (nrow(X)!=nrow(Z)) stop("Dimensions of Z and X do not match")
      #}
  
      #if (class(Z)!= "matrix") stop("Z is not a matrix")
      #if (nrow(Z)!=length(U)) stop("Dimensions of Z and U do not match")
      #if (nrow(Z)!=length(Delta)) stop("Dimensions of Z and Delta do not match")


	# ====================================== #
	# Construct Kernel Matrix 
	# ====================================== #
	

		Kij = kernel


	nn = length(U); Gi = matrix(rnorm(npert*nn),nrow=nn)
	
	tl = U[Delta==1] # failure time
	stl = sort(tl); nl = length(tl) 
	if (is.null(X)) V = rep(1,nn) else V = as.vector(exp(X%*%gamma))
	pi0hat.stl = sum.I(stl,"<=", U, V) ## nl X 1
	dLam.stl <- 1/pi0hat.stl; Lam.stl <- cumsum(dLam.stl) ## nl X 1   
	# For each subject, tmpind is the number of the failure time is smaller than the time
    tmpind <- sum.I(U,">=",stl)+1; Lam.U <- c(0,Lam.stl)[tmpind] # sum{I(Xi>=stl)dLam.stl} n X 1
    Mhati <- Delta - Lam.U*V ## n X 1
    

    
    if (!is.null(X)){
    	q0 = ncol(X)
    	# S^(1), S^(1) S'^(1), S^(2)
    	pi1hat.stl <- sum.I(stl,"<=",U,X*V)
    	pi1hat2.stl <- matrix(0,length(stl),q0^2)
    	for (k in 1:length(stl)) {pi1hat2.stl[k,] <- as.vector(pi1hat.stl[k,]%*%t(pi1hat.stl[k,]))}
    	X2.mat <- matrix(0,nn,q0^2)
    	for (i in 1:nn){ X2.mat[i,] <- as.vector(X[i,]%*%t(X[i,]))}
    	pi2hat.stl = sum.I(stl,"<=",U,X2.mat*V)
    	Ahat = matrix(apply(pi2hat.stl/pi0hat.stl-pi1hat2.stl/pi0hat.stl^2,2,sum),q0,q0,byrow=T)/nn

	   #-----------------------------------------------
	   # change in V3 - save invAhat
       invAhat <- solve(Ahat)
	   #-----------------------------------------------
	   	    	
       W.gamma = W.gamma.t1 = W.gamma.t2 = W.gamma.t3 = matrix(0,nn,q0) 
       W.gamma.t1 = X*Mhati
       W.gamma.t2[Delta==1,] = pi1hat.stl[rank(tl),]/pi0hat.stl[rank(tl)]
       int.t3 = apply(pi1hat.stl/(pi0hat.stl)^2,2,cumsum)
       tmpind <- sum.I(U,">=",stl)+1
       W.gamma.t3 = rbind(rep(0,q0),int.t3)[tmpind,]*V
       #-----------------------------------------------
	   # change in V3 - save invAhat
	   W.gamma <- t( invAhat%*%t(W.gamma.t1-W.gamma.t2+W.gamma.t3) )
       #W.gamma = t(solve(Ahat)%*%t(W.gamma.t1-W.gamma.t2+W.gamma.t3))   
       #-----------------------------------------------
       X.tilde <- X*(Lam.U*V) 
     	}

	#-----------------------------------------------     	
    # change in V3 - moved this block below 	
    term1.ii <- c(0,cumsum(dLam.stl/pi0hat.stl))[tmpind]
    term2.ij <- pmin(matrix(U,nrow=nn,ncol=nn),matrix(U,nrow=nn,ncol=nn,byrow=T))
    term2.ij <- matrix(sum.I(c(term2.ij),">=",stl,dLam.stl/pi0hat.stl),ncol=nn)
    Mhati.Mhatl = Mhati%*%t(Mhati)
	#----------------------------------------------- 
	   
    score = 0; wptb = rep(0,npert)
    term1 <- sum(diag(Kij)*(diag(Mhati.Mhatl)- (Lam.U*V - term1.ii*V^2)))
    term2 <- sum((Kij-diag(diag(Kij)))*(term2.ij*(V%*%t(V))+Mhati.Mhatl))
    score <- term1+term2
        
    Keig<-eigen(Kij,symmetric=T); Keig.value<-pmax(Keig$values,0); Keig.vec <- Keig$vectors # n times m
    ri.mat1 = Keig.vec*Mhati
    omega.stl = sum.I(stl,"<=",U,Keig.vec*V)/pi0hat.stl #nt x m matrix 
    ri.mat2 = Keig.vec*0 # n x m matrix
    ri.mat2[Delta==1,] = omega.stl[rank(tl),] 
    ri.mat2 = ri.mat2 - sum.I(U, ">=", stl, omega.stl*dLam.stl)*V # n \times m     
    ri.mat = ri.mat1 - ri.mat2
    if (!is.null(X)){
    	ri.mat2.3 = t(pi1hat.stl%*%t(W.gamma)*dLam.stl)%*%omega.stl/nn
       ri.mat3 = W.gamma%*%t(X.tilde)%*%Keig.vec/nn
       ri.mat = ri.mat + ri.mat2.3 - ri.mat3
    }
    out = VTM(sqrt(Keig.value),nn)*ri.mat
    wptb = apply((t(out)%*%Gi)^2,2,sum)
 	
 	## p value obtained from perturbation
 	wptb.c = wptb - mean(wptb)
 	p.value.ptb = mean((wptb.c-score)>0)
 	
 	## p value obtained from chi-square approximation
 	wptb.mu = mean(wptb); wptb.var = var(wptb)
 	chisq.scale = wptb.var/wptb.mu/2
 	chisq.df = wptb.mu/chisq.scale
 	p.value.chisq = 1-pchisq((score+wptb.mu)/chisq.scale,df=chisq.df)
 	
 	return(list(n=nn,nsnps=ncol(Z), shat=score, sptb = wptb, pvalue.ptb=p.value.ptb,pvalue.Q=p.value.chisq, df.chisq=chisq.df))     
}

sum.I <- function(yy,FUN,Yi,Vi=NULL)
{
    if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
    # for each distinct ordered failure time t[j], number of Xi < t[j]
    pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')    
    if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
    if (!is.null(Vi)) {
       ## if FUN contains '=', tmpind is the order of decending
        if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
        ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
        Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
        return(rbind(0,Vi)[pos+1,])
    } else return(pos)
}

Kern.FUN.Lin <- function(zz,weights=NULL)
{
    zz = as.matrix(zz);
    zz%*%t(zz)
}

# Kern.FUN.IBS modified April 28, 2009
Kern.FUN.IBS <- function(zz,weights=NULL)
  {
    # computes the average IBS between two individuals 
    # sum over IBS at all markers and divide by total no. of markers
    # only non-missing markers used (missing value coded as NA)
    # missing data is not included in both the numerator and denominator
    # max possible value is 1 if IBS=2 at ALL non-missing markers
    temp <- zz
    Ina <- 1*(is.na(temp)==F)
    zz[is.na(zz)] <- -9
    I0 = 1*(zz==0); I1 = 1*(zz==1); I2 = 1*(zz==2); 
    maxsnps <- Ina%*%t(Ina)
    (I0%*%t(I0)+I1%*%t(I1)+I2%*%t(I2)+0.5*I1%*%t(I0+I2)+0.5*(I0+I2)%*%t(I1))/max(maxsnps,0.00001)
  }

# Kern.FUN.weightedIBS modified June 28, 2010
Kern.FUN.weightedIBS <- function(zz, weights)
  {
    # computes the average IBS between two individuals 
    # weights is a vector of length p*n (not a matrix)
    # if there are S SNPs and n samples and we are weighting inversely by sqrt(MAF),
    # weights is of the form c(rep(MAF1, nsamples), rep(MAF2, nsamples), ... , rep(MAFp, nsamples))
    # sum over weighted IBS at all markers and divide by total no. of markers
    # input matrix CANNOT contain missing genotypes
    # max possible value is 1 if IBS=2 at ALL markers

    I0 = 1*(zz==0); I1 = 1*(zz==1); I2 = 1*(zz==2); 
    I0 = I0*(1/weights)^0.25; I1 = I1*(1/weights)^0.25; I2 = I2*(1/weights)^0.25;
    I0%*%t(I0)+I1%*%t(I1)+I2%*%t(I2)+0.5*I1%*%t(I0+I2)+0.5*(I0+I2)%*%t(I1)
  }


# Date: June 5, 2009
mafcall <- function(x){
	min(sum(x, na.rm=T), 2*sum(is.na(x)==F)-sum(x, na.rm=T))/(2*sum(is.na(x)==F)) 
}

# Date: June 5, 2009
checkpolymorphic <- function(x){
	length(unique(x))>1
}

VTM<-function(vc, dm){
    matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}
