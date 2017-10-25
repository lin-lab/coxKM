# Last Edited: Oct 2015

coxKM <- function(Z=NULL, U, Delta, X = NULL, gamma = NULL, kernel = "linear", weights = NULL, npert = 10^4, 
			npert.check=TRUE, npert.upper=10^8, npert.threshold=50,
			impute.method = "fixed", is_check_genotype=TRUE, is_dosage=FALSE, missing_cutoff=0.15, SetID=NULL){




	  #-------------------------------------
          # check dimensions when kernel matrix is not supplied
	  # if kernel matrix is supplied, Z is not required

	  if(class(kernel) != "matrix"){ 
	  	
	  	
	  #-------------------------------------
	  # only change in V0.3 is to move this under the "if(class(kernel) != "matrix"){ " statement
      # check that weights are supplied when kernel="IBS.weighted" or "linear.weighted"
	  if(is.null(weights)==TRUE & kernel=="IBS.weighted") stop("For weighted kernels, weights must be specified.")
      if(is.null(weights)==TRUE & kernel=="linear.weighted") stop("For weighted kernels, weights must be specified.")
	  #-------------------------------------
	  	
	  	if(class(Z)!= "matrix") stop("Z is not a matrix.")
        	if(nrow(Z)!=length(U)) stop("Dimensions of Z and U do not match")
        	if(nrow(Z)!=length(Delta)) stop("Dimensions of Z and Delta do not match")
	  

	 	if(!is.null(X)){
	   	   if(class(X)!= "matrix") stop("X is not a matrix.")
        	   if(nrow(X)!= nrow(Z)) stop("Dimensions of X and Z don't match.")
	  	}


	  	if(is.null(weights)==FALSE){
        		if(ncol(Z)!= length(weights)) stop("Dimensions of Z and weights don't match.")
			if(sum(weights<0)!=0) stop("weights have to be non-negative and cannot be missing.")
	  	}
	  }



	  #-------------------------------------
          # check dimensions when kernel matrix is supplied
	  # if kernel matrix is supplied, Z is not required

	  if(class(kernel) == "matrix"){ 
	  	if(nrow(kernel)!= ncol(kernel)) stop("kernel is not a square matrix.")
        	if(nrow(kernel)!=length(Delta)) stop("Dimensions of U and Delta do not match")
        	if(length(U)!=length(Delta)) stop("Dimensions of U and Delta do not match")
	  

	 	 if(!is.null(X)){
	   	   if(class(X)!= "matrix") stop("X is not a matrix.")
        		if(nrow(X)!= length(U)) stop("Dimensions of X and U don't match.")
	  	}
	  }



	  #-------------------------------------
	  # check for missing values in X, U, Delta

	  if(is.null(X)==FALSE){
        	if(sum(is.na(X))!= 0) stop("X cannot have any missing values.")
	  }

	  if(sum(is.na(U))!= 0) stop("U cannot have any missing values.")
	  if(sum(is.na(Delta))!= 0) stop("Delta cannot have any missing values.")



	  #-------------------------------------
	  # check that X doesn't have intercept
	  if(is.null(X)==FALSE){

		if(ncol(X)==1){
			if(checkpolymorphic(X)==FALSE) stop("X should not include intercept and must have more than one level.")	
		}else{
        		if(sum(apply(X, 2, checkpolymorphic))!= ncol(X)) stop("X should not include intercept and must have more than one level.")

		}
	  }



	#-------------------------------------
	# check Z and impute when kernel matrix is not supplied

	if(class(kernel) != "matrix"){ 
 	if(is_dosage ==TRUE){
		impute.method="fixed"
	}

	if(is_check_genotype==TRUE | is_dosage==TRUE){
		Z.out <- coxSKAT_MAIN_Check_Z(Z=Z, SetID=SetID, weights=weights, impute.method=impute.method,  missing_cutoff=missing_cutoff)
		Z <- as.matrix(Z.out$Z.test)
		weights <- as.vector(Z.out$weights.Z.test)
	}
	
	if(is.null(Z)==TRUE){ 

		if(is.null(SetID)){
			msg <- sprintf("The Z matrix has no SNPs." )
		} else {
			msg <- sprintf("In %s, the Z matrix has no SNPs.", SetID )
		}

		warning(msg,call.=FALSE)
	      return(list(p.value=NA, Q=NA, n.marker.test=NA, n.indiv=NA, df= NA))
	}



	}



	#-------------------------------------
	# Get the kernel matrix

	if (class(kernel) == "matrix") {
		kernel.matrix <- kernel
	}else{
		kernel.matrix <- lskmTest.GetKernel(Z, kernel, weights, n=nrow(Z), m=ncol(Z)) # n x n, has to be a matrix
	}



	#-------------------------------------
	# Run coxSKAT

	if(npert.check==FALSE){
		if(npert<=10^4){
			coxSKAT.out <- KMTest.surv(Z=NULL, U=U, Delta=Delta, X=X, gamma=gamma, kernel=kernel.matrix, npert=npert)
			return(list(p.value=max(0.5/npert, coxSKAT.out$pvalue.ptb), Q=coxSKAT.out$shat, n.marker.test=ncol(Z), 
				n.indiv=coxSKAT.out$n, df= coxSKAT.out$df.chisq))
		} else{
			n.iterations <- ceiling(npert/10^4)
			p.iterated <- c()
			for(ttt in n.iterations){
				ptemp <- KMTest.surv(Z=NULL, U=U, Delta=Delta, X=X, gamma=gamma, kernel=kernel.matrix, npert=10^4)
				p.iterated <- c(p.iterated, ptemp$pvalue.ptb)
			}
			return(list(p.value=max(0.5/(10^4*ceiling(npert/10^4)), mean(p.iterated)), Q=ptemp$shat, n.marker.test=ncol(Z),
                        	n.indiv=ptemp$n, df= ptemp$df.chisq))

	
		}	
	}else{

	       coxSKAT.out <- KMTest.surv(Z=NULL, U=U, Delta=Delta, X=X, gamma=gamma, kernel=kernel.matrix, npert=npert)
		if(coxSKAT.out$pvalue.ptb<=npert.threshold/npert){

			n.iterations <- ceiling(npert.upper/10^4)
                        p.iterated <- c()
                        for(ttt in n.iterations){
                                ptemp <- KMTest.surv(Z=NULL, U=U, Delta=Delta, X=X, gamma=gamma, kernel=kernel.matrix, npert=10^4)
                                p.iterated <- c(p.iterated, ptemp$pvalue.ptb)
                        }

                 return(list(p.value=max(0.5/(n.iterations*10^4),mean(p.iterated)), Q=coxSKAT.out$shat, n.marker.test=ncol(Z),
                                n.indiv=coxSKAT.out$n, df= coxSKAT.out$df.chisq))


		}else{
			 return(list(p.value=max(0.5/npert, coxSKAT.out$pvalue.ptb), Q=coxSKAT.out$shat, n.marker.test=ncol(Z),
                                n.indiv=coxSKAT.out$n, df= coxSKAT.out$df.chisq))
		}
	}



}














#     modified from iSKAT_MAIN_Check_Z()
#	Check the Z, and do imputation
#     iSKAT_MAIN_Check_Z() modified significantly from SKAT_MAIN_Check_Z from V0.78
#
coxSKAT_MAIN_Check_Z <- function(Z, SetID, weights=NULL, impute.method,  missing_cutoff){

	# check.Z.error = 0 : no snps removed, but some snps possibly imputed
	# check.Z.error = 1 : all snps removed, returns NULL matrix for Z
	# check.Z.error = 2 : some snps removed, remainder snps may have been imputed

	check.Z.error <- 0
	n <- nrow(Z)
	##############################################
	# Recode Missing to NA

	IDX_MISS <- union(which(is.na(Z)), which(Z == 9))
	if(length(IDX_MISS) > 0){
		Z[IDX_MISS] <- NA
	} 

	###################################################
	# Check missing rates and exclude any SNPs with missing rate > missing_cutoff
	# Also exclude non-polymorphic SNPs
	m <- ncol(Z)
	ID_INCLUDE_SNP <- NULL
	for(i in 1:m){
		missing.ratio <- length(which(is.na(Z[,i])))/n
		sd1 <- sd(Z[,i], na.rm=TRUE)
		if(missing.ratio < missing_cutoff && sd1 > 0){
			ID_INCLUDE_SNP <- c(ID_INCLUDE_SNP,i)
		}
	}
	
	if(length(ID_INCLUDE_SNP) == 0){

		if(is.null(SetID)){
			msg <- sprintf("ALL SNPs have either high missing rates or no-variation. ")
		} else {
			msg <- sprintf("In %s, ALL SNPs have either high missing rates or no-variation. ",SetID )
		}
		warning(msg, call.=FALSE)
		
		re <- list(Z.test=NULL, weights.Z.test=NULL, check.Z.error=1) 
		  

	} else{

 		if(m - length(ID_INCLUDE_SNP) > 0 ){

			if(is.null(SetID)){
				msg <- sprintf("%d SNPs with either high missing rates or no-variation are excluded!", m - length(ID_INCLUDE_SNP))
			} else {
				msg <- sprintf("In %s, %d SNPs with either high missing rates or no-variation are excluded!",SetID, m - length(ID_INCLUDE_SNP) )
			}

			warning(msg, call.=FALSE)	
			Z <- as.matrix(Z[,ID_INCLUDE_SNP])
			if(is.null(weights)==FALSE) weights <- weights[ID_INCLUDE_SNP]
			check.Z.error <- 2
			IDX_MISS <- which(is.na(Z))
		}
	
	

		##########################################
		# Missing Imputation
		if(length(IDX_MISS) > 0){
			if(is.null(SetID)){
				msg <- sprintf("The missing genotype rate is %f. Imputation is applied.", (length(IDX_MISS))/length(Z) )
			} else {
				msg <- sprintf("In %s, the missing genotype rate is %f. Imputation is applied.", SetID, (length(IDX_MISS))/length(Z) )
			}

			warning(msg,call.=FALSE)
			Z <- Impute_coxSKAT(Z,impute.method)
		} 
		re <- list(Z.test=Z, weights.Z.test=weights, check.Z.error=check.Z.error)
	
	}

	return(re)
}




# copied without modifcation from SKAT V0.78, renamed Impute_coxSKAT()
# Simple Imputation
# Z : an n x p genotype matrix with n samples and p SNPs
# Missing has to be NA: a missing genotype value.

Impute_coxSKAT <-function(Z, impute.method){
	
	p <- dim(Z)[2]

	if(impute.method =="random"){
		for(i in 1:p){
			IDX <- which(is.na(Z[,i]))
			if(length(IDX) > 0){
				maf1 <- mean(Z[-IDX,i])/2
				Z[IDX,i] <- rbinom(length(IDX),2,maf1)
			}
		}
	} else if(impute.method =="fixed"){
		for(i in 1:p){
			IDX<-which(is.na(Z[,i]))
			if(length(IDX) > 0){
				maf1 <- mean(Z[-IDX,i])/2
				Z[IDX,i] <- 2 * maf1
			}
		}
	} else {
		stop("Error: Imputation method shoud be either \"fixed\" or \"random\"! ")
	}

	return(Z)
}
