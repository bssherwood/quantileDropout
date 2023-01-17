library(quantreg)
library(xtable)
library(splines)

### Function for generating errors for responses
genError <- function(type, n, m, SIGMA){
	ss <- eigen(SIGMA)
	ss <- ss$vectors %*% diag(ss$values^(.5)) %*% t(ss$vectors)
	if( type == "N" )
		matrix(rnorm(n*m), n,m) %*% ss
	else if( type == "t3" )
		matrix(rt(n*m, 3), n,m) %*% ss
	else if( type=="assymetric" )
		matrix(rchisq(n*m, 3), n,m) %*% ss
}

####################################################
### Simulation Settings
error = "Hetero"#"Hetero" # "N" # "t3" # 
rhos = .75#c(.75, .5, .25)#.75#c(.5,.75)#
ns = c(100, 500,1000,5000) # ns = c(300, 500) # 
ms = c(2,3)# 5)
taus = .5#c(.5, .7) # taus = c(.3, .5, .7) # 
nsim = 100
sigma2 = 2^2 # 1 # .7^2 # 
### Coefficients for model : y = XB + error


outTotal <- NULL
for(ycoef in c(-1)){

	B = rep(1,4)#c( 2, 2, 2, 2 )
	#if( error == "N" )
	#	B[1] <- .5
	#bmod <- c(.25,1,3)
	#xmissCoef <- c(-2,0,2)
	#ymissMult <- c(

	### Coefficients for log-odds of dropout: P(R_ij=1|R_ij-1=1)
	# Coefs: Int, x1, x2, x3, x4, z1, z2, y1 (y2, y3, ...)  
	B2 = c(    4,  0,  0,  0,  0,  0,  0, rep(0,0), ycoef)   
	B3 = c(    4,  0,  0,  0,  0,  0,  0, rep(0,1), ycoef) 
	B4 = c(    4,  0,  0,  0,  0,  0,  0, rep(0,2), ycoef) 
	B5 = c(    4,  0,  0,  0,  0,  0,  0, rep(0,3), ycoef) 
	####################################################

	# SIGMA <- matrix(rho, m,m)
	# SIGMA <- sigma2*SIGMA^(abs(row(SIGMA)-col(SIGMA)))

	Bdrop <- lapply( 2:max(ms), function(xx) get(paste("B",xx,sep="")) )
	names(Bdrop) <- paste("B", 2:max(ms), sep="")
	BTrue <- B



	for( n in ns ){
		set.seed(4)

		### Generate covariates for all settings
		x1 <- runif(n, 0,1)
		x2 <- rnorm(n)
		x3 <- rnorm(n)
		x4 <- rnorm(n)
		z1 <- runif(n,  0,1)
		z2 <- runif(n, -1,1)

		X <- cbind(x1,x2,x3,x4)
		lin.coefs <- 1:ncol(X)
		Z <- cbind(z1,z2)
		Zbs <- cbind(1, bs(z1, degree=3, Boundary.knots=c(0,1)),bs(z2, degree=3, Boundary.knots=c(-1,1)) )
		XZbs <- cbind(X, Zbs)
		XB <- X%*%B
		Zeffect <- sin(2*pi*z1) + z2^3
		Effect <- matrix(XB, n, max(ms)) + matrix(Zeffect, n, max(ms))

		X.dropout <- cbind(1, x1, x2, x3, x4, z1, z2)

		### Initialize tables for saving Bias, MSE, and gMSE, ConfIntWidth, ConfIntCoverage
		Bias <- matrix(0, 0, length(B)+5)
		colnames(Bias) <- c("tau", "m", "rho", "n","ycoef", paste("B",1:length(B), sep=""))
		MSE <- Bias
		colnames(MSE) <- c("tau", "m", "rho", "n","ycoef", paste("B_mse",1:length(B), sep=""))
		gMSE <- matrix(0,0, 6)
		colnames(gMSE) <- c("tau", "m", "rho", "n","ycoef", "gMSE")

		Dropout <- matrix(0,0, 1+max(ms))
		colnames(Dropout) <- c("rho", "n", paste("Time",2:max(ms)))

		ConfIntWidth <- Bias
		colnames(ConfIntWidth) <- c("tau", "m", "rho", "n","ycoef", paste("B_IntervalWidth",1:length(B), sep=""))

		ConfIntCoverage <- Bias
		colnames(ConfIntCoverage) <- c("tau", "m", "rho", "n","ycoef", paste("B_IntervalCoverage",1:length(B), sep=""))



		for( rho in rhos ){
			### Generate Y for all values of m
			SIGMA <- matrix(rho, max(ms),max(ms))
			SIGMA <- sigma2*SIGMA^(abs(row(SIGMA)-col(SIGMA)))

			if( error == "N" ){
				Ys <- lapply( 1:nsim, function(xx) Effect + genError("N", n, max(ms), SIGMA) )
			} else { # error == "Hetero"
				Ys <- lapply( 1:nsim, function(xx) Effect + x1*genError("N", n, max(ms), SIGMA) )
			}

			### Generate probabilities and dropout
			XY.dropout <- lapply( Ys, function(xx) cbind(X.dropout, xx) )
			PTrue <- lapply( 1:nsim, function(xx) matrix(1, n,1) )
			R <- PTrue
			for( j in 1:length(Bdrop) ){
				odds <- lapply(XY.dropout, function(xx) exp(xx[,which(Bdrop[[j]]!=0)] %*% Bdrop[[j]][which(Bdrop[[j]]!=0)]))
				PTrue <- mapply(function(xx,yy) cbind(xx, yy/(1+yy)), PTrue, odds, SIMPLIFY=FALSE)
				R <- mapply(function(xx,yy) cbind(xx, xx[,j]*rbinom(n,1,yy[,j+1])), R, PTrue, SIMPLIFY=FALSE)
			}

			## Delete these to help with memory issues
			rm(PTrue)
			rm(odds)

			### Estimate Probabilities
			Pest <- lapply(1:nsim, function(xx) rep(1,n))
			Pest.old <- Pest
			index <- lapply(R, function(xx) which(xx==1))
			rowIndex <- lapply( R, function(xx) c(row(xx)*xx)[c(xx)!=0] )
			for( j in 1:length(Bdrop) ){ # j = 1 corresponds to Time Point 2 (length(Bdrop) = m-1)
				last.time.1 <- lapply(R, function(xx) which(xx[,j]==1))
				this.time.1 <- mapply(function(xx,yy) which(xx[yy,j+1]==1), R,last.time.1, SIMPLIFY=FALSE)

				XY.temp <- mapply(function(xx, yy) xx[yy, which(Bdrop[[j]]!=0)], XY.dropout, last.time.1, SIMPLIFY=FALSE )
				R.temp <- mapply(function(xx,yy) xx[yy,j+1], R, last.time.1, SIMPLIFY=FALSE)

				Pest.old <- mapply(function(xx,yy,zz,ww) zz[ww]*glm.fit(xx,yy, family=binomial())$fitted.values[ww], XY.temp, R.temp, Pest.old, this.time.1, SIMPLIFY=FALSE)
				Pest <- mapply(function(xx,yy) c(xx,yy), Pest, Pest.old, SIMPLIFY=FALSE)
			}

			## Delete these to help with memory issues
			rm(XY.dropout)
			rm(last.time.1)
			rm(this.time.1)
			rm(XY.temp)
			rm(R.temp)
			rm(Pest.old)

			drops <- sapply(R, function(xx) colMeans(xx[,-1]))
			mins <- round(100*apply(drops, 1, min), 0)
			maxs <- round(100*apply(drops, 1, max), 0)
			means <- round(100*rowMeans(drops), 0)
			Dropout <- rbind(Dropout, 
				c(rho, n, apply(cbind(mins, maxs, means), 1, function(xx) paste(xx[1],"-",xx[2]," (",xx[3],")", sep=""))))

			rm(drops)


			for( m in ms ){
				### Only use relevant Ys
				Rm <- lapply( R, function(xx) sum(xx[,1:m]) ) # Get total number of observations
				Ym <- mapply(function(xx,yy,zz) c(xx)[yy][1:zz], Ys, index, Rm, SIMPLIFY=FALSE)
				XZbsm <- mapply(function(xx,yy) XZbs[xx[1:yy],], rowIndex, Rm, SIMPLIFY=FALSE)
				Zbsm <- mapply(function(xx,yy) Zbs[xx[1:yy],], rowIndex, Rm, SIMPLIFY=FALSE)
				ZEffectm <- mapply(function(xx,yy) Zeffect[xx[1:yy]], rowIndex, Rm, SIMPLIFY=FALSE)
				Pestm <- mapply(function(xx,yy) xx[1:yy], Pest, Rm, SIMPLIFY=FALSE)

				### Oracle data
				Yo <- lapply(Ys, function(xx) c(xx[,1:m]))
				XZbso <- XZbs[rep(1:n,m),]
				Zbso <- Zbs[rep(1:n,m),]
				Zeffecto <- Zeffect[rep(1:n,m)]

				## Delete these to help with memory issues
				rm(Rm)



				for( tau in taus ){

					cat("tau =", tau, "  m =", m, "  rho =", rho, "  n =", n, "\n")
					cat("error is",error, " sigma2 =", sigma2, "\n\n")

					if( error=="Hetero" )
						BTrue[1] <- B[1] + qnorm(tau, mean=0, sd=sqrt(SIGMA[1]))

					settings <- c(tau=tau, m=m, rho=rho, n=n,ycoef=ycoef)

					### Oracle Model
					coefs <- mapply(function(xx) rq.fit( XZbso, xx, tau=tau )$coefficients, Yo, SIMPLIFY=FALSE)

					bias <- sapply(coefs, function(xx) xx[lin.coefs] - BTrue)
					ghat <- mapply(function(xx,yy) xx%*%yy[-lin.coefs], Zbsm, coefs)
					
					rm(coefs)

					mse <- apply(bias, 2, function(xx) xx^2)
					Bias <- rbind( Bias, Oracle=c(settings, abs(rowMeans(bias))) )
					MSE  <-  rbind( MSE, Oracle=c(settings, rowMeans(mse)) )
					
					rm(bias)
					rm(mse)

					gmse <- mapply(function(xx,yy) mean((xx-mean(xx) - yy)^2), ghat, ZEffectm)
					gMSE <- rbind( gMSE, Oracle=c(settings, mean(gmse)) )
					rm(ghat)
					rm(gmse)

					bootN <- mapply(function(xx) boot.rq( XZbso, xx, tau=tau, R=1000,bsmethod="xy" )$B[,1:length(B)], Yo, SIMPLIFY=FALSE)
					CI <- lapply(bootN, function(xx) apply(xx, 2, quantile, probs=c(0.05, 0.95)) )
					width <- sapply(CI, function(xx) apply(xx, 2, diff) )
					coverage <- sapply(CI, function(xx) 1*(xx[1,]<=BTrue & BTrue<=xx[2,]) )

					ConfIntWidth <- rbind( ConfIntWidth, Oracle=c(settings, rowMeans(width)) )
					ConfIntCoverage <- rbind( ConfIntCoverage, Oracle=c(settings, rowMeans(coverage)) )
					rm(bootN)
					rm(CI)
					rm(width)
					rm(coverage)



					### Naive Model
					coefs <- mapply(function(xx,yy) rq.fit( xx, yy, tau=tau )$coefficients, XZbsm, Ym, SIMPLIFY=FALSE)

					bias <- sapply(coefs, function(xx) xx[lin.coefs] - BTrue)
					ghat <- mapply(function(xx,yy) xx%*%yy[-lin.coefs], Zbsm, coefs)
					
					rm(coefs)

					mse <- apply(bias, 2, function(xx) xx^2)
					Bias <- rbind( Bias, Naive=c(settings, abs(rowMeans(bias))) )
					MSE  <-  rbind( MSE, Naive=c(settings, rowMeans(mse)) )
					
					rm(bias)
					rm(mse)

					gmse <- mapply(function(xx,yy) mean((xx-mean(xx) - yy)^2), ghat, ZEffectm)
					gMSE <- rbind( gMSE, Naive=c(settings, mean(gmse)) )
					rm(ghat)
					rm(gmse)

					bootN <- mapply(function(xx,yy) boot.rq( xx, yy, tau=tau, R=1000,bsmethod="xy" )$B[,1:length(B)], XZbsm, Ym, SIMPLIFY=FALSE)
					CI <- lapply(bootN, function(xx) apply(xx, 2, quantile, probs=c(0.05, 0.95)) )
					width <- sapply(CI, function(xx) apply(xx, 2, diff) )
					coverage <- sapply(CI, function(xx) 1*(xx[1,]<=BTrue & BTrue<=xx[2,]) )

					ConfIntWidth <- rbind( ConfIntWidth, Naive=c(settings, rowMeans(width)) )
					ConfIntCoverage <- rbind( ConfIntCoverage, Naive=c(settings, rowMeans(coverage)) )
					rm(bootN)
					rm(CI)
					rm(width)
					rm(coverage)



					### IPW Model
					coefs <- mapply(function(xx,yy,zz) rq.fit( pmin(1/zz, 20)*xx, pmin(1/zz, 20)*yy, tau=tau )$coefficients, XZbsm, Ym, Pestm, SIMPLIFY=FALSE)

					bias <- sapply(coefs, function(xx) xx[lin.coefs] - BTrue)
					ghat <- mapply(function(xx,yy) xx%*%yy[-lin.coefs], Zbsm, coefs)
					
					rm(coefs)

					mse <- apply(bias, 2, function(xx) xx^2)
					Bias <- rbind( Bias, IPW=c(settings, abs(rowMeans(bias))) )
					MSE  <-  rbind( MSE, IPW=c(settings, rowMeans(mse)) )
					
					rm(bias)
					rm(mse)

					gmse <- mapply(function(xx,yy) mean((xx-mean(xx) - yy)^2), ghat, ZEffectm)
					gMSE <- rbind( gMSE, IPW=c(settings, mean(gmse)) )
					rm(ghat)
					rm(gmse)

					bootN <- mapply(function(xx,yy,zz) boot.rq( pmin(1/zz, 20)*xx, pmin(1/zz, 20)*yy, tau=tau, R=1000,bsmethod="xy" )$B[,1:length(B)], XZbsm, Ym, Pestm, SIMPLIFY=FALSE)
					CI <- lapply(bootN, function(xx) apply(xx, 2, quantile, probs=c(0.05, 0.95)) )
					width <- sapply(CI, function(xx) apply(xx, 2, diff) )
					coverage <- sapply(CI, function(xx) 1*(xx[1,]<=BTrue & BTrue<=xx[2,]) )

					ConfIntWidth <- rbind( ConfIntWidth, IPW=c(settings, rowMeans(width)) )
					ConfIntCoverage <- rbind( ConfIntCoverage, IPW=c(settings, rowMeans(coverage)) )
					rm(bootN)
					rm(CI)
					rm(width)
					rm(coverage)


	# sd( colSums(mse) )/sqrt(nsim)
	# sd( gmse )/sqrt(nsim)

				} # Close tau loop
			} # Close m loop
		} # Close rho loop


		## Print results for a single value of n
		out <- data.frame( Method=row.names(Bias),
				cbind(Bias, MSE=MSE[,-1:-5], gMSE=gMSE[,-1:-5],ConfIntWidth[,-1:-5],ConfIntCoverage[,-1:-5]), row.names=NULL )
		# print(xtable( out , 
		# 	digits=c(1,1,1, 0, 2, 0, rep(3, ncol(Bias)+2-4)) ), include.rownames=FALSE)
		# cat("\n\n")

		## Print dropout rates for a single value of n (should be similar for all n)
		# print(xtable(Dropout), include.rownames=FALSE)
		# cat("\n\n\n")

		saveRDS( list(n=n, ms=ms, rhos=rhos, taus=taus, BTrue=BTrue, error=error, 
						sigma2=sigma2, results=out, Dropout=Dropout), 
						file=paste("goodHeteroxyBig5_n",n,"_2022.RDS", sep="") )

		cat("Finished simulations for all n =", n, "\n\n")

	} # Close n loop
	outTotal <- rbind(outTotal, out)
}

saveRDS(outTotal, "goodHeteroxyBig.RDS")
#out1 <- readRDS("Nn500var4bInter-2_2022.RDS")$results
#out2 <- readRDS("Nn500var4bInter2_2022.RDS")$results
#out3 <- readRDS("Nn500var4bInter0.1_2022.RDS")$results

#outTotal <- rbind(out1,out2,out3)

#saveRDS(outTotal,"checkingMissing.RDS")







