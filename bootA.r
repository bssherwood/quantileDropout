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
	else if( type=="asymetric" )
		exp( matrix(rnorm(n*m), n,m) %*% ss) - exp(SIGMA[1]/2)
}

####################################################
### Simulation Settings
nboot = 2#100
error = "asymetric" # "N" # "t3" # 
rhos = .75#c(.75, .5, .25)
ns = 100#c(100, 500,1000,5000)#,10000) # ns = c(300, 500) # 
ms = c(2, 3)
taus = c(.5,.7)#c(.5, .7) # taus = c(.3, .5, .7) # 
nsim = 100
sigma2 = 2^2 # 1 # .7^2 # 
### Coefficients for model : y = XB + error
B = rep(1,4)#c( 0, 2, 1, 2 )

### Coefficients for log-odds of dropout: P(R_ij=1|R_ij-1=1)
# Coefs: Int, x1, x2, x3, x4, z1, z2, y1 (y2, y3, ...)  
B2 = c(    4,  0,  0,  0,  0,  0,  0, rep(0,0), -1)   
B3 = c(    4,  0,  0,  0,  0,  0,  0, rep(0,1), -1) 
B4 = c(    4,  0,  0,  0,  0,  0,  0, rep(0,2), -1) 
B5 = c(    4,  0,  0,  0,  0,  0,  0, rep(0,3), -1)  
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
	lin.coefs <- 1:ncol(X)+7
	Z <- cbind(z1,z2)
	Zbs <- cbind(1, bs(z1, degree=3, Boundary.knots=c(0,1)),bs(z2, degree=3, Boundary.knots=c(-1,1)) )
	XZbs <- cbind(Zbs,X)
	XB <- X%*%B
	Zeffect <- sin(2*pi*z1) + z2^3
	Effect <- matrix(XB, n, max(ms)) + matrix(Zeffect, n, max(ms))

	X.dropout <- cbind(1, x1, x2, x3, x4, z1, z2)

	### Initialize tables for saving Bias, MSE, and gMSE, ConfIntWidth, ConfIntCoverage
	Bias <- matrix(0, 0, length(B)+4)
	colnames(Bias) <- c("tau", "m", "rho", "n", paste("B",1:length(B), sep=""))
	MSE <- Bias
	colnames(MSE) <- c("tau", "m", "rho", "n", paste("B_mse",1:length(B), sep=""))
	gMSE <- matrix(0,0, 5)
	colnames(gMSE) <- c("tau", "m", "rho", "n", "gMSE")

	Dropout <- matrix(0,0, 1+max(ms))
	colnames(Dropout) <- c("rho", "n", paste("Time",2:max(ms)))

	ConfIntWidth <- Bias
	colnames(ConfIntWidth) <- c("tau", "m", "rho", "n", paste("B_IntervalWidth",1:length(B), sep=""))

	ConfIntCoverage <- Bias
	colnames(ConfIntCoverage) <- c("tau", "m", "rho", "n", paste("B_IntervalCoverage",1:length(B), sep=""))



	for( rho in rhos ){
		### Generate Y for all values of m
		SIGMA <- matrix(rho, max(ms),max(ms))
		SIGMA <- sigma2*SIGMA^(abs(row(SIGMA)-col(SIGMA)))

		if( error != "Hetero" ){
			Ys <- lapply( 1:nsim, function(xx) Effect + genError(error, n, max(ms), SIGMA) )
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
		# rm(XY.dropout)
		rm(last.time.1)
		rm(this.time.1)
		rm(XY.temp)
		rm(R.temp)
		rm(Pest.old)


		### Prepare for a bootstrap where we re-sample 
		## the entire observation vector for each patient
		## and estimate the probabilities
		# Strategy is to make a list of lists where outer layer of lists corresponds
		# to each simulation round 
		# and each inner layer corresponds to each bootstrap for that simulation round
		boot.index <- mapply(function(xx) lapply(1:nboot, function(yy) sample(1:n,n,replace=TRUE)), 1:nsim, SIMPLIFY=FALSE)
		R.boot <- mapply( function(xx,yy) lapply(yy, function(aa) xx[aa,]),
					R, boot.index, SIMPLIFY=FALSE )
		XY.dropout.boot <- mapply( function(xx,yy) lapply(yy, function(aa) xx[aa,]),
					XY.dropout, boot.index, SIMPLIFY=FALSE )
		rm(XY.dropout)

		Pest.boot <- lapply( 1:nsim, function(xx) 
						lapply(1:nboot, function(yy) rep(1,n)) )
		Pest.old.boot <- Pest.boot

		# I got lost in all the m/lapply so I switched to loops
		# The list structure still holds
		# This is probably less efficient
		index.boot <- 
			rowIndex.boot <- 
			last.time.1.boot <- 
			this.time.1.boot <- 
			XY.temp.boot <- 
			R.temp.boot <- 
				lapply(1:nsim, function(xx) lapply(1:nboot, function(aa) 1))

		for( jj in 1:nsim ){
			for( kk in 1:nboot ){
				index.boot[[jj]][[kk]] <- which(R.boot[[jj]][[kk]]==1)
				# index.boot <- lapply(R, function(xx) which(xx==1))
			
				rowIndex.boot[[jj]][[kk]] <- c( row( R.boot[[jj]][[kk]] )*R.boot[[jj]][[kk]] )[ c(R.boot[[jj]][[kk]])!=0 ]
				# rowIndex <- lapply( R, function(xx) c(row(xx)*xx)[c(xx)!=0] )
			
				for( j in 1:length(Bdrop) ){ # j = 1 corresponds to Time Point 2 (length(Bdrop) = m-1)
					last.time.1.boot[[jj]][[kk]] <- which( R.boot[[jj]][[kk]][,j]==1 )
					# last.time.1 <- lapply(R, function(xx) which(xx[,j]==1))
			
					this.time.1.boot[[jj]][[kk]] <- which( R.boot[[jj]][[kk]][ last.time.1.boot[[jj]][[kk]], j+1]==1 )
					# this.time.1 <- mapply(function(xx,yy) which(xx[yy,j+1]==1), R,last.time.1, SIMPLIFY=FALSE)

					XY.temp.boot[[jj]][[kk]] <- (XY.dropout.boot[[jj]][[kk]])[last.time.1.boot[[jj]][[kk]], which(Bdrop[[j]]!=0)] 
					# XY.temp <- mapply(function(xx, yy) xx[yy, which(Bdrop[[j]]!=0)], XY.dropout, last.time.1, SIMPLIFY=FALSE )

					R.temp.boot[[jj]][[kk]] <- (R.boot[[jj]][[kk]])[last.time.1.boot[[jj]][[kk]],j+1]
					# R.temp <- mapply(function(xx,yy) xx[yy,j+1], R, last.time.1, SIMPLIFY=FALSE)

					Pest.old.boot[[jj]][[kk]] <- (Pest.old.boot[[jj]][[kk]])[this.time.1.boot[[jj]][[kk]]]*glm.fit(XY.temp.boot[[jj]][[kk]],R.temp.boot[[jj]][[kk]], family=binomial())$fitted.values[this.time.1.boot[[jj]][[kk]]]
					# Pest.old <- mapply(function(xx,yy,zz,ww) zz[ww]*glm.fit(xx,yy, family=binomial())$fitted.values[ww], XY.temp, R.temp, Pest.old, this.time.1, SIMPLIFY=FALSE)

					Pest.boot[[jj]][[kk]] <- c(Pest.boot[[jj]][[kk]], Pest.old.boot[[jj]][[kk]])
					# Pest <- mapply(function(xx,yy) c(xx,yy), Pest, Pest.old, SIMPLIFY=FALSE)
				}
			}
		}

		## Delete these to help with memory issues
		rm(XY.dropout.boot)
		rm(last.time.1.boot)
		rm(this.time.1.boot)
		rm(XY.temp.boot)
		rm(R.temp.boot)
		rm(Pest.old.boot)


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

			### Now do the bootstrap part
			Rm.boot <- 
				Ym.boot <- 
				XZbsm.boot <- 
				Zbsm.boot <- 
				ZEffectm.boot <- 
				Pestm.boot <- 
				Ys.boot <- 
				XZbs.boot <- 
				Zbs.boot <- 
				Zeffect.boot <- 
				Yo.boot <- 
				XZbso.boot <- 
				Zbso.boot <- 
				Zeffecto.boot <- 
					lapply(1:nsim, function(xx) lapply(1:nboot, function(aa) 1))
			for( jj in 1:nsim ){
			for( kk in 1:nboot ){
				Ys.boot[[jj]][[kk]] <- (Ys[[jj]])[boot.index[[jj]][[kk]], ]
				XZbs.boot[[jj]][[kk]] <- XZbs[boot.index[[jj]][[kk]], ]
				Zbs.boot[[jj]][[kk]] <- Zbs[boot.index[[jj]][[kk]], ]
				Zeffect.boot[[jj]][[kk]] <- Zeffect[ boot.index[[jj]][[kk]] ]

				Rm.boot[[jj]][[kk]] <- sum( (R.boot[[jj]][[kk]])[,1:m] ) # Get total number of observations
				Ym.boot[[jj]][[kk]] <- c(Ys.boot[[jj]][[kk]])[index.boot[[jj]][[kk]]][ 1:Rm.boot[[jj]][[kk]] ]
				XZbsm.boot[[jj]][[kk]] <- (XZbs.boot[[jj]][[kk]])[(rowIndex.boot[[jj]][[kk]])[1:Rm.boot[[jj]][[kk]]],]
				Zbsm.boot[[jj]][[kk]] <- (Zbs.boot[[jj]][[kk]])[(rowIndex.boot[[jj]][[kk]])[1:Rm.boot[[jj]][[kk]]],]
				ZEffectm.boot[[jj]][[kk]] <- (Zeffect.boot[[jj]][[kk]])[(rowIndex.boot[[jj]][[kk]])[1:Rm.boot[[jj]][[kk]]]]
				Pestm.boot[[jj]][[kk]] <- (Pest.boot[[jj]][[kk]])[1:Rm.boot[[jj]][[kk]]]

				Yo.boot[[jj]][[kk]] <- c(Ys.boot[[jj]][[kk]][,1:m])
				XZbso.boot[[jj]][[kk]] <- XZbs.boot[[jj]][[kk]][rep(1:n,m),]
				Zbso.boot[[jj]][[kk]] <- Zbs.boot[[jj]][[kk]][rep(1:n,m),]
				Zeffecto.boot[[jj]][[kk]] <- Zeffect.boot[[jj]][[kk]][rep(1:n,m)]


			}
			}

			rm(Rm.boot)


			for( tau in taus ){

				cat("tau =", tau, "  m =", m, "  rho =", rho, "  n =", n, "\n")
				cat("error is",error, " sigma2 =", sigma2, "\n\n")

				if( error=="Hetero" )
					BTrue[1] <- B[1] + qnorm(tau, mean=0, sd=sqrt(SIGMA[1]))

				settings <- c(tau=tau, m=m, rho=rho, n=n)

				### Oracle Model
				coefs <- mapply(function(xx) rq.fit( XZbso, xx, tau=tau)$coefficients, Yo, SIMPLIFY=FALSE)

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

				bootN <- lapply(1:nsim, function(xx) 1)
				for( jj in 1:nsim ){
					bootN[[jj]] <- t(mapply( function(xx,yy) rq.fit( xx, yy, tau=tau)$coefficients[lin.coefs], XZbso.boot[[jj]], Yo.boot[[jj]] ))
				}

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
				coefs <- mapply(function(xx,yy) rq.fit( xx, yy, tau=tau)$coefficients, XZbsm, Ym, SIMPLIFY=FALSE)

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

				bootN <- lapply(1:nsim, function(xx) 1)
				for( jj in 1:nsim ){
					bootN[[jj]] <- t(mapply( function(xx,yy) rq.fit( xx, yy, tau=tau)$coefficients[lin.coefs], XZbsm.boot[[jj]], Ym.boot[[jj]] ))
				}
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
				coefs <- mapply(function(xx,yy,zz){
									rq.fit( pmin(1/zz, 20)*xx, pmin(1/zz, 20)*yy , tau=tau)$coefficients}, XZbsm, Ym, Pestm, SIMPLIFY=FALSE)

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

				bootN <- lapply(1:nsim, function(xx) 1)
				for( jj in 1:nsim ){
					bootN[[jj]] <- t(mapply( function(xx,yy,zz){
									rq.fit( pmin(1/zz, 20)*xx, pmin(1/zz, 20)*yy , tau=tau)$coefficients[lin.coefs]}, XZbsm.boot[[jj]], Ym.boot[[jj]], Pestm.boot[[jj]] ))
				}
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
			cbind(Bias, MSE=MSE[,-1:-4], gMSE=gMSE[,-1:-4],ConfIntWidth[,-1:-4],ConfIntCoverage[,-1:-4]), row.names=NULL )
	# print(xtable( out , 
	# 	digits=c(1,1,1, 0, 2, 0, rep(3, ncol(Bias)+2-4)) ), include.rownames=FALSE)
	# cat("\n\n")

	## Print dropout rates for a single value of n (should be similar for all n)
	# print(xtable(Dropout), include.rownames=FALSE)
	# cat("\n\n\n")

	saveRDS( list(n=n, ms=ms, rhos=rhos, taus=taus, BTrue=BTrue, error=error, 
					sigma2=sigma2, results=out, Dropout=Dropout), 
					file=paste(error,"n",n,"var",sigma2,"xyTrue.RDS", sep="") )

	cat("Finished simulations for all n =", n, "\n\n")

}

