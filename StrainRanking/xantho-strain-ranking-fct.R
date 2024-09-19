
## R code contributing to the analyses using StrainRanking and presented in Bandhari et al. (prep., Climatic fluctuations modulate eco-evolutionary processes associated with pathogen dynamics and continue to fuel bacterial spot epidemics in tomato).
## Underlying functions used in the main code
## Version 2024-09-17
## Author: Samuel Soubeyrand (INRAE, BioSP, 84914 Avignon, France; <samuel.soubeyrand@inrae.fr>)


## Extract and format data for using them in the StrainRanking analyses
## Genetic frequencies at time1 and time2 are averaged via a weighted mean with weights
## prop and (1-prop), respectively
## Covariates values at time1 and time2 are averaged via a weighted mean with weights
## prop.covar and (1-prop.covar), respectively
format.data=function(x,time1,time2,prop,prop.covar){
	id=unique(x$Farm_id)
	SC.col=8:15
	COVAR.col=c(3,(max(SC.col)+1):ncol(x))
	d=NULL
	for(i in id){
		Sev1=x$Severity[x$Farm_id==i & x$Year==time1[1] & x$Time==time1[2]]
		Sev2=x$Severity[x$Farm_id==i & x$Year==time2[1] & x$Time==time2[2]]
		if(length(Sev1)>0 & length(Sev2)>0){
			Freq1=x[x$Farm_id==i & x$Year==time1[1] & x$Time==time1[2],SC.col]
			Freq2=x[x$Farm_id==i & x$Year==time2[1] & x$Time==time2[2],SC.col]
			Freq=prop*Freq1+(1-prop)*Freq2
			Freq=Freq/sum(Freq)
			Covar1=x[x$Farm_id==i & x$Year==time1[1] & x$Time==time1[2],COVAR.col]
			Covar2=x[x$Farm_id==i & x$Year==time2[1] & x$Time==time2[2],COVAR.col]
			Covar=prop.covar*Covar1+(1-prop.covar)*Covar2
			if(sum(is.na(Freq)==0)){
				if(!is.na(Sev1+Sev2) & sum(Freq)>0){
					d$demographic.coord=rbind(d$demographic.coord,c(x$longKM[x$Farm_id==i][1],
						x$latKM[x$Farm_id==i][1]))
					d$demographic.measures=rbind(d$demographic.measures,c(Sev1,Sev2))
					d$genetic.frequencies=rbind(d$genetic.frequencies,Freq)
					d$covar=rbind(d$covar,Covar)
				}
			}
		}
	}
	d$genetic.coord=d$demographic.coord
	return(d)
}


## Data used for merging sub-datasets corresponding, for example, to different study years
## noise.sd>0 allows to noise the locations with a Gaussian distribution
## longlag>0 or latlag>0 allows to translate the locations with the provided lag value in the
##   aim of "separating in space" the sub-datasets
merge.data=function(xlist,noise.sd=0,longlag=0,latlag=0){
	x=list(NULL,NULL,NULL,NULL,NULL)
	names(x)=names(xlist[[1]])
	for(i in 1:length(xlist)){
		u=xlist[[i]]$demographic.coord
		if(noise.sd>0 | longlag>0 | latlag>0){
			u[,1]=rnorm(nrow(u),u[,1],noise.sd)+longlag*(i-1)
			u[,2]=rnorm(nrow(u),u[,2],noise.sd)+latlag*(i-1)
		}
		for(j in 1:length(x)){
			v=xlist[[i]][[j]]
			if((noise.sd>0 | longlag>0 | latlag>0) & grepl("coord",names(xlist[[i]])[j])){
				v=u
			}
			x[[j]]=rbind(x[[j]],v)
		}
	}
	return(x)
}


## Function used for accounting of an eventual thresholding of the growth variable
lm.s=function(Z,p,Zlim){
	par0=lm(Z ~ -1 + p)$coef
	f=function(par){
		mean((Z-pmax(Zlim[1],pmin(Zlim[2],p%*%par)))^2)
	}
	opt=optim(par0,f)
	opt$coefficients=opt$par
	opt$fitted.values=pmax(Zlim[1],pmin(Zlim[2],p%*%opt$par))
	return(opt)
}


## Modification of the ranking.strains() function in StrainRanking, updated for including
## thresholding of the growth variable (new argument Zlim)
ranking.strains.s=function (DGobject, bw, nb.mcsimul, plots = FALSE, kernel.type = "Quadratic",
	Zlim=c(-Inf,Inf) ) {
    if (!(kernel.type %in% c("Linear", "Quadratic", "Power3", 
        "Power4"))) {
        kernel.type = "Quadratic"
        print("WARNING: specified kernel.type is unknown. Default Quadratic kernel is used.")
    }
    geneticData = DGobject@genetic
    x = DGobject@demographic[, 1:2]
    Z = DGobject@demographic[, 3]
    nbStrains = ncol(geneticData) - 2
    N = nrow(x)
    pSEstimes = StrainRanking:::.estimation.pS(geneticData, x, bw, kernel.type)
#    regression <- lm(Z ~ -1 + pSEstimes)
    regression <- lm.s(Z,pSEstimes,Zlim)
    z = regression$coef
    zStar = matrix(0, nb.mcsimul, length(z))
    for (i in 1:nb.mcsimul) {
        permut = (1:N)[!is.na(pSEstimes[, 1])]
        pSEstimesStar = pSEstimes
        pSEstimesStar[permut, ] = pSEstimesStar[sample(permut), 
            ]
#        regressionStar = lm(Z ~ -1 + pSEstimesStar)
        regressionStar = lm.s(Z,pSEstimesStar,Zlim)
        zStar[i, ] = regressionStar$coefficients
    }
    pvals = NULL
    pvals.name = NULL
    for (i in 1:(ncol(zStar) - 1)) {
        for (j in (i + 1):ncol(zStar)) {
            pvals = c(pvals, mean(zStar[, j] - zStar[, i] >= 
                z[j] - z[i]))
            pvals.name = c(pvals.name, paste("(", j, "-", i, 
                ")>0", sep = ""))
            pvals = c(pvals, mean(zStar[, j] - zStar[, i] <= 
                z[j] - z[i]))
            pvals.name = c(pvals.name, paste("(", j, "-", i, 
                ")<0", sep = ""))
            pvals = c(pvals, mean(abs(zStar[, j] - zStar[, i]) >= 
                abs(z[j] - z[i])))
            pvals.name = c(pvals.name, paste("(", j, "-", i, 
                ")neq0", sep = ""))
        }
    }
    names(pvals) = pvals.name
    if (plots) {
        nrow.plot = ceiling((2 + nbStrains)/3)
        par(mfrow = c(nrow.plot, 3), mar = c(2, 2, 1, 0.2))
        plot(x, cex = (Z - min(Z))/(max(Z) - min(Z)) * 3, pch = 19, 
            asp = 1, axes = F, xlab = "", ylab = "", main = "Growth variable")
        plot(x, cex = 1, pch = 1, asp = 1, axes = F, xlab = "", 
            ylab = "", main = "Sampling sites")
        points(geneticData[, 1:2], pch = 19)
        for (j in 1:nbStrains) {
            zseq = seq(min(c(z, zStar), na.rm = TRUE), max(c(z, 
                zStar), na.rm = TRUE), l = 20)
            hist(zStar[, j], xlim = range(zseq), breaks = zseq, 
                xlab = "", ylab = "", main = paste("Coef. strain ", 
                  j))
            abline(v = z[j], lwd = 2, lty = "dashed", col = 2)
        }
    }
    names(z) = paste("strain", 1:nbStrains, sep = "")
    colnames(zStar) = paste("strain", 1:nbStrains, sep = "")
    return(list(permutation.estimates = zStar, estimates = z, 
        p.values = pvals))
}


## Function for estimating fitness uncertainty (under a Gaussian assumption of the error) and
## testing whether the fitness is larger or lower than 0 (the original ranking.strains() function
## in StrainRanking was only comparing fitness in a pairwise manner)
estim.fitness.strains=function (DGobject, bw, nb.mcsimul, plots = FALSE, kernel.type = "Quadratic",
	option=list(xlim="same",coverage=0.95,zero=TRUE),hyp="permut",Zlim=c(-Inf,Inf)){
    if (!(kernel.type %in% c("Linear", "Quadratic", "Power3", 
        "Power4"))) {
        kernel.type = "Quadratic"
        print("WARNING: specified kernel.type is unknown. Default Quadratic kernel is used.")
    }
    geneticData = DGobject@genetic
    x = DGobject@demographic[, 1:2]
    Z = DGobject@demographic[, 3]
    nbStrains = ncol(geneticData) - 2
    N = nrow(x)
    pSEstimes = StrainRanking:::.estimation.pS(geneticData, x, bw, kernel.type)
#    regression <- lm(Z ~ -1 + pSEstimes)
    regression <- lm.s(Z,pSEstimes,Zlim)
    z=regression$coefficients
    zStar = matrix(0, nb.mcsimul, length(z))
    zStar.gaussian = matrix(0, nb.mcsimul, length(z))
   	regression0 <- lm(Z ~ -1 + pSEstimes)
    sigmahat=sigma(regression0)
    pred=regression$fitted.values
    for (i in 1:nb.mcsimul) {
    	## permutation
   	   		permut = (1:N)[!is.na(pSEstimes[, 1])]
   	   		pSEstimesStar = pSEstimes
   	   		pSEstimesStar[permut, ] = pSEstimesStar[sample(permut), ]
#   	   		regressionStar = lm(Z ~ -1 + pSEstimesStar)
   	   		regressionStar = lm.s(Z,pSEstimesStar,Zlim)
   	   		zStar[i, ] = regressionStar$coefficients
   	   	## gaussian hypothesis
		    ZStar=Z
		    ZStar[regression$na]=NA
	   	 	ZStar[!is.na(ZStar)]=pred+rnorm(length(pred),0,sigmahat)
#   		 	regressionStar <- lm(ZStar ~ -1 + pSEstimes)
   		 	regressionStar <- lm.s(ZStar,pSEstimes,Zlim)
   	   		zStar.gaussian[i, ] = regressionStar$coefficients
    }
    if (plots) {
    	nrow.plot = ceiling((nbStrains)/3)
        par(mfrow = c(nrow.plot, 3), mar = c(2, 2, 1, 0.2))
        for (j in 1:nbStrains) {
            zseq = seq(min(c(z, zStar, zStar.gaussian), na.rm = TRUE), max(c(z, 
                zStar, zStar.gaussian), na.rm = TRUE), l = 20)
            if(option$xlim=="same" | sum(!is.na(zStar[, j]))==0){ 
        		hist(zStar[, j], xlim = range(zseq), breaks = zseq, 
	                xlab = "", ylab = "", main = paste("Coef. strain ", j))
            } else { 
        		hist(zStar[, j], xlab = "", ylab = "", main = paste("Coef. strain ", j),
        			xlim=range(c(z[j], zStar[,j], zStar.gaussian[,j]),na.rm=TRUE))
            }
            abline(v = z[j], lwd = 2, lty = "dashed", col = 2)
            if(!is.null(option$coverage) & sum(!is.na(zStar[, j]))>0){
            	alpha=1-option$coverage
            	abline(v = quantile(zStar[, j],c(alpha/2,1-alpha/2)), 
            			lwd = 1, lty = "dashed", col = 3)
            	abline(v = quantile(zStar.gaussian[, j],c(alpha/2,1-alpha/2)), 
            			lwd = 1, lty = "dashed", col = 4)
            }
            if(option$zero){
            	abline(v =  0, lwd = 2, lty = "dashed", col = 1)
            }
        }
    }
    
    pvals = pvals.gaussian = NULL
    pvals.name = NULL
    for (j in 1:ncol(zStar)) {
    	## permutation
            pvals = c(pvals, mean(zStar[, j] >= z[j]))
            pvals.name = c(pvals.name, paste("(", j, ")>0", sep = ""))
            pvals = c(pvals, mean(zStar[, j] <= z[j]))
            pvals.name = c(pvals.name, paste("(", j, ")<0", sep = ""))
		## gaussian hypothesis
        	pvals.gaussian = c(pvals.gaussian, mean(zStar.gaussian[, j] <= 0))
        	pvals.gaussian = c(pvals.gaussian, mean(zStar.gaussian[, j] >= 0))
    }
    names(pvals) = names(pvals.gaussian) = pvals.name
    
    names(z) = paste("strain", 1:nbStrains, sep = "")
    colnames(zStar) = paste("strain", 1:nbStrains, sep = "")
    colnames(zStar.gaussian) = paste("strain", 1:nbStrains, sep = "")

    return(list(distribution.estimates = zStar, distribution.estimates.gaussian = zStar.gaussian, 			p.values=pvals, p.values.gaussian=pvals.gaussian, estimates = z))
}

## Function used for formating data in the exploration of pathogen-pathogen interaction
interaction=function(genetic,threshold=0){
	n=ncol(genetic)
	g=cbind(genetic,matrix(0,nrow(genetic),n*(n-1)/2))
	k=0
	for(i in 1:n){
		g[genetic[,i] > threshold & rowSums(cbind(genetic[,-i]) > threshold)>=1,i]=0
		if(i<n){
			for(j in (i+1):n){
				k=k+1
				filter=(genetic[,i]>threshold & genetic[,j]>threshold)
				if(sum(filter)>0){
					g[filter,n+k]=rowSums(rbind(genetic[filter,c(i,j)]))
				}
				colnames(g)[n+k]=paste(colnames(genetic)[c(i,j)],collapse="x")
			}
		}
	}
	filter=rowSums(genetic)<rowSums(g)
	u=rbind(g[filter,(n+1):ncol(g)])
	g[filter,(n+1):ncol(g)]=u/rowSums(u)*
		(rowSums(rbind(genetic[filter,]))-rowSums(rbind(g[filter,1:n])))
	g=g[,colSums(g)>0]
	return(g)
}

## Function used for formating data in the exploration of pathogen-environment interaction
environment=function(strain,covar,type,name=NULL){
	if(type=="median"){ threshold=median(covar) }
	if(type=="median.specific"){ threshold=median(covar[strain>0]) }
	if(type=="mean"){ threshold=mean(covar) }
	if(type=="mean.specific"){ threshold=mean(covar[strain>0]) }
	s1=s0=strain
	s1[covar<threshold]=0
	s0[covar>=threshold]=0
	return(cbind(s1,s0))
}

