
## R code contributing to the analyses using StrainRanking and presented in Bandhari et al. (prep., Climatic fluctuations modulate eco-evolutionary processes associated with pathogen dynamics and continue to fuel bacterial spot epidemics in tomato).
## Main analysis code, requiring to load a set of underlying functions contained in the R file xantho-strain-ranking-fct.R
## Version 2024-09-17
## Author: Samuel Soubeyrand (INRAE, BioSP, 84914 Avignon, France; <samuel.soubeyrand@inrae.fr>)


########################################################################################
#### LOAD PACKAGES AND FUNCTIONS

library(OpenStreetMap) ## openmap, openproj
library(sp) ## SpatialPoints
library(StrainRanking)

source("xantho-strain-ranking-fct.R")


########################################################################################
#### LOAD AND FORMAT DATA 

xantho=read.table("Strain_Ranking_Data_env.csv",sep=";",header=TRUE,dec=",")
names(xantho)[4]="lat"
xantho=xantho[,!apply(xantho,2,function(u) sum(!is.na(u))==0)]

## Project data in a metric system
Blat=range(xantho$lat,na.rm=TRUE)+c(-1,1)
Blong=range(xantho$long,na.rm=TRUE)+c(-1,1)

map0KM=openmap(c(Blat[2],Blong[1]),c(Blat[1],Blong[2]))
coords=SpatialPoints(coords=cbind(xantho$long,xantho$lat),
	proj4string=CRS("+proj=longlat +ellps=WGS84"))
coordsKM=spTransform(coords,CRS=CRS(map0KM$tiles[[1]][[3]]@projargs)) #"+init=EPSG:3857"

xantho$longKM=coordsKM@coords[,1]
xantho$latKM=coordsKM@coords[,2]

## Display sites in a map
plot(map0KM)
points(coordsKM)


########################################################################################
#### RUN STRAINRANKING

## Run a loop for accounting for different values of the tuning parameters
for(p in c(0.1,0.5,0.9)){
	print(p)
	p.covar=p
	## p: Weigth of the Mid sample in the genetic frequencies (between 0 and 1)
	## p.covar: Weigth of the Mid values of the covariates in the analysis of the 
	##   pathogen-environment interactions (between 0 and 1)
	## saturation of the growth function, before the log transformation, the growth ranges 
	##   between -8 and +8
	saturation=c(-log(8),log(8))

	#### FORMAT DATA FOR THEIR UTILIZATION IN STRAINRANKING

	## Dataset formed for each year
	x1=format.data(xantho,time1=c(2020,"Mid"),time2=c(2020,"End"),prop=p,prop.covar=p.covar)
	x2=format.data(xantho,time1=c(2021,"Mid"),time2=c(2021,"End"),prop=p,prop.covar=p.covar)
	x3=format.data(xantho,time1=c(2022,"Mid"),time2=c(2022,"End"),prop=p,prop.covar=p.covar)

	## Merged annual datasets "separating in space" fields with same coordinates (with 
	##   argument noise.sd) and samples of successive years (with argument longlag)
	## The separation in space is necessary for avoiding that StrainRanking combines data
	##   from different locations or years
	x=merge.data(list(x1,x2,x3),noise.sd=10^2,longlag=2*10^6)
	x$covar$Y2020=0+(x$covar$Year==2020)
	x$covar$Y2021=0+(x$covar$Year==2021)
	x$covar$Y2022=0+(x$covar$Year==2022)
	plot(x$genetic.coord,asp=1)


	## Display fields and years that are taken into account
	## (some farms are not included because of non paired samples Mid-End)
	temp=xantho[as.numeric(rownames(x$genetic.frequencies)),c("Farm_id","Year")]
	temp
	temp[order(temp[,1]),]
	nrow(temp)

	#### APPLICATION OF STRAINRANKING TO EXPLORE MAIN EFFECTS
	#### -> MATERIAL FOR FIGURE 3A IN THE MANUSCRIPT (AS WELL AS SUPPL. FIG. S7 AND S8)

	## Set bandwidth for smoothed estimates of SC proportions (in meters)
	b=1 # the value b = 1 meter leads StrainRanking to not combine data at different locations

	## Check what this bandwidth correspond to
	DIST=dist(x$genetic.coord)
	hist(DIST/10^3,breaks=1000,xlim=c(0,1600),xlab="Distance (km)")
	abline(v=b/10^3,col=2)
	sort(DIST/10^3)[1:10]
	## Test that none of the distance are lower than b
	sum(DIST<b)==0

	## Rank the 8 Sequence Clusters (SC) and estimate their fitness
	DGdata=DGobj.rawdata(demographic.coord=x$demographic.coord,genetic.coord=x$genetic.coord,
		demographic.measures=x$demographic.measures,genetic.frequencies=x$genetic.frequencies)
	colnames(DGdata@genetic)[-(1:2)]=paste("SC",1:8,sep="")
	#### Rank SC (comparison of SC fitness)
	RS=ranking.strains.s(DGobject=DGdata, bw=b, nb.mcsimul=10^4, plots=TRUE,
		kernel.type="Quadratic",Zlim=saturation)
	RS$p.values
	cat("Correspondance:",paste(1:length(RS$estimates),
		colnames(DGdata@genetic[,-(1:2)]),sep="->"),"\n")
	#### Estimate fitness of each SC
	EFS=estim.fitness.strains(DGobject=DGdata, bw=b, nb.mcsimul=10^4, plots=TRUE,
    	kernel.type="Quadratic",option=list(xlim="adaptive",coverage=0.95,zero=TRUE),
    	Zlim=saturation)
	EFS$estimates
	EFS$p.values
	## Confidence interval of fitness under a Gaussian hypothesis
	apply(EFS$distribution.estimates.gaussian,2,quantile,c(0.025,0.975),na.rm=TRUE)

	## Rank the 3 dominant SC (SC3, SC4 and SC6) and the remaining merged SC that are 
	##   at low frequency and estimate their fitness
	## merge low freq strains
	DGdata.subresol=DGdata
	weight=colSums(DGdata@genetic[,-(1:2)])
	weight
	DGdata.subresol@genetic=cbind(DGdata.subresol@genetic[,1:2],
		DGdata.subresol@genetic[,-(1:2)][,weight>=1],
		strainLowFreq=rowSums(DGdata.subresol@genetic[,-(1:2)][,weight<1]))
	colnames(DGdata.subresol@genetic)[-(1:2)]=paste("SC",c(3,4,6,"LowFreq"),sep="")
	## weights of SC3, SC4, SC6 and SCLowFreq
	colSums(DGdata.subresol@genetic[,-(1:2)])
	#### Rank SC (comparison of SC fitness)
	RS=ranking.strains.s(DGobject=DGdata.subresol, bw=b, nb.mcsimul=10^4, plots=TRUE,
    	kernel.type="Quadratic",Zlim=saturation)
	RS$p.values
	cat("Correspondance:",paste(1:length(RS$estimates),
		colnames(DGdata.subresol@genetic[,-(1:2)]),sep="->"),"\n")
	#### Evaluation of each strain
	EFS=estim.fitness.strains(DGobject=DGdata.subresol, bw=b, nb.mcsimul=10^4, plots=TRUE,
    	kernel.type="Quadratic",option=list(xlim="adaptive",coverage=0.95,zero=TRUE),
    	Zlim=saturation)
	EFS$estimates
	EFS$p.values
	## Confidence interval of fitness under a Gaussian hypothesis
	apply(EFS$distribution.estimates.gaussian,2,quantile,c(0.025,0.975),na.rm=TRUE)

	## Store results of this section
	est=EFS$estimates
	pvals=c(RS$p.values,EFS$p.values)
	ic=apply(EFS$distribution.estimates.gaussian,2,quantile,c(0.025,0.975),na.rm=TRUE)
	EST=rbind(c(est,rep(NA,6-length(est))))
	PVALS=rbind(c(pvals,rep(NA,57-length(pvals))))
	IC=cbind(ic,matrix(NA,2,6-ncol(ic)))
	NAMES=rbind(c("SC3","SC4","SC6","SClowfreq",rep(NA,6-length(est))))
	rownames(PVALS)=rownames(EST)=rownames(NAMES)="Main effects"

	#### APPLICATION OF STRAINRANKING TO EXPLORE PATHOGEN-PATHOGEN INTERACTION
	#### -> MATERIAL FOR FIGURE 3B IN THE MANUSCRIPT (AS WELL AS SUPPL. FIG. S7 AND S8)

	## Rank "the interaction SC3xSC4: sum of SC3 and SC4 frequencies when SC3>20% & SC4>20%",
	##   SC3 when "SC3<=20% or SC4<=20%", SC4 when "SC3<=20% or SC4<=20%", SC6 and 
	##   the remaining merged low freq SC, and estimate their fitness
	interaction.threshold=0.2 ## minimal proportion for determining cases where interaction occurs
	DGdata.interaction=DGdata.subresol
	DGdata.interaction@genetic=cbind(DGdata.interaction@genetic[,1:2],
		interaction(DGdata.interaction@genetic[,2+c(1,2)],threshold=interaction.threshold),
		DGdata.interaction@genetic[,2+c(3,4)])
	colnames(DGdata.interaction@genetic)[-(1:2)]=paste("SC",c(3,4,"3xSC4",6,"LowFreq"),sep="")
	## weights of SC3, SC4, SC3xSC4, SC6 and SCLowFreq
	colSums(DGdata.interaction@genetic[,-(1:2)])
	#### Rank SC (comparison of SC fitness)
	RS=ranking.strains.s(DGobject=DGdata.interaction, bw=b, nb.mcsimul=10^4, plots=TRUE,
    	kernel.type="Quadratic",Zlim=saturation)
	RS$p.values
	cat("Correspondance:",paste(1:length(RS$estimates),
		colnames(DGdata.interaction@genetic[,-(1:2)]),sep="->"),"\n")
	#### Evaluation of each strain
	EFS=estim.fitness.strains(DGobject=DGdata.interaction, bw=b, nb.mcsimul=10^4, plots=TRUE,
	    kernel.type="Quadratic",option=list(xlim="adaptive",coverage=0.95,zero=TRUE),
    	Zlim=saturation)
	EFS$estimates
	EFS$p.values
	## Confidence interval of fitness under a Gaussian hypothesis
	apply(EFS$distribution.estimates.gaussian,2,quantile,c(0.025,0.975),na.rm=TRUE)

	## Store results of this section
	est=EFS$estimates
	pvals=c(RS$p.values,EFS$p.values)
	ic=apply(EFS$distribution.estimates.gaussian,2,quantile,c(0.025,0.975),na.rm=TRUE)
	EST=rbind(EST,c(est,rep(NA,6-length(est))))
	PVALS=rbind(PVALS,c(pvals,rep(NA,57-length(pvals))))
	IC=rbind(IC,cbind(ic,matrix(NA,2,6-ncol(ic))))
	NAMES=rbind(NAMES,c("SC3","SC4","SC3xSC4","SC6","SCLowFreq",rep(NA,6-length(est))))
	rownames(PVALS)[nrow(PVALS)]=rownames(EST)[nrow(EST)]=rownames(NAMES)[nrow(NAMES)]="Interaction"

	#### APPLICATION OF STRAINRANKING TO EXPLORE PATHOGEN-ENVIRONMENT INTERACTION
	#### -> MATERIAL FOR FIGURE 3C IN THE MANUSCRIPT (AS WELL AS SUPPL. FIG. S7 AND S8)

	## Rank "SC3 & covariate > mean value specific to SC3", 
	## "SC3 & covariate <= mean value specific to SC3", 
	## "SC4 & covariate > mean value specific to SC4", 
	## "SC4 & covariate <= mean value specific to SC4",
	## SC6 and the remaining merged low freq SC, and estimate their fitness
	est=NULL
	pvals=NULL
	ic=NULL
	for(k in 2:ncol(x$covar)){
		print(names(x$covar)[k])
		DGdata.env=DGdata.subresol
		colSums(DGdata.env@genetic[,-(1:2)])
		DGdata.env@genetic=cbind(DGdata.env@genetic[,1:2],
			environment(DGdata.env@genetic[,3],x$covar[,k],"mean.specific"),
			environment(DGdata.env@genetic[,4],x$covar[,k],"mean.specific"),
			DGdata.env@genetic[,5:6])
		colnames(DGdata.env@genetic)[-(1:2)]=paste("SC",c("3covar+","3covar-","4covar+","4covar-",
			"6","LowFreq"),sep="")
		colSums(DGdata.env@genetic[,-(1:2)])
		#### Comparison between strains
		RS=ranking.strains.s(DGobject=DGdata.env, bw=b, nb.mcsimul=10^4, plots=TRUE,
		    kernel.type="Quadratic",Zlim=saturation)
		#### Evaluation of each strain
		EFS=estim.fitness.strains(DGobject=DGdata.env, bw=b, nb.mcsimul=10^4, plots=TRUE,
		    kernel.type="Quadratic",option=list(xlim="adaptive",coverage=0.95,zero=TRUE),
	    	Zlim=saturation)
		est=rbind(est,EFS$estimates)
		ic=rbind(ic,apply(EFS$distribution.estimates.gaussian,2,quantile,c(0.025,0.975),na.rm=TRUE))
		pvals=rbind(pvals,c(RS$p.values,EFS$p.values))
	}
	rownames(pvals)=rownames(est)=names(x$covar)[-1]
	head(pvals)	
	colnames(DGdata.env@genetic)[-(1:2)]
	apply(pvals,2,function(u) mean(u<0.05,na.rm=TRUE))
	pvals[apply(pvals[,c(1:3,28:30)],1,function(u) sum(u<0.10,na.rm=TRUE)>0),]
	par(mfrow=c(3,3))
	apply(pvals,2,hist,breaks=20)

	## Store results of this section
	filter=apply(pvals[,c(1:3,28:30)],1,function(u) sum(u<0.10,na.rm=TRUE)>0)
	EST=rbind(EST,est[filter,])
	PVALS=rbind(PVALS,pvals[filter,])
	IC=rbind(IC,ic[as.numeric(sapply(1:nrow(est),function(u) c(2*(u-1)+1:2))[,filter]),])
	NAMES=rbind(NAMES,paste("SC",c("3covar+","3covar-","4covar+","4covar-","6","LowFreq"),sep=""))

	## Display fitness values
	## -> Fitness estimates in Figures 3, S7 and S8
	COL=rbind(c(rgb(0.7,0,0.8),rgb(0,0.7,0.8),"red","grey",NA,NA),
		c(rgb(0.7,0,0.8),rgb(0,0.7,0.8),"green","red","grey",NA),
		c(rgb(c(0.5,0.9),0,0.8),rgb(0,c(0.5,0.9),0.8),"red","grey"))
	pdf(paste("fitness_pctfreqMid=",p*100,"_pctcovarMid=",p.covar*100,".pdf",sep=""),
		width=7,height=10)
	par(mfrow=c(1,1),mar=c(4,0,0,0))
	plot(0,0,col="white",xlim=range(ic)+c(-1.5,0),ylim=c(nrow(EST),0)*2,xlab="Fitness",ylab="",
		axes=FALSE)
	axis(1,at=unique(round(c(IC,EST))))
	abline(v=0,col=1,lty="dashed")
	j=0
	for(i in c(2+order(abs(EST[-(1:2),4]-EST[-(1:2),3])),2,1)){
		j=j+1
		text(min(IC,na.rm=TRUE),2*nrow(EST)+1-2*j+ncol(IC)/2/ncol(IC),rownames(EST)[i],pos=2,cex=0.9)
		rank=(1:ncol(EST))[!is.na(EST[i,])] #order(EST[i,])
		points(EST[i,rank],2*nrow(EST)+1-2*j+(1:max(rank))/max(rank),col=COL[min(i,3),rank],pch=19)
		apply(rbind(1:max(rank),IC[(i-1)*2+1:2,rank]),2,function(u) 	
			lines(u[2:3],2*nrow(EST)+1-2*j+rep(u[1],2)/max(rank),col=COL[min(i,3),rank[u[1]]]))
	}	
	dev.off()

	## Display growth curves
	## -> Epidemic curves in Figures 3, S7 and S8
	pdf(paste("growthcurves_pctfreqMid=",p*100,"_pctcovarMid=",p.covar*100,".pdf",sep=""),
		width=10,height=10)
	par(mfrow=c(3,3))
	sev.max=exp(saturation[2])
	for(i in 1:nrow(EST)){
		plot(0:sev.max,pmax(0,pmin(sev.max,exp(EST[i,1])*(1+0:sev.max)-1)),asp=1,type="l",
			xlim=c(0,sev.max),ylim=c(0,sev.max),col=COL[min(i,3),1],
			xlab="Mid severity",ylab="End severity",main=rownames(EST)[i])
		sapply((1:ncol(EST))[!is.na(EST[i,])], function(j) lines(0:sev.max,
			pmax(0,pmin(sev.max,exp(EST[i,j])*(1+0:sev.max)-1)),col=COL[min(i,3),j]))
		abline(0,1,lty="dashed")
		points(x$demographic.measures)
	}
	dev.off()

	## Save results
	saveRDS(list(data=x,
		parameters=list(p=p,p.covar=p.covar,saturation=saturation,b=b),
		output=list(EST=EST,PVALS=PVALS,IC=IC,NAMES=NAMES,COL=COL)),
		paste("results_pctfreqMid=",p*100,"_pctcovarMid=",p.covar*100,".rds",sep=""))

}


