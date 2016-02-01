
#set working directory

	dir.in <- paste(dir.prj,"1_Data_SWInput",sep="/")
	dir.sw <- paste(dir.in,"swrun",sep="/")
	dir.weather <- paste(dir.prj,'1_Data_SWInput/treatments/LookupWeatherFolder',sep="/")

#site names
			   
	dir.sites <- paste('Weather/Weather_',sites,sep="")
	files.in.sites <- paste('files_v31',sites,'in',sep=".")
			   
#get input data

	swIn <- weatherList <- list()
	
	for(j in 1:length(sites)){

	  #create the input object for Soilwat runs
		swIn[[j]] <- sw_inputDataFromFiles(dir=dir.sw, files.in=files.in.sites[j]) 
		
    #get the years for the runs, defined in "/1_Data_SWInput/swrun/Input/years.in"
		years <- c(Rsoilwat31::swYears_StartYear(swIn[[j]]),Rsoilwat31::swYears_EndYear(swIn[[j]]))


		#This imports all the weather data
    weatherList[[j]] <- getWeatherData_folders(LookupWeatherFolder=dir.sw, 
								Label=dir.sites[[j]], 
    							filebasename="weath", startYear=years[1], endYear=years[2])

	}
	
	#set transpiration regions
		TRdepth <- cbind(1,3)
			colnames(TRdepth) <- c('ndx','layer')
		
		for(j in 1:length(sites)){
			Rsoilwat31::swSite_TranspirationRegions(swIn[[j]]) <- TRdepth
			}

	
	#import predictions for each soil moisture sensor
	
			
		#get sensor layer numbers
		
			sens.Lyr <- read.csv(paste(dir.in,   #matrix (in list form) with layer numbers for each sensor depth
									   "fielddata",
									   "FieldSensors_MappedTo_SoilWatLayers.csv",sep='/'),
								 header=TRUE)
								 
		#get observations of soil water content

			obs <- list()
			
			for(j in 1:length(sites)){
				obs[[j]] <- read.csv(paste(dir.in,
										   'fielddata',
										   sites[j],
										   '3_Field_Output',
										   paste(sites[j],'SoilWater.csv',sep="_"),sep="/"),
									 header=TRUE)
									 
				obs[[j]] <- obs[[j]][,c(1,2,grep('VWC',colnames(obs[[j]])))]
				obs[[j]] <- cbind(obs[[j]][,'Date'],
								  year=substr(obs[[j]][,'Date'],1,4),
								  day = obs[[j]][,'doy'],
								  obs[[j]][,-1])
				obs[[j]][-1,'day'] <-  cumsum(abs(diff(obs[[1]][, "doy"])-1))+obs[[j]][-1,'day']
				
				obs[[j]] <- obs[[j]][obs[[j]][,'year'] %in% years[1]:years[2],]
				
				}
				
			#columns of a given matrix within obs with soil moisture data
      obs.col <- list()
			
			for(j in 1:length(sites))
				obs.col[[j]] <- grep('VWC',colnames(obs[[j]]))

##Priors and Proposal Information

	priors <- read.csv(paste(dir.in,paste('priors',runname,vers,'csv',sep='.'),sep='/'),header=TRUE)

	#composition for each site
	comp.prior.mean <- get.priors('composition','all')
		names(comp.prior.mean) <- c("grass", "shrub", "tree","forb","bareground")
  
  #prior for soils textures is taken as the mean across the profile
  soils.prior.mean <- apply(swSoils_Layers(swIn[[j]])[,c('sand_frac','clay_frac')],2,mean)
	soils.prior.V <- diag(rep(0.0001,2))

	#deep drainage
    drain.prior.mean <- get.priors('drain','all')[1]
	  drain.prior.var  <- get.priors('drain','all')[2]
	  drain.prior.lo   <- get.priors('drain','all')[3]
	  drain.prior.hi   <- get.priors('drain','all')[4]
	
  #evaporation for each site
	  alpha.prior.mean <- get.priors('evap','all')[c(1,5)]
	  alpha.prior.var  <- get.priors('evap','all')[c(2,6)]
	  alpha.prior.lo   <- get.priors('evap','all')[c(3,7)]
	  alpha.prior.hi   <- get.priors('evap','all')[c(4,8)]
	
	#critical soil water potential
	#see WoodlandPars.xlsx for sources of prior information
		crit.tree.prior.mean  <- get.priors('crit','tree')[1]
		crit.tree.prior.var   <- get.priors('crit','tree')[2]
		crit.tree.prior.lo    <- get.priors('crit','tree')[3]
		crit.tree.prior.hi    <- get.priors('crit','tree')[4]

	#based on expert opinion?
		crit.shrub.prior.mean  <- get.priors('crit','shrub')[1]
		crit.shrub.prior.var   <- get.priors('crit','shrub')[2]
		crit.shrub.prior.lo    <- get.priors('crit','shrub')[3]
		crit.shrub.prior.hi    <- get.priors('crit','shrub')[4]

		crit.grass.prior.mean  <- get.priors('crit','grass')[1]
		crit.grass.prior.var   <- get.priors('crit','grass')[2]
		crit.grass.prior.lo    <- get.priors('crit','grass')[3]
		crit.grass.prior.hi    <- get.priors('crit','grass')[4]

		crit.forb.prior.mean  <- get.priors('crit','forb')[1]
		crit.forb.prior.var   <- get.priors('crit','forb')[2]
		crit.forb.prior.lo    <- get.priors('crit','forb')[3]
		crit.forb.prior.hi    <- get.priors('crit','forb')[4]

	#phenology -- Litter, Biomass, Live_pct, and LAI_Conv
		tmp <- cbind(get.priors('Litter','tree'),get.priors('Biomass','tree'),
				     get.priors('Live_pct','tree'),get.priors('LAI_conv','tree'))
		phen.tree.prior.mean  <-      tmp[1,]
		phen.tree.prior.var   <- diag(tmp[2,])
		phen.tree.prior.lo    <-      tmp[3,]
		phen.tree.prior.hi    <-      tmp[4,]

	#based on expert opinion?
		tmp <- cbind(get.priors('Litter','shrub'),get.priors('Biomass','shrub'),
				     get.priors('Live_pct','shrub'),get.priors('LAI_conv','shrub'))
		phen.shrub.prior.mean  <-      tmp[1,]
		phen.shrub.prior.var   <- diag(tmp[2,])
		phen.shrub.prior.lo    <-      tmp[3,]
		phen.shrub.prior.hi    <-      tmp[4,]

		tmp <- cbind(get.priors('Litter','grass'),get.priors('Biomass','grass'),
				     get.priors('Live_pct','grass'),get.priors('LAI_conv','grass'))
		phen.grass.prior.mean  <-      tmp[1,]
		phen.grass.prior.var   <- diag(tmp[2,])
		phen.grass.prior.lo    <-      tmp[3,]
		phen.grass.prior.hi    <-      tmp[4,]

		tmp <- cbind(get.priors('Litter','forb'),get.priors('Biomass','forb'),
				     get.priors('Live_pct','forb'),get.priors('LAI_conv','forb'))
		phen.forb.prior.mean  <-      tmp[1,]
		phen.forb.prior.var   <- diag(tmp[2,])
		phen.forb.prior.lo    <-      tmp[3,]
		phen.forb.prior.hi    <-      tmp[4,]

		rm(tmp)
		
	#transpiration coefficients model priors
	#see Schenk et al. 2002; d50 and c
		tmp <- cbind(get.priors('d50','tree'),get.priors('c','tree'))
		
		beta.tree.prior.mean  <- tmp[1,]
		beta.tree.prior.var   <- tmp[2,]
		beta.tree.prior.lo    <- tmp[3,]
		beta.tree.prior.hi    <- tmp[4,]

		tmp <- cbind(get.priors('d50','shrub'),get.priors('c','shrub'))
		
		beta.shrub.prior.mean  <- tmp[1,]
		beta.shrub.prior.var   <- tmp[2,]
		beta.shrub.prior.lo    <- tmp[3,]
		beta.shrub.prior.hi    <- tmp[4,]

		tmp <- cbind(get.priors('d50','grass'),get.priors('c','grass'))
		
		beta.grass.prior.mean  <- tmp[1,]
		beta.grass.prior.var   <- tmp[2,]
		beta.grass.prior.lo    <- tmp[3,]
		beta.grass.prior.hi    <- tmp[4,]

		tmp <- cbind(get.priors('d50','forb'),get.priors('c','forb'))
		
		beta.forb.prior.mean  <- tmp[1,]
		beta.forb.prior.var   <- tmp[2,]
		beta.forb.prior.lo    <- tmp[3,]
		beta.forb.prior.hi    <- tmp[4,]

		rm(tmp)
				
#Initialize proposal objects
	#prior covariances  and bounds for phenology proporsals
		phen.tree.prior.V <- diag(diag(phen.tree.prior.var)[c(4,rep(c(1,2,3),length(sites)))])
		phen.tree.lo  <- phen.tree.prior.lo[c(4,rep(1:3,times=length(sites)))]
		phen.tree.hi  <- phen.tree.prior.hi[c(4,rep(1:3,times=length(sites)))]

		phen.shrub.prior.V <- diag(diag(phen.shrub.prior.var)[c(4,rep(c(1,2,3),length(sites)))])
		phen.shrub.lo  <- phen.shrub.prior.lo[c(4,rep(1:3,times=length(sites)))]
		phen.shrub.hi  <- phen.shrub.prior.hi[c(4,rep(1:3,times=length(sites)))]

		phen.grass.prior.V <- diag(diag(phen.grass.prior.var)[c(4,rep(c(1,2,3),length(sites)))])
		phen.grass.lo  <- phen.grass.prior.lo[c(4,rep(1:3,times=length(sites)))]
		phen.grass.hi  <- phen.grass.prior.hi[c(4,rep(1:3,times=length(sites)))]

		phen.forb.prior.V <- diag(diag(phen.forb.prior.var)[c(4,rep(c(1,2,3),length(sites)))])
		phen.forb.lo  <- phen.forb.prior.lo[c(4,rep(1:3,times=length(sites)))]
		phen.forb.hi  <- phen.forb.prior.hi[c(4,rep(1:3,times=length(sites)))]
		

	#prior covariances and bounds for depth proposals-- evaporation
  	alpha.prior.V <- diag(rep(alpha.prior.var,length(sites)))
	  alpha.lo  <- rep(alpha.prior.lo,length(sites))
	  alpha.hi  <- rep(alpha.prior.hi,length(sites))
	
	#prior covariances and bounds for depth proposals -- transpiration
			beta.tree.prior.V <- diag(rep(beta.tree.prior.var,length(sites)))
			beta.tree.lo  <- rep(beta.tree.prior.lo,length(sites))
			beta.tree.hi  <- rep(beta.tree.prior.hi,length(sites))
			beta.shrub.prior.V <- diag(rep(beta.shrub.prior.var,length(sites)))
			beta.shrub.lo  <- rep(beta.shrub.prior.lo,length(sites))
			beta.shrub.hi  <- rep(beta.shrub.prior.hi,length(sites))
			beta.grass.prior.V <- diag(rep(beta.grass.prior.var,length(sites)))
			beta.grass.lo  <- rep(beta.grass.prior.lo,length(sites))
			beta.grass.hi  <- rep(beta.grass.prior.hi,length(sites))
			beta.forb.prior.V <- diag(rep(beta.forb.prior.var,length(sites)))
			beta.forb.lo  <- rep(beta.forb.prior.lo,length(sites))
			beta.forb.hi  <- rep(beta.forb.prior.hi,length(sites))

	#prior covariances and bounds for critical soil water potential
	
			crit.tree.prior.V <- matrix(crit.tree.prior.var,1,1)
			crit.tree.lo  <- rep(crit.tree.prior.lo,length(sites))
			crit.tree.hi  <- rep(crit.tree.prior.hi,length(sites))
			crit.shrub.prior.V <- matrix(crit.shrub.prior.var,1,1)
			crit.shrub.lo  <- rep(crit.shrub.prior.lo,length(sites))
			crit.shrub.hi  <- rep(crit.shrub.prior.hi,length(sites))
			crit.grass.prior.V <- matrix(crit.grass.prior.var,1,1)
			crit.grass.lo  <- rep(crit.grass.prior.lo,length(sites))
			crit.grass.hi  <- rep(crit.grass.prior.hi,length(sites))
			crit.forb.prior.V <- matrix(crit.forb.prior.var,1,1)
			crit.forb.lo  <- rep(crit.forb.prior.lo,length(sites))
			crit.forb.hi  <- rep(crit.forb.prior.hi,length(sites))
			
			if(length(sites)>1){
				crit.tree.prior.V <- diag(rep(crit.tree.prior.var,length(sites)))
				crit.shrub.prior.V <- diag(rep(crit.shrub.prior.var,length(sites)))
				crit.grass.prior.V <- diag(rep(crit.grass.prior.var,length(sites)))
				crit.forb.prior.V <- diag(rep(crit.forb.prior.var,length(sites)))
			}

	#prior covariances and bounds for deep drainage
	
	drain.prior.V <- matrix(drain.prior.var,1,1)
	drain.lo  <- rep(drain.prior.lo,length(sites))
	drain.hi  <- rep(drain.prior.hi,length(sites))
	
	if(length(sites)>1){
	  drain.prior.V <- diag(rep(drain.prior.var,length(sites)))
	}
	
	#process error
		sig <- 0.1   #initial value of the variance
			s1 <- 1000 #parameter for the inverse gamma variance prior
			s2 <- sig * (s1 - 1) #parameter for the inverse gamma variance prior
			
		rho  <- 0.0     #initial value of the autocorrelation
		rhoV <- 0.0001  #proposal variance for metropolis algorithm
		rho.lo <- -0.75 #minimum autocorrelation
		rho.hi <- 0.75  #maximum autocorrelation
		
		rho.prior.mean <- 0    #prior mean autocorrelation
		rho.prior.var  <- 0.01 #prior variance autocorrelation (small values mean stronger prior)
		
		
##set initial values
	
	comp <- (comp.prior.mean/sum(comp.prior.mean))[rep(1:5,length(sites))]
		names(comp) <- paste('comp',rep(c('grass','shrub','tree','forb','bareground'),
										times = length(sites)),
									rep(sites,each=5),sep='.')
  
  for(j in 1:length(sites)){
  
    tmp.soils <- swSoils_Layers(swIn[[j]])[,names(soils.prior.mean)]
      tmp.layers <- c(which(diff(tmp.soils[,1])!=0 | diff(tmp.soils[,1]!=0)),nrow(tmp.soils))
      tmp.ind    <- findInterval(1:nrow(tmp.soils),tmp.layers+1)+1
      tmp.soils <- as.vector(tmp.soils[tmp.layers,])
    
    names(tmp.soils) <- paste(rep(names(soils.prior.mean),each=length(tmp.soils)/2),
                                 paste('depth',rep(1:(length(tmp.soils)/2),times=2),sep=""),
                                 rep(sites[j],length(tmp.soils)),sep=".")
    
    
    if(j == 1){
      
      soils <- tmp.soils
      soils.layers <- rep(tmp.layers,times=2)
      soils.ind <- tmp.ind
      
    }

    if(j > 1){
      
      soils <- c(soils,tmp.soils)
      soils.layers <- c(soils.layers,rep(tmp.layers,times=2))
      soils.ind <- c(soils.ind,tmp.ind)
      
    }
    
  }

  soils.obs <- soils #assume that the values in 
  
  #ranges of acceptable soil texture values for each layer. Basically, here we assume that the 
  #actual soil texture proportions cannot differ by more than 0.25 from the inputs in 
  #"/1_Data_SWInput/swrun/Input/soils"
  soils.lo <- soils.obs - .25
  soils.hi <- soils.obs + .25
  
  soils.lo[soils.lo<.01] <- .01
  soils.hi[soils.hi>.99] <- .99
  
  soils.Vpar <- diag(soils.Vpar,length(soils))
  
	alpha <- rep(alpha.prior.mean,length(sites))
	  names(alpha) <- paste(c('evap.max','evap.rate'),sites,sep='.')

  beta.tree <- rep(beta.tree.prior.mean,length(sites))
		names(beta.tree) <- paste(c('d50.tree','c.tree'),sites,sep='.')
	beta.shrub <- rep(beta.shrub.prior.mean,length(sites))
		names(beta.shrub) <- paste(c('d50.shrub','c.shrub'),sites,sep='.')
	beta.forb <- rep(beta.forb.prior.mean,length(sites))
		names(beta.forb) <- paste(c('d50.forb','c.forb'),sites,sep='.')
	beta.grass <- rep(beta.grass.prior.mean,length(sites))
		names(beta.grass) <- paste(c('d50.grass','c.grass'),sites,sep='.')

	crit.tree <- rep(crit.tree.prior.mean,length(sites))
		names(crit.tree) <- paste(c('crit.tree'),sites,sep='.')
	crit.shrub <- rep(crit.shrub.prior.mean,length(sites))
		names(crit.shrub) <- paste(c('crit.shrub'),sites,sep='.')
	crit.forb <- rep(crit.forb.prior.mean,length(sites))
		names(crit.forb) <- paste(c('crit.forb'),sites,sep='.')
	crit.grass <- rep(crit.grass.prior.mean,length(sites))
		names(crit.grass) <- paste(c('crit.grass'),sites,sep='.')

	drain <- rep(drain.prior.mean,length(sites))
	  names(drain) <- paste(c('drain'),sites,sep='.')
	
	phen.tree <- phen.tree.prior.mean[c(4,rep(c(1,2,3),length(sites)))]
		names(phen.tree) <- 
			c('LAI_conv.tree.ALL',
	          apply(expand.grid(c('Litter.tree','Biomass.tree','Live_pct.tree'),sites),
	                1,paste,collapse="."))
	phen.shrub <- phen.shrub.prior.mean[c(4,rep(c(1,2,3),length(sites)))]
		names(phen.shrub) <- 
			c('LAI_conv.shrub.ALL',
	          apply(expand.grid(c('Litter.shrub','Biomass.shrub','Live_pct.shrub'),sites),
	                1,paste,collapse="."))
	phen.forb <- phen.forb.prior.mean[c(4,rep(c(1,2,3),length(sites)))]
		names(phen.forb) <- 
			c('LAI_conv.forb.ALL',
	          apply(expand.grid(c('Litter.forb','Biomass.forb','Live_pct.forb'),sites),
	                1,paste,collapse="."))
	phen.grass <- phen.grass.prior.mean[c(4,rep(c(1,2,3),length(sites)))]
		names(phen.grass) <- 
			c('LAI_conv.grass.ALL',
	          apply(expand.grid(c('Litter.grass','Biomass.grass','Live_pct.grass'),sites),
	                1,paste,collapse="."))


##proposals -- what parameters to propose

	#identify part of model to allow to vary
	prop.names <- unique(as.character(priors[priors[,'include']==1,'par.names']))
				   
	var.names  <- list(phenology.tree = names(phen.tree),
	                   phenology.shrub = names(phen.shrub),
	                   phenology.forb = names(phen.forb),
	                   phenology.grass = names(phen.grass),
	                   ETcoef.tree = names(beta.tree),
	                   ETcoef.shrub = names(beta.shrub),
	                   ETcoef.forb = names(beta.forb),
	                   ETcoef.grass = names(beta.grass),
	                   crit.tree = names(crit.tree),
	                   crit.shrub = names(crit.shrub),
	                   crit.forb = names(crit.forb),
	                   crit.grass = names(crit.grass),
	                   composition = names(comp),
                     evap = names(alpha),
                     soils = names(soils),
                     drain = names(drain))
	gibbs.names <- c('pTgibbs','pSgibbs','pFgibbs','pGgibbs',
					 'bTgibbs','bSgibbs','bFgibbs','bGgibbs',
					 'cTgibbs','cSgibbs','cFgibbs','cGgibbs',
					 'cgibbs','agibbs','sgibbs','dgibbs')
		gibbs.names <- c(gibbs.names[names(var.names) %in% prop.names],'vgibbs')
		
	    var.names <- var.names[names(var.names) %in% prop.names]
	    var.names <- var.names[match(prop.names,names(var.names))]
	    
	
	#get list of functional groups
	groups <- c('tree','shrub','forb','grass')
		groups <- groups[groups %in% unlist(strsplit(prop.names,'\\.'))]
		
	

