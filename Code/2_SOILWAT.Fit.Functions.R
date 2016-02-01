
#Additional functions used for Soilwat fitting

	


tnorm <- function(n,lo,hi,mu,sig){   
  
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig)
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == Inf]  <- lo[z == Inf]
  z[z == -Inf] <- hi[z == -Inf]
  z
}


tnorm.mvt <- function(avec,muvec,smat,lo,hi){   
  
  # truncated multvariate normal
  # muvec is the vector of means
  # smat is the covariance matrix 
  
  # smat <- nearPD(smat)$mat
  
  # avec    <- muvec
  avec[1] <- tnorm(1,lo[1],hi[1],muvec[1],sqrt(smat[1,1]))
  
  for(k in 2:length(muvec)){
    skk <- smat[-k,-k]
    #  skk <- as.matrix(nearPD(smat[-k,-k])$mat)
    
    testv <- try(chol(skk),T)
    if(inherits(testv,'try-error')){
      
      avec[k] <- 0
      next
      
    }
    
    piece1 <- smat[-k,k] %*% chol2inv(testv)
    muk <- muvec[k] + piece1 %*% (avec[-k] - muvec[-k])
    sgk <- as.numeric(smat[k,k] - piece1 %*% smat[k,-k])
    if(sgk < .000000001)sgk <- .000000001
    avec[k] <- tnorm(1,lo[k],hi[k],muk,sqrt(sgk))
  }
  avec
}

get.priors <- function(param,group) 
	as.numeric(unlist(strsplit(as.character(priors[priors[,'parameter']==param &
												   priors[,'group']==group,'pars']),' ')))

T_depth.Jackson1996 <- function(d, params){ #d = depth, params = c(beta); NO CURRENTLY IN USE
	
	dseq   <- seq(1,max(d),by=1)
	
	r.dseq <- 1 - params[1]^dseq
		r.d    <- r.dseq[d]
		
		r.d[-1] <- diff(r.d)
		r.d <- r.d/sum(r.d)
		
	list(r.d = r.d, dseq = dseq, r.dseq = r.dseq)
	
	}
	
T_depth.Zeng2001 <- function(d, params){ #d = depth, params = c(a, b); NO CURRENTLY IN USE

	dseq   <- seq(1,max(d),by=1)
	
	r.dseq <- 1 - 0.5 * (exp(-params[1]*dseq) + exp(-params[2]*dseq))
		r.d    <- r.dseq[d]

		r.d[-1] <- diff(r.d)
		r.d <- r.d/sum(r.d)
		
	list(r.d = r.d, dseq = dseq, r.dseq = r.dseq)

	}
	
T_depth.Schenk2002 <- function(d, params){ #d = depth, params = c(d50, c)

	dseq   <- seq(1,max(d),by=1) #dseq = 1 cm depths from 1 to the maximum depth
	
	r.dseq <- (1 + (dseq / params[1])^params[2])^-1 #r.dseq = cumulative transpiration proportions
		r.d    <- r.dseq[d]

		r.d[-1] <- diff(r.d)
		r.d <- r.d/sum(r.d) # differences and standardized r.d, which is the proportion of transpiration from that 1-cm layer
		
	list(r.d = r.d, dseq = dseq, r.dseq = r.dseq) 

	}
	
E_depth_Wythers1999 <- function(d, sand, clay, params){

	# soil texture influence: Wythers KR, Lauenroth WK, Paruelo JM (1999) 
	#                               Bare-Soil Evaporation Under Semiarid Field Conditions. 
	#                               Soil Science Society of America Journal, 63, 1341-1349. 
	
	#sand.mean and clay.mean is mean sand and clay in relavent layers
	
	# max = 15 cm: Torres EA, Calera A (2010) Bare soil evaporation under high 
	#                evaporation demand: a proposed modification to the FAO-56 model. 
	#                Hydrological Sciences Journal-Journal Des Sciences Hydrologiques, 55, 303-315.
	
	ld <- d             #place depths in cm into ld
	ld0 <- c(0,d[-length(d)])
  
  lw <- c(0,diff(ld)) #get widths by taking the differences of the vector ld
	
	bsEvap.depth.min <- min(ld)    #evaporation must impact the first soil layer
	bsEvap.depth.max <- params[1]  #max depth equals Dmax (a parameter)
			
	temp <- 4.1984+0.6695*sand^2+168.7603*clay^2	#calculate potential max depth of evaporation based on mean clay and sand of profile
	
	stopifnot(bsEvap.depth.min < bsEvap.depth.max) 

	temp <- matrix(data=c(temp, 
						  rep(bsEvap.depth.min, times=length(temp)), 
						  rep(bsEvap.depth.max, times=length(temp))), 
				   ncol=3, byrow=FALSE) 

	bsEvap.depth <- apply(temp, 
						  MARGIN=1, 
						  FUN=function(x) min(c(x[3], max(x[1:2], na.rm=TRUE)), na.rm=TRUE)) 
  
  
  bsEvap.coef <- exp(1-params[2]*ld0/bsEvap.depth)/exp(1) - exp(1-params[2]*ld/bsEvap.depth)/exp(1)
	  bsEvap.coef <- bsEvap.coef / sum(bsEvap.coef)
	
	bsEvap.coef[d < params[1]] <- 0
	  bsEvap.coef <- bsEvap.coef / sum(bsEvap.coef)
	
  return(bsEvap.coef)
  
	}
	
#calculate density using parallel cluster

#for Soilwat model parameters
calc.dens.params <- function(j,i){
				
			#get timesteps for which there are data after the first year	
				keep <- which(is.finite(obs[[j]][,i]))
				if(gibbs.thin) keep <- keep[keep %in% obs.thin]
        keep <- keep[keep>365]
			
			#initialize variables to hold summed log densities
				p.old <- p.new <- 0
				
			#continue only if there is at least one observation after the first year
				if(length(keep)>0) {
				
			#get inverse Cholesky decompositon of covariance matrix; used by 
				rooti <- backsolve(chol(sig*rmat(nrow(obs[[j]]),rho)[keep,keep]),
								   diag(length(keep)))
			
			#calculate density given old SOILWAT parameters
				p.old <- sum(bayesm::lndMvn(    x = obs[[j]][keep,i],#contribution of site j
						 			   mu = pred.vwc.old[[j]][keep,i],
						 			  rooti = rooti),
						 			  #lower = rep(0,length(keep)),
						 			  #upper = rep(1,length(keep)),
						 			   # log = TRUE),
						 			    na.rm=TRUE)

			#calculate density given new SOILWAT parameters
				p.new <- sum(bayesm::lndMvn(    x = obs[[j]][keep,i],#contribution of site j
						 			   mu = pred.vwc.new[[j]][keep,i],
						 			  rooti = rooti),
						 			  #lower = rep(0,length(keep)),
						 			  #upper = rep(1,length(keep)),
						 		    	#log = TRUE),
						 		    	na.rm=TRUE)
				}
			
			#output
			list(p.old = p.old, p.new = p.new)
			}

#for error terms
calc.dens.rho <- function(j,i){ 
				
			#get timesteps for which there are data after the first year	
				keep <- which(is.finite(obs[[j]][,i]))
				if(gibbs.thin) keep <- keep[keep %in% obs.thin]
				keep <- keep[keep>365]
			
			#initialize variables to hold summed log densities
				p.old <- p.new <- 0
				
			#continue only if there is at least one observation after the first year
				if(length(keep)>0) {
				
			#get inverse Cholesky decompositon of covariance matrix
				rooti <- backsolve(chol(sig*rmat(nrow(obs[[j]]),rho)[keep,keep]),
								   diag(length(keep)))
			
			#calculate density given old SOILWAT parameters
				p.old <- sum(bayesm::lndMvn(    x = obs[[j]][keep,i],#contribution of site j
						 			   mu = pred.vwc.old[[j]][keep,i],
						 			  rooti = rooti),
						 			  #lower = rep(0,length(keep)),
						 			  #upper = rep(1,length(keep)),
						 			   # log = TRUE),
						 			    na.rm=TRUE)

			#propose new sigma
				
				
			#get inverse Cholesky decompositon of covariance matrix
				rooti <- backsolve(chol(sig*rmat(nrow(obs[[j]]),rho.prop)[keep,keep]),
								   diag(length(keep)))

			#calculate density given new SOILWAT parameters
				p.new <- sum(bayesm::lndMvn(    x = obs[[j]][keep,i],#contribution of site j
						 			   mu = pred.vwc.old[[j]][keep,i],
						 			  rooti = rooti),
						 			  #lower = rep(0,length(keep)),
						 			  #upper = rep(1,length(keep)),
						 		    	#log = TRUE),
						 		    	na.rm=TRUE)
				}
			
			
			#output
			list(p.old = p.old, p.new = p.new)
			}

#inverse of rho matrix
inv.rmat <- function(n,rho=0.9){
	rinv <- diag((1 + rho^2),n,n)
	rinv[1,1] <- 1
	rinv[n,n] <- 1
	rinv[row(rinv) == (col(rinv)-1)] <- -rho
	rinv[row(rinv) == (col(rinv)+1)] <- -rho
	return(rinv)
	}		
	
#calculate R from covariance matrix
rmat  <- function(n, rho = 0.9) {
	mat <- diag(rep(1,n))
    mat <- rho^abs(row(mat)-col(mat))
    ((1 - rho^2)^-1)*mat
	}

#execute soilwat based on Input objsect	
	get.pred <- function(Input,file.in=TRUE){
		swOut <- list()
	
		for(j in 1:length(sites)){
			
			f.in <- ""
			if(file.in) f.in <- files.in.sites[j]
		
			swOut[[j]] <- Rsoilwat31::sw_exec(inputData=Input[[j]], 
						     weatherList=weatherList[[j]], 
						     dir="", 
						     files.in=f.in, 
				    		 echo=FALSE, 
					    	 quiet=FALSE)
		}
		return(swOut)
	}

#function extracts parameters for phenology for a given function group GR
get.phenology <- function(GR,Input){

	#existing vegetation parameter values
  old <- get(paste('phen',GR,sep='.'))
  
  #propose new vegeattion parameter values using a truncated normal distribution
    new <- tnorm.mvt(avec = old,
	                  muvec = old,
	                   smat = get(paste('phen',GR,'var',sep='.')),
	                     lo = get(paste('phen',GR,'lo',sep='.')),
	                     hi = get(paste('phen',GR,'hi',sep='.'))) 
	
	#obejct to hold both new and old parameter values
    pars <- rbind(as.vector(old),as.vector(new))
	
	
	phen.prop <- list()

	for(j in 1:length(sites)){

		if(GR == 'tree')  phen.prop[[j]] <- swProd_MonProd_tree(Input[[j]])
		if(GR == 'shrub')  phen.prop[[j]] <- swProd_MonProd_shrub(Input[[j]])
		if(GR == 'forb')  phen.prop[[j]] <- swProd_MonProd_forb(Input[[j]])
		if(GR == 'grass')  phen.prop[[j]] <- swProd_MonProd_grass(Input[[j]])
    
		#preserves phenology and rescales the seasonal pattern to a new maximum
				
		phen.prop[[j]][,4] <- new[1] * (phen.prop[[j]][,4]/max(phen.prop[[j]][,4])) 
						
		for(i in 1:3){
				
			phen.prop[[j]][,i] <- new[i+1 + 3*(j-1)] * (phen.prop[[j]][,i]/max(phen.prop[[j]][,i]))
				
			}
		}
				    
	list(pars = pars, phen.prop = phen.prop)
	
	}
	
#function extracts parameters for evaporation
get.Evap <- function(Input){
  
  
  old <- alpha
  new <- tnorm.mvt(avec = old,
                  muvec = old,
                   smat = alpha.var,
                     lo = alpha.lo,
                     hi = alpha.hi) 
  
  pars <- rbind(as.vector(old),as.vector(new))
  
  evap.prop <- list()
  
  for(j in 1:length(sites)){
    
    d  <- swSoils_Layers(Input[[j]])[,1]
    sand  <- swSoils_Layers(Input[[j]])[,"sand_frac"]
    clay  <- swSoils_Layers(Input[[j]])[,"clay_frac"]
    
    evap.prop[[j]] <- E_depth_Wythers1999(d,sand,clay,new)
  }
  
  list(pars = pars, evap.prop = evap.prop) ##return here##
  
  
}


#function extracts parameters for phenology for a given functional group GR
get.Trans <- function(GR,Input){
	
	
	old <- get(paste('beta',GR,sep="."))
	new <- tnorm.mvt(avec = old,
	                muvec = old,
	                 smat = get(paste('beta',GR,'var',sep='.')),
	                   lo = get(paste('beta',GR,'lo',sep='.')),
	                   hi = get(paste('beta',GR,'hi',sep='.'))) 
	
	pars <- rbind(as.vector(old),as.vector(new))
	
	ET.prop <- list()

	for(j in 1:length(sites)){
				
		d  <- swSoils_Layers(Input[[j]])[,1]
											
		ET.prop[[j]] <- T_depth.Schenk2002(d,new[2*(j-1)+1:2])
		}

	list(pars = pars, ET.prop = ET.prop) ##return here##
	
	
	}

#function extracts parameters for critical soil water potential

get.crit <- function(GR, Input){

	old <- get(paste('crit',GR,sep='.'))
	
  if(length(old)>1)
    new <- tnorm.mvt(avec = old,
	                  muvec = old,
	                   smat = get(paste('crit',GR,'var',sep='.')),
	                     lo = get(paste('crit',GR,'lo',sep='.')),
	                     hi = get(paste('crit',GR,'hi',sep='.'))) 
	
  if(length(old)==1)
    new <- tnorm(n=1,
               lo = get(paste('crit',GR,'lo',sep='.')),
               hi = get(paste('crit',GR,'hi',sep='.')),
               mu = old,
              sig = sqrt(get(paste('crit',GR,'var',sep='.')))) 
	
	pars <- rbind(as.vector(old),as.vector(new))

	return(pars)
	
	}
	
#function extracts and proposes parameters for deep drainage

get.drain <- function(Input){
  
  old <- drain
  
  if(length(old)>1)
    new <- tnorm.mvt(avec = old,
                     muvec = old,
                     smat = get(paste('drain','var',sep='.')),
                     lo = get(paste('drain','lo',sep='.')),
                     hi = get(paste('drain','hi',sep='.'))) 
  
  if(length(old)==1)
    new <- tnorm(n=1,
                 lo = get(paste('drain','lo',sep='.')),
                 hi = get(paste('drain','hi',sep='.')),
                 mu = old,
                 sig = sqrt(get(paste('drain','var',sep='.')))) 
  
  pars <- rbind(as.vector(old),as.vector(new))
  
  return(pars)
  
}

	params.prop <- function(Input){
	
		prop.par <- vector('list',length=length(prop.names))
		names(prop.par) <- prop.names

		if('drain' %in% prop.names){
		  
		  tmp <- get.drain(Input)
		  
		  for(j in 1:length(sites))
		    Input <- model.set(Input = Input, set.names = "drain", set.par = tmp[2,j], j = j)	
		  
		  #save parameters into prop.par
		  prop.par[[which(prop.names == 'drain')]] <- tmp
		  colnames(prop.par[[which(prop.names == "drain")]]) <- 
		    var.names[[which(prop.names == "drain")]]
		  rownames(prop.par[[which(prop.names == "drain")]]) <- c('old','new')
		  
		  #delete temporary object
		  rm(tmp)
		  
		}
		
		for(GR in groups){

			if(paste('crit',GR,sep='.') %in% prop.names){
			
				tmp <- get.crit(GR,Input)
				
				for(j in 1:length(sites))
				  Input <- model.set(Input = Input, set.names = paste('crit',GR,sep='.'), set.par = tmp[2], j = j)	
				
			#save parameters into prop.par
				prop.par[[which(prop.names == paste('crit',GR,sep='.'))]] <- tmp
					colnames(prop.par[[which(prop.names == paste('crit',GR,sep='.'))]]) <- 
						var.names[[which(prop.names == paste('crit',GR,sep='.'))]]
					rownames(prop.par[[which(prop.names == paste('crit',GR,sep='.'))]]) <- c('old','new')

			#delete temporary object
				rm(tmp)
			
			}

			
			if(paste('phenology',GR,sep='.') %in% prop.names){
				
				tmp <- get.phenology(GR,Input)
								
				for(j in 1:length(sites))
				  Input <- model.set(Input = Input, set.names = paste('phenology',GR,sep='.'), set.par = tmp$phen.prop[[j]], j = j)	
				
			#save parameters into prop.par
				prop.par[[which(prop.names == paste('phenology',GR,sep='.'))]] <- tmp$pars
					colnames(prop.par[[which(prop.names == paste('phenology',GR,sep='.'))]]) <- 
						var.names[[which(prop.names == paste('phenology',GR,sep='.'))]]
					rownames(prop.par[[which(prop.names == paste('phenology',GR,sep='.'))]]) <- c('old','new')

			#delete temporary object
				rm(tmp)
			}

		if(paste('ETcoef',GR,sep='.') %in% prop.names){ ##needs to be based on some sort of exponential decay
									 ## perhaps we predict the root distribution as a parametric function
		
			tmp <- get.Trans(GR,Input)

      for(j in 1:length(sites))
			  Input <- model.set(Input = Input, set.names = paste('ETcoef',GR,sep='.'), set.par = tmp$ET.prop[[j]]$r.d/sum(tmp$ET.prop[[j]]$r.d), j = j)	
			
			#save parameters into prop.par
				prop.par[[which(prop.names == paste('ETcoef',GR,sep='.'))]] <- tmp$pars
					colnames(prop.par[[which(prop.names == paste('ETcoef',GR,sep='.'))]]) <- 
						var.names[[which(prop.names == paste('ETcoef',GR,sep='.'))]]
					rownames(prop.par[[which(prop.names == paste('ETcoef',GR,sep='.'))]]) <- c('old','new')

			#delete temporary object
				rm(tmp)
							
			}
		}
		
		if('evap' %in% prop.names){
		  
		  tmp <- get.Evap(Input)
		  
		  for(j in 1:length(sites))
        Input <- model.set(Input = Input, set.names = "evap", set.par = tmp$evap.prop[[j]], j = j)	
		  
		  #save parameters into prop.par
		  prop.par[[which(prop.names == 'evap')]] <- tmp$pars
		    colnames(prop.par[[which(prop.names == 'evap')]]) <- 
		    var.names[[which(prop.names == 'evap')]]
		  rownames(prop.par[[which(prop.names == 'evap')]]) <- c('old','new')
		  
		  #delete temporary object
		  rm(tmp)
		  
		}
    
		list(Input = Input, prop.par = prop.par)
		}
		
	comp.prop <- function(Input){
	
		prop.par <- vector('list',length=length(prop.names))
		names(prop.par) <- prop.names

		if('composition' %in% prop.names){
			
			old <- new <- comp
			for(j in 1:length(sites)) {
				
				prop.keep <- which(comp.prior.mean>0)[which(comp.prior.mean>0)<5]
				
				new[prop.keep + (j-1)*5] <- 
					rtmvnorm(n = 1,
				    		mean = old[prop.keep],
				    		sigma = comp.Vpar[prop.keep + (j-1)*5,prop.keep + (j-1)*5],
				    		lower = rep(0.01,length(prop.keep)),
				    		upper = rep(0.99,length(prop.keep))) 
				
				if(sum(new[prop.keep + (j-1)*5])>1) new[prop.keep + (j-1)*5] <- 0.99*new[prop.keep + (j-1)*5]/sum(new[prop.keep + (j-1)*5]) #ensures vegetation cannot sum to more than 0.99
				    	
				new[j*5] <- 1 - sum(new[prop.keep + (j-1)*5])
								
				}

			for(j in 1:length(sites)){
				
				tmp <- new[grep(sites[j],names(old))]/sum(new[grep(sites[j],names(old))])
				names(tmp) <- c('Grasses','Shrubs','Trees','Forbs','Bare Ground')
				
        for(j in 1:length(sites)) 
          Input <- model.set(Input = Input, set.names = "composition", set.par = tmp, j = j)	
			
			}
		
			prop.par[[which(prop.names == 'composition')]] <- rbind(old,new)
				colnames(prop.par[[which(prop.names == 'composition')]]) <- 
					var.names[[which(prop.names == 'composition')]]
				rownames(prop.par[[which(prop.names == 'composition')]]) <- c('old','new')
				
			rm(list = c('old','new'))
		
		}

			

		list(Input = Input, prop.par = prop.par)
		}
		
  soils.prop <- function(Input){
  
    prop.par <- vector('list',length=length(prop.names))
    names(prop.par) <- prop.names
  
    if('soils' %in% prop.names){
    
      old <- new <- soils
      
      for(j in 1:length(sites)) {
        new.tmp <- new[grep(sites[j],names(new))]
        old.tmp <- old[grep(sites[j],names(old))]
          
        for(dj in 1:(length(old.tmp)/2)){
          new[(length(old.tmp)/2)*c(0,1)+dj] <- 
            tnorm.mvt(avec  = old.tmp[(length(old.tmp)/2)*c(0,1)+dj],
                      muvec = old.tmp[(length(old.tmp)/2)*c(0,1)+dj],
                       smat = soils.Vpar[(length(old.tmp)/2)*c(0,1)+dj,(length(old.tmp)/2)*c(0,1)+dj],
                         lo = soils.lo[(length(old.tmp)/2)*c(0,1)+dj],
                         hi = soils.hi[(length(old.tmp)/2)*c(0,1)+dj]) 
      
          if(sum(new[(length(old.tmp)/2)*c(0,1)+dj]>1)) 
            new[(length(old.tmp)/2)*c(0,1)+dj] <- 0.99*new[(length(old.tmp)/2)*c(0,1)+dj] /
                                                   sum(new[(length(old.tmp)/2)*c(0,1)+dj]) #ensures vegetation cannot sum to more than 0.99
        }
    }
    
    for(j in 1:length(sites)){
      new.tmp <- new[grep(sites[j],names(new))]
      tmp <- swSoils_Layers(Input[[j]])[,c('sand_frac','clay_frac')]
      tmp[,'sand_frac'] <- new.tmp[grep("sand",names(new))][soils.ind] #soils.ind needs to be fixed for multiple sites
      tmp[,'clay_frac'] <- new.tmp[grep("clay",names(new))][soils.ind]
           
      for(j in 1:length(sites)) 
        Input <- model.set(Input = Input, set.names = "soils", set.par = tmp, j = j)  
      
    }
    
    ##need to get rid of the list format so everything is vectors
    prop.par[[which(prop.names == 'soils')]] <- rbind(old,new)
    colnames(prop.par[[which(prop.names == 'soils')]]) <- 
      var.names[[which(prop.names == 'soils')]]
    rownames(prop.par[[which(prop.names == 'soils')]]) <- c('old','new')
    
    rm(list = c('old','new'))
    
  }
  
  
  
  list(Input = Input, prop.par = prop.par)
}

#function to initialize runs based on inputs. Should be similar to some of the param functions

model.set <- function(Input,set.names,set.par,j){ #Input is a list of swIn objects, 
                                                       #set.names is vector of names,
                                                       #set.par are the values, and
							                                         #j is the site

		if('drain' %in% set.names)
		  swSite_DrainageCoefficient(Input[[j]])  <- set.par
		
		for(GR in groups){
      if(paste('crit',GR,sep='.') %in% set.names){
			
			if(GR == 'tree')  swProd_CritSoilWaterPotential(Input[[j]])['Trees']   <- set.par
			if(GR == 'shrub') swProd_CritSoilWaterPotential(Input[[j]])['Shrubs']  <- set.par
			if(GR == 'grass') swProd_CritSoilWaterPotential(Input[[j]])['Grasses'] <- set.par
			if(GR == 'forb')  swProd_CritSoilWaterPotential(Input[[j]])['Forbs']   <- set.par

		  }
			
	  	if(paste('phenology',GR,sep='.') %in% set.names){
				
  		if(GR=='tree')  swProd_MonProd_tree(Input[[j]])  <- set.par
						
      if(GR=='shrub') swProd_MonProd_shrub(Input[[j]]) <-  set.par
					
			if(GR=='grass') swProd_MonProd_grass(Input[[j]]) <-  set.par
					
			if(GR=='forb')  swProd_MonProd_forb(Input[[j]])  <-  set.par
			}				
			

		  	if(paste('ETcoef',GR,sep='.') %in% set.names){ ##needs to be based on some sort of exponential decay
									 ## perhaps we predict the root distribution as a parametric function
		
				if(GR == 'tree')  swSoils_Layers(Input[[j]])[,'transpTree_frac']  <- set.par

				if(GR == 'shrub') swSoils_Layers(Input[[j]])[,'transpShrub_frac']  <- set.par

				if(GR == 'grass') swSoils_Layers(Input[[j]])[,'transpGrass_frac'] <- set.par

				if(GR == 'forb')  swSoils_Layers(Input[[j]])[,'transpForb_frac']  <- set.par

					
  			}
		}
		
		if('evap' %in% set.names)
      swSoils_Layers(Input[[j]])[,'EvapBareSoil_frac']  <- set.par
		  
		if("composition" %in% set.names)
		  swProd_Composition(Input[[j]]) <- set.par
    
		if('soils' %in% set.names)
		  swSoils_Layers(Input[[j]])[,c('sand_frac','clay_frac')] <- set.par
		
		

		return(Input)
		}
		
model.init <- function(set.names,Input, gg){
  set.par <- vector('list',length=length(set.names))
  names(set.par) <- set.names
  
  if("drain" %in% set.names)
    set.par$drain <- dgibbs[gg]
  
  if("crit.shrub" %in% set.names){
    set.par$crit.shrub <- cSgibbs[gg]
      names(set.par$crit.shrub) <- "Shrubs"
  }
  
  if("crit.grass" %in% set.names){
    set.par$crit.grass <- cGgibbs[gg]
      names(set.par$crit.grass) <- "Grasses"
  }

  if("phenology.shrub" %in% set.names){
    phen <- pSgibbs[gg,]
    phen.mat <- swProd_MonProd_shrub(Input[[j]])
  
    for(phen.col in 1:4){
      phen.mat[,phen.col] <- phen[c(2,3,4,1)[phen.col]] * 
        phen.mat[,phen.col] / max(phen.mat[,phen.col])
      }
    set.par$phenology.shrub <- phen.mat
  }

  if("phenology.grass" %in% set.names){
    phen <- pGgibbs[gg,]
    phen.mat <- swProd_MonProd_grass(Input[[j]])
  
    for(phen.col in 1:4){
      phen.mat[,phen.col] <- phen[c(2,3,4,1)[phen.col]] * 
        phen.mat[,phen.col] / max(phen.mat[,phen.col])
      }
    set.par$phenology.grass <- phen.mat
  }

  if("ETcoef.shrub" %in% set.names){
    d  <- swSoils_Layers(Input[[j]])[,1]
    ET <- T_depth.Schenk2002(d,bSgibbs[gg,])$r.d
    ET <- ET / sum(ET)
    set.par$ETcoef.shrub <- ET
  }
  
  if("ETcoef.grass" %in% set.names){
    d  <- swSoils_Layers(Input[[j]])[,1]
    ET <- T_depth.Schenk2002(d,bGgibbs[gg,])$r.d
    ET <- ET / sum(ET)
    set.par$ETcoef.grass <- ET
  }
  
  if("evap" %in% set.names){
    d     <- swSoils_Layers(Input[[j]])[,"depth_cm"]
    sand  <- swSoils_Layers(Input[[j]])[,"sand_frac"]
    clay  <- swSoils_Layers(Input[[j]])[,"clay_frac"]
    set.par$evap  <- E_depth_Wythers1999(d,sand,clay,agibbs[gg,])
  }
  
  if("soils" %in% set.names){
    tmp <- swSoils_Layers(Input[[j]])[,c('sand_frac','clay_frac')]
    tmp[,'sand_frac'] <- sgibbs[gg,grep("sand",colnames(sgibbs))][soils.ind] #soils.ind needs to be fixed for multiple sites
    tmp[,'clay_frac'] <- sgibbs[gg,grep("clay",colnames(sgibbs))][soils.ind]
    set.par$soils <- tmp
    rm(tmp)
  }
                              
  if("composition" %in% set.names){
    set.par$composition <- cgibbs[gg,]
    names(set.par$composition) <- c('Grasses','Shrubs','Trees','Forbs','Bare Ground')
  }
  
  return(set.par)
}

#new likelihood function for independent residuals

calc.dens.params.Ind <- function(j,i){
  
  #get timesteps for which there are data after the first year	
  keep <- which(is.finite(obs[[j]][,i]) & 
                !is.na(pred.vwc.old[[j]][,i]) &
                !is.na(pred.vwc.new[[j]][,i]))
  if(gibbs.thin) keep <- keep[keep %in% obs.thin]
  keep <- keep[keep>365]
  
  #initialize variables to hold summed log densities
  p.old <- p.new <- 0
  
  #continue only if there is at least one observation after the first year
  if(length(keep)>0) {
    
    
    #calculate density given old SOILWAT parameters
    p.old <- sum(dnorm(x = obs[[j]][keep,i],#contribution of site j
                       mean = pred.vwc.old[[j]][keep,i],
                       sd = sqrt(sig),
                       log = TRUE))
    
    #calculate density given new SOILWAT parameters
    p.new <- sum(dnorm(x = obs[[j]][keep,i],#contribution of site j
                       mean = pred.vwc.new[[j]][keep,i],
                       sd = sqrt(sig),
                       log = TRUE))
  }
  
  #output
  list(p.old = p.old, p.new = p.new)
}





