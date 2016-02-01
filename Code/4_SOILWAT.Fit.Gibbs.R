
#alter to use model.set function

		cl <- makeCluster(4)
	 	registerDoParallel(cl)
		
#WORKFLOW  --  Embed within Gibbs Sampler
	

	bTgibbs <- matrix(NA,nrow=ng,ncol=length(beta.tree))
		colnames(bTgibbs) <- names(beta.tree)
	pTgibbs <- matrix(NA,nrow=ng,ncol=length(phen.tree))
		colnames(pTgibbs) <- names(phen.tree)
	cTgibbs <- matrix(NA,nrow=ng,ncol=length(crit.tree))
		colnames(cTgibbs) <- names(crit.tree)

	bSgibbs <- matrix(NA,nrow=ng,ncol=length(beta.shrub))
		colnames(bSgibbs) <- names(beta.shrub)
	pSgibbs <- matrix(NA,nrow=ng,ncol=length(phen.shrub))
		colnames(pSgibbs) <- names(phen.shrub)
	cSgibbs <- matrix(NA,nrow=ng,ncol=length(crit.shrub))
		colnames(cSgibbs) <- names(crit.shrub)

	bGgibbs <- matrix(NA,nrow=ng,ncol=length(beta.grass))
		colnames(bGgibbs) <- names(beta.grass)
	pGgibbs <- matrix(NA,nrow=ng,ncol=length(phen.grass))
		colnames(pGgibbs) <- names(phen.grass)
	cGgibbs <- matrix(NA,nrow=ng,ncol=length(crit.grass))
		colnames(cGgibbs) <- names(crit.grass)

	bFgibbs <- matrix(NA,nrow=ng,ncol=length(beta.forb))
		colnames(bFgibbs) <- names(beta.forb)
	pFgibbs <- matrix(NA,nrow=ng,ncol=length(phen.forb))
		colnames(pFgibbs) <- names(phen.forb)
	cFgibbs <- matrix(NA,nrow=ng,ncol=length(crit.forb))
		colnames(cFgibbs) <- names(crit.forb)

	dgibbs <- matrix(NA,nrow=ng,ncol=length(drain))
		colnames(cFgibbs) <- names(drain)
		
  cgibbs <- matrix(NA,nrow=ng,ncol=length(comp))
		colnames(cgibbs) <- names(comp)

	sgibbs <- matrix(NA,nrow=ng,ncol=length(soils))
		colnames(sgibbs) <- names(soils)
		
	agibbs <- matrix(NA,nrow=ng,ncol=length(alpha))
		colnames(agibbs) <- names(alpha)
		
	vgibbs <- matrix(NA,nrow=ng,ncol=2)
		colnames(vgibbs) <- c('sig','rho')
	
	
	acc <- racc <- cacc <- sacc <- track <- 0 #acceptance counters
	
	
	check <- seq(check.step,ng,by=check.step)
	
	adapt <- seq(adapt.step*2,ng,by=adapt.step)
		adapt <- adapt[adapt<=10000]

	for(g in 1:ng){

   # oi <- readline()
    #if(oi == "y") stop
    
		track <- track + 1

	#(1) Generate predictions of volumetris soil water content for current parameters

			
			tmp <- get.pred(swIn.old,file.in=FALSE)
			
			swOut.old <- tmp
						
			pred.vwc.old <- obs
				
				for(j in 1:length(sites)) {
				
					pred.vwc.old[[j]][,-(1:4)] <- NA
				
					for(i in 5:ncol(obs[[j]])){
						pred.vwc.old[[j]][,i] <- swOut.old[[j]]@VWCBULK@Day[,2+sens.Lyr[sens.Lyr[,'Label']==sites[j],i-1]]

					}
					
				}
				
			
	#(2) Propose new parameter values		

			tmp <- params.prop(swIn.old)
			
      params[-c(grep('composition',prop.names),grep('soils',prop.names))] <- 
        tmp$prop.par[-c(grep('composition',prop.names),grep('soils',prop.names))]
			swIn.new <- tmp$Input
					

	#(3) Generate predictions of volumetris soil water content for proposed parameters
		
			tmp <- get.pred(swIn.new,file.in=FALSE)
			
			swOut.new <- tmp
			
			pred.vwc.new <- obs
				
				for(j in 1:length(sites)) {
				
					pred.vwc.new[[j]][,-(1:4)] <- NA
				
					for(i in 5:ncol(obs[[j]])){
						pred.vwc.new[[j]][,i] <- swOut.new[[j]]@VWCBULK@Day[,2+sens.Lyr[sens.Lyr[,'Label']==sites[j],i-1]]

					}
				}

		
	#(4) Acceptance Step
	
		p.old <- p.new <- 0
			
			obslist <- cbind(1,obs.col[[1]])
			if(length(obs.col)>1) 
				for(j in 2:length(obs.col))
					obslist <- rbind(obslist,cbind(j,obs.col[[j]]))
			
			
			
		#get densities for each sensor based on old and new parameters
			out <- foreach(k=1:nrow(obslist)) %dopar% {
			  if(!est.rho) 	tmp <- calc.dens.params.Ind(obslist[k,1],obslist[k,2])
			  if(est.rho) 	tmp <- calc.dens.params(obslist[k,1],obslist[k,2])
        return(tmp)
			}
					
		#sum densities
			p.old <- p.new <- 0
			for(ij in 1:length(out)){
			
				p.old <- p.old + out[[ij]]$p.old
				p.new <- p.new + out[[ij]]$p.new
				
			
				}
		
		#add prior effect to densities
		####need to make this work for all
				
					
			if('phenology.tree' %in% prop.names)
				p.old <- p.old + 
					lndMvn(params$phenology.tree['old',], 
							phen.tree.prior.mean[c(4,rep(c(1,2,3),length(sites)))],
							SparseM::backsolve(chol(phen.tree.prior.V),
									  diag(length(params$phenology.tree['old',]))))			#phenology.tree

			if('phenology.shrub' %in% prop.names)
				p.old <- p.old + 
					lndMvn(params$phenology.shrub['old',], 
							phen.shrub.prior.mean[c(4,rep(c(1,2,3),length(sites)))],
							SparseM::backsolve(chol(phen.shrub.prior.V),
									  diag(length(params$phenology.shrub['old',]))))			#phenology.shrub

			if('phenology.grass' %in% prop.names)
				p.old <- p.old + 
					lndMvn(params$phenology.grass['old',], 
							phen.grass.prior.mean[c(4,rep(c(1,2,3),length(sites)))],
							SparseM::backsolve(chol(phen.grass.prior.V),
									  diag(length(params$phenology.grass['old',])))) 		#phenology.grass

			if('phenology.forb' %in% prop.names)
				p.old <- p.old + 
					lndMvn(params$phenology.forb['old',], 
							phen.forb.prior.mean[c(4,rep(c(1,2,3),length(sites)))],
							SparseM::backsolve(chol(phen.forb.prior.V),
									  diag(length(params$phenology.forb['old',]))))			#phenology.forb
	 				
	    if('evap' %in% prop.names)
	      p.old <- p.old + 
	        lndMvn(params$evap['old',], 
	            rep(alpha.prior.mean,length(sites)),
	            SparseM::backsolve(chol(alpha.prior.V),
	                            diag(length(params$evap['old',])))) 			#ETcoef.tree
	
      if('ETcoef.tree' %in% prop.names)
				p.old <- p.old + 
					lndMvn(params$ETcoef.tree['old',], 
							rep(beta.tree.prior.mean,length(sites)),
							SparseM::backsolve(chol(beta.tree.prior.V),
									  diag(length(params$ETcoef.tree['old',])))) 			#ETcoef.tree

			if('ETcoef.shrub' %in% prop.names)
				p.old <- p.old + 
					lndMvn(params$ETcoef.shrub['old',], 
							rep(beta.shrub.prior.mean,length(sites)),
							SparseM::backsolve(chol(beta.shrub.prior.V),
									  diag(length(params$ETcoef.shrub['old',])))) 			#ETcoef.shrub

			if('ETcoef.grass' %in% prop.names)
				p.old <- p.old + 
					lndMvn(params$ETcoef.grass['old',], 
							rep(beta.grass.prior.mean,length(sites)),
							SparseM::backsolve(chol(beta.grass.prior.V),
									  diag(length(params$ETcoef.grass['old',])))) 			#ETcoef.grass

			if('ETcoef.forb' %in% prop.names)
				p.old <- p.old + 
					lndMvn(params$ETcoef.forb['old',], 
							rep(beta.forb.prior.mean,length(sites)),
							SparseM::backsolve(chol(beta.forb.prior.V),
									  diag(length(params$ETcoef.forb['old',]))))			#ETcoef.forb

			if('crit.tree' %in% prop.names)
				p.old <- p.old + 
					lndMvn(params$crit.tree['old',], 
							rep(crit.tree.prior.mean,length(sites)),
							SparseM::backsolve(chol(beta.tree.prior.V),
									  diag(length(params$crit.tree['old',])))) 			#crit.tree

			if('crit.shrub' %in% prop.names)
				p.old <- p.old + 
					lndMvn(params$crit.shrub['old',], 
							rep(crit.shrub.prior.mean,length(sites)),
							SparseM::backsolve(chol(crit.shrub.prior.V),
									  diag(length(params$crit.shrub['old',])))) 			#crit.shrub

			if('crit.grass' %in% prop.names)
				p.old <- p.old + 
					lndMvn(params$crit.grass['old',], 
							rep(crit.grass.prior.mean,length(sites)),
							SparseM::backsolve(chol(crit.grass.prior.V),
									  diag(length(params$crit.grass['old',])))) 			#crit.grass

			if('crit.forb' %in% prop.names)
				p.old <- p.old + 
					lndMvn(params$crit.forb['old',], 
							rep(crit.forb.prior.mean,length(sites)),
							SparseM::backsolve(chol(crit.forb.prior.V),
									  diag(length(params$crit.forb['old',]))))			#crit.forb
					
	    if('drain' %in% prop.names)
	      p.old <- p.old + 
	        lndMvn(params$drain['old',], 
	            rep(drain.prior.mean,length(sites)),
	            SparseM::backsolve(chol(drain.prior.V),
	                   diag(length(params$drain['old',]))))			    #drain
	
			
			if('phenology.tree' %in% prop.names)
				p.new <- p.new + 			
					lndMvn(params$phenology.tree['new',], 
							phen.tree.prior.mean[c(4,rep(c(1,2,3),length(sites)))],
							SparseM::backsolve(chol(phen.tree.prior.V),
									  diag(length(params$phenology.tree['new',]))))			#phenology.tree

			if('phenology.shrub' %in% prop.names)
				p.new <- p.new + 			
					lndMvn(params$phenology.shrub['new',], 
							phen.shrub.prior.mean[c(4,rep(c(1,2,3),length(sites)))],
							SparseM::backsolve(chol(phen.shrub.prior.V),
									  diag(length(params$phenology.shrub['new',]))))			#phenology.shrub

			if('phenology.grass' %in% prop.names)
				p.new <- p.new + 			
					lndMvn(params$phenology.grass['new',], 
							phen.grass.prior.mean[c(4,rep(c(1,2,3),length(sites)))],
							SparseM::backsolve(chol(phen.grass.prior.V),
									  diag(length(params$phenology.grass['new',])))) 		#phenology.grass

			if('phenology.forb' %in% prop.names)
				p.new <- p.new + 			
					lndMvn(params$phenology.forb['new',], 
							phen.forb.prior.mean[c(4,rep(c(1,2,3),length(sites)))],
							SparseM::backsolve(chol(phen.forb.prior.V),
									  diag(length(params$phenology.forb['new',]))))			#phenology.forb

	    if('evap' %in% prop.names)
	      p.new <- p.new + 
	        lndMvn(params$evap['new',], 
	            rep(alpha.prior.mean,length(sites)),
	            SparseM::backsolve(chol(alpha.prior.V),
	                            diag(length(params$evap['new',])))) 			#ETcoef.tree
	
  
      if('ETcoef.tree' %in% prop.names)
				p.new <- p.new + 
					lndMvn(params$ETcoef.tree['new',], 
							rep(beta.tree.prior.mean,length(sites)),
							SparseM::backsolve(chol(beta.tree.prior.V),
									  diag(length(params$ETcoef.tree['new',])))) 			#ETcoef.tree

			if('ETcoef.shrub' %in% prop.names)
				p.new <- p.new + 
					lndMvn(params$ETcoef.shrub['new',], 
							rep(beta.shrub.prior.mean,length(sites)),
							SparseM::backsolve(chol(beta.shrub.prior.V),
									  diag(length(params$ETcoef.shrub['new',])))) 			#ETcoef.shrub

			if('ETcoef.grass' %in% prop.names)
				p.new <- p.new + 
					lndMvn(params$ETcoef.grass['new',], 
							rep(beta.grass.prior.mean,length(sites)),
							SparseM::backsolve(chol(beta.grass.prior.V),
									  diag(length(params$ETcoef.grass['new',])))) 			#ETcoef.grass

			if('ETcoef.forb' %in% prop.names)
				p.new <- p.new + 
					lndMvn(params$ETcoef.forb['new',], 
							rep(beta.forb.prior.mean,length(sites)),
							SparseM::backsolve(chol(beta.forb.prior.V),
									  diag(length(params$ETcoef.forb['new',]))))			#ETcoef.forb
			

			if('crit.tree' %in% prop.names)
				p.new <- p.new + 
					lndMvn(params$crit.tree['new',], 
							rep(crit.tree.prior.mean,length(sites)),
							SparseM::backsolve(chol(crit.tree.prior.V),
									  diag(length(params$crit.tree['new',])))) 			#crit.tree

			if('crit.shrub' %in% prop.names)
				p.new <- p.new + 
					lndMvn(params$crit.shrub['new',], 
							rep(crit.shrub.prior.mean,length(sites)),
							SparseM::backsolve(chol(crit.shrub.prior.V),
									  diag(length(params$crit.shrub['new',])))) 			#crit.shrub

			if('crit.grass' %in% prop.names)
				p.new <- p.new + 
					lndMvn(params$crit.grass['new',], 
							rep(crit.grass.prior.mean,length(sites)),
							SparseM::backsolve(chol(crit.grass.prior.V),
									  diag(length(params$crit.grass['new',])))) 			#crit.grass

			if('crit.forb' %in% prop.names)
				p.new <- p.new + 
					lndMvn(params$crit.forb['new',], 
							rep(crit.forb.prior.mean,length(sites)),
							SparseM::backsolve(chol(crit.forb.prior.V),
									  diag(length(params$crit.forb['new',]))))			#crit.forb
	
      if('drain' %in% prop.names)
        p.new <- p.new + 
	        lndMvn(params$drain['new',], 
	            rep(drain.prior.mean,length(sites)),
	            SparseM::backsolve(chol(drain.prior.V),
	                  diag(length(params$drain['new',]))))			#drain
	
		#metropolis step
		
			a <- exp(p.new - p.old)
			z <- runif(1,0,1)
			if(a > z){
				
				if('crit.tree' %in% prop.names) crit.tree  <- params$crit.tree['new',]
				if('ETcoef.tree' %in% prop.names) beta.tree  <- params$ETcoef.tree['new',]
				if('phenology.tree' %in% prop.names) phen.tree  <- params$phenology.tree['new',]
				if('crit.shrub' %in% prop.names) crit.shrub <- params$crit.shrub['new',]
				if('ETcoef.shrub' %in% prop.names) beta.shrub  <- params$ETcoef.shrub['new',]
				if('phenology.shrub' %in% prop.names) phen.shrub <- params$phenology.shrub['new',]
				if('crit.grass' %in% prop.names) crit.grass <- params$crit.grass['new',]
				if('ETcoef.grass' %in% prop.names) beta.grass  <- params$ETcoef.grass['new',]
				if('phenology.grass' %in% prop.names) phen.grass <- params$phenology.grass['new',]
				if('crit.forb' %in% prop.names) crit.forb  <- params$crit.forb['new',]
				if('ETcoef.forb' %in% prop.names) beta.forb  <- params$ETcoef.forb['new',]
				if('phenology.forb' %in% prop.names) phen.forb  <- params$phenology.forb['new',]
				if('evap' %in% prop.names) alpha  <- params$evap['new',]
				if('drain' %in% prop.names) drain  <- params$drain['new',]
				
				swIn.old  <- swIn.new
				swOut.old <- swOut.new
				
				acc <- acc+1
				
        }


if('composition' %in% prop.names){
	#(1) Generate predictions of volumetris soil water content for current parameters

			#tmp <- get.pred(swIn.old,file.in=FALSE)
			
			#swOut.old <- tmp
						
			pred.vwc.old <- obs
				
				for(j in 1:length(sites)) {
				
					pred.vwc.old[[j]][,-(1:4)] <- NA
				
					for(i in 5:ncol(obs[[j]])){
						pred.vwc.old[[j]][,i] <- swOut.old[[j]]@VWCBULK@Day[,2+sens.Lyr[sens.Lyr[,'Label']==sites[j],i-1]]

					}
				}
				
			
	#(2) Propose new parameter values		

			tmp <- comp.prop(swIn.old)
			
			params$composition   <- tmp$prop.par$composition
			swIn.new <- tmp$Input
					

	#(3) Generate predictions of volumetris soil water content for proposed parameters
		
			tmp <- get.pred(swIn.new,file.in=FALSE)
			
			swOut.new <- tmp
			
			pred.vwc.new <- obs
				
				for(j in 1:length(sites)) {
				
					pred.vwc.new[[j]][,-(1:4)] <- NA
				
					for(i in 5:ncol(obs[[j]])){
						pred.vwc.new[[j]][,i] <- swOut.new[[j]]@VWCBULK@Day[,2+sens.Lyr[sens.Lyr[,'Label']==sites[j],i-1]]

            
					}
				}

		
	#(4) Acceptance Step
	
		p.old <- p.new <- 0
			
			obslist <- cbind(1,obs.col[[1]])
			if(length(obs.col)>1) 
				for(j in 2:length(obs.col))
					obslist <- rbind(obslist,cbind(j,obs.col[[j]]))
			
			
			
		#get densities for each sensor based on old and new parameters
			out <- foreach(k=1:nrow(obslist)) %dopar% {
			    if(!est.rho) 	tmp <- calc.dens.params.Ind(obslist[k,1],obslist[k,2])
			    if(est.rho) 	tmp <- calc.dens.params(obslist[k,1],obslist[k,2])
          return(tmp)
			}
					
		#sum densities
			p.old <- p.new <- 0
			for(ij in 1:length(out)){
			
				p.old <- p.old + out[[ij]]$p.old
				p.new <- p.new + out[[ij]]$p.new
				
			
				}
		
		#add prior effect to densities
		####need to make this work for all
		
			pcomp.old <- pcomp.new <- 0

			for(j in 1:length(sites)){
			
				comp.keep <- which(comp!=0)
			
				pcomp.old <- pcomp.old + log(MCMCpack::ddirichlet(params$composition['old',1:5 + (j-1)*5][comp.keep],comp.prior.mean[comp.keep]))
				pcomp.new <- pcomp.new + log(MCMCpack::ddirichlet(params$composition['new',1:5 + (j-1)*5][comp.keep],comp.prior.mean[comp.keep]))
			
				}
		
			p.old <- p.old + pcomp.old
					
			p.new <- p.new + pcomp.new
			
				
		#metropolis step
		
			a <- exp(p.new - p.old)
			z <- runif(1,0,1)
			if(a > z){
				
				 comp  <- params$composition['new',]
								
				swIn.old  <- swIn.new
				swOut.old <- swOut.new
				
				cacc <- cacc+1
				
				}

		
		}


if('soils' %in% prop.names){
  #(1) Generate predictions of volumetris soil water content for current parameters
  
  #tmp <- get.pred(swIn.old,file.in=FALSE)
  
  #swOut.old <- tmp
  
  pred.vwc.old <- obs
  
  for(j in 1:length(sites)) {
    
    pred.vwc.old[[j]][,-(1:4)] <- NA
    
    for(i in 5:ncol(obs[[j]])){
      pred.vwc.old[[j]][,i] <- swOut.old[[j]]@VWCBULK@Day[,2+sens.Lyr[sens.Lyr[,'Label']==sites[j],i-1]]

    }
  }
  
  
  #(2) Propose new parameter values		
  
  tmp <- soils.prop(swIn.old)
  
  params$soils   <- tmp$prop.par$soils
  swIn.new <- tmp$Input
  
  
  #(3) Generate predictions of volumetris soil water content for proposed parameters
  
  tmp <- get.pred(swIn.new,file.in=FALSE)
  
  swOut.new <- tmp
  
  pred.vwc.new <- obs
  
  for(j in 1:length(sites)) {
    
    pred.vwc.new[[j]][,-(1:4)] <- NA
    
    for(i in 5:ncol(obs[[j]])){
      pred.vwc.new[[j]][,i] <- swOut.new[[j]]@VWCBULK@Day[,2+sens.Lyr[sens.Lyr[,'Label']==sites[j],i-1]]

      
    }
  }
  
  
  #(4) Acceptance Step
  
  p.old <- p.new <- 0
  
  obslist <- cbind(1,obs.col[[1]])
  if(length(obs.col)>1) 
    for(j in 2:length(obs.col))
      obslist <- rbind(obslist,cbind(j,obs.col[[j]]))
  
  
  
  #get densities for each sensor based on old and new parameters
  out <- foreach(k=1:nrow(obslist)) %dopar% {
    if(!est.rho) tmp <- calc.dens.params.Ind(obslist[k,1],obslist[k,2])
    if(est.rho) tmp <- calc.dens.params(obslist[k,1],obslist[k,2])
    return(tmp)
  }
  
  #sum densities
  p.old <- p.new <- 0
  for(ij in 1:length(out)){
    
    p.old <- p.old + out[[ij]]$p.old
    p.new <- p.new + out[[ij]]$p.new
    
    
  }
  
  #add prior effect to densities
  ####need to make this work for all
  
  psoils.old <- psoils.new <- 0
  
  for(j in 1:length(sites)){
    
    soils.keep <- grep(sites[j],names(soils))
    
    for(dj in 1:(length(soils.keep)/2)){
      dkeep <- soils.keep[length(soils.keep)/2 * c(0,1) + dj]
      psoils.old <- psoils.old + log(MCMCpack::ddirichlet(c(params$soils['old',dkeep],1-sum(params$soils['old',dkeep])),
                                                          c(soils.prior.mean,1-sum(soils.prior.mean))))
      psoils.new <- psoils.new + log(MCMCpack::ddirichlet(c(params$soils['new',dkeep],1-sum(params$soils['new',dkeep])),
                                                          c(soils.prior.mean,1-sum(soils.prior.mean))))
      
    }
  }
  
  p.old <- p.old + psoils.old
  
  p.new <- p.new + psoils.new
  
  #add observations and a variance
  
  p.old <- p.old + sum(dmvnorm(x = soils.obs,#contribution of site j
                            mean = params$soils['old',],
                           sigma = diag(.01,length(soils)),
                                  #lower = rep(0,length(soils)),
                                  #upper = rep(1,length(soils)),
                             log = TRUE),
                       na.rm=TRUE)

  p.new <- p.new + sum(dmvnorm(x = soils.obs,#contribution of site j
                               mean = params$soils['new',],
                               sigma = diag(.01,length(soils)),
                               #lower = rep(0,length(soils)),
                               #upper = rep(1,length(soils)),
                               log = TRUE),
                       na.rm=TRUE)
  
  #metropolis step
  
  a <- exp(p.new - p.old)
  z <- runif(1,0,1)
  if(a > z){
    
    soils  <- params$soils['new',]
    
    swIn.old  <- swIn.new
    swOut.old <- swOut.new
    
    sacc <- sacc+1
    
  }
  
  
}

  #update variance
			
				#get new predictions
       pred.vwc <- obs

				for(j in 1:length(sites)) {
				
					pred.vwc[[j]][,-(1:4)] <- NA
				
					for(i in 5:ncol(obs[[j]])){
						pred.vwc[[j]][,i] <- swOut.old[[j]]@VWCBULK@Day[,2+sens.Lyr[sens.Lyr[,'Label']==sites[j],i-1]]
					}
				
				}

#ensure that no NA values enter into the calculations
  keep.sig <- which(is.finite(obs[[j]][,i]) & 
                  !is.na(pred.vwc[[j]][,i]))
  if(gibbs.thin) keep.sig <- keep.sig[keep.sig %in% obs.thin]
  keep.sig <- keep.sig[keep.sig>365]


if(!est.rho){
  
  T <- 0
  res <- 0
  
  for(ij in 1:nrow(obslist)) {
    tmp <- (as.numeric(obs[[obslist[ij,1]]][cbind(keep.sig,obslist[ij,2])]) - 
              as.numeric(pred.vwc[[obslist[ij,1]]][,obslist[ij,2]]))^2
    
    tmp <- tmp[is.finite(tmp)]
    
    
    T <- T + length(tmp)
    res <- res + sum(tmp)
    
  }
  
  sig <- 1/rgamma(1,s1 + .5*T,s2 + .5*(res))
  
}

#get intermediate forms
      if(est.rho){
      T <- rep(0,nrow(obslist))
			res <- invR <- list()
				for(ij in 1:nrow(obslist)) {
					res[[ij]] <-      obs[[obslist[ij,1]]][,obslist[ij,2]] - 
						                pred.vwc[[obslist[ij,1]]][,obslist[ij,2]]
					
						         
					invR[[ij]]   <- inv.rmat(length(res[[ij]]),rho)
					
					T[ij] <- length(which(is.finite(res[[ij]])))
			
					invR[[ij]] <- invR[[ij]][is.finite(res[[ij]]),is.finite(res[[ij]])]
					res[[ij]] <- res[[ij]][is.finite(res[[ij]])]
					
				}
			
			sumT <- sum(T)
			sumR <- 0
				for(ij in 1:length(res)){
					if(length(res[[ij]])==0) next 
					
					#get the inverse of the matrix R[[ij]] using the SparseM package
					
					sumR <- sumR + t(res[[ij]]) %*% invR[[ij]]  %*% res[[ij]]
					}
			
			sig <- 1/rgamma(1,s1 + .5*sumT,
							 s2 + .5*sum(sumR))
			
			#density for rho

	
			rho.prop <- as.vector(rtmvnorm(n = 1,
									 mean = rho,
									 sigma = rhoV,
									 lower = rho.lo,
									 upper = rho.hi))

			
			out <- foreach(k=1:nrow(obslist)) %dopar% {
						calc.dens.rho(obslist[k,1],obslist[k,2])
					}
				
			#sum densities
			p.old <- p.new <- 0
			for(ij in 1:length(out)){
			
				p.old <- p.old + out[[ij]]$p.old
				p.new <- p.new + out[[ij]]$p.new
				
			
				}
				
			p.old <- p.old + dtmvnorm.marginal(xn = rho,
									  mean = rho.prior.mean,
									  sigma = rho.prior.var,
									  lower = rho.lo,
									  upper = rho.hi,log=TRUE)
					
			p.new <- p.new + dtmvnorm.marginal(xn = rho.prop,
									  mean = rho.prior.mean,
									  sigma = rho.prior.var,
									  lower = rho.lo,
									  upper = rho.hi,log=TRUE)
					
		#metropolis step
		
			a <- exp(p.new - p.old)
			z <- runif(1,0,1)
			if(a > z){
				
				rho <- rho.prop
				
				racc <- racc+1
				
				}
				}
		
		print('===========================')
		print(paste('step = ',g))
		#print(beta.tree)
		#print(phen.tree)
		print(paste('sigma =',round(sig,4)))
		#print(paste('rho =',round(rho,4)))
		print(paste('param acceptances = ',round(acc/track,4)))
		print(paste('comp acceptances = ',round(cacc/track,4)))
    print(paste('soils acceptances = ',round(sacc/track,4)))
    print(paste('rho acceptances = ',round(racc/track,4)))
		print('===========================')

    if('evap' %in% prop.names)              agibbs[g,] <- alpha

    if('crit.tree' %in% prop.names)        cTgibbs[g] <- crit.tree
		if('ETcoef.tree' %in% prop.names)     bTgibbs[g,] <- beta.tree
		if('phenology.tree' %in% prop.names)  pTgibbs[g,] <- phen.tree

		if('crit.shrub' %in% prop.names)       cSgibbs[g] <- crit.shrub
		if('ETcoef.shrub' %in% prop.names)    bSgibbs[g,] <- beta.shrub
		if('phenology.shrub' %in% prop.names) pSgibbs[g,] <- phen.shrub

		if('crit.grass' %in% prop.names)       cGgibbs[g] <- crit.grass
		if('ETcoef.grass' %in% prop.names)    bGgibbs[g,] <- beta.grass
		if('phenology.grass' %in% prop.names) pGgibbs[g,] <- phen.grass

		if('crit.forb' %in% prop.names)        cFgibbs[g] <- crit.forb
		if('ETcoef.forb' %in% prop.names)     bFgibbs[g,] <- beta.forb
		if('phenology.forb' %in% prop.names)  pFgibbs[g,] <- phen.forb
		if('composition' %in% prop.names)     cgibbs[g,]  <- comp
		
    if('drain' %in% prop.names)           dgibbs[g,]  <- drain

    if('soils' %in% prop.names)           sgibbs[g,]  <- soils

		vgibbs[g,] <- c(sig,rho)
		
		
		if(g %in% adapt){
		
				if(acc/track > .4){															
					#get a proportional scalar to modify the proposal covariance matrix
					# based on Haario et al. 2001. An adaptive Metropolis algorithm. Bernoulli 7: 223-242.
					#          Gelamn et al. 1996. Efficient Metropolis jumping rules. In J. M. Bernardo, J. O. Berger, A. F. David, and A. F. M. Smith (eds.). Bayesian Statistics V, pp. 599-608. Oxford: Oxford University Press.
					prop.scale <- 2.38^2 / (3 + 2 * (1 + length(sites) * (3 + 2 + 1)))
					
					if('phenology.tree' %in% prop.names) phen.tree.var   <- prop.scale*var(pTgibbs[(g-adapt.step+1):(g-1),])*diag(ncol(pTgibbs))
					if('phenology.shrub' %in% prop.names) phen.shrub.var <- prop.scale*var(pSgibbs[(g-adapt.step+1):(g-1),])*diag(ncol(pSgibbs))
					if('phenology.grass' %in% prop.names) phen.grass.var <- prop.scale*var(pGgibbs[(g-adapt.step+1):(g-1),])*diag(ncol(pGgibbs))
					if('phenology.forb' %in% prop.names) phen.forb.var   <- prop.scale*var(pFgibbs[(g-adapt.step+1):(g-1),])*diag(ncol(pFgibbs))

					if('evap' %in% prop.names) alpha.var   <- prop.scale*var(agibbs[(g-adapt.step+1):(g-1),])*diag(ncol(agibbs))
					
          			  	if('ETcoef.tree' %in% prop.names) beta.tree.var   <- prop.scale*var(bTgibbs[(g-adapt.step+1):(g-1),])*diag(ncol(bTgibbs))
					if('ETcoef.shrub' %in% prop.names) beta.shrub.var <- prop.scale*var(bSgibbs[(g-adapt.step+1):(g-1),])*diag(ncol(bSgibbs))
					if('ETcoef.grass' %in% prop.names) beta.grass.var <- prop.scale*var(bGgibbs[(g-adapt.step+1):(g-1),])*diag(ncol(bGgibbs))
					if('ETcoef.forb' %in% prop.names) beta.forb.var   <- prop.scale*var(bFgibbs[(g-adapt.step+1):(g-1),])*diag(ncol(bFgibbs))

					if('crit.tree' %in% prop.names) crit.tree.var   <- prop.scale*var(cTgibbs[(g-adapt.step+1):(g-1),])*diag(ncol(cTgibbs))
					if('crit.shrub' %in% prop.names) crit.shrub.var <- prop.scale*var(cSgibbs[(g-adapt.step+1):(g-1),])*diag(ncol(cSgibbs))
					if('crit.grass' %in% prop.names) crit.grass.var <- prop.scale*var(cGgibbs[(g-adapt.step+1):(g-1),])*diag(ncol(cGgibbs))
					if('crit.forb' %in% prop.names) crit.forb.var   <- prop.scale*var(cFgibbs[(g-adapt.step+1):(g-1),])*diag(ncol(cFgibbs))
					
          				if('drain' %in% prop.names) drain.var   <- prop.scale*var(dgibbs[(g-adapt.step+1):(g-1),])*diag(ncol(dgibbs))

				}
        
				if(acc/track  < .20){
          
				  if('phenology.tree' %in% prop.names) phen.tree.var   <- phen.tree.var * .5
				  if('phenology.shrub' %in% prop.names) phen.shrub.var <- phen.shrub.var * .5
				  if('phenology.grass' %in% prop.names) phen.grass.var <- phen.grass.var * .5
				  if('phenology.forb' %in% prop.names) phen.forb.var   <- phen.forb.var * .5
          
				  if('evap' %in% prop.names) alpha.var   <- alpha.var * .5
				  
				  if('ETcoef.tree' %in% prop.names) beta.tree.var   <-beta.tree.var * .5
				  if('ETcoef.shrub' %in% prop.names) beta.shrub.var <- beta.shrub.var * .5
				  if('ETcoef.grass' %in% prop.names) beta.grass.var <- beta.grass.var * .5
				  if('ETcoef.forb' %in% prop.names) beta.forb.var   <- beta.forb.var * .5
				  
				  if('crit.tree' %in% prop.names) crit.tree.var   <- crit.tree.var * .5
				  if('crit.shrub' %in% prop.names) crit.shrub.var <- crit.shrub.var * .5
				  if('crit.grass' %in% prop.names) crit.grass.var <- crit.grass.var * .5
				  if('crit.forb' %in% prop.names) crit.forb.var   <- crit.forb.var * .5
				  if('drain' %in% prop.names) drain.var   <- drain.var * .5
          
				}
        
					if('composition' %in% prop.names) {
					  if(cacc/track  < .2) comp.Vpar <- comp.Vpar * .5
            if(cacc/track  >.4) comp.Vpar <- var(cgibbs[(g-adapt.step+1):(g-1),])
					}	
				
				if('soils' %in% prop.names) {
				 				  
         			  prop.scale <- 2.38^2 / (length(soils))
				  if(sacc/track  < .2) soils.Vpar <- soils.Vpar * .5
				  if(sacc/track  >.4) soils.Vpar <- prop.scale*var(sgibbs[(g-adapt.step+1):(g-1),])*diag(ncol(sgibbs))
				}	
				
				rhoV <- var(vgibbs[(g-adapt.step+1):(g-1),2])
				
			acc <- 0
			racc <- 0
			cacc <- 0
			sacc <- 0
			track <- 0
			
		}
		
		if(g %in% check){
						
			}
			#plot gibbs chains as the loop progresses
			
			if(TRUE){
			
				par(mfrow=c(4,4),mar=c(4,4,2,.5))
				
          plot.keep <- max(1,g-1000):g
        
					plot(plot.keep,bTgibbs[plot.keep,1],type='l',main='Critical Threshold',
						 ylim=range(cbind(cTgibbs[plot.keep,(1:length(sites))],
						 				  cGgibbs[plot.keep,(1:length(sites))],
						 				  cSgibbs[plot.keep,(1:length(sites))],
						 				  cFgibbs[plot.keep,(1:length(sites))]),na.rm=TRUE))
						lines(plot.keep,cSgibbs[plot.keep,1],lty=2)
						lines(plot.keep,cGgibbs[plot.keep,1],lty=1,lwd=2)
						lines(plot.keep,cFgibbs[plot.keep,1],lty=2,lwd=2)
					
					if(length(sites)>1){
						for(j in 2:length(sites)) {
							lines(plot.keep,cTgibbs[plot.keep,j],col=j)
							lines(plot.keep,cSgibbs[plot.keep,j],col=j,lty=2)
							lines(plot.keep,cGgibbs[plot.keep,j],col=j,lty=1,lwd=2)
							lines(plot.keep,cFgibbs[plot.keep,j],col=j,lty=2,lwd=2)
							}
						}

					plot(plot.keep,bTgibbs[plot.keep,1],type='l',main='d50',
						 ylim=range(cbind(bTgibbs[plot.keep,(1:length(sites))*2-1],
						 				  bGgibbs[plot.keep,(1:length(sites))*2-1],
						 				  bSgibbs[plot.keep,(1:length(sites))*2-1],
						 				  bFgibbs[plot.keep,(1:length(sites))*2-1]),na.rm=TRUE))
						lines(plot.keep,bSgibbs[plot.keep,1],lty=2)
						lines(plot.keep,bGgibbs[plot.keep,1],lty=1,lwd=2)
						lines(plot.keep,bFgibbs[plot.keep,1],lty=2,lwd=2)
					
					if(length(sites)>1){
						for(j in 2:length(sites)) {
							lines(plot.keep,bTgibbs[plot.keep,(j - 1)*2 + 1],col=j)
							lines(plot.keep,bSgibbs[plot.keep,(j - 1)*2 + 1],col=j,lty=2)
							lines(plot.keep,bGgibbs[plot.keep,(j - 1)*2 + 1],col=j,lty=1,lwd=2)
							lines(plot.keep,bFgibbs[plot.keep,(j - 1)*2 + 1],col=j,lty=2,lwd=2)
							}
						}
						

					plot(plot.keep,bTgibbs[plot.keep,2],type='l',main='c',
						 ylim=range(cbind(bTgibbs[plot.keep,(1:length(sites))*2],
						 				  bGgibbs[plot.keep,(1:length(sites))*2],
						 				  bSgibbs[plot.keep,(1:length(sites))*2],
						 				  bFgibbs[plot.keep,(1:length(sites))*2]),na.rm=TRUE))
						lines(plot.keep,bSgibbs[plot.keep,2],lty=2)
						lines(plot.keep,bGgibbs[plot.keep,2],lty=1,lwd=2)
						lines(plot.keep,bFgibbs[plot.keep,2],lty=2,lwd=2)
					
					if(length(sites)>1){
						for(j in 2:length(sites)) {
							lines(bTgibbs[plot.keep,plot.keep,(j - 1)*2 + 2],col=j)
							lines(bSgibbs[plot.keep,plot.keep,(j - 1)*2 + 2],col=j,lty=2)
							lines(bGgibbs[plot.keep,plot.keep,(j - 1)*2 + 2],col=j,lty=1,lwd=2)
							lines(bFgibbs[plot.keep,plot.keep,(j - 1)*2 + 2],col=j,lty=2,lwd=2)
							}
						}
						
					legend('bottomleft',legend=sites,text.col=1:length(sites),cex=.8,bg='white')
					
				plot(plot.keep,agibbs[plot.keep,1],type='l',main='evap.max',
				     ylim=range(agibbs[plot.keep,(1:length(sites))*2-1],na.rm=TRUE))
				
				if(length(sites)>1){
				  for(j in 2:length(sites)) {
				    lines(plot.keep,agibbs[plot.keep,(j - 1)*2 + 1],col=j)
				  }
				}				
				plot(plot.keep,agibbs[plot.keep,2],type='l',main='evap.rate',
				     ylim=range(agibbs[plot.keep,(1:length(sites))*2],na.rm=TRUE))
				
				if(length(sites)>1){
				  for(j in 2:length(sites)) {
				    lines(plot.keep,agibbs[plot.keep,(j - 1)*2 + 2],col=j)
				  }
				}
				
        plot(plot.keep,pTgibbs[plot.keep,1],type='l',main='LAI_conv',log='y',
						 ylim=c(min(min(pTgibbs[plot.keep,1],na.rm=TRUE),
								   min(pGgibbs[plot.keep,1],na.rm=TRUE),
								   min(pSgibbs[plot.keep,1],na.rm=TRUE),
								   min(pFgibbs[plot.keep,1],na.rm=TRUE)),
							  max(max(pTgibbs[plot.keep,1],na.rm=TRUE),
								   max(pGgibbs[plot.keep,1],na.rm=TRUE),
								   max(pSgibbs[plot.keep,1],na.rm=TRUE),
								   max(pFgibbs[plot.keep,1],na.rm=TRUE))))
						lines(plot.keep,pSgibbs[plot.keep,1],lty=2)
						lines(plot.keep,pGgibbs[plot.keep,1],lty=1,lwd=2)
						lines(plot.keep,pFgibbs[plot.keep,1],lty=2,lwd=2)
					legend('topleft',legend=c('Grass','Shrub','Tree','Forb'),lty=c(1,2,1,2),lwd=c(2,1,1,2))
						

					plot(plot.keep,pTgibbs[plot.keep,2],type='l',log='y',
						ylim=c(min(min(pTgibbs[plot.keep,3*(1:length(sites)-1)+2],na.rm=TRUE),
								   min(pGgibbs[plot.keep,3*(1:length(sites)-1)+2],na.rm=TRUE),
								   min(pSgibbs[plot.keep,3*(1:length(sites)-1)+2],na.rm=TRUE),
								   min(pFgibbs[plot.keep,3*(1:length(sites)-1)+2],na.rm=TRUE)),
							  max(max(pTgibbs[plot.keep,3*(1:length(sites)-1)+2],na.rm=TRUE),
								   max(pGgibbs[plot.keep,3*(1:length(sites)-1)+2],na.rm=TRUE),
								   max(pSgibbs[plot.keep,3*(1:length(sites)-1)+2],na.rm=TRUE),
								   max(pFgibbs[plot.keep,3*(1:length(sites)-1)+2],na.rm=TRUE))),
						main='Litter')
						
						lines(plot.keep,pSgibbs[plot.keep,2],lty=2)
						lines(plot.keep,pGgibbs[plot.keep,2],lty=1,lwd=2)
						lines(plot.keep,pFgibbs[plot.keep,2],lty=2,lwd=2)
					
					if(length(sites)>1){
						for(j in 2:length(sites)) {
							lines(plot.keep,pTgibbs[plot.keep,3*(j-1)+2],col=j)
							lines(plot.keep,pSgibbs[plot.keep,3*(j-1)+2],col=j,lty=2)
							lines(plot.keep,pGgibbs[plot.keep,3*(j-1)+2],col=j,lty=1,lwd=2)
							lines(plot.keep,pFgibbs[plot.keep,3*(j-1)+2],col=j,lty=2,lwd=2)
							}
						}
						
					plot(plot.keep,pTgibbs[plot.keep,3],type='l',log='y',
						ylim=c(min(min(pTgibbs[plot.keep,3*(1:length(sites)-1)+3],na.rm=TRUE),
								   min(pGgibbs[plot.keep,3*(1:length(sites)-1)+3],na.rm=TRUE),
								   min(pSgibbs[plot.keep,3*(1:length(sites)-1)+3],na.rm=TRUE),
								   min(pFgibbs[plot.keep,3*(1:length(sites)-1)+3],na.rm=TRUE)),
							  max(max(pTgibbs[plot.keep,3*(1:length(sites)-1)+3],na.rm=TRUE),
								   max(pGgibbs[plot.keep,3*(1:length(sites)-1)+3],na.rm=TRUE),
								   max(pSgibbs[plot.keep,3*(1:length(sites)-1)+3],na.rm=TRUE),
								   max(pFgibbs[plot.keep,3*(1:length(sites)-1)+3],na.rm=TRUE))),
						main='Biomass')
						
						lines(plot.keep,pSgibbs[plot.keep,3],lty=2)
						lines(plot.keep,pGgibbs[plot.keep,3],lty=1,lwd=2)
						lines(plot.keep,pFgibbs[plot.keep,3],lty=2,lwd=2)
					
					if(length(sites)>1){
						for(j in 2:length(sites)) {
							lines(plot.keep,pTgibbs[plot.keep,3*(j-1)+3],col=j)
							lines(plot.keep,pSgibbs[plot.keep,3*(j-1)+3],col=j,lty=2)
							lines(plot.keep,pGgibbs[plot.keep,3*(j-1)+3],col=j,lty=1,lwd=2)
							lines(plot.keep,pFgibbs[plot.keep,3*(j-1)+3],col=j,lty=2,lwd=2)
							}
						}
						
					plot(plot.keep,pTgibbs[plot.keep,4],type='l',log='y',
						ylim=c(min(min(pTgibbs[plot.keep,3*(1:length(sites)-1)+4],na.rm=TRUE),
								   min(pGgibbs[plot.keep,3*(1:length(sites)-1)+4],na.rm=TRUE),
								   min(pSgibbs[plot.keep,3*(1:length(sites)-1)+4],na.rm=TRUE),
								   min(pFgibbs[plot.keep,3*(1:length(sites)-1)+4],na.rm=TRUE)),
							  max(max(pTgibbs[plot.keep,3*(1:length(sites)-1)+4],na.rm=TRUE),
								   max(pGgibbs[plot.keep,3*(1:length(sites)-1)+4],na.rm=TRUE),
								   max(pSgibbs[plot.keep,3*(1:length(sites)-1)+4],na.rm=TRUE),
								   max(pFgibbs[plot.keep,3*(1:length(sites)-1)+4],na.rm=TRUE))),
						main='Live_pct')
						
						lines(plot.keep,pSgibbs[plot.keep,4],lty=2)
						lines(plot.keep,pGgibbs[plot.keep,4],lty=1,lwd=2)
						lines(plot.keep,pFgibbs[plot.keep,4],lty=2,lwd=2)

					if(length(sites)>1){
						for(j in 2:length(sites)) {
							lines(plot.keep,pTgibbs[plot.keep,3*(j-1)+4],col=j)
							lines(plot.keep,pSgibbs[plot.keep,3*(j-1)+4],col=j,lty=2)
							lines(plot.keep,pGgibbs[plot.keep,3*(j-1)+4],col=j,lty=1,lwd=2)
							lines(plot.keep,pFgibbs[plot.keep,3*(j-1)+4],col=j,lty=2,lwd=2)
							}
						}
						
			if('composition' %in% prop.names){
				for(j in 1:length(sites)){
					plot(plot.keep,cgibbs[plot.keep,1 + (j-1) * 5],main='composition',type='l',ylim=c(0,1))
						for(k in 2:5) lines(plot.keep,cgibbs[plot.keep,k + (j-1) * 5],type='l',col=k)
				}
				legend('topleft',legend=c('Grass','Shrub','Tree','Forb','BareGround'),text.col=1:5)
			}

					plot(plot.keep,vgibbs[plot.keep,1],main='sigma^2',type='l',log="y")
					plot(plot.keep,vgibbs[plot.keep,2],main='rho',type='l')
          plot(plot.keep,dgibbs[plot.keep,1],main='drainage',type='l')
				}


		}
		
		
stopCluster(cl)

save.image(paste(dir.prj,paste(runname,vers,'RData',sep="."),sep=""))
