###############################
## date last modified: 2014.02.10
## author: David M. Bell (2014.02.01 - present)
## TO DO:	check soil texture files
##			add composition
##			add phenology and ETcoef for grasses and shrubs
##			add critical soil water content
##			add evaporation coefficients
###############################

#clear workspace
  
  try(stopCluster(cl))
  rm(list = ls())

#load R libraries
  library('mvtnorm')
  library('tmvtnorm')
  library('foreach')
  library('MCMCpack')
  library('bayesm')
  library('SparseM')
  library('doParallel')
  library(Rsoilwat31)
  packageVersion("Rsoilwat31")


#set name of run and version of priors to use
  runname <- 'Bayes_ReynoldsCreek' #This is the name of the folder from which all data will be extracted
	vers <- 'SG.20150210' #vers of run, to help distonguish between different runs
	
#Toggles  
  gibbs.thin <- FALSE # if true, thin the observations for the purposes of fitting by 1/4

#provide directory for code, data, etc.
  dir.model <- 'E:/Soilwat'
	dir.prj <- paste(dir.model,runname,sep='/')


	#set working directory
	setwd(dir.prj)

  #choose site(s)
	sites <- '2029_SCAN_ReynoldsCreek' #Should match the name of the field data folder for the site
                                     #located in the SW Input folder
	
	est.rho <- FALSE #toggle to determine whether or not autocorrelative variances are estimated
                  #this is currently pretty slow ... might be able to accelerate it in the future
	
	ng <- 20000       #number of Gibbs steps
	check.step <- 100  
	adapt.step <- 100  #number of steps between adaptations

  obs.thin <- seq(1,ng,by=1) #sub select steps of the Gibbs sampler to retain for the purposes of 
                             #reducing autocorrelation in the chain of estimates and reducing output 
                             #data size

  #set values for metropolis algorithm proposal variances

    #evaporation parameter for Wythers et al. 1999
    alpha.var <- c(.1,0.001)
      alpha.var <- diag(rep(alpha.var, length(sites)))

    #transpiration parameters for Schenk and Jackson 2002
    beta.tree.var <- c(.01,0.00001)
		beta.shrub.var <- c(.01,0.00001)
		beta.grass.var <- c(.01,0.00001)
		beta.forb.var <- c(.01,0.00001)
		
			beta.tree.var <- diag(rep(beta.tree.var, length(sites)))
			beta.shrub.var <- diag(rep(beta.shrub.var, length(sites)))
			beta.grass.var <- diag(rep(beta.grass.var, length(sites)))
			beta.forb.var <- diag(rep(beta.forb.var, length(sites)))

    #critical soil water potential parameter
		crit.tree.var  <- matrix(.00001,1,1)
		crit.shrub.var <- matrix(.00001,1,1)
		crit.grass.var <- matrix(.00001,1,1)
		crit.forb.var  <- matrix(.00001,1,1)
		
		if(length(sites)>1){
			crit.tree.var <- diag(rep(crit.tree.var, length(sites)))
			crit.shrub.var <- diag(rep(crit.shrub.var, length(sites)))
			crit.grass.var <- diag(rep(crit.grass.var, length(sites)))
			crit.forb.var <- diag(rep(crit.forb.var, length(sites)))
		}

    #deep drainage parameter
    drain.var  <- matrix(.0000001,1,1)

    if(length(sites)>1){
      drain.var <- diag(rep(drain.var, length(sites)))
    }

    #phenological parameters
		phen.tree.var  <- diag(c(1,rep(c(1,1,0.00001),times=length(sites))))
		phen.shrub.var <- diag(c(1,rep(c(1,1,0.00001),times=length(sites))))
		phen.grass.var <- diag(c(1,rep(c(1,1,0.00001),times=length(sites))))
		phen.forb.var  <- diag(c(1,rep(c(1,1,0.00001),times=length(sites))))

		#composition parameter
    comp.Vpar <- diag(rep(.00001,5*length(sites)))
    
    #soil parameter
    soils.Vpar <- 0.00001

  #import functions
	source(paste(dir.model,'2_SOILWAT.Fit.Functions.r',sep='/'))

  #initialize model fitting
	source(paste(dir.model,'3_SOILWAT.Fit.Initialize.R',sep='/'))
	

	#initialize swIN, the initial Soilwat run
		tmp <-  params.prop(swIn)
			swIn.old <- tmp$Input
			params <- tmp$prop.par

    tmp <-  comp.prop(swIn.old)
      swIn.old  <- tmp$Input
      params$composition <- tmp$prop.par$composition

    tmp <-  soils.prop(swIn.old)
      swIn.old <- tmp$Input
      params$soils <- tmp$prop.par$soils

  #take a look at the precip data  
    prec <- swIn[[1]]@weatherHistory[[1]]@data[,4]
      for(j in 2:length(swIn[[1]]@weatherHistory))
	      prec <- c(prec,swIn[[1]]@weatherHistory[[j]]@data[,4])

    plot(prec,type='l',col='blue',pch=20)

  #run the model

  	source(paste(dir.model,'4_SOILWAT.Fit.Gibbs.R',sep='/'))
	
#plotting code for exploring variance estimation issues
if(FALSE){
  
  g1 <- 92
  g2 <- 80
  par(mfrow=c(1,1))
  plot(obs[[1]][,5],type="l",ylim=c(0,1))
  lines(pred.vwc.save[[g1]][[1]][,5],col="red")
  lines(pred.vwc.save[[g2]][[1]][,5],col="blue")
  
}
	