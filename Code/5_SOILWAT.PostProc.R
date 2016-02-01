#SOILWAT post processing, including saving output as .csv files, plotting, and 
#statistic calculation

	load('E:/SOILWAT/Bayes_ReynoldsCreek/Bayes_ReynoldsCreekBayes_ReynoldsCreek.SG.20150126.RData')
#
	ng <- g-1
	burnin <- 5000
	post.keep <- burnin:ng

  library('mvtnorm')
  library('tmvtnorm')
  library('foreach')
  library('MCMCpack')
  library('bayesm')
  library('SparseM')
  library('doParallel')
  library(Rsoilwat31)
  packageVersion("Rsoilwat31")

#set up output directory

	dir.run <- paste(dir.prj,'3_Runs',sep='/')

	setwd(dir.run)

	if(!('Bayes_Output' %in% list.files())) dir.create('Bayes_Output')

	dir.out <- paste(dir.run,'Bayes_Output',sep='/')
	
	setwd(dir.out)

#save workspace

	if(!('data' %in% list.files())) dir.create('data')
	
	dir.out.data <- paste(dir.out,'data',sep='/')
	
	setwd(dir.out.data)

		save.image(paste(runname,vers,'RData',sep='.'))	

	setwd(dir.out)
	
	rm(list = ls(pattern = 'gibbs')[-which(ls(pattern = 'gibbs') %in% c('gibbs.names',gibbs.names))])

#calculate and save statistics

	if(!('tables' %in% list.files())) dir.create('tables')
	
	dir.out.tables <- paste(dir.out,'tables',sep='/')
	
	setwd(dir.out.tables)

		#list of Gibbs matrices and prop.names
			rm(list = c('gibbs.list','gibbs.prop','gibbs.names'))
			
			gibbs.list <- gibbs.prop <- ls(pattern='gibbs')
				if("agibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "agibbs"] <- "evap"

				if("cgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "cgibbs"] <- "composition"
				if("dgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "dgibbs"] <- "drain"
				if("sgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "sgibbs"] <- "soils"
				if("vgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "vgibbs"] <- "variance"

				if("bFgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "bFgibbs"] <- "ETcoef.forb"
				if("bGgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "bGgibbs"] <- "ETcoef.grass"
				if("bSgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "bSgibbs"] <- "ETcoef.shrub"
				if("bTgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "bTgibbs"] <- "ETcoef.tree"

				if("cFgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "cFgibbs"] <- "crit.forb"
				if("cGgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "cGgibbs"] <- "crit.grass"
				if("cSgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "cSgibbs"] <- "crit.shrub"
				if("cTgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "cTgibbs"] <- "crit.tree"

				if("pFgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "pFgibbs"] <- "phenology.forb"
				if("pGgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "pGgibbs"] <- "phenology.grass"
				if("pSgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "pSgibbs"] <- "phenology.shrub"
				if("pTgibbs" %in% gibbs.prop) gibbs.prop[gibbs.prop == "pTgibbs"] <- "phenology.tree"

		#means, sds, and credible intervals
		
			stats <- vector(mode = 'list', length = length(gibbs.list))
				names(stats) <- gibbs.prop
		
			for(j in 1:length(gibbs.list)){
			
				tmp <- matrix(NA,nrow=4,ncol=ncol(get(gibbs.list[j])))
					rownames(tmp) <- c('mean','sd','2.5%','97.5%')
					colnames(tmp) <- colnames(get(gibbs.list[j]))
				
				if(ncol(get(gibbs.list[j]))>1){
					tmp['mean',]  <- apply(get(gibbs.list[j])[post.keep,],2,mean)
					tmp['sd',]    <- apply(get(gibbs.list[j])[post.keep,],2,sd)
					tmp['2.5%',]  <- apply(get(gibbs.list[j])[post.keep,],2,quantile,.025)
					tmp['97.5%',] <- apply(get(gibbs.list[j])[post.keep,],2,quantile,.975)
					}
				if(ncol(get(gibbs.list[j]))==1){
					tmp['mean',]  <- mean(get(gibbs.list[j])[post.keep,])
					tmp['sd',]    <- sd(get(gibbs.list[j])[post.keep,])
					tmp['2.5%',]  <- quantile(get(gibbs.list[j])[post.keep,],.025)
					tmp['97.5%',] <- quantile(get(gibbs.list[j])[post.keep,],.975)
					}
					
				stats[[j]] <- tmp
				
				write.csv(tmp,paste('stats',gibbs.prop[j],'csv',sep='.'))
			
			}
			
		rm(tmp)
		
		#save tables
	tmp <- round(cor(cbind(agibbs,dgibbs,cSgibbs,cGgibbs,bSgibbs,bGgibbs,pSgibbs,pGgibbs,cgibbs[,c(1,2,5)])[post.keep,]),2)
	write.csv(tmp,'covariation.csv')
			

	setwd(dir.out)

#plotting

	if(!('figures' %in% list.files())) dir.create('figures')
	
	dir.out.figures <- paste(dir.out,'figures',sep='/')
	
	setwd(dir.out.figures)

	#transpiration functions
		
	ET.vec <- grep('ETcoef',gibbs.prop)
		
	if(length(ET.vec)>0){
		
		bgibbs <- get(gibbs.list[ET.vec[1]])[post.keep,]
		
		depths <- swSoils_Layers(swIn[[1]])[,1]
		
		dseq <- 1:max(depths)
		
		R.D <- R.DSEQ <- list()
		
		r.d <- r.dseq <- matrix(NA,nrow=nrow(bgibbs),ncol=length(dseq))
		
		for(g in 1:nrow(bgibbs)){
			tmp        <- T_depth.Schenk2002(d = 1:max(depths), params = bgibbs[g,])
			r.d[g,]    <- tmp$r.d/sum(tmp$r.d)
			r.dseq[g,] <- tmp$r.dseq/max(tmp$r.dseq)
			rm(tmp)
		}
		
		R.D[[1]]    <- r.d
		R.DSEQ[[1]] <- r.dseq

		if(length(ET.vec)>1){
			
			for(j in 2:length(ET.vec)){
				bgibbs <- get(gibbs.list[ET.vec[j]])[post.keep,]
			
				for(g in 1:nrow(bgibbs)){
					tmp <- T_depth.Schenk2002(d = 1:max(depths), params = bgibbs[g,])
					r.d[g,] <- tmp$r.d/sum(tmp$r.d)
					r.dseq[g,] <- tmp$r.dseq/max(tmp$r.dseq)
					rm(tmp)
				}
				
				R.D[[j]]    <- r.d
				R.DSEQ[[j]] <- r.dseq
			}
		}		
		
		rm(list = c('r.d','r.dseq'))
			
		jpeg('ET.Responses.jpg',width=6,height=3,units="in",res=300)
			
		par(mfrow=c(1,2),mar=c(4,4,1,1))
		
		plot(apply(R.DSEQ[[1]],2,mean),-dseq,type='l',xlab='cumulative proportion of transpiration', 
			ylab='depth (cm)')	
		lines(apply(R.DSEQ[[1]],2,quantile,.025),-dseq,lty=2)	
		lines(apply(R.DSEQ[[1]],2,quantile,.975),-dseq,lty=2)	


		if(length(ET.vec)>1){
			for(j in 2:length(ET.vec)){

				lines(apply(R.DSEQ[[j]],2,mean),-dseq,lty=2,col=j)	
				lines(apply(R.DSEQ[[j]],2,quantile,.025),-dseq,lty=2,col=j)	
				lines(apply(R.DSEQ[[j]],2,quantile,.975),-dseq,lty=2,col=j)	

			}					
		}
			
			legend('bottomleft',legend=gibbs.prop[ET.vec],text.col=1:length(ET.vec))
	
		xlim <- range(apply(R.D[[1]],2,quantile,c(.025,.975)))
		if(length(ET.vec)>1) {
			
			for(j in 2:length(ET.vec)) 
				xlim <- range(c(xlim,range(apply(R.D[[j]],2,quantile,c(.025,.975)))))
			
		}
		
		
		plot(apply(R.D[[1]],2,mean ),-dseq,type='l',xlab='proportion of transpiration', 
			ylab='depth (cm)',xlim=xlim)	
		lines(apply(R.D[[1]],2,quantile,.025),-dseq,lty=2)	
		lines(apply(R.D[[1]],2,quantile,.975),-dseq,lty=2)	
		
			
		if(length(ET.vec)>1){
			for(j in 2:length(ET.vec)){

				lines(apply(R.D[[j]],2,mean),-dseq,lty=2,col=j)	
				lines(apply(R.D[[j]],2,quantile,.025),-dseq,lty=2,col=j)	
				lines(apply(R.D[[j]],2,quantile,.975),-dseq,lty=2,col=j)	

			}					
		}
		
		dev.off()
		
	}

  
	rm(list=c("bgibbs","R.D","invR","R.DSEQ"))

  ####
    panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r = round((cor(x, y)),digits)
         txt <- format(c(r, 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) cex <- 2.5#/strwidth(txt)
         text(0.5, 0.5, txt, cex = max(cex*sqrt(abs(r)),.1))
     }
####
jpeg('pairs.jpg',width=12,height=12,units='in',res=300)

samp <- cbind(agibbs,dgibbs,bSgibbs,bGgibbs,pSgibbs,pGgibbs,cSgibbs,cGgibbs,cgibbs[,c(1,2,5)])[sample(post.keep,500,replace=FALSE),]
	colnames(samp) <- c('evap.max','evap.rate','drain','d50\nShrub','c\nShrub','d50\nGrass','c\nGrass',
	  		   'LAI_conv\nShrub','Litter\nShrub','Biomass\nShrub','Live.pct\nShrub',
	  		   'LAI_conv\nGrass','Litter\nGrass','Biomass\nGrass','Live.pct\nGrass',
	  		   'crit\nshrub','crit\ngrass','grass','shrub','bareground')
pairs(samp,lower.panel=panel.smooth, upper.panel=panel.cor,cex.labels=1.3)
dev.off()

png('CWP.png',width=3,height=3,units="in",res=300)
par(mfrow=c(1,1),mar=c(3,5,1,1))
boxplot(cbind(grass=cGgibbs[post.keep,],shrub=cSgibbs[post.keep,]),ylab='Critical Water Potential')
dev.off()

	setwd(dir.out)

#data fit
	if(!('fit' %in% list.files())) dir.create('fit')
	
	dir.out.fit <- paste(dir.out,'fit',sep='/')
	
	setwd(dir.out.fit)

	#time-series

	
	#extract data from soilwat runs
	
	j <- 1
	
    samp <- sample(min(post.keep):max(post.keep),2000, replace = FALSE)
		
		#variables to predict
		pred.vars <- c(paste(rep("VWCBULK",5),1:5,sep="_"))
	  et.vars <- c(paste(rep("TRANSP",length(depths)),depths,sep="_"))
	
		#matrix to hold predictions
    pred.mat <- PV.mat <- array(NA, dim = c(length(pred.vars), length(obs[[j]][,1]), length(samp)),
                          dimnames = list(vars = pred.vars, time = 1:length(obs[[j]][,1]),samp = samp))
	  etGrass.mat <- etShrub.mat <- array(NA, dim = c(length(et.vars), length(obs[[j]][,1]), length(samp)),
	                                      dimnames = list(vars = et.vars, time = 1:length(obs[[j]][,1]),samp = samp))
	
    swIn.gibbs <- swIn.new
	
    scalar <- 10000
  
    for(i in 1:length(samp)){
      #Gibbs step
      gg <- samp[i] 
      
      #use all variables allowed to vary in model
      set.names <- prop.names

      for(nn in 1:length(set.names)){
        set.par <- model.init(set.names = set.names[nn],Input = swIn.gibbs, gg = gg)[[1]]
        
        swIn.gibbs <- model.set(swIn.gibbs, set.names[nn], set.par, j = 1)
      }
		
		  #run model
        gibbs.run <- get.pred(swIn.gibbs,file.in=FALSE)
      
      #extract predictions
		  
        if(length(grep("VWCBULK",pred.vars))>0){
          
          dd <- as.numeric(sens.Lyr[sens.Lyr[,'Label']==sites[j],4:8])
          
          for(k in 1:length(grep("VWCBULK",pred.vars))){
            
            if(is.na(dd[k])) next #if no censor, ignore depth
            
            pred.mat[k,,i] <- as.integer(trunc(scalar * gibbs.run[[j]]@VWCBULK@Day[,2+dd[k]]))
            
            PV.mat[k,,i] <- as.integer(trunc(scalar * rnorm(length(pred.mat[k,,i]), pred.mat[k,,i]/scalar, sqrt(vgibbs[gg,1]))))
            
            }

          for(k in 1:length(depths)){
            
            etGrass.mat[k,,i] <- as.integer(trunc(scalar * gibbs.run[[j]]@TRANSP@Day[,2+length(depths)*4+k]))
            etShrub.mat[k,,i] <- as.integer(trunc(scalar * gibbs.run[[j]]@TRANSP@Day[,2+length(depths)*2+k]))
            
           
          }
          
          }
      
		  }
	
  for(i in 1:5){
    jpeg(paste("fit",depths[dd[i]], "cm","jpg",sep="."),width=6,height=2.5,units='in',res=1000)
    par(mfrow=c(1,1),mar=c(4,4,1,1))
    plot(obs[[j]][,"day"],apply(pred.mat[i,,],1,mean)/scalar,type="l",ylim=c(0,.5),xlab="Day", ylab = paste("Bulk VWC at",depths[dd[i]], "cm"))
	    lines(obs[[j]][,"day"],pmax(0,apply(pred.mat[i,,],1,quantile,.025)/scalar),lty=2)
	    lines(obs[[j]][,"day"],pmax(0,apply(pred.mat[i,,],1,quantile,.975)/scalar),lty=2)
	    
      lines(obs[[j]][,"day"],pmax(0,apply(PV.mat[i,,],1,quantile,0.16)/scalar),col="gray40")
	    lines(obs[[j]][,"day"],pmax(0,apply(PV.mat[i,,],1,quantile,0.84)/scalar),col="gray40")
	  
      lines(obs[[j]][,"day"],obs[[j]][,4+i],col='red')
#	    lines(obs[[j]][,"day"],pred.vwc.new[[1]][,i+4],type='l',col='blue')
    dev.off()
  }
  
	
	jpeg("fit.timeseries.jpg",width=6,height=8,units='in',res=1000)
	par(mfrow=c(5,1),mar=c(1,4,1,1))
	for(i in 1:5){
	  plot(obs[[j]][,"day"],apply(pred.mat[i,,],1,mean)/scalar,type="l",ylim=c(0,.5),xlab="Day", ylab = paste("Bulk VWC at",depths[dd[i]], "cm"))
	  lines(obs[[j]][,"day"],pmax(0,apply(pred.mat[i,,],1,quantile,.025)/scalar),lty=2)
	  lines(obs[[j]][,"day"],pmax(0,apply(pred.mat[i,,],1,quantile,.975)/scalar),lty=2)
	  
	  lines(obs[[j]][,"day"],pmax(0,apply(PV.mat[i,,],1,quantile,0.16)/scalar),col="gray40")
	  lines(obs[[j]][,"day"],pmax(0,apply(PV.mat[i,,],1,quantile,0.84)/scalar),col="gray40")
	  
	  lines(obs[[j]][,"day"],obs[[j]][,4+i],col='red')
	  #	    lines(obs[[j]][,"day"],pred.vwc.new[[1]][,i+4],type='l',col='blue')
	}
	dev.off()
	
  
  ##need to calculate total transpiration and transpiration by 20 cm bins
  
  et.all <- function() etGrass.mat + etShrub.mat
  week <- rep(1:1000,each=7)[1:length(obs[[j]][,"day"])]

  et.total <- tapply(apply(apply(et.all(),c(1,2),mean),2,sum),week,sum)/scalar
  et.hi <- et.lo <- et.total*0

  tmp <- apply(et.all(),c(2,3),sum)
  for(w in 1:max(week)){
    wtmp <- apply(tmp[week==w,],2,sum)
    et.lo[w] <- quantile(wtmp,.025)/scalar
    et.hi[w] <- quantile(wtmp,.975)/scalar
    
  }
  
  rm(list = c("tmp","wtmp"))
  
  
	dbins <- seq(0,165,by=33)
  
  et.dbin <- matrix(NA, nrow = max(week), ncol = length(dbins) - 1)

	et.tmp <- apply(et.all(),c(1,2),mean)
	
  for(j in 1:(length(dbins)-1)){
    
    maxd <- max(which(depths < dbins[j + 1]))
    mind <- min(which(depths >= dbins[j]))
    
    et.dbin[,j] <- tapply(apply(et.tmp[mind:maxd,],2,sum)/scalar + 
                          et.tmp[maxd+1,]/scalar* (dbins[j+1] - depths[maxd])/ diff(depths)[maxd],
                          week,sum)
  }
  
jpeg("ProportionTrans.Depth.jpg",width=8,height=6,units="in",res=100)
  par(mfrow=c(2,1),mar=c(1,4,1,1))
	
  plot(et.total,type="l",ylab="Mean Total Transpiration (cm/week)")
	  lines(et.lo,lty=2)
	  lines(et.hi,lty=2)
	
  plot(et.dbin[,1]/et.total,type="l",ylim=c(0,1),ylab="proportion of transpiration")
  for(j in 2:5) lines(et.dbin[,j]/et.total,col=j)
  legend("right",legend=paste(dbins[-length(dbins)],"cm to",dbins[-1],"cm"),text.col=1:5,
         bty="n",cex=.75)
  dev.off()
  
  
  plot(tapply(apply(et.all()[1,,],1,mean),week,sum)/scalar*c(depths[1],diff(depths))[1],
       type="l",ylab="Weekly transpiration (mm?)")
  for(d in 2:length(depths))
    lines(tapply(apply(et.all()[d,,],1,mean),week,sum)/scalar*c(depths[1],diff(depths))[d],
          col=d)
	
library(ade4)
  
  s.ind <- grep("sand",colnames(sgibbs))
  c.ind <- grep("clay",colnames(sgibbs))
  
  for(d in 1:length(s.ind)){
  
    mat <- as.data.frame(cbind(sand = sgibbs[samp,s.ind[d]],
                               clay = sgibbs[samp,c.ind[d]],
                               silt = 1,
                              depth = d))
    mat[,3] <- mat[,3] - apply(mat[,1:2],1,sum)

    if(d == 1) MAT <- mat  
    if(d > 1)  MAT <- rbind(MAT,mat)
  }
  
  jpeg("soil.textures.jpg", width=6,height=6,units="in",res=1000)
  Tplt <- triangle.plot(ta = MAT[,-4],
                       scale = FALSE,
               labeltriangle = TRUE,
               show.position = FALSE,
                        min3 = c(0,0,0),
                        max3 = c(1,1,1),
                      cpoint = 0)
  points(Tplt, col = MAT[,"depth"],pch=20,cex=.5)
  
  d1 <- c(0,depths[-length(depths)])
  d2 <- depths
  
  legend("topright",legend=paste(tapply(d1,soils.ind,min), "cm to",tapply(d2,soils.ind,max),"cm"),
         pch=20,col=1:length(depths),bty="n")
  
  dev.off()
  
  setwd(dir.out)
	
	#pred vs. obs

