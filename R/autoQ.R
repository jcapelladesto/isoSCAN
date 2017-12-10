#' Read files and quantify each isotopologue
#' 
#' Read each file and quantify both Area and Maxo of each metabolite isotopologues, by evaluating the RT window of each metabolite to minimize background noise
#' @param SampleFiles List of files to be processed (see \emph{base::list.files()})
#' @param formulaTable Data.frame including CompoundName, Formula, mz, RT and NumAtoms (see \link{buildFormulaTable})
#' @param RTwin Retention time window (in seconds) to extract the metabolite isotopologues
#' @param IsoMass Stable isotope mass difference, for C13: 1.003355
#' @param ScanE Absolute error in which to search the Maxo scan
#' @param IntE Absolute error in which to integrate the area of the isotopologue. Should be larger than \emph{ScanE} and smaller than 0.5
#' @param fit.p Max linear model p-value for a peak to have its area quantified, otherwise only Maxo is returned
#' @return A data.frame including the Area and Maxo for each file and metabolite isotopologues. If peak was missing, values are NA.
#' @export

autoQ <- 
	function(SampleFiles=NULL,formulaTable=NULL, IsoMass = 1.003355, RTwin = 5, toplot = F, ScanE = 0.1, IntE= 0.4, fit.p = 0.05)
{
	integrations_files <- lapply(SampleFiles,function(fn){
		cat(fn)
		# read raw data
		aa <- mzR::openMSfile(fn)
		h <- mzR::header(aa)
		p <- mzR::peaks(aa)
		h <- lapply(1:length(p),function(x){
			rt <- h[which(h$acquisitionNum==x),"retentionTime"]
			res <- cbind(p[[x]],rep(x=rt,times=nrow(p[[x]])),rep(x=x,times=nrow(p[[x]])))
			return(res)
		})
		h <- do.call("rbind",h)
		h <- as.data.frame(h)
		colnames(h) <- c("basePeakMZ","basePeakIntensity","retentionTime","ScanNum")
		fn  <- gsub("\\..*$","",fn)
		fn  <- gsub(pattern="[[:punct:]]",replacement="_",fn)
		# process data
		if(toplot){pdf(file=paste0("./QQQ_IntegrationPlots_",fn,".pdf"))}
		autoInt <- lapply(1:nrow(formulaTable), function(i){ #1:nrow(formulaTable)
			vr <- formulaTable[i,]
			message(paste0(c("Processing: ",as.character(vr$CompoundName))))
			vrMZ  <-  vr$mz
			vrRT <- vr$RT
			vrRTRan <- c(vrRT-RTwin,vrRT+RTwin)
			vrMZRan <- round(sapply(0:vr$NumAtoms,function(x) vrMZ+IsoMass*x),2)
			eic <- h[which(h$retentionTime>vrRTRan[1] &	h$retentionTime<vrRTRan[2]) ,]
			
			## evaluate and recalculate RT for each metabolite isotopologues
			suppressWarnings(RTWin_mz <- lapply(vrMZRan, function (j){ 
				#j <- vrMZRan[1]
				jRT <- NA
				mzdiff <- abs(eic$basePeakMZ-j)
				if (min(mzdiff)<ScanE){ # 
					mzScan <- which(mzdiff==min(mzdiff)) #Method 1 Scan, min m/z diff = center of peak
					toInt <- eic[mzScan,]
					# 1 max
					maxScan <- which(toInt$basePeakIntensity==max(toInt$basePeakIntensity))[1]
					maxScanRT <- toInt$retentionTime[maxScan]
					maxScanI <- toInt$basePeakIntensity[maxScan]
					# 1st minimum
					minScan <- which(toInt$basePeakIntensity==min(toInt$basePeakIntensity))[1]
					minScanRT <- toInt$retentionTime[minScan]
					minScanI <- toInt$basePeakIntensity[minScan]
# 					# median intensity in the EIC: Background in the area
# 					medI <- median(toInt$basePeakIntensity)
					timeWindowTest <- F
					# if the intensity distribution is not normal, highly probable a peak
					if(shapiro.test(toInt$basePeakIntensity)$p.value < 0.01) {timeWindowTest <- T}
					
 					#plot(toInt$retentionTime,toInt$basePeakIntensity);Sys.sleep(0.2)
					
					while(timeWindowTest){
						if(minScanRT<maxScanRT){ # look for second min: left side 
							tempi <- seq(from=maxScan,to=nrow(toInt))
							R_L <- T
						}else{ # right side
							tempi <- seq(from=1,to=maxScan)	
							R_L <- F
						}
						if(length(tempi)<5){ # not enought scans at both sides of "peak"   ## consider delete
							timeWindowTest <- T
						}else{
						temp <- toInt[tempi,]
						# 2nd minimum
						minScan2 <- which(temp$basePeakIntensity==min(temp$basePeakIntensity))[1]
						minScanI2 <- temp$basePeakIntensity[minScan2]
						if(R_L){
							minScan2 <- maxScan+minScan2-1
							minScanRT2 <- toInt$retentionTime[minScan2]		  	
						}else{
							minScanRT2 <- toInt$retentionTime[minScan2]	 #?	
						}
						# Peak intensity cutoff
						if(abs(minScan-minScan2)<=10){ 
						# enough scans to define a peak?
							jRT <- NA
							timeWindowTest <- T  ##*	
						}else{
							medI <- mean(toInt$basePeakIntensity[-(minScan:minScan2)])
							mNoise <- mean(c(minScanI,minScanI2)) # fit a lm?
							if(is.na(medI)){ medI <- mNoise	}
							evalSN <- ceiling((maxScanI-mNoise)/(medI-mNoise)) #use only this?
							if (mNoise >= medI){ evalSN <- (maxScanI-mNoise)/medI }
							if(evalSN>=3){
								jRT <- sort(c(minScanRT,minScanRT2))
# 							plot(toInt$retentionTime,toInt$basePeakIntensity);Sys.sleep(0.2)
								# do peak fit in RT set by user abs(minScanRT-minScanRT2)>RTwin
								if(vrRT>jRT[1]&vrRT<jRT[2]){ 
									timeWindowTest <- F
								}
							}else{ # if SN is too bad
								timeWindowTest <- T 
							}
						}
						}

						if(timeWindowTest){
							#bad max, delete Scan, recalculate max and min
							toInt <- toInt[-maxScan,]
							maxScan <- which(toInt$basePeakIntensity==max(toInt$basePeakIntensity))[1]
							maxScanRT <- toInt$retentionTime[maxScan]
							maxScanI <- toInt$basePeakIntensity[maxScan]
							
							minScan <- which(toInt$basePeakIntensity==min(toInt$basePeakIntensity))[1]
							minScanRT <- toInt$retentionTime[minScan]
							minScanI <- toInt$basePeakIntensity[minScan]
						}
						# add flag IF PEAK IS TOO WIDE?
						if(maxScanI==0|nrow(toInt)<=10){ #cannot define a peak 
							jRT <- NA
							timeWindowTest <- F
						}
					}
				}
				return(jRT)
			}))
			
			vrMZRani <- sapply(RTWin_mz,function(x) is.na(x)[1])
			RTWin_mz <- RTWin_mz[!vrMZRani] #  RTWin_mz <- RTWin_mz[which(sapply(RTWin_mz,function(x) is.na(x)[1])==F)]
			
			if(all(is.na(RTWin_mz))){
				print("Not found")
				return()
}
			
			# new RT and window are calculated
			RTweights <- sapply(RTWin_mz,function(x) abs(mean(x) - vrRT) )
			RTweights <- (1/(10^RTweights))^2
			RTweights <- RTweights/sum(RTweights)
			vrRTRan <- sapply(c(1,2),function(x){
				weighted.mean(sapply(RTWin_mz,function(y){y[x]}),w=RTweights)
			})
			vrRT <- mean(vrRTRan)
			eic <- h[which(h$retentionTime>vrRTRan[1] &	h$retentionTime<vrRTRan[2]),]
			
			# Now Integrate
			whoInt <- lapply(1:length(vrMZRan), function (j){
				if(vrMZRani[j]){
					return(c(NA,1))
				}else{
				mz <- vrMZRan[j]
				mzdiff <- abs(eic$basePeakMZ-mz)
				mzScan <- which(mzdiff==min(mzdiff)) 
				toInt <- eic[mzScan,]
#  				plot(toInt$retentionTime,toInt$basePeakIntensity,main=c(i,j));Sys.sleep(0.2)
				# define noise prior to integration
				maxidx <- which(toInt$basePeakIntensity==max(toInt$basePeakIntensity))[1]
				maxScanI <- toInt$basePeakIntensity[maxidx]
				minScanI  <- c(min(toInt[1:maxidx,"basePeakIntensity"]),
							    min(toInt[maxidx:nrow(toInt),"basePeakIntensity"]))
				minidx <- c(which(toInt$basePeakIntensity==minScanI[1]),which(toInt$basePeakIntensity==minScanI[2]))
				minScanI <- mean(minScanI)
				delidx <- (minidx[1]:minidx[2])
				if(length(delidx)==nrow(toInt)){
					medI <- minScanI
				}else{
					medI <- mean(toInt$basePeakIntensity[-delidx])
				}

				if(toplot){
					plot(toInt$retentionTime,toInt$basePeakIntensity,main=c(i,mz))
					Sys.sleep(0.25)
				}
				# Evaluate peak
				if(nrow(toInt)<10){
					return(c(NA,1))
				}else{
					maxidx <- which(toInt$basePeakIntensity==max(toInt$basePeakIntensity))[1]
# 					toInt <- toInt[-maxidx,]
# 					maxidx <- which(toInt$basePeakIntensity==max(toInt$basePeakIntensity))[1]
					minScanI  <- c(min(toInt[1:maxidx,"basePeakIntensity"]),
											 min(toInt[maxidx:nrow(toInt),"basePeakIntensity"]))
					cor_idx <- list(1:(maxidx-1),(maxidx+1):nrow(toInt))
					c_m <- sapply(cor_idx,length)
					cor_idx <- list(cor_idx[[1]][(1+c_m[[1]]-min(c_m)):c_m[[1]]],
													cor_idx[[2]][min(c_m):1])
					c_m <- sum(sapply(cor_idx,length))
					lmr <- suppressWarnings(lm(x~y,data=data.frame("x"=toInt$basePeakIntensity[cor_idx[[1]]],
        															                 "y"=toInt$basePeakIntensity[cor_idx[[2]]])))
				
				if(is.na(lmr$coefficients[["y"]])|c_m<=4){
					lmr_pval <- 1
				}else{
					sum_lmr <- summary(lmr)
					lmr_pval <- sum_lmr$coefficients["y","Pr(>|t|)"]
					if(is.na(lmr_pval)){lmr_pval <- 1}
				}
				# try rerunning lm fit cleaning "outliers"
				if(lmr_pval>fit.p & lmr_pval<0.1){
					torm <- round(log10(c_m/2))
					if(torm>0){
					lmr_pval <- min(sapply(1:torm,function(x){ 
						torm_x <- order(abs(lmr$residuals),decreasing=T)[x]
						cor_idx <- lapply(cor_idx,function(y){
							if(torm_x %in% y){
								y <- y[-which(y==torm_x)]
							}else{
								y <- y[1:(length(y)-1)]
							}
							return(y)
						})
# 						cor_idx[[1]] <- cor_idx[[1]][-torm_x]
# 						cor_idx[[2]] <- cor_idx[[2]][-torm_x]
						lmr <- suppressWarnings(lm(x~y,data=data.frame("x"=toInt$basePeakIntensity[cor_idx[[1]]],
																													 "y"=toInt$basePeakIntensity[cor_idx[[2]]])))
						sum_lmr <- summary(lmr)
						lmr_pval <- sum_lmr$coefficients["y","Pr(>|t|)"]
						return(lmr_pval)
					}))
					}
				}
				}
				}
				return(c(medI,lmr_pval))
				
				})
		
			IntAll <- sum(sapply(whoInt,function(x) (x[2]<0.05)),na.rm=T) > (0.5*sum(!sapply(whoInt,function(x) is.na(x[1]))))
			#IntAll <- sum((whoInt<0.05),na.rm=T)>(0.5*sum(!is.na(whoInt)))
			
			Integrate <- sapply(1:length(vrMZRan),function(j){
				res <- c(NA,NA)
				if(is.na(whoInt[[j]][1]) && !IntAll){
					return(res)
				}else{
					mz <- vrMZRan[j]
					mzdiff <- abs(eic$basePeakMZ-mz)
					# Now get the area
					mzScan <- which(abs(eic$basePeakMZ-mz)<IntE)
					toInt <- eic[mzScan,]
					maxScanI <- max(toInt$basePeakIntensity)
					# 					Ilim <- median(toInt$basePeakIntensity)
					Ilim <- whoInt[[j]][1] # get the mean noise intensity
					if(is.na(Ilim)){
						Ilim <- median(toInt$basePeakIntensity)*1.2 # estimate noise with median
					}else{
						Ilim <- 1.2*sqrt(maxScanI-Ilim)+Ilim # function is ok but idk why
						# 					Ilim <- Ilim*1.2 # cut a bit above the noise level  Ilim*(1+(1/log10(maxScanI-Ilim)))	
					}
					toInt$basePeakIntensity <- toInt$basePeakIntensity-Ilim 
					toInt <- toInt[which(toInt$basePeakIntensity>0),]
#  					plot(toInt$retentionTime,toInt$basePeakIntensity)
					if(nrow(toInt)>10){
						if(whoInt[[j]][2]<fit.p || IntAll){
							res <- c(sum(toInt$basePeakIntensity),max(toInt$basePeakIntensity)) 
						}else{
							res <- c(NA,max(toInt$basePeakIntensity))
						}
					}
				}
				return(res)
			})
			
			# output format
			intRes <- data.frame("CompoundName"=rep(vr$CompoundName,times=vr$NumAtoms+1),
													 "Isotopologue"=paste0("M+",0:vr$NumAtoms),
													 "Area"=Integrate[1,],
													 "Maxo"=Integrate[2,])
			return(intRes)
		})
		if(toplot){dev.off()}
		return(autoInt)
	})
	names(integrations_files) <- SampleFiles
	
	check_int <- sapply(integrations_files,sapply,length)
	check_int_c <- apply(check_int,1,function(x) all(x==0))
	if (any(check_int_c)){ 	
		print(c("Compound/s lost:", as.character(formulaTable$CompoundName[which(check_int_c)])))
}
	check_int_compounds <- c(1:nrow(check_int))[!check_int_c]
	check_int_f <- apply(check_int,2,function(x) all(x==0))
	if (any(check_int_f)){
		print(c("Files with 0 compounds quantifiable:", as.character(SampleFiles[which(check_int_f)])))
}
	check_int_files <- c(1:ncol(check_int))[!check_int_f]
		
	# put results in order
	FinalInt <- lapply(check_int_compounds,function(i){
		Metdata <- lapply(check_int_files, function(j){
			x <- SampleFiles[j]
			m <- formulaTable$CompoundName[i]
			mdata <- integrations_files[[x]][[i]]
			if(is.null(mdata)){
				return()
			}else{
				if(j==1){
					return(mdata)
				}else{
					return(mdata[,3:4])
				}
			}
		})
		check_null <- sapply(Metdata,is.null)
		check_null_i <- which(check_null)
		if(any(check_null)){
			for (z in check_null_i){
				if(z==1){
					Metdata[[1]] <- data.frame("CompoundName"=rep(formulaTable$CompoundName[i],times=formulaTable$NumAtoms[i]+1),
																		 "Isotopologue"=paste0("M+",0:(formulaTable$NumAtoms[i])),
																		 "Area"=rep(NA,times=formulaTable$NumAtoms[i]+1),
																		 "Maxo"=rep(NA,times=formulaTable$NumAtoms[i]+1))
				}else{
					Metdata[[z]] <- data.frame("Area"=rep(NA,times=formulaTable$NumAtoms[i]+1),
																		 "Maxo"=rep(NA,times=formulaTable$NumAtoms[i]+1))
				}
			}
		}
		Metdata <- do.call("cbind",Metdata)
		colnames(Metdata)[3:ncol(Metdata)] <- paste(colnames(Metdata)[3:ncol(Metdata)],
																								rep(gsub(pattern=".mzXML",replacement="",SampleFiles)
																										,each=2), sep="_")
		return(Metdata)
	})
	FinalInt <- do.call("rbind",FinalInt) 
	
	FinalInt$Isotopologue <- factor(FinalInt$Isotopologue,levels=paste0("M+",0:(max(formulaTable$NumAtoms)+1)))
	FinalInt$CompoundName <- factor(FinalInt$CompoundName,levels=unique(FinalInt$CompoundName))
	return(FinalInt)
}