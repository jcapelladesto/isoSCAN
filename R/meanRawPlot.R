#' Read all scans raw data, calculate average spectrum for all files and plot
#' 
#' Raw data plotting average spectra of all SampleFiles and compound for quality control pourposes
#' @param SampleFiles List of files to be plotted (see \emph{base::list.files()})
#' @param formulaTable Data.frame including CompoundName, Formula, mz, RT and NumAtoms (see \link{buildFormulaTable})
#' @param RTwin Retention time window (in seconds) to extract the metabolite isotopologues
#' @param topdf File directory where the plots will be saved in pdf format. One pdf file per mzML file in Sampledir
#' @return Plots per each file and compound including all isotopologues
#' @note To plot files separately, use \link{rawPlot}
#' @export

meanRawPlot <- function(SampleFiles=NULL,formulaTable=NULL,RTwin=5,topdf=NULL){
		try(dev.off(),silent=T)
message("Reading raw data")
rawdat <- lapply(SampleFiles,function(fn){
	# add progress information?
	aa <- openMSfile(fn)
	h <- header(aa)
	p <- peaks(aa)
	h <- lapply(1:length(p),function(x){
		rt <- h[which(h$acquisitionNum==x),"retentionTime"]
		res <- cbind(p[[x]],rep(x=rt,times=nrow(p[[x]])),rep(x=x,times=nrow(p[[x]])))
		return(res)
	})
	h <- do.call("rbind",h)
	h <- as.data.frame(h)
	colnames(h) <- c("basePeakMZ","basePeakIntensity","retentionTime","ScanNum")
	return(h)
})
message("Plotting spectra")
if(!is.null(topdf)) {pdf(file=topdf)}
for (i in 1:nrow(formulaTable)){ 
	vr <- formulaTable[i,]
	vrMZ  <-  vr$mz
	vrRT <- vr$RT
	vrRTRan <- c(vrRT-RTwin,vrRT+RTwin)
	vrMZRan <- c(vrMZ-1,vrMZ+vr$NumAtoms+1) 
	Metdata <- lapply(1:length(SampleFiles),function(j){ # read all files and get data.frame  
		msdata <- rawdat[[j]]
		msdata <- msdata[which(msdata $retentionTime>vrRTRan[1] &
													 	msdata $retentionTime<vrRTRan[2] &
													 	msdata $basePeakMZ>vrMZRan[1] & 
													 	msdata $basePeakMZ<vrMZRan[2]	),
										 c("retentionTime","basePeakMZ","basePeakIntensity") ]
		msdata$sample <- j
		return(msdata)
	})

	Metdata <- do.call("rbind",Metdata)
	Metdata <- as.data.frame(Metdata)
 	
	utimes <- sort(unique(Metdata$retentionTime))	 # unique rt, some need corrections, small errors
 	sorted_times <- sort(Metdata$retentionTime)
	round_times <- round(sorted_times,digits=2) # all rt, rounded, no error
	ttimes <- sapply(utimes,function(rt){ length(which(Metdata$retentionTime==rt )) } )
	rtimes <- round(utimes,digits=2) # round utimes, no corrections
 	rrtimes <- sapply(rtimes,function(rt){ length(which(round_times==rt )) } )
	timesdf <- data.frame(utimes,rtimes,ttimes,rrtimes)
	max_rrtimes <- max(rrtimes) # number of ideal mz scans per rt point

	# loop over the rts so each rt point has ideal mz scan number
	timesList <- list()
	ignore_idx <- c()
	for(m in 1:nrow(timesdf)){
		if(!(m%in%ignore_idx)){
			return_vals <- T
			tempdf <- timesdf[m,]
			rtimes_idx <- which(timesdf$rtimes==tempdf$rtimes)
			while(sum(timesdf$ttimes[rtimes_idx])<max_rrtimes){ # some wont work with rounded so use this
				rtimes_idx <- c(rtimes_idx,rtimes_idx[length(rtimes_idx)]+1)
				if(is.na(sum(timesdf$ttimes[rtimes_idx]))){ # if we dont have ideal scan number in that rt, lose it
					return_vals <- F
					break
				}
			}
			if(return_vals){
				temptimes <- timesdf$utimes[rtimes_idx]
				timesList <- append(timesList,list(temptimes))
				ignore_idx <- unique(append(ignore_idx,rtimes_idx))
			}
		}
	}

	newtimes <- sapply(timesList,mean) # mean rt corrected
	# rtwin to calculate mean intensity
	timesList <- lapply(timesList,function(x) return(c(min(x),max(x)))) 
	Metdata2 <- sapply(timesList,function(rt){  # calculate mean intensity
		sapply(sort(unique(Metdata$basePeakMZ)),function(mz){
			idx <- which(Metdata$basePeakMZ==mz 
									 & Metdata$retentionTime>=rt[1] 
									 & Metdata$retentionTime<rt[2])
			if(rt[1]==rt[2]){idx <- which(Metdata$basePeakMZ==mz 
																		& Metdata$retentionTime==rt[2])}
			return(mean(Metdata$basePeakIntensity[idx]) )
		})
	})
	
	Metdata2 <- as.data.frame(Metdata2) # matrix : mz by rt intensity values
	colnames(Metdata2) <- newtimes
	rownames(Metdata2) <- sort(unique(Metdata$basePeakMZ))
	nacol <- which(apply(Metdata2,2,function(x) all(is.na(x))))
	if(length(nacol)>0){	Metdata2 <- Metdata2[,-nacol];	newtimes <- colnames(Metdata2)}
	
	Metdata2 <- lapply(colnames(Metdata2),function(colN){ # my melt style apply
		vals <- Metdata2[,colN]
		res <- data.frame("basePeakIntensity"=as.numeric(vals),
							 "retentionTime"=as.numeric(rep(colN,times=length(vals))),
							 "basePeakMZ"=as.numeric(rownames(Metdata2)))
		return(res)
	})
	Metdata2 <- do.call("rbind",Metdata2)
	narow <- which(is.na(Metdata2$basePeakIntensity))
	if(length(narow)>0){	Metdata2 <- Metdata2[-narow,]}
	
	myPlot <- ggplot(data=Metdata2,
									 aes(x=basePeakMZ, y= retentionTime, color = basePeakIntensity))+
		ggtitle(vr$CompoundName)+
		geom_hline(yintercept=vrRT,color="red")+
		scale_colour_gradientn(colours = terrain.colors(5))+
		geom_point(size=5)+ylim(vrRTRan)
	
	tryCatch({print(myPlot);Sys.sleep(0)}, error=function(e){
		myPlot <- ggplot(data=Metdata2,
										 aes(x=basePeakMZ, y= retentionTime))+ggtitle(vr$CompoundName)+
			geom_hline(yintercept=vrRT,color="red")+
			geom_point(size=5)+ylim(vrRTRan)
		
		print(myPlot);Sys.sleep(0)
	})
}
if(!is.null(topdf)) {	dev.off()}
}
