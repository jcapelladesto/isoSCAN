#' Read all scans raw data and plot
#' 
#' Raw data plotting per each file and compound for quality control pourposes
#' @param SampleFiles List of files to be plotted (see \emph{base::list.files()})
#' @param formulaTable Data.frame including CompoundName, Formula, mz, RT and NumAtoms (see \link{buildFormulaTable})
#' @param RTwin Retention time window (in seconds) to extract the metabolite isotopologues
#' @param topdf File directory where the plots will be saved in pdf format. One pdf file per mzML file in Sampledir
#' @return Plots per each file and compound including all isotopologues
#' @note It is recommended to run it for a few files only, plotting a high number of files may take time.
#' @export

rawPlot <- 
	function(SampleFiles=NULL,formulaTable=NULL,RTwin=5,topdf=NULL,plotall=F){
		try(dev.off(),silent=T)
		for (fn in SampleFiles){
			# add progress information?
			oMS <- mzR::openMSfile(fn)
			h <- mzR::header(oMS)
			p <- mzR::peaks(oMS)
			h <- lapply(1:length(p),function(x){
				rt <- h$retentionTime[x]
				res <- cbind(p[[x]],rep(x=rt,times=nrow(p[[x]])))
				return(res)
			})
			h <- do.call("rbind",h)
			h <- as.data.frame(h)
			colnames(h) <- c("basePeakMZ","basePeakIntensity","retentionTime")
			fn  <- rmfileExt(fn,ext.regexp = "\\.mz(x)?ML")
			if(!is.null(topdf)) {pdf(file=paste0(topdf,fn,".pdf"))}
			for (i in 1:nrow(formulaTable)){
			  fTi <- formulaTable[i,]
			  message(paste0(c("Processing: ",as.character(fTi$CompoundName))))
			  fTiMZ  <-  fTi$mz
			  fTiRT <- fTi$RT
			  fTiRTRan <- c(fTiRT-RTwin,fTiRT+RTwin)
			  fTiFormula <- as.character(fTi$Formula)
			  natom <- fTi$NumAtoms
			  fTiMZRan <- c(fTiMZ-1,fTiMZ+(natom+1))
				
				if (plotall){
					myPlot <- ggplot(data=h,aes(x=basePeakMZ, y= retentionTime, color = log10(basePeakIntensity)))+
						ggtitle("All scans")+
						scale_colour_gradientn(colours = terrain.colors(5))+
						geom_point()+
					  theme(legend.title = NULL)
				
				}
				
				msdata <- h[which(h$retentionTime>fTiRTRan[1] &
														h$retentionTime<fTiRTRan[2] &
														h$basePeakMZ>fTiMZRan[1] & 
														h$basePeakMZ<fTiMZRan[2]	),c("retentionTime","basePeakMZ","basePeakIntensity") ]
				
				myPlot <- ggplot(data=msdata,
												 aes(x=basePeakMZ, y= retentionTime, color = basePeakIntensity))+
				  ggtitle(fTi$CompoundName)+
					geom_hline(yintercept=fTiRT,color="red")+
					scale_colour_gradientn(colours = terrain.colors(5))+
					geom_point(size=4)+
				  ylim(fTiRTRan)
						
				tryCatch({print(myPlot);Sys.sleep(0)}, error=function(e){
					myPlot <- ggplot(data=msdata,
													 aes(x=basePeakMZ, y= retentionTime))+
					  ggtitle(fTi$CompoundName)+
						geom_hline(yintercept=fTiRT,color="red")+
						geom_point(size=4)+
					  ylim(fTiRTRan)
					
					print(myPlot);Sys.sleep(0)
				})
				plotall <- F
			}
			if(!is.null(topdf)) {	dev.off()}
		}
		# return()
	}
