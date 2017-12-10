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
			fn  <- gsub("\\..*$","",fn)
			fn  <- gsub(pattern="[[:punct:]]",replacement="_",fn)
			if(!is.null(topdf)) {pdf(file=paste0(topdf,fn,".pdf"))}
			for (i in 1:nrow(formulaTable)){
				vr <- formulaTable[i,]
				vrMZ  <-  vr$mz
				vrRT <- vr$RT
				vrRTRan <- c(vrRT-RTwin,vrRT+RTwin)
				vrMZRan <- c(vrMZ-1,vrMZ+vr$NumAtoms+1)
				
				if (plotall){
					myPlot <- ggplot(data=h,aes(x=basePeakMZ, y= retentionTime, color = basePeakIntensity))+
						ggtitle("All scans")+
						scale_colour_gradientn(colours = terrain.colors(5))+
						geom_point()
				}
				
				msdata <- h[which(h$retentionTime>vrRTRan[1] &
														h$retentionTime<vrRTRan[2] &
														h$basePeakMZ>vrMZRan[1] & 
														h$basePeakMZ<vrMZRan[2]	),c("retentionTime","basePeakMZ","basePeakIntensity") ]
				
				myPlot <- ggplot(data=msdata,
												 aes(x=basePeakMZ, y= retentionTime, color = basePeakIntensity))+ggtitle(vr$CompoundName)+
					geom_hline(yintercept=vrRT,color="red")+
					scale_colour_gradientn(colours = terrain.colors(5))+
					geom_point(size=5)+ylim(vrRTRan)
						
				tryCatch({print(myPlot);Sys.sleep(0)}, error=function(e){
					myPlot <- ggplot(data=msdata,
													 aes(x=basePeakMZ, y= retentionTime))+ggtitle(vr$CompoundName)+
						geom_hline(yintercept=vrRT,color="red")+
						geom_point(size=5)+ylim(vrRTRan)
					
					print(myPlot);Sys.sleep(0)
				})
				plotall <- F
			}
			if(!is.null(topdf)) {	dev.off()}
		}
		return()
	}
