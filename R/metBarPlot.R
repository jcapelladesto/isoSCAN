#' Creates a PDF file with barplots per each compound
#' 
#' Creates a PDF file with barplots (Mean +- SD) per each compound and isotopologues using ggplot2
#' @param x The data.frame result from \link{autoQ} or \link{QTransform}
#' @param topdf File directory where the plot should be saved. If \code{NULL} plots will be sent to R plotting device
#' @param val.to.plot Values to use for plotting. Either "Area" or "Maxo". If values were transformed with \code{Qtransform}, indicate the value used for transformation
#' @param groups Sample groups. A vector of length equal to number of samples that will match the samples order
#' @param ylabel Label for y axis in the plot. If not indicated, val.to.plot value will be used as default.
#' @param ... other arguments that will be passed to \emph{grDevices::pdf}
#' @return Barplots per each compound including deviation in each isotopologue.
#' @export

# label for they yaxis , delete label decision alg
metBarPlot <- 
	function(x=NULL, topdf=NULL, val.to.plot = "Area", groups = c(), ylabel = NULL, ...)
{
	try(dev.off(),silent=T)
	FinalInt <- x
	sleepval <- 0.1
	if(is.null(ylabel)){
		warning("'ylabel' value is NULL. Using 'val.to.plot' value as default label for y axis.")
		ylabel <- val.to.plot 
	}
	if(!is.null(topdf)){pdf(topdf,...);sleepval <- sleepval-0.1}
	plotapply <- lapply(levels(FinalInt$CompoundName),function(x){
		cat(c("Plotting compoundName:",x,"\n"))
		selcol2 <- c(1:2,grep(val.to.plot,colnames(FinalInt)))
		pos <- which(FinalInt$CompoundName==x)
		plotdata <- FinalInt[pos,selcol2]
		
		suppressMessages(resmat <- melt(plotdata,id.vars=c("CompoundName","Isotopologue")))
    resmat$variable <- rep(groups,each=nrow(plotdata))
		resmat$variable <- as.factor(resmat$variable)
		if(all(is.na(resmat$value))){
			message("Not quantified")
			return()
		}
		Myplot <- ggplot(data=resmat,
										 aes(x=Isotopologue,y=value, colour = Isotopologue, 
										 		fill=Isotopologue,group=variable))+
			stat_summary(fun.data="myFunc.errorbars",geom="bar")+
			stat_summary(fun.data ="myFunc.errorbars",
									 geom = "errorbar", width = 0.5,size = 0.7)+
			facet_grid(.~variable)+theme(axis.title.x=element_blank(),
															 axis.text.x=element_blank(),
															 axis.ticks.x=element_blank(),
															 panel.grid.minor=element_blank(),
															 panel.grid.minor.x=element_blank(),
															 panel.grid.major.x=element_blank())+
			ggtitle(x)+ylab(ylabel)
		

		if (all(resmat$value<=1)) {
			try(Myplot <- Myplot+scale_y_continuous(breaks=seq(0,max(resmat$value,na.rm=T),0.1)),silent=T)
		}
		suppressWarnings(print(Myplot));Sys.sleep(sleepval)

		return()
	})
	if(!is.null(topdf)){dev.off()}
}