#' Mean +- Standard Deviation
#' 
#' Error plotting for metBarPlot
#' @export
myFunc.errorbars  <-  function(x) {
	x <- na.omit(x)
	data.frame(y = mean(x),
						 ymin = mean(x)-sd(x),
						 ymax = mean(x)+sd(x))
}