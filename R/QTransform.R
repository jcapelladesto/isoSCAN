#' Area or Maxo value transformation 
#' 
#' Creates a data.frame per each compound and isotopologues calculating required transformation, prior to export or plotting
#' @param x The data.frame result from \link{autoQ}
#' @param val.to.use Type of quantification to transform. Either "Area" or "Maxo"
#' @param val.trans Set it to "P" for labelling percentages or "R" for relative quantification to the Monoisotopic peak(for isotopic pattern observation)
#' @return  A data.frame with transformed values. Can be input to \link{metBarPlot}.
#' @export

QTransform <- function(x=NULL,  val.to.use = "area", val.trans = c("P","R")){
	if(length(val.trans)==2){stop("Please select an option for 'value.trans': 'P' for percentages or 'R' for relative quantification")}	
	FinalInt <- x
	Transapply <- lapply(unique(FinalInt$CompoundName),function(x){ 
		selcol2 <- grep(val.to.use,colnames(FinalInt))
		pos <- which(FinalInt$CompoundName==x)
		y <- FinalInt[pos, 1:2,drop=F]
		resmat <- apply(FinalInt[pos, selcol2,drop=F], 2, function(y) {
			if(val.trans=="P"){res <- sum(y,na.rm=T)/100}
			if(val.trans=="R"){res  <- y[1] }
			res <- (y/res)
			return(res)
		})
			if(all(is.na(resmat))){
			warning(call.=F,paste(c(x," has no ",val.to.use," values.")))
		}
		# resmat <- cbind(FinalInt[pos,1:2],resmat)
		resmat <- data.frame(resmat)
		if(nrow(resmat)>nrow(y)){resmat <- t(resmat)}
		colnames(resmat) <- colnames(FinalInt)[selcol2]
		resmat <- cbind(y, resmat)
		return(resmat)
	})
	Transapply <- do.call("rbind",Transapply)
	Transapply<- as.data.frame(Transapply)
	rownames(Transapply) <- NULL
	return(Transapply)
}