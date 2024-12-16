#' Corrects natural abundance
#' 
#' Corrects natural abundance given the abundance predicted within autoQ, using only M+0 isotopologue.
#' This abundance is previously calculated by enviPat functions. 
#' @param autoQres The data.frame result from \link{autoQ}
#' @return A data frame with a row per each isotopologue and compound with corrected abundances per sample.
#' @export
#' 

subNaturalM0 <- function(autoQres){
	ppmi <- grep("ppm",colnames(autoQres))
	newres <- lapply(unique(autoQres$CompoundName),function(x){
		pos <- which(autoQres$CompoundName==x)
		resi <- autoQres[pos,]
		abu <- resi$abundance
		abu <- (abu/100)
		abu[1] <- 0
		valcols <- 5:ncol(resi)
		valcols <- valcols[-which(valcols%in%ppmi)]
		for(j in valcols){
			# y <- resi[,j]
			# y <- y-(y[1]*abu)
			# y[y<0] <- 0
			# resi[,j] <- y
			y <- resi[, j]
			y[is.na(y)] <- 0
			y <- (y - (y[1] * abu))
			y[y < 0] <- 0
			resi[, j] <- y
		}
		return(resi)
	})
	newres <- do.call("rbind",newres)
	
	return(newres) 
}


#' Corrects natural abundance
#' 
#' Corrects natural abundance given the abundance predicted within autoQ, using all quantified isotopologues.
#' This abundance is previously calculated by enviPat functions. 
#' @param autoQres The data.frame result from \link{autoQ}
#' @return A data frame with a row per each isotopologue and compound with corrected abundances per sample.
#' @export
#' 

subNaturalAll <- function(autoQres){
	ppmi <- grep("ppm",colnames(autoQres))
	newres <- lapply(unique(autoQres$CompoundName),function(x){
		pos <- which(autoQres$CompoundName==x)
		resi <- autoQres[pos,]
		isos <- resi$Isotopologue
		isos <- gsub("M.","",isos)
		nisos <- as.numeric(isos)
		abu <- resi$abundance/100
		valcols <- 5:ncol(resi)
		valcols <- valcols[-which(valcols%in%ppmi)]
		for(j in valcols){
			y <- resi[,j]
			for(n in unique(isos)){
				idxi <- which(isos==n)
				maxi <- which.max(abu[idxi])
				topiso <- idxi[maxi]
				excl <- unique(c(topiso,which(nisos<n)))
				newy <- rep(0,length(isos))
				if(length(excl)<length(y)){
					abuy <- abu/abu[topiso]
					newy[-excl] <- (y[topiso]*abuy[-excl])
					newy[newy<0] <- 0
				}
				if(any(is.na(y))){ newy[is.na(y)] <- NA}
				y <- y-newy
				y[y<0] <- 0
			}
			resi[,j] <- y
		}
		return(resi)
	})
	newres <- do.call("rbind",newres)
	return(newres)
}
