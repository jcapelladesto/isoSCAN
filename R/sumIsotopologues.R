#' Sums equivalent isotopologues abundances
#' 
#' Sums abundances of ions equivalent to each isotopologue (M+x). Necessary to execute metBarPlot function 
#' after running autoQ function on high-resolution data.
#' @param autoQres The data.frame result from \link{autoQ}
#' @return A data frame with a row per each isotopologue and compound with summed abundances per sample.
#' @export

sumIsotopologues <- function(autoQres){
  ppmi <- grep("ppm",colnames(autoQres))
  if(length(ppmi)!=0){
  	autoQres <- autoQres[,-ppmi]
  	warning("ppm values were removed for Isotopologue sum, make sure you did perform natural abundance substraction if necessary")
  }
  autoQres$m.z <- NULL
  autoQres$abundance <- NULL
  newres <- lapply(unique(autoQres$CompoundName),function(x){
    pos <- which(autoQres$CompoundName==x)
    resi <- autoQres[pos,]
    temp <- resi[1,]
    isoi <- lapply(unique(resi$Isotopologue),function(i){
      isoidx <- which(resi$Isotopologue==i)
      temp$Isotopologue <- i
      for(j in 3:ncol(resi)){
        y <- sum(resi[isoidx,j],na.rm=T)
        temp[,j] <- y
      }
      return(temp)
    })
    isoi <- do.call("rbind",isoi)
    if(any(isoi==0)){isoi[isoi==0] <- NA}
    return(isoi)
  })
  newres <- do.call("rbind",newres)

 return(newres) 
}
