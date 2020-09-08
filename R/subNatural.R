#' Corrects natural abundance
#' 
#' Corrects natural abundance given the abundance predicted within autoQ. This abundance is calculated by enviPat functions. 
#' @param autoQres The data.frame result from \link{autoQ}
#' @return A data frame with a row per each isotopologue and compound with corrected abundances per sample.
#' @export
#' 
#' 
#' https://bioconductor.org/packages/devel/bioc/vignettes/IsoCorrectoRGUI/inst/doc/IsoCorrectoRGUI.html#why-perform-correction-for-natural-stable-isotope-abundance-and-tracer-purity

subNatural <- function(autoQres){
  ppmi <- grep("ppm",colnames(autoQres))
  newres <- lapply(levels(autoQres$CompoundName),function(x){
    pos <- which(autoQres$CompoundName==x)
    resi <- autoQres[pos,]
    abu <- resi$abundance
    abu <- (abu/100)
    abu[1] <- 0
    # sumabu <- abu
    # labu <- length(sumabu)
    # if(labu>3){
    #   for(i in 3:labu){
    #     sumabu[i] <- sumabu[i]+prod(abu[3:i])
    #   }
    # }
    # abu <- sumabu
    valcols <- 5:ncol(resi)
    valcols <- valcols[-which(valcols%in%ppmi)]
    for(j in valcols){
      y <- resi[,j]
      y <- y-(y[1]*abu)
      y[y<0] <- 0
      resi[,j] <- y
    }
    return(resi)
  })
  newres <- do.call("rbind",newres)
  
  return(newres) 
}



