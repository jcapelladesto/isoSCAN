#' Calculate Fractional Contribution (FrC)
#' 
#' Fractional contribution quantifies the contribution of a labeled nutrient to the metabolite of interest.
#' @param finalres Dataframe resulting from \emph{autoQ} function
#' @param val.to.use Type of quantification to transform. Either "Area" or "Maxo"
#' @param formulaTable Dataframe including CompoundName, Formula, mz, RT and NumAtoms (used in  \link{autoQ})
#' @return A data.frame containing a FrC value for every sample and compound.
#' @references doi: 10.1016/j.copbio.2015.02.003
#' @export
#' 

calcFrC <- function(FinalInt, val.to.use, formulaTable){
  frcres <- lapply(unique(FinalInt$CompoundName), function(y) {
    selcol2 <- grep(val.to.use, colnames(FinalInt))
    x <- finalres[which(finalres$CompoundName==y),]
    natoms <- x$Isotopologue
    natoms <- as.numeric(gsub("M.","",natoms))
    res <- sapply(selcol2, function(j){
      vals <- x[,j]
      ss <- sum(sapply(1:length(natoms), function(i) natoms[i]*vals[i]),na.rm=T)
      total <- formulaTable$NumAtoms[which(formulaTable$CompoundName==y)]
      total <- total*sum(vals,na.rm=T)
      ss/total
    })
    res <- t(data.frame(res))
    colnms <- colnames(FinalInt)[selcol2]
    colnms <- gsub(val.to.use,"FrC",colnms)
    CompoundName <- as.character(x$CompoundName[1])
    suppressWarnings(res <- data.frame(CompoundName,res))
    colnames(res)[2:ncol(res)] <- colnms
    rownames(res) <- NULL
    return(res)
  })
  frcres <- do.call("rbind",frcres)  
}