#' Simplify data.frame of quantified isotopologues 
#' 
#' Simplify data.frame of quantified isotopologues by clustering isotopologue ions within a certain ppm mass difference. 
#' These differences are caused by variation in recorded ion m/z through that may not match from sample to sample.
#' @param finalres Dataframe resulting from \emph{autoQ} function from high-resolution data.
#' @param ppm Relative mass tolerance to cluster isotopologue masses
#' @return A data.frame with a reduced number of rows in which similar mass isotopologues have been collapsed into a single row.
#' @export
simplifybyPpm <- function(finalres,ppm=NULL){
  p <- grep("ppm",colnames(finalres[[1]]))
  finalres <- lapply(finalres,function(res){
    if(!is.null(res)){
      i <- 1
      while(i<nrow(res)){
        z <- res$m.z[i]
        if(is.null(ppm)){
          ppmz <- as.numeric(res[i,p])
          ppmz <- max(ppmz,na.rm=T)*1.05
        }else{
          ppmz <- ppm
        }
        ppmd <- z*ppmz/1e6
        mzdiff <- abs(res$m.z-z)
        idx <- which(ppmd>mzdiff)
        idx <- idx[-which(idx==i)]
        if(length(idx)!=0 ){
          dfa <- res[i,]
          for (j in 1:length(idx)){
            dfrepair <- res[idx[j],]
            dfa[] <-Map(function(x,y) {x[is.na(x)] <- y[is.na(x)]; x}, dfa, dfrepair)
          }
          dfa$m.z <- mean(res$m.z[c(i,idx)])
          dfa$abundance <- mean(res$abundance[c(i,idx)])
          dfa$Isotopologue <- max(res$Isotopologue[c(i,idx)])
          res <- res[-idx,]
          res[i,] <- dfa
          i <- 1 #reset?
        }
        i <- i+1
      }
    }
    return(res)
  })
  return(finalres)
}

simplifybyIso <- function(finalres){
  p <- grep("ppm",colnames(finalres[[1]]))
  finalres <- lapply(finalres,function(res){
    if(!is.null(res)){
      i <- 1
      while(i<nrow(res)){
        z <- res$m.z[i]
        ppmz <- as.numeric(res[i,p])
        ppmz <- max(ppmz,na.rm=T)*1.25
        ppmd <- z*ppmz/1e6
        mzdiff <- abs(res$m.z-z)
        idx <- which(ppmd>mzdiff)
        idx <- idx[-which(idx==i)]
        if(length(idx)!=0 & all(res$Isotopologue[idx]==res$Isotopologue[i])){
          dfa <- res[i,]
          for (j in 1:length(idx)){
            dfrepair <- res[idx[j],]
            dfa[] <-Map(function(x,y) {x[is.na(x)] <- y[is.na(x)]; x}, dfa, dfrepair)
          }
          dfa$m.z <- mean(res$m.z[c(i,idx)])
          dfa$abundance <- mean(res$abundance[c(i,idx)])
          res <- res[-idx,]
          res[i,] <- dfa
          i <- 1 #reset?
        }
        i <- i+1
      }
    }
    return(res)
  })
  return(finalres)
}