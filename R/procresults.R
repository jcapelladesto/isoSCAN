rmOverlap <- function(isotab){
  maxi <- ncol(isotab)-1
  maxv <- isotab[,maxi]
  if(any(duplicated(maxv))){
    i <- 1
    while(i<nrow(isotab)){
      maxv <- isotab[,maxi]
      idx <- which(maxv==maxv[i])
      if(length(idx)>1){
        abunds <- isotab[idx,2]
        if(length(idx)>1){
          idxrm <- which(abunds!=max(abunds))
          isotab <- isotab[-idx[idxrm],]
          idx <- idx[-idxrm]
        }
        if(length(idx)==1){
          i <- 1
        }
      }
      i <- i+1
    }
  }
  return(isotab)
}

checkRows <- function(ppm_df,ppm_values){
  rws <- sapply(ppm_df,nrow)
  ppm_values <- rep(ppm_values,2)
  
  if(any(rws==0)){
    rws0 <- which(rws==0)
    ppm_df <- ppm_df[-rws0]
    rws <- rws[-rws0]
    ppm_values <- ppm_values[-rws0]
  }
  
  if(length(ppm_df)==0){
    return(data.frame())
  }
  
  if(length(ppm_df)==1){
    return(ppm_df[[1]]) 
  }
  
  a <- 2:length(ppm_df)
  b <- 1:(length(ppm_df)-1)
  rwsdiff <- c(1,rws[a]/rws[b])
  
  nas <- sapply(ppm_df,function(x) sum(is.na(x)))
  
  maxosum <- sapply(ppm_df,function(x) sum(x[,grep("maxo",colnames(x))],na.rm=T) )
  mdiff <- c(1,maxosum[a]/maxosum[b])
  torm <- is.nan(mdiff)
  if(any(torm)){
    mdiff[which(torm!=0)] <- 1
  }
  torm <- is.na(mdiff)+is.infinite(mdiff)+(mdiff==0)+is.nan(mdiff)
  torm[is.na(torm)] <- 1
  suppressWarnings(if(any(torm)){
    mdiff[which(torm!=0)] <- 1
  })
  
  areasum <- sapply(ppm_df,function(x) sum(x[,grep("area",colnames(x))],na.rm=T) )
  adiff <- c(1,areasum[a]/areasum[b])
  torm <- is.nan(adiff)
  if(any(torm)){
    adiff[which(torm!=0)] <- 1
  }
  torm <- is.na(adiff)+is.infinite(adiff)+(adiff==0)
  if(any(torm)){
    adiff[which(torm!=0)] <- 1
  }
  
  ppmerr <- sapply(ppm_df,function(x) sum(x$ppm,na.rm=T)) 
  ppmerr <- (ppmerr/ppm_values)/rws # added rws factor
  
  abundiff <- (adiff/mdiff)
  
  m0found <- sapply(ppm_df,function(x) any(x$abundance==100))
  
  dfidx <- which(
    ppmerr<quantile(ppmerr,0.85)  &
      rws>=quantile(rws,0.25)   &
      abundiff>=quantile(abundiff,0.25)  ) # was 0.33
  
  if(length(dfidx)>1 & any(m0found)){
    m0found <- m0found[dfidx]
    dfidx <- dfidx[which(m0found)]
  }
  
  if(length(dfidx)>1){
    nas <- nas[dfidx]
    dfidx <- dfidx[which(nas==min(nas))]
  }
  if(length(dfidx)>1){
    abundiff <- (adiff/mdiff)[dfidx]
    dfidx <- dfidx[which(abundiff==min(abundiff))]
  }
  if(length(dfidx)>1){
    areasum<- areasum[dfidx]
    dfidx <- dfidx[which(areasum==max(areasum))]
  }
  if(length(dfidx)>1){
    ppmerr <- ppm_values[dfidx/((dfidx%/%length(ppm_values)+1))]
    dfidx <- dfidx[which(ppmerr==min(ppmerr))][1]
  }
  
  if(length(dfidx)==0){
    return(data.frame())
  }else{
    dfres <- ppm_df[[dfidx]]
    natest <- all(is.na(dfres[,-(1:4)]))
    if(natest){
      return(data.frame())
    }else{
      return(dfres)        
    }
  }
}

createNAdf <- function(res,HR){
  
  if(HR){
    temp <- lapply(res,function(r) { if(nrow(r)>0) return(r[,1:3]) else return() })
    temp <- do.call("rbind",temp)
    ir <- data.frame("m.z" =  sort(unique(temp$m.z)) )
    ir$abundance <- sapply(ir$m.z, function(iso) mean(temp$abundance[which(temp$m.z==iso)] ))
    ir$Isotopologue <-
      sapply(ir$m.z, function(iso) names(sort(table(temp$Isotopologue[which(temp$m.z==iso)]))[1]) )
    ir$ppm <- NA 
  }else{
    temp <- lapply(res,function(r) { if(nrow(r)>0) return(r[,1:2]) else return() })
    temp <- do.call("rbind",temp)
    ir <- data.frame("m.z" =  sort(unique(temp$m.z)) )
    ir$Isotopologue <-
      sapply(ir$m.z, function(iso) names(sort(table(temp$Isotopologue[which(temp$m.z==iso)]))[1]) )
  }
  
  ir$maxo <- NA
  ir$area <- NA
  
  abu <- ir[,3]
  if(any(duplicated(abu))){
    i <- 1
    while(i<nrow(ir)){
      abu <- ir[,3]
      iso <- ir[,2]
      idx <- which(abu==abu[i] & iso==iso[i])
      if(length(idx)>1){
        ir <- ir[-idx[-1], ]
      }
      i <- i+1
    }
  }
  rownames(ir) <- NULL
  return(ir)
}

prepRes <- function(res,ir,rows,SampleFiles,HR){
  
  mr <- nrow(ir)
  idx <- which(rows==0)
  if(length(idx)>0){
    for(j in idx){
      sampname <- rmfileExt(SampleFiles[j],"\\.mz.*")
      temp <- ir
      colnames(temp)[(ncol(temp)-1):ncol(temp)] <- paste(
        colnames(temp)[(ncol(temp)-1):ncol(temp)],sampname,sep="_")
      res[[j]] <- temp
    }
  }
  
  for(r in 1:length(res)){
    temp1 <- res[[r]]
    if(rows[r]!=mr){
      # temp2 <- ir[which(!(ir$Isotopologue%in%temp1$Isotopologue)),,drop=F]
      temp2 <- ir
      colnames(temp2) <- colnames(temp1)
      temp2 <- temp2[-sapply(temp1$m.z,function(x) which.min(abs(temp2$m.z-x))),]
      temp1 <- rbind(temp1,temp2)
      temp1 <- temp1[order(as.character(temp1$m.z)),,drop=F]
    }
    if(r!=1){
      if(HR){
        temp1 <- temp1[,-c(1:3),drop=F]
      }else{
        temp1 <- temp1[,-c(1:2),drop=F]
      }
    }
    res[[r]] <- temp1
  }
  return(res)
}