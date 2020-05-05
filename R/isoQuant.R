extractEIC <- function(myfile,h,myscans,mzmin,mzmax){
  aa <- mzR::openMSfile(myfile)
  p <- mzR::peaks(aa,scans = myscans )
  eic <- lapply(1:length(p),function(x){
    # rt <- h[which(h$acquisitionNum==myscans[x]),"retentionTime"]
    rt <- h[myscans[x],"retentionTime"]
    mp <- p[[x]]
    mpidx <- which(mp[,1]>mzmin &  mp[,1]<mzmax)
    mp <- mp[mpidx,,drop=F]
    res <- cbind(mp,rep(x=rt,times=nrow(mp)),rep(x=x,times=nrow(mp)))
    return(res)
  })
  eic <- do.call("rbind",eic)
  eic <- as.data.frame(eic)
  colnames(eic) <- c("basePeakMZ","basePeakIntensity","retentionTime","ScanNum")
  return(eic)
}

searchEIC.HR <- function(targetmz,eic,ppm){
  patppmi <- lapply(targetmz$m.z,function(x){ ppmDiff(eic$basePeakMZ,x) })
  patfoundi <- lapply(1:length(patppmi),function(x){ 
    ppmi <- targetmz$ppm[x]
    x <- patppmi[[x]]
    if(ppmi<ppm){ppm <- ppmi}
    idx <-  which(x<ppm)
    return(idx)
  })
  patppmi <- sapply(1:length(patppmi),function(x){ 
    ppmd <-  patppmi[[x]]
    idx <- patfoundi[[x]]
    median(ppmd[idx])
  })
  if(any(is.na(patppmi ))){
    patfoundi <- patfoundi[-which(is.na(patppmi))]
    targetmz <- targetmz[-which(is.na(patppmi)),]
    patppmi <- patppmi[-which(is.na(patppmi))]
  }
  targetmz$ppm <- patppmi
  
  return(list("A"=targetmz,"B"=patfoundi))
}


searchEIC.LR <- function(targetmz,eic,mzerror){
  patppmi <- lapply(targetmz$m.z,function(x){ abs(eic$basePeakMZ-x) })
  patfoundi <- lapply(patppmi,function(x){ 
    idx <- which(x<mzerror)
    return(idx)
  })
  patppmi <- sapply(patfoundi,length)
  
  if(any(patppmi==0)){
    patfoundi <- patfoundi[-which(patppmi==0)]
    targetmz <- targetmz[-which(patppmi==0),]
  }
  
  return(list("A"=targetmz,"B"=patfoundi))
}

procEIC <- function(fTi,targetmz, patfoundi, eic, minwidth, maxwidth, minscans, SNR, fit.p, s){
  
  if(nrow(targetmz)>0){
    isotab <- findIsos(fTi,targetmz, patfoundi, eic, minwidth, maxwidth, minscans, SNR, fit.p, s) 
  }else{
    return(data.frame())
  }
  if(is.null(isotab)){
    return(data.frame())
  }
  
  sampname <- rmfileExt(SampleFiles[s],"\\.mz.*")
  
  colnames(isotab)[(ncol(isotab)-1):ncol(isotab)] <- paste(
    colnames(isotab)[(ncol(isotab)-1):ncol(isotab)],sampname,sep="_")
  colnames(isotab)[1] <- "m.z"
  
  isotab <- rmOverlap(isotab)
  
  return(isotab)
}


isoQuant.HR <- function(SampleFiles, formulaTable, SNR, minscans , RTwin, fit.p,
                        maxppm, minwidth, maxwidth, isotopes, labelatom,
                        resolution, HR){
  
  # preload sample headers
  hdlist <- lapply(1:length(SampleFiles),function(s){
    oMS <- mzR::openMSfile(SampleFiles[s])
    h <- mzR::header(oMS)
    h <- h[,c("retentionTime","acquisitionNum")]
    cbind(h,sample=s)
  })
  
  finalres <-lapply(1:nrow(formulaTable), function(i){
    
    fTi <- formulaTable[i,]
    message(paste0(c("Processing: ",as.character(fTi$CompoundName))))
    fTiMZ  <-  fTi$mz
    fTiRT <- fTi$RT
    fTiRTRan <- c(fTiRT-RTwin,fTiRT+RTwin)
    fTiFormula <- as.character(fTi$Formula)
    natom <- fTi$NumAtoms
    
    atom <- isotopes$element[which(isotopes$isotope==labelatom)]
    atom <- atom[-grep("\\[",atom)]
    isomass <- isotopes$mass[which(isotopes$isotope==labelatom & isotopes$element==atom)]
    atmass <- isotopes$mass[which(isotopes$element==atom & isotopes$isotope!=labelatom)]
    
    massdiff <- isomass-atmass
    mzmin <- fTiMZ-1
    mzmax <- ceiling(fTiMZ+((natom+1)*massdiff))
    
    if(length(resolution)==1){
      # fixed Res
      Res <- as.numeric(resolution)
    }else{
      # variable Res
      Res <- calc.Res(resolution[1],resolution[2],fTiMZ)
    }
    
    targetisos <- iso.generator(fTiFormula,labelatom,1:natom,isotopes)
    
    targetpat <- getTargetPat(targetisos,isotopes,mzmax,labelatom)
    
    targetmz.theo <- doTheoPat(targetpat,labelatom)
    
    targetmz.conv <- doConvPat(targetpat,Res,labelatom)
    
    ppm_values <- seq(1,maxppm,0.5)
    
    res <- lapply(1:length(SampleFiles),function(s){
      
      h <- hdlist[[s]]
      # myscans <- h$acquisitionNum[which(h$retentionTime>fTiRTRan[1] & h$retentionTime<fTiRTRan[2])]
      myscans <- which(h$retentionTime>fTiRTRan[1] & h$retentionTime<fTiRTRan[2])
      eic <- extractEIC(SampleFiles[s],h,myscans,mzmin,mzmax)
      
      ppm_df <- lapply(list(targetmz.conv,targetmz.theo),function(x){
        lapply(ppm_values, function(ppm){
          eicL <- searchEIC.HR(x,eic,ppm)
          targetmz <- eicL$A
          patfoundi <- eicL$B
          procEIC(fTi,targetmz, patfoundi, eic, minwidth, maxwidth, minscans, SNR, fit.p, s)
          
        })
      })
      ppm_df <- do.call("c",ppm_df)
      
      r <- checkRows(ppm_df,ppm_values)
      rownames(r) <- NULL
      return(r)
    })
    
    rows <- sapply(res,nrow)
    
    if(all(rows==0)){
      message("Compound lost")
      return()
    }
    ir <- createNAdf(res,HR)
    
    res <- prepRes(res,ir,rows,SampleFiles,HR)
    res <- do.call("cbind", res)
    
    torm <- which(apply(res[,grep("maxo|area",colnames(res))],1,function(x) all(is.na(x))))
    
    if(length(torm)>0){
      res <- res[-torm,]
    }
    
    if(nrow(res)==0){ # remove?
      return()
    }
    
    res <- cbind("CompoundName"=fTi$CompoundName,res)
    return(res)
  })
  
  finalres <- simplifybyIso(finalres)
  
  finalres <- simplifybyPpm(finalres)
  
  finalres <- do.call("rbind", finalres)
  
  rownames(finalres) <- NULL
  return(finalres)
}


isoQuant.LR <- function(SampleFiles, formulaTable, SNR, minscans , RTwin, fit.p,
                        mzerror, massdiff, minwidth, maxwidth, HR){
  
  # preload sample headers
  hdlist <- lapply(1:length(SampleFiles),function(s){
    oMS <- mzR::openMSfile(SampleFiles[s])
    h <- mzR::header(oMS)
    h <- h[,c("retentionTime","acquisitionNum")]
    cbind(h,sample=s)
  })
  
  finalres <-lapply(1:nrow(formulaTable), function(i){
    
    fTi <- formulaTable[i,]
    message(paste0(c("Processing: ",as.character(fTi$CompoundName))))
    fTiMZ  <-  fTi$mz
    fTiRT <- fTi$RT
    fTiRTRan <- c(fTiRT-RTwin,fTiRT+RTwin)
    natom <- fTi$NumAtoms
    targetisos <- c(0:natom)
    targetmz <- fTiMZ+(targetisos*massdiff)
    
    targetmz <- data.frame("m.z"=targetmz,"Isotopologue"=targetisos)
    
    mzmin <- floor(fTiMZ-1)
    mzmax <- ceiling(max(targetmz))+1
    
    res <- lapply(1:length(SampleFiles),function(s){
      h <- hdlist[[s]]
      # myscans <- h$acquisitionNum[which(h$retentionTime>fTiRTRan[1] & h$retentionTime<fTiRTRan[2])]
      # myscans <- h$acquisitionNum[which(h$retentionTime>fTiRTRan[1] & h$retentionTime<fTiRTRan[2])]
      myscans <- which(h$retentionTime>fTiRTRan[1] & h$retentionTime<fTiRTRan[2])
      eic <- extractEIC(SampleFiles[s],h,myscans,mzmin,mzmax)
      # Isolabel
      eicL <- searchEIC.LR(targetmz,eic,mzerror)
      targetmz <- eicL$A
      patfoundi <- eicL$B
      # sapply(patfoundi,function(i) plot(eic[i,3:2]))
      
      r <- procEIC(fTi,targetmz, patfoundi, eic, minwidth, maxwidth, minscans, SNR, fit.p, s)
      
      return(r)
    })
    
    rows <- sapply(res,nrow)
    
    if(all(rows==0)){
      message("Compound lost")
      return()
    }
    ir <- createNAdf(res,HR)
    
    res <- prepRes(res,ir,rows,SampleFiles,HR)
    res <- do.call("cbind", res)
    
    torm <- which(apply(res[,grep("maxo|area",colnames(res))],1,function(x) all(is.na(x))))
    
    if(length(torm)>0){
      res <- res[-torm,]
    }
    
    if(nrow(res)==0){ # remove?
      return()
    }
    
    res <- cbind("CompoundName"=fTi$CompoundName,res)
    return(res)
  })
  
  finalres <- do.call("rbind", finalres)
  rownames(finalres) <- NULL
  return(finalres)
}