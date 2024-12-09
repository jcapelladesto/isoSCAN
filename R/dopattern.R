sumIsoLabel <- function(pattern.list,labelatom){ # sums columns 13C 
  
  for(i in 1:length(pattern.list)){
    temp <- pattern.list[[i]]
    tidx <- grep(labelatom,colnames(temp))
    tidx2 <- tidx
    if(i!=1){
      if(length(tidx)==1){
      	tempe <- pattern.list[[1]]
      	nname <- setdiff(colnames(tempe), colnames(temp))
      	temp <- cbind(temp, 0)
      	colnames(temp)[ncol(temp)] <- nname
      	temp <- temp[, colnames(tempe), drop = F]
      	tidx <- grep(labelatom, colnames(temp))
      	tidx2 <- tidx
      }else{
      	tempi <- temp[, tidx[1], drop=F]
      	temp[, tidx[2]] <- temp[, tidx[2]] + tempi
      	temp <- temp[, (-1) * (tidx[1]), drop=F]
      	tidx2 <- (tidx[length(tidx)] - 1)
      	abuprev <- pattern.list[[i - 1]][, 1:2, drop=F]
      	mzdiff <- abs(temp[1, 1] - abuprev[, 1])
      	abuprev <- as.numeric(abuprev[which(mzdiff == 
      																				min(mzdiff)), 2])/100
      	temp[, 2] <- temp[, 2] * abuprev
      }
    	temp <- temp[which(temp[, tidx2] %in% c(i - 1, i)),, drop = F]
    }
    # temp <- temp[which(temp[,tidx2]%in%c(i-1,i)),] #,i+1
    pattern.list[[i]] <- temp
  }
  
  return(pattern.list)
}

removeDuplicates <- function(isopat.matrix){
  
  temp <- sapply(1:nrow(isopat.matrix),function(i){
    x <- isopat.matrix[i,(grep("abundance",colnames(isopat.matrix))+1):ncol(isopat.matrix)]
    paste(x,collapse="")
  })
  
  if(any(duplicated(temp))){
    i <- 1
    while(i<nrow(isopat.matrix)){
      idx <- which(temp==temp[i])
      if(length(idx)>1){
        idx1 <- idx[1]
        idx2 <- idx[2:length(idx)]
        isopat.matrix[idx1,2] <- max(isopat.matrix[idx,2]) # this was out before
        # isopat.matrix[idx1,2] <- sum(isopat.matrix[idx,2]) 
        isopat.matrix[idx1,1] <- mean(isopat.matrix[idx,1])
        isopat.matrix <- isopat.matrix[-idx2,]
        temp <- temp[-idx2]
      }
      i <- i+1
    }
  }
  
  temp <- isopat.matrix[,2]
  # isopat.matrix[,2] <- 10 # improve detection of valleys in envelopes for low nat abund
  # isopat.matrix[which(temp>0.1),2] <- 25 # improve detection of valleys in envelopes for low nat abund
  # isopat.matrix[which(temp>1),2] <- 50 # improve detection of valleys in envelopes for low nat abund
  # isopat.matrix[which(temp>5),2] <- 100 # improve detection of valleys in envelopes for low nat abund
  isopat.matrix <- isopat.matrix[order(isopat.matrix[,1]),]
  return(isopat.matrix)
}

getTargetPat <- function(targetisos,isotopes,maxmz,labelatom,thr=thr){
  targetpat <- enviPat::isopattern(isotopes=isotopes,
                                   chemforms = targetisos,
                                   charge = 1,verbose = F,threshold = thr)
  
  for (i in 1:length(targetpat)) {
  	temp <- targetpat[[i]]
  	temp <- temp[which(temp[, 1] < maxmz), ,drop=F]
  	targetpat[[i]] <- temp
  }
  rm(temp)
  
  targetpat <- sumIsoLabel(targetpat,labelatom) 
  
  return(targetpat)
}

doTheoPat <- function(targetpat,labelatom){
  targetmz <- data.frame(do.call("rbind",targetpat))
  targetmz <- targetmz[order(targetmz$m.z),]
  
  targetmz <- removeDuplicates(targetmz)
  
  lbli <- grep(labelatom,colnames(targetmz))
  lbli <- targetmz[,lbli]
  targetmz <- targetmz[,-c(3:ncol(targetmz))]
  targetmz$Isotopologue <- lbli
  
  targetmz$ppm <- sapply(1:length(targetmz$m.z),function(z){
    m <- targetmz$m.z[z]
    diff <- ppmDiff(targetmz$m.z[-z],m)  
    diff[which(diff==min(diff))]*0.90
  })
  
  return(targetmz)
}

doConvPat <- function(targetpat,Res,labelatom){
  
  targetall <- do.call("rbind",targetpat)
  targetall <- removeDuplicates(targetall)
  
  envires <-  enviPat::envelope(
    list(targetall),
    ppm = F, dmz = "get", frac = 1/4, env =  "Gaussian",
    resolution = Res, plotit =F, verbose = F)
  
  envicent <-  enviPat::vdetect(
    envires,
    detect="centroid",plotit=F,
    verbose=F)
  
  targetmz <- data.frame(envicent[[1]])
  targetall <- data.frame(targetall)
  
  lbli <- grep(labelatom,colnames(targetall))
  
  ppmi <- sapply(1:nrow(targetmz),function(z){
    m <- targetmz$m.z[z]
    diff <- ppmDiff(targetmz$m.z[-z],m)
    diff[which(diff==min(diff))]*0.90
  })
  
  check_isos <- floor(targetmz$m.z)
  check_isos <- (all(table(check_isos)==1) + all(diff(check_isos)==1))==2
  
  if (check_isos) {
  	targetmz$Isotopologue <- seq(0, nrow(targetmz) - 1)
  } else {
  	targetisos <- lapply(1:nrow(targetmz), function(z) {
  		target <- targetmz[z, ]
  		ppmz <- ppmi[z]
  		z <- target$m.z
  		ppmd <- z * ppmz/1e+06
  		mzdiff <- abs(targetall$m.z - z)
  		res <- targetall[which(ppmd > mzdiff), ]
  		if (nrow(res) > 1) {
  			res <- res[which(res$abundance >= (max(res$abundance)*0.95)),  ]
  			if(nrow(res)>1){
  				res <- c(sum(res[,2]),max(res[,lbli]))
  			}else{
  				res <- res[1,lbli]
  			}
  			
  		} else {
  			res <- res[1, lbli]
  		}
  		return(res)
  	})
  	targetmz$Isotopologue <- 0
  	for(j in seq(length(targetisos))){
  		jj <- targetisos[[j]]
  		if(length(jj)==1){
  			targetmz$Isotopologue[j] <- jj
  		}else{
  			targetmz$Isotopologue[j] <- jj[2]
  			targetmz$abundance[j] <- jj[1]
  		}
  	}
  }
  
  targetmz$ppm <- ppmi
  return(targetmz)
}

