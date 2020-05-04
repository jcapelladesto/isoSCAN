calc.Res <- function(R,m,mz){ R*sqrt(m/mz)}

ppmDiff <- function(a,b){ abs(a-b)*1e6/b }

rmfileExt <- function(x,ext.regexp){
  x <- basename(x)
  x <- gsub(ext.regexp,"",x)
  return(x)
}

decompose.formula <- function(formula){
  rgm <- gregexpr("[A-Z]",formula)[[1]]
  rgm <- as.numeric(rgm)
  r <- sapply(1:length(rgm),function(i){
    j <- rgm[i]
    if(rgm[i+1]%in%rgm){
      ele <-  substring(formula,j,j)
      num <- substring(formula,j+1,rgm[i+1]-1)
      if(is.na(num)){num <- substring(formula,j+1,nchar(formula))}
    }else{
      if(grepl("[0-9]",substring(formula,j+1,j+1))){
        ele <- substring(formula,j,j)
        num <- substring(formula,j+1,rgm[i+1]-1)
        if(is.na(num)){num <- substring(formula,j+1,nchar(formula))}
      }else{
        ele <- substring(formula,j,j+1)
        num <- substring(formula,j+2,rgm[i+1]-1)
        if(is.na(num)){num <- substring(formula,j+2,nchar(formula))}
      }
      
    }
    if(num==""){num <- 1}
    return(c(ele,num))
  })
  Element <- r[1,]
  Number <- r[2,]
  r <- data.frame(Element,Number,stringsAsFactors = F)  
  # colnames(r) <- c("Element","Number")
  rownames(r) <- NULL
  return(r)
}

iso.generator <- function(formula,atom,nisos,isotopes){
  
  isos <- isotopes[which(isotopes$isotope==atom),]
  isos <- isos[grep("\\[",isos$element),]
  
  isoel <- gsub("[0-9]","",atom)
  
  decform <- decompose.formula(formula)
  isoidx <- which(decform$Element==isoel)
  isoeln <- as.numeric(decform$Number[isoidx])
  decform <- decform[-isoidx,]
  
  newiso <- isoeln-nisos
  newiso <- sapply(newiso,function(x) { 
    if(x>0){
      paste0(isoel,x) 
    }else{
      ""
    }
  })
  
  varvec <- paste0(paste0(isos$element,nisos),newiso)
  fixed <- paste(sapply(1:nrow(decform),function(x) paste0(decform[x,],collapse="")),collapse="")
  
  res <- paste0(varvec,fixed)
  
  res <- c(formula,res)
  return(res)
}
