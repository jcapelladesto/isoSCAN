find.LocalMax <- function(rt,int, minwidth, maxwidth,minscans){
  res <- NA
  mxw2 <- maxwidth/2
  miw2 <- minwidth/2
  ### find max
  # maxScan <- which(int==max(int))[1]
  # maxScanRT <- rt[maxScan]
  # maxScanI <- int[maxScan]
  # 
  # ### find 1st min within X scans/secs?
  # # LEFT 
  # idx_l <- which(rt>(maxScanRT-mxw2) & rt<(maxScanRT-miw2) )
  # 
  # ### find 2nd min within X scans/secs?
  # # RIGHT
  # idx_r <- which(rt<(maxScanRT+(mxw2)) & rt>(maxScanRT+(miw2)) )

  validMax <- T
  while(validMax){
    maxScan <- which(int==max(int))[1]
    maxScanRT <- rt[maxScan]
    maxScanI <- int[maxScan]

    idx_l <- which(rt>(maxScanRT-mxw2) & rt<(maxScanRT-miw2) )

    idx_r <- which(rt<(maxScanRT+(mxw2)) & rt>(maxScanRT+(miw2)) )

    if(length(idx_l)==0 | length(idx_r)==0){
      rt <- rt[-maxScan]
      int <- int[-maxScan]
    }else{
      validMax <- F
    }
    if(length(rt)==0){
      validMax <- F
    }

  }
  
  # aa <- c(idx_l,maxScan,idx_r); plot(rt[aa],int[aa])
  
  if( (length(idx_l)>0) & (length(idx_r)>0) ){
    
    rt_l <- rt[idx_l]
    int_l <- int[idx_l]
    min_l <- which(int_l==min(int_l))[1]
    
    rt_r <- rt[idx_r]
    int_r <- int[idx_r]
    min_r <- which(int_r==min(int_r))[1]
    
    ## define noise using mins
    noise <- mean(c(int_r[c(min_r-1,min_r,min_r+1)],
                    int_l[c(min_l-1,min_l,min_l+1)]),na.rm=T)
    # maxScan
    pkscns <- c(idx_l[min_l]:idx_r[min_r])
    pk <- int[pkscns]
    
    ## return RTmaxo
    maxoRT <- maxScanRT
    ## calculate SNR
    SNR <- maxScanI/noise
    # calculate maxo
    maxo <- maxScanI
    ## define RTrange of peak 
    RTmin <- rt_l[min_l]
    RTmax <- rt_r[min_r]
    # calculate area (sum)
    area <- pk-noise
    area <- sum(area[which(area>0)])
    
    # report: peakid RTrange RTmaxo SNR maxo area
    if(length(pkscns)>=minscans){
      res <- list(pkscns,data.frame(RTmin, RTmax, maxoRT, SNR, maxo, area))
      
      # fit curve or linear model ??
    }else{
      res <- NA
    }
  }
  return(res)
}


# find.All.LocalMax <- function(rt,int,minwidth,maxwidth,minscans){
#   
#   res <- find.LocalMax(rt,int,minwidth,maxwidth,minscans)
#   
#   if(all(is.na(res))){
#     df <- NA
#     
#   }else{
#     
#     pkscns <- res[[1]]
#     df <- res[[2]]
#     
#     input <- data.frame(rt,int) 
#     
#     input <- list(input[1:pkscns[1],],
#                   input[pkscns[length(pkscns)]:nrow(input),])
#     
#     in_nrow <- sapply(input,nrow)
#     
#     new_check <- (in_nrow>minscans)
#     
#     doWhile <- any(new_check)
#     # if(doWhile){pkscns <- list(c(pkscns))}
#     
#     while(doWhile){ # while check?
#       # print("start While")
#       doWhile <- F
#       # print(doWhile)
#       
#       input <- input[which(new_check)]
#       
#       resList <- lapply(input,function(x) find.LocalMax(x$rt,x$int,minwidth,maxwidth,minscans))
#       
#       natest <- sapply(1:length(resList),function(r) is.na(resList[[r]][1] ))
#       
#       if(any(!natest)){
#         resList <- resList[which(!natest)]
#         pkscns <- list()
#         for(i in 1:length(resList)){
#           res <- resList[[i]]
#           
#           df <- rbind(df,res[[2]])
#           
#           pkscns <- c(pkscns,c(res[1]))
#         }
#         input <- lapply(1:length(pkscns),function(i){
#           p <- pkscns[[i]]
#           inp <- input[[i]]
#           list(inp[1:p[1],], inp[p[length(p)]:nrow(inp),])
#         })
#         input <- unlist(input,recursive =F)
#         
#         in_nrow <- sapply(input,nrow)
#         new_check <- (in_nrow>minscans)
#         doWhile <- any(new_check)
#       }
#       # print(doWhile)
#       # print("END")
#     }
#     
#     return(df)
#   }
# }
# 

find.All.LocalMax <- function(rt,int,minwidth,maxwidth,minscans){
  
  res <- find.LocalMax(rt,int,minwidth,maxwidth,minscans)
  
  if(all(is.na(res))){
    df <- NA
    
  }else{
    
    pkscns <- res[[1]]
    df <- res[[2]]
    
    rt <- rt[-pkscns]
    int <- int[-pkscns]
    
    doWhile <- length(rt)>minscans
    
    while(doWhile){ # while check?
      # print("start While")
      # doWhile <- F
      # print(doWhile)
      
      res <- find.LocalMax(rt,int,minwidth,maxwidth,minscans)
      
      natest <- is.na(res[1])
      
      if(!(natest)){
        
        ndf <- res[[2]]
        rtcheck <- any(ndf$RTmin>=df$RTmin & ndf$RTmax<=df$RTmax)
        maxcheck <- any(ndf$maxoRT==df$maxoRT & ndf$maxo==df$maxo)
        
        if(rtcheck & maxcheck){
          doWhile <- F
        }else{
          df <- rbind(df,ndf)
          pkscns <- res[[1]]
          
          rt <- rt[-pkscns]
          int <- int[-pkscns]
          
          # plot(rt,int)
        }
        
      }else{
        doWhile <- F
      }
    }
    
    return(df)
  }
}


peakSim <- function(int,minscans,fit.p){
  #return True or False 
  
  int <- int[which(int>0)] # 
  
  if(length(int)<minscans){return(c(-1))}
  
  maxidx <- which(int==max(int))[1]
  
  cor_idx <- list(1:(maxidx-1),(maxidx+1):length(int))
  c_m <- sapply(cor_idx,length)
  cor_idx <- list(cor_idx[[1]][(1+c_m[[1]]-min(c_m)):c_m[[1]]],
                  cor_idx[[2]][min(c_m):1])
  
  c_m <- sum(sapply(cor_idx,length))
  lmr <- suppressWarnings(lm(x~y,data=data.frame("x"=log10(int[cor_idx[[1]]]),
                                                 "y"=log10(int[cor_idx[[2]]])
  )))
  
  if(is.na(lmr$coefficients[["y"]])){
    lmr_pval <- 1
  }else{
    sum_lmr <- summary(lmr)
    lmr_pval <- sum_lmr$coefficients["y","Pr(>|t|)"]
    if(is.na(lmr_pval)){lmr_pval <- 1}
  }
  # try rerunning lm fit cleaning "outliers"
  if(lmr_pval>fit.p & lmr_pval<0.1){
    torm <- round(log10(c_m/2))
    if(torm>0){
      lmr_pval <- min(sapply(1:torm,function(x){ 
        torm_x <- order(abs(lmr$residuals),decreasing=T)[x]
        cor_idx <- lapply(cor_idx,function(y){
          if(torm_x %in% y){
            y <- y[-which(y==torm_x)]
          }else{
            y <- y[1:(length(y)-1)]
          }
          return(y)
        })
        lmr <- suppressWarnings(lm(x~y,data=data.frame("x"=int[cor_idx[[1]]],
                                                       "y"=int[cor_idx[[2]]])))
        sum_lmr <- summary(lmr)
        lmr_pval <- sum_lmr$coefficients["y","Pr(>|t|)"]
        return(lmr_pval)
      }))
    }
  }
  return(lmr_pval)
}

findIsos <- function(fTi, targetmz,patfoundi,eic,minwidth,maxwidth,minscans,SNR,fit.p,
                     SampleFiles,s) {
  
  res <- lapply(1:length(patfoundi), function(i) {
    iix <- patfoundi[[i]]
    if (length(iix) >= minscans) {
      
      isoeic <- eic[iix, ]
      rt <- isoeic$retentionTime
      int <- isoeic$basePeakIntensity

      res <- find.All.LocalMax(rt, int, minwidth, maxwidth, minscans)
      
      if(is.na(res[[1]][1])){return()}
      res$patfoundi <- i

      
      check_sim <- sapply(1:nrow(res),function(r) 
        peakSim(int[which(rt>=res$RTmin[r] & rt<=res$RTmax[r])],
                minscans,fit.p))
      
      if (any(check_sim==c(-1))) {
        res <- res[-which(check_sim==c(-1)), ]
        check_sim <- check_sim[-which(check_sim==c(-1))]
      }
      
      if (any(check_sim>=fit.p)) {
        res$area[which(check_sim>=fit.p)] <- NA
      }
      
      res <- res[which(res$SNR > SNR), ]
      
      if(nrow(res)==0){return()}

      if(nrow(res) > 1) { # A or B?
        check_RTp <-
          apply(res, 1, function(x) {
            x["RTmin"] < fTi$RT & x["RTmax"] > fTi$RT
          })
        if (any(!check_RTp)) {
          res <- res[-which(!check_RTp), ]
        }
      }
      
      if(nrow(res) > 1) { # A or B?
        check_RTdev <- abs(res$maxoRT-fTi$RT)/fTi$RT
        res <- res[which(check_RTdev==min(check_RTdev)), ]
      }
      
      
      if(nrow(res)==1) {
        targetmz2 <- targetmz[i, ]
        targetmz2 <- cbind(targetmz2, res[, c("maxo", "area")])
      }else{
        if(nrow(res)>0) {
          print("Too many peaks for")
          print(targetmz[i])
          print(SampleFiles[s])
          return()
        }
      }
      
    } else{
      return()
    }
    
  })
  res <- do.call("rbind", res)
  return(res)
}