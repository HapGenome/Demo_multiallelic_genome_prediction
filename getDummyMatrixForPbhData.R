library(tidyverse)
library(caret)

#### convertToDummy ####
# this function takes the condensed data and lineNames as input and then output dummy matrix

#step1: load condensedStr data 

#strPbhData <- readRDS("/mnt/data/strCondensed_ancestralHaploData_combined_LSH_PBH_20201201_2.rds")
#strPbhData <- readRDS("/domino/datasets/pedigreeBasedHaplotypes/strCondensed_ancestralHaploData_20201201.rds")


# supporting functions

convertBack2PBH <- function(cndStr){
  # cndStr is a condense string from condensStr
  # eg. A::4;T::3;C::1;A::2;T::2 => AAAATTTCAATT
  tmp <- unlist(lapply(unlist(strsplit(cndStr,";")),
                       function(x){y <- unlist(strsplit(x,"::")); z <- as.numeric(y[2]);
                       names(z) <- y[1];return(z)}))
  tmp1 <- data.frame(name=names(tmp),repTimes = as.numeric(tmp),stringsAsFactors = F)
  str1 <- as.character(unlist(apply(tmp1,1,function(x){return(rep(x[1],x[2]))})))
  return(str1)
}


getDummyMatrixForPbhData <- function(subPbh,maf = 0.01,numOfBinsPerRun = 6000,totalBins = 17116,groupNA=TRUE){

#getDummyMatrixForPbhData <- function(lineNames,strPbhData,maf = 0.01,numOfBinsPerRun = 6000,totalBins = 17116,groupNA=TRUE){
  # @parameters
  # 1) lineNames is a vector includes line names  
  # 2) strPbhData is the condensed P-BHs data
  # 3) maf is the minimum allele frequency allowed. Haplotypes lower than this will be set as missing (NA)
  # 4) numOfBinsPerRun is used to set 17116 bins into many trunks. When converting dummy variables, too many bins included will exceed the limit and casue error
  # 5) totalBins is the total number of 0.1 cM genetic bins 
  
  # @ information
  # The default setting worked well for 4255 NA females and the whole process took about 7 mins
  # the output file size is about 3 Gb (4255 lines x 88804 hapBins)
  # Please set the RAM according, Double the output size. 
  
  #allLines <- lineNames[lineNames %in% strPbhData$lineName]
  #print(paste0("preChecking .... ",length(allLines),"/",length(lineNames), " lines are in the P-BHs"))
  #print(paste0("preChecking .... ",length(lineNames)-length(allLines),"/",length(lineNames), " lines are NOT in the P-BHs"))
  
  # step1 : deCompress the condesed data
  #print(paste("step1 .... deCondense PBH for",length(allLines),"lines"))
  #subPbh <- sapply(strPbhData[allLines,"pbh"],convertBack2PBH)
  #colnames(subPbh) <- allLines
  #subPbh <- t(subPbh) # flip the dataset
  
  # step2: 
  print(paste("Step2 ...... Set rare haplotypes as missing ...... "))
  haploOut <- setNAtoRareHaplos(subPbh,maf,totalBins)
  newOutHomo <- haploOut[["newOutHomo"]]
  newOutHet <- haploOut[["newOutHet"]]
  binAlleleCount = apply(newOutHomo,2,function(x){length(table(x))})
  newOutHomo = newOutHomo[,which(binAlleleCount > 1)] # the bin has only 1 count is fixed, no need to continue
  newOutHet = newOutHet[,which(binAlleleCount > 1)]
  rm(haploOut,subPbh)
  
  # step3
  print(paste("Step3 ...... Converting common homozygous haplotypes to dummy matrix ...... "))
  dummyHomo <- convertHomoToDummmy(newOutHomo,numOfBinsPerRun,ncol(newOutHomo),groupNA=TRUE)
  
  # step4 
  print(paste("Step4 ...... Add heterozygous haplotypes to dummy matrix ...... "))
  
  dummyHomo <- AddHetToDummy(newOutHet,dummyHomo,ncol(newOutHet))
  
  print(paste("All Done !!! "))
  
  return(dummyHomo)
  
}


setNAtoRareHaplos <- function(subPbh,maf,totalBins){
  # set rare haplotypes to NA, based on the maf
  # took about 20 seconds for 4255 lines x 17116 bins
  # output a list with two datasets, one homo and one Het data
  
  numOfLines <- nrow(subPbh)
  minCounts <- ceiling(numOfLines * maf)
  newOutHomo <- subPbh
  newOutHet <- subPbh
  
  for(i in 1:ncol(subPbh)){
    #if(i == 1 | i == ncol(subPbh)){print(paste("Running",i,Sys.time()))}
    haplos <- as.character(subPbh[,i])
    counts <- sort(table(haplos),decreasing = T)
    commonHaplos <- names(counts[which(counts >= minCounts & !str_detect(names(counts),"&"))])
    indexHomo <- which(haplos %in% commonHaplos)
    
    if(length(indexHomo) > 0){
      newOutHomo[-indexHomo,i] <- "NA"
    }else{ # no common haplo
      newOutHomo[,i] <- "NA"
      newOutHet[,i] <- "NA"
      next;
    }
    
    hetHaplos <- names(counts[str_detect(names(counts),"&")])
    commonHetHaplos <- as.character()
    for(hetHaplo in hetHaplos){
      if(all(unlist(strsplit(hetHaplo,"&")) %in% commonHaplos)){
        commonHetHaplos  <- c(commonHetHaplos,hetHaplo)
      }
    }
    
    indexHet <- which(haplos %in% commonHetHaplos)
    if(length(indexHet) > 0){
      newOutHet[-indexHet,i] <- "NA"
    }else{ # no common het
      newOutHet[,i] <- "NA"
    }
    
    
  }
  
  colnames(newOutHomo) <- paste0("HB",1:totalBins,"__")
  return(list(newOutHomo = newOutHomo, newOutHet = newOutHet))
  
}

convertHomoToDummmy <- function(newOutHomo,numOfBinsPerRun,totalBins,groupNA){
  ### dummy variables for homozygous haplotype
  ### 6000 per run, otherwise exceed the limit
  ### ~1 min for 4255 NA femal line x 17116 bins 
  dummyHomo <- matrix(, nrow = nrow(newOutHomo), ncol = 0)
  numOfHaplos <- apply(newOutHomo,2,function(x){return(length(unique(x)))})
  if(length(which(numOfHaplos == 1)) > 0){
    newOutHomo <- newOutHomo[,-(which(numOfHaplos == 1))]
  }
  
  numOfRun <- ceiling(min(totalBins,ncol(newOutHomo))/numOfBinsPerRun)
  for( run in 1:numOfRun){
    st <- (run-1)*numOfBinsPerRun + 1
    ed <- min(run*numOfBinsPerRun,min(totalBins,ncol(newOutHomo)))
    print(paste("    ... Processing bins:",st,"-",ed,";"))
    test <- dummyVars("~.",data=newOutHomo[,st:ed],sep=NULL) %>% predict(newOutHomo[,st:ed])
    dummyHomo <- cbind(dummyHomo,test)
  }
  dummyHomo <- dummyHomo * 2
  allBinNames <- colnames(dummyHomo)
  
  if(!groupNA){
    allBinNames <- allBinNames[!str_detect(allBinNames,"__NA$")]
  }
  
  dummyHomo <- dummyHomo[,allBinNames]
  return(dummyHomo)
}

AddHetToDummy <- function(newOutHet,dummyHomo,totalBins){
  ### dummy variables for heterozygous haplotype
  ### run bin by bin; parse het into 2 homos and add into dummyHomo
  ### took ~6 mins for 4255 NA females x 17116 bins
  
  for(bin in 1:totalBins){
    if(bin == 1 | bin %% 5000 == 0 | bin == totalBins){print(paste0("    .... processed: ",bin," bin(s);"))}
    tmp <- newOutHet[,bin]
    names(tmp) = row.names(newOutHet)
    tmp <- tmp[!tmp == "NA"]
    tmp <- tmp[str_detect(tmp,"&")]
    for(hetHaplo in unique(tmp)){
      subHetLines <- names(which(tmp == hetHaplo))
      haploBins <- paste0("HB",bin,"__",unlist(strsplit(hetHaplo,"&")))
      dummyHomo[subHetLines,haploBins] <- 1
    }
  }
  return(dummyHomo)
  
}

getDummy_H_comb <- function(lineData,dummyMatrix){
  # consider male and femlae together
  timestamp()
  dummyM_h <- matrix(0,nrow=nrow(lineData),ncol=ncol(dummyMatrix),
                     dimnames = list(lineData$LINE_NAME,colnames(dummyMatrix)))
  
  for(i in 1:nrow(lineData) ){
    #if((i-1) %% 1000 ==0){print(paste(i,Sys.time()))}
    p1 <- lineData$P1[i]
    p2 <- lineData$P2[i]
    dummyM_h[i,] <- as.numeric(dummyMatrix[p1,])/2 + as.numeric(dummyMatrix[p2,])/2
  }
  return(dummyM_h)
}


getDummy_H_sex <- function(lineData,dummyMatrix_F,dummyMatrix_M){
  #consider male and female separately
  timestamp()
  dummyM_h <- matrix(0,nrow=nrow(lineData),ncol=ncol(dummyMatrix_F) + ncol(dummyMatrix_M),
                     dimnames = list(lineData$LINE_NAME,c(colnames(dummyMatrix_F),colnames(dummyMatrix_M))))
  
  for(i in 1:nrow(lineData) ){
    #if((i-1) %% 1000 ==0){print(paste(i,Sys.time()))}
    p1 <- lineData$P1[i]
    p2 <- lineData$P2[i]
    dummyM_h[i,] <- c(as.numeric(dummyMatrix_F[p1,])/2,as.numeric(dummyMatrix_M[p2,])/2)
  }
  return(dummyM_h)
}



# # To get the converted dummy matrix
# #Depending on your downstream analysis, you may like to remove duplicates from the output. 
# 
# lines <- c("B73","01DKD2","FIDA240","90DJD28","CM105","PA2121","1067-1","01DHD10","FBLL","ABCD")
# dummyMatrix <- getDummyMatrixForPbhData(lines,strPbhData)
# dummyMatrix <- dummyMatrix[,!duplicated(t(dummyMatrix))]
# 
# #res <- res[,!duplicated(t(res))] # Remmember to remove duplicates if needed
# 
# 
# 
# 
# 
# # phenotype data
# lineData <- readRDS("/mnt/data/lineData_8429.rds")
# lineData <- lineData[,c("Name","Gender","Origin","YLD_BE_GCA","yearCoded")]
# sex <- "F"
# fLines <- subset(lineData,Gender == sex & !is.na(YLD_BE_GCA))$Name
# 
# 
# dummyMatrix <- getDummyMatrixForPbhData(fLines,strPbhData)
# dummyMatrix <- dummyMatrix[,!duplicated(t(dummyMatrix))]
# saveRDS(dummyMatrix,"/mnt/BLUP_data/GWS/data/dummyMatrix_blupData_pbh_F.rds")