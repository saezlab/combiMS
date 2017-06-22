# ######################################################################################################
# NORMALIZATION PIPELINE FOR COMBIMS
#
# IS STILL A BIG TO DO
#
#
#
# ######################################################################################################

# Load libraries
library(reshape2)
library(CellNOptR)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Find data files
file_list = list.files('../../data/phosphos_raw_clean/')

# ######################################################################################################
# FUNCTION FOR COMBINATION OF BOTH TIMEPOINTS
combine_time_points = function(file_name){}

# ######################################################################################################
# FUNCTION FOR HILL-NORMALIZATION
normalise_with_Hill=function(file_name,HillCoef=2){
  
  ## Read MIDAS file, convert into CNOlist, and extract data columns
  M = readMIDAS(paste0("../../data/phosphos_merged/", file_name), verbose=FALSE)
  M2 = makeCNOlist(M, subfield = TRUE, verbose=FALSE)
  data = M2$valueSignals[[2]]
  row.names(data) = M2$namesStimuli
  colnames(data) = M2$namesSignals
  
  # remove BSA and PE (internal controls)
  data = data[,-c(18,19)]
  
  #************* step1: X= log2 (rawx)
  data_log_2=log2(data)
  
  #************* step2: compute Fold Change, since we have log2 we just substract
  data_fc = t(apply(data_log_2, 1, FUN=function(x){return(x-data_log_2[1,])}))

  # ************ step3: for each protein, find out which distribuitions are skewed, and to which side
  # calculate t-test and keep pvalue. for each protein, a patient has a distribution
  # prepare a matrix to hold the pvalue for each patient, that this protein is skewed to the right
  numProt=length(colnames(data_fc))
  pvalueGreater=vector()
  pvalueLess=vector()
  # to the right?  (i.e. for that prot, if significant the patient has a signfivantly positive fold change)
  for (i in 1:numProt){
    pvalueGreater[i]=t.test(data_fc[-1,i],mu=0,alternative="greater")$p.value
  }
  positiveProt=which(pvalueGreater<0.05)
  # to the left?  (i.e. for that prot, if significant the patient has a signfivantly negative fold change)
  for (i in 1:numProt){
    pvalueLess[i]=t.test(data_fc[-1,i],mu=0,alternative="less")$p.value
  }
  negativeProt=which(pvalueLess<0.05)
  # proteins that fail both tests are not significantly pos or neg phosphorylated
  for (i in 1:numProt){
    noiseProt=intersect(which(pvalueGreater>0.05),which(pvalueLess>0.05))
  }
  
  #************** step4: correct pos and neg distributions (non-significant still need to be corrected)
  # in positive distributions: leave positives unchanged. between 0 and -1, replace by 0. below -1 replace by NA
  numExp=dim(data_fc)[1]
  if (length(positiveProt)>0){
    
    for(i in 1:numExp){
      for (j in 1:length(positiveProt)){
        if (data_fc[i,positiveProt[j]]>(-1) && data_fc[i,positiveProt[j]]<0 && !is.na(data_fc[i,positiveProt[j]])){
          data_fc[i,positiveProt[j]]=0
        }
        
        if(data_fc[i,positiveProt[j]]<(-1) && !is.na(data_fc[i,positiveProt[j]])){
          data_fc[i,positiveProt[j]]=NA
        }
      }
    }
    
  }
  # in negative distributions: between 0 and 1, replace by 0. above 1 replace by NA. 
  # Last, turn remaining negatives into positives.
  if (length(negativeProt)>0){
    print("there's negative distributions")
    for(i in 1:numExp){
      for (j in 1:length(negativeProt)){
        if (data_fc[i,negativeProt[j]]>0 && data_fc[i,negativeProt[j]]<1 && !is.na(data_fc[i,negativeProt[j]])){
          data_fc[i,negativeProt[j]]=0
        } else if(data_fc[i,negativeProt[j]]>1 && !is.na(data_fc[i,negativeProt[j]])){
          data_fc[i,negativeProt[j]]=NA
        } else if (data_fc[i,negativeProt[j]] && !is.na(data_fc[i,negativeProt[j]])){
          data_fc[i,negativeProt[j]]=abs(data_fc[i,negativeProt[j]]) 
        }
      }
    }
    
  }
  
  # ************* step 5: correct non-significant distributions (not done before because it crashes normalization by 0)  
  if(length(noiseProt)>0){
    for(i in 1:numExp){
      for (j in 1:length(noiseProt)){
        
        if(data_fc[i,noiseProt[j]]>(-1) && data_fc[i,noiseProt[j]]<1 && !is.na(data_fc[i,noiseProt[j]])){
          data_fc[,noiseProt[j]]=0} else if ( (data_fc[i,noiseProt[j]]<(-1) || data_fc[i,noiseProt[j]]>1) && !is.na(data_fc[i,positiveProt[j]])){
            data_fc[,noiseProt[j]]=NA
          }
        
      }
    }
  }
  
  
  
  # ************* step 6: Hill function on FCz, we are using as EC50 the median of each signal
  # taula2 = 1 / (1 + (EC50 / taula2)^HillCoef)
  taula3=data_fc
  indexAllButNoise=setdiff(seq(1:ncol(data_fc)), noiseProt)
  taula3[,indexAllButNoise]=apply(data_fc[,indexAllButNoise], 2, function(column){
    1 / (1 + (median(column,na.rm=TRUE) / column) ^ HillCoef)
  })
  
  #EC50 <- colMedians(as.matri(taula2))
  #taula2 = 1 / (1 + (EC50 / taula2)^HillCoef)
  
  
  
  
  # ************* step7: shift negatives
  if (length(negativeProt)>0){
    for (i in 1:length(negativeProt)){
      taula3[,negativeProt[i]] = 1- abs(taula3[,negativeProt[i]])
    }
  }
  
  #************* replace DV values in midas with normalized values
  M2$namesSignals = colnames(taula3)
  M2$valueSignals[[1]] = t(matrix(rep(taula3[1,],21), ncol=21))
  M2$valueSignals[[2]] = taula3
  M2$valueVariances[[1]] = M2$valueVariances[[1]][,-c(18, 19)]
  M2$valueVariances[[2]] = M2$valueVariances[[2]][,-c(18, 19)]
  
  #************* save normalized values for this patient
  # writeMIDAS(M2, filename='~/Desktop/test.csv')
  
  return(taula)
  
}

# ######################################################################################################
# FUNCTION FOR FURTHER PROCESSING
process_normalized = function(file_name){}