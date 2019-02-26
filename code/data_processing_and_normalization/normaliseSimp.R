normaliseSimp=function(taula,HillCoef=2){
   
   library(ggplot2)
   library(reshape2)
   
   
   #************* create a subset of the matrix containing only the DV values
   taula2=taula[,grep('DV',colnames(taula))]
   # rename columns for improved plotting by removing "DV"
   for (i in 1:length(names(taula2))){
      colnames(taula2)[i]=strsplit(names(taula2[i]),":")[[1]][2]                     
   }
   
   #************* step1: X= log2 (rawx)
   taula2=log2(taula2)
   
   #************* step2: compute Fold Change, since we have log2 we just substract
   row1=taula2[1,]
   for (i in 1:dim(taula2)[1]){
      taula2[i,] = taula2[i,] - row1
   }
   
   #************* plot log-FC distribution
   meltPatient=melt(taula2)
   names(meltPatient)=c('protein','logFC')
   
   ggplot(meltPatient, aes(logFC, colour=protein)) + 
      geom_density() +
      #geom_bar(position='dodge') +
      # facet_grid(protein~., scales='free') +
      geom_vline(xintercept=0, colour='grey') +
      theme_bw()+ggtitle("patient IB015, 21 experiments x protein")
   
   # ************ for each protein, find out which distribuitions are skewed, and to which side
   # calculate t-test and keep pvalue. for each protein, a patient has a distribution
   # prepare a matrix to hold the pvalue for each patient, that this protein is skewed to the right
   numProt=length(colnames(taula2))
   pvalueGreater=vector()
   pvalueLess=vector()
   # to the right?  (i.e. for that prot, if significant the patient has a signfivantly positive fold change)
   for (i in 1:numProt){
      pvalueGreater[i]=t.test(taula2[i],mu=0,alternative="greater")$p.value
   }
   positiveProt=which(pvalueGreater<0.05)
   # to the left?  (i.e. for that prot, if significant the patient has a signfivantly negative fold change)
   for (i in 1:numProt){
      pvalueLess[i]=t.test(taula2[i],mu=0,alternative="less")$p.value
   }
   negativeProt=which(pvalueLess<0.05)
   # proteins that fail both tests are not significantly pos or neg phosphorylated
   
   

   
    noiseProt=intersect(which(pvalueGreater>=0.05),which(pvalueLess>=0.05))    

   
   #************** plot proteins that are negatively and not significantly phospho
   # plot proteins that are not significantly phosphorylated
   nameProts=names(taula2)[noiseProt]
   subsetNoise=meltPatient[which(meltPatient$protein %in% nameProts),]
   ggplot(subsetNoise, aes(logFC, colour=protein)) + 
      geom_density() +
      #geom_bar(position='dodge') +
      # facet_grid(protein~., scales='free') +
      geom_vline(xintercept=0, colour='grey') +
      theme_bw()+ggtitle("patient IB015, unphosphorylated proteins")
   # plot proteins that are negatively significantly phosphorylated
   nameProts=names(taula2)[negativeProt]
   subsetNegative=meltPatient[which(meltPatient$protein %in% nameProts),]
   ggplot(subsetNegative, aes(logFC, colour=protein)) + 
      geom_density() +
      #geom_bar(position='dodge') +
      # facet_grid(protein~., scales='free') +
      geom_vline(xintercept=0, colour='grey') +
      theme_bw()+ggtitle("patient IB015, negatively phosphorylated proteins")
   # plot proteins that are positively significantly phosphorylated
   nameProts=names(taula2)[positiveProt]
   subsetPositive=meltPatient[which(meltPatient$protein %in% nameProts),]
   ggplot(subsetPositive, aes(logFC, colour=protein)) + 
      geom_density() +
      #geom_bar(position='dodge') +
      # facet_grid(protein~., scales='free') +
      geom_vline(xintercept=0, colour='grey') +
      theme_bw()+ggtitle("patient IB015, positively phosphorylated proteins")
   
   #************** step3: correct pos and neg distributions (non-significant still need to be corrected)
   # in positive distributions: leave positives unchanged. between 0 and -1, replace by 0. below -1 replace by NA
   numExp=dim(taula2)[1]
   if (length(positiveProt)>0){
      
      for(i in 1:numExp){
         for (j in 1:length(positiveProt)){
            if (taula2[i,positiveProt[j]]>=(-1) && taula2[i,positiveProt[j]]<0 && !is.na(taula2[i,positiveProt[j]])){  

               taula2[i,positiveProt[j]]=0
            }
            
            if(taula2[i,positiveProt[j]]<(-1) && !is.na(taula2[i,positiveProt[j]])){
               taula2[i,positiveProt[j]]=NA
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
      
             if (taula2[i,negativeProt[j]]>0 && taula2[i,negativeProt[j]]<=1 && !is.na(taula2[i,negativeProt[j]])){   
               
            taula2[i,negativeProt[j]]=0
            } else if(taula2[i,negativeProt[j]]>1 && !is.na(taula2[i,negativeProt[j]])){
               taula2[i,negativeProt[j]]=NA
            } else if (taula2[i,negativeProt[j]]!=0 && !is.na(taula2[i,negativeProt[j]])){                           

               taula2[i,negativeProt[j]]=abs(taula2[i,negativeProt[j]]) 
            }
         }
      }
      
   }
   
   # ************* step 4: correct non-significant distributions (not done before because it crashes normalization by 0)  
   if(length(noiseProt)>0){
      for(i in 1:numExp){
         for (j in 1:length(noiseProt)){
            
            if(taula2[i,noiseProt[j]]>=(-1) && taula2[i,noiseProt[j]]<=1 && !is.na(taula2[i,noiseProt[j]])){                  

               taula2[,noiseProt[j]]=0} else if ( (taula2[i,noiseProt[j]]<(-1) || taula2[i,noiseProt[j]]>1) && !is.na(taula2[i,positiveProt[j]])){
                  taula2[,noiseProt[j]]=NA
               }
            
         }
      }
   }
   
   
   
   # ************* step 5: Hill function on FCz, we are using as EC50 the median of each signal
   # taula2 = 1 / (1 + (EC50 / taula2)^HillCoef)
   taula3=taula2
   indexAllButNoise=setdiff(seq(1:ncol(taula2)), noiseProt)
   taula3[,indexAllButNoise]=apply(taula2[,indexAllButNoise], 2, function(column){
      1 / (1 + (median(column,na.rm=TRUE) / column) ^ HillCoef)
   })
   
   #EC50 <- colMedians(as.matri(taula2))
   #taula2 = 1 / (1 + (EC50 / taula2)^HillCoef)
   
   
   
   
   #   #************* step6: shift negatives
   if (length(negativeProt)>0){
      for (i in 1:length(negativeProt)){
         taula3[,negativeProt[i]] = 1- abs(taula3[,negativeProt[i]])
      }
   }
   
   #************* plot normalised distributions
   
   normalizedPatient=melt(data.frame(taula3))
   names(normalizedPatient)=c('protein','logFC')
   
   # plot proteins positively phosphpo prots after normalization
   nameProts=names(taula2)[positiveProt]
   posNorma=normalizedPatient[which(normalizedPatient$protein %in% nameProts),]
   ggplot(posNorma, aes(logFC, colour=protein)) + 
      geom_density() +
      #geom_bar(position='dodge') +
      # facet_grid(protein~., scales='free') +
      geom_vline(xintercept=0, colour='grey') +
      theme_bw()+ggtitle("patient IB015, normalized positively phosphorylated proteins")
   
   # plot proteins negatively phosphpo prots after normalization
   nameProts=names(taula2)[negativeProt]
   negNorma=normalizedPatient[which(normalizedPatient$protein %in% nameProts),]
   ggplot(negNorma, aes(logFC, colour=protein)) + 
      geom_density() +
      #geom_bar(position='dodge') +
      # facet_grid(protein~., scales='free') +
      geom_vline(xintercept=0, colour='grey') +
      theme_bw()+ggtitle("patient IB015, normalized negatively phosphorylated proteins")
   
   
   # plot proteins not significantly phosphpo after normalization
   nameProts=names(taula2)[noiseProt]
   noiseNorma=normalizedPatient[which(normalizedPatient$protein %in% nameProts),]
   ggplot(noiseNorma, aes(logFC, fill=protein)) + 
      #geom_density() +
      geom_bar(position='dodge') +
      # facet_grid(protein~., scales='free') +
      geom_vline(xintercept=0, colour='grey') +
      theme_bw()+ggtitle("patient IB015, normalized unphosphorylated proteins")
   
   
   
   # plot all normalized proteins
   ggplot(normalizedPatient, aes(logFC, colour=protein)) + 
      geom_density() +
      #geom_bar(position='dodge') +
      # facet_grid(protein~., scales='free') +
      geom_vline(xintercept=0, colour='grey') +
      theme_bw()+ggtitle("normalized patient IB015, 21 experiments x protein")
   
   #************* replace DV values in midas with normalized values
   taula[,grep('^DV',colnames(taula))]=taula3
   
   return(taula)
   
}

