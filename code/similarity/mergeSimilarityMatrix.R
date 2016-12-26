similarity_folder="/Users/marti/Documents/R/combiMS/cluster/similarity/median/"
pat_pairs=read.table("/Users/marti/Documents/R/combiMS/modelling/combinationMatrix.txt")
cat("pairs of patients",length(pat_pairs$V1),"\n")
sim_files=list.files(similarity_folder,pattern="*.RData",full.names=FALSE)

sim_files_2=lapply(sim_files,function(x){
  unlist(strsplit(x,"\\."))[[1]]
})

sim_files_3=unlist(sim_files_2)
length(sim_files_3)
which(!as.character(pat_pairs$V1) %in% sim_files_3)
pat_pairs$V1[which(!as.character(pat_pairs$V1) %in% sim_files_3)]
# missing CH029_UZ023 CH051_KI025 IB009_IB028

#************************************ 
#generate empty matrix with patient names
#************************************ 

data_folder="/Users/marti/Documents/ebi/combiMS/data/phosphosMergedAbsMax/processed/normalized/secondRoundProcessedMidas/"
patients=list.files(data_folder,pattern="*.csv",full.names=FALSE)

patients=unlist(lapply(patients,function(x){
  unlist(strsplit(x,"\\."))[[1]]
}))

similarityMatrix=matrix(nrow=169,ncol=169)
rownames(similarityMatrix)=patients
colnames(similarityMatrix)=patients

#************************************ 
#****************** fill matrix
#************************************ 
successfulTries=0
errorLoading=F
for(i in 1:length(pat_pairs$V1)){
  
  result=tryCatch(
    {load(paste0(similarity_folder,as.character(pat_pairs$V1[i]),".RData")); successfulTries=successfulTries+1; errorLoading=F},
    error=function(e) {cat("error!!!",as.character(pat_pairs$V1[i]),"\n"); errorLoading=T; return(errorLoading)}
    )


  pat1=strsplit(as.character(pat_pairs$V1[i]),"\\_")[[1]][1]
  pat2=strsplit(as.character(pat_pairs$V1[i]),"\\_")[[1]][2]
  positionRow=grep(pat1,rownames(similarityMatrix))
  positionCol=grep(pat2,colnames(similarityMatrix))

  if (result == T) {
    similarityMatrix[positionRow,positionCol]=NA
    similarityMatrix[positionCol,positionRow]=NA
    
  } else {
    
    score_similarity=distanceModels
    similarityMatrix[positionRow,positionCol]=score_similarity
    similarityMatrix[positionCol,positionRow]=score_similarity
  
  }
}

save(similarityMatrix,file='/Users/marti/Documents/R/combiMS/cluster/analysis/similarityMatrix.RData')

#************************************ 
#****************** test missing files
#************************************ 
pat1="CH029"
pat2="UZ023"
positionRow=grep(pat1,rownames(similarityMatrix))
positionCol=grep(pat2,colnames(similarityMatrix))
similarityMatrix[positionRow,positionCol]
similarityMatrix[positionCol,positionRow]

#************************************ 
#****************** test non missing
#************************************ 
pat1="UZ038"
pat2="UZ039"
positionRow=grep(pat1,rownames(similarityMatrix))
positionCol=grep(pat2,colnames(similarityMatrix))
similarityMatrix[positionRow,positionCol]
similarityMatrix[positionCol,positionRow]
load(paste0(similarity_folder,pat1,"_",pat2,".RData"))
distanceModels




