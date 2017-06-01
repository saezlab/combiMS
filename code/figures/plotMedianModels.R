# QC of combiMS models after fixing median strategy to achieve a single model. 
# When optimization did not converge due to excessive solutions within given relative tolerance of best model,
# optimisation was repeated until convergence was achieved 10 times. Then the resulting models were merged using mergeModels.R. 
#
# Here, we analyse and discard artefacts such as model performance being affected by model size, or by number of solutions grouped.
# The script also produces figure S2 on model QC. 
# This script must be run after mergeModels.R -or requires some vars to be loaded-.
# Marti Bernardo-Faura. June 2015. 

# Adapted for directory structure of the GitHub project by Jakob Wirbel in 2017

library(ggplot2)
library(ggrepel)
library(corrplot)
library(CellNOptR)

# *************** set working directory for relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../utils/newJaccard.R")
source("../utils/calculateScore.R")

# calculateScore.R was benchmarked by editing computeScoreT1multiple.R, where the difference with computeScoreT1.R was to return not just the score but also MSE and num measurements
#source("/Users/marti/Documents/r/combiMS/computeScoreT1multiple.R")
#source("/Users/marti/Documents/r/combiMS/getFitmultipe.R")

#***********************************************************************
# *************** as a reminder, this is the metadata about the models
#***********************************************************************

load("../../files/model_merging/statsModels.RData")
names(total_patientDF)
# [1] "name"               "num_nws"            "avg_score"          "total_nws"          "best_score"         "merged_model_score" "nws_new_reltol"     "convergence"        "best_of_bests"      "run_id"  
total_patientDF[which(total_patientDF$name=="IB030"),]
# name num_nws avg_score total_nws best_score merged_model_score nws_new_reltol convergence best_of_bests run_id
# 671 IB030  281182 0.2041323   4437921  0.1989447          0.1859097           2838           3     0.1859091      1
# 672 IB030  570384 0.2382348   4437921  0.2356783          0.1859097           2838           3     0.1859091      2
# 673 IB030  150570 0.2056265   4437921  0.1859091          0.1859097           2838           3     0.1859091      3
# 674 IB030   38201 0.2296705   4437921  0.1859091          0.1859097           2838           3     0.1859091      4
# 675 IB030  444667 0.2356595   4437921  0.2324375          0.1859097           2838           3     0.1859091      5
# 676 IB030  586940 0.2100592   4437921  0.2070877          0.1859097           2838           3     0.1859091      6
# 677 IB030  610442 0.2378742   4437921  0.2356772          0.1859097           2838           3     0.1859091      7
# 678 IB030  691897 0.2378967   4437921  0.2356767          0.1859097           2838           3     0.1859091      8
# 679 IB030  552197 0.2097010   4437921  0.2070866          0.1859097           2838           3     0.1859091      9
# 680 IB030  511441 0.2101008   4437921  0.2085658          0.1859097           2838           3     0.1859091     10

#and these are the median models in allMedianNetworks loaded by mergeModels.R
load("../../files/median_models/allMedianModels.RData")
numInteractions=length(allMedianNetworks[1,])
numPatients=169

#***********************************************************************
# *************** calculate num edges in each median model
#***********************************************************************
edgesKept=apply(allMedianNetworks,1, function(x){
  length(which(x==1))
})

#***********************************************************************
# *************** calculate mse to investigate whether it is linked to num edges (mse does not account for num edges, as opposed to score)
#***********************************************************************
data_folder='../../data/phosphos_processed/'
patientData=list.files(data_folder,pattern="*.csv",full.names=FALSE)
model_path='../../files/model/combiMSplaneCUT.sif'
fileName=patientData[1]

midas=CNOlist(paste(data_folder,fileName,sep=""))
model=readSIF(model_path)  
sprintf("*********Before preprocessing:  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))

# model=preprocessing(midas,model,expansion=FALSE)
# numInteractions=length(model$reacID)
# sprintf("*********After compressing: %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))

#midas MUST be that of right patient
pb = txtProgressBar(min=0, max=numPatients)
medianNetworkMSE=list()
for (index in 1:numPatients){
  fileName=patientData[index]
  midas=CNOlist(paste(data_folder,fileName,sep=""))
  setTxtProgressBar(pb, (pb$getVal()+1))
  # cat("calculating MSE of median network, data is",fileName,"\n")    
  medianNetworkMSE[[index]]=calculateScore(model=model,CNOlist = midas,bString=allMedianNetworks[index,])
}

MSEs = unlist(medianNetworkMSE)
# MSEs=unlist(lapply(medianNetworkMSE,function(x){
  # (x$deviationPen)/(x$numMeasurements)
# }))

#*************************************************************
# *************** create simpler metadata structure, rejecting info about all runs, keeping only best
#***********************************************************************
patient_names=unique(total_patientDF$name)
mini_DF=data.frame(name=character(),total_nws=vector(),merged_model_score=vector(),nws_new_reltol=vector(),convergence=character(),best_of_bests=vector())
# mini_DF=data.frame()
# names(mini_DF)=names(total_patientDF[c(1,4,6,7,8,9)])
for (i in 1:length(patient_names)){
  positions_patient=which(total_patientDF$name==patient_names[i])
  mini_DF=rbind(mini_DF, total_patientDF[positions_patient[1],c("name","total_nws","merged_model_score","nws_new_reltol","convergence","best_of_bests")])
}
mini_DF$edges=edgesKept
mini_DF$MSE=MSEs
# *************** Save intermediate results
save(mini_DF,file="../../files/model_merging/mini_DF.RData")

#***********************************************************************
# *************** plot distributions
#***********************************************************************
#merged
ggplot(data=mini_DF, aes(reorder(name, merged_model_score), merged_model_score,fill=convergence))+geom_bar(stat="identity") + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#best
ggplot(data=mini_DF, aes(reorder(name, merged_model_score), best_of_bests,fill=convergence))+geom_bar(stat="identity") + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#***********************************************************************
# *************** plot linear model best error vs merged model error
#***********************************************************************
a=ggplot(data=mini_DF,aes(best_of_bests,merged_model_score))+geom_point(shape=1)+geom_smooth(method=lm)
# and label those patients that are from diagonal
a +  geom_text(aes(label=ifelse(abs(best_of_bests-merged_model_score)>0.01,as.character(name),'')),hjust=0,vjust=1) + theme_bw()
# scatter plot to color by convergence
b=ggplot(data=mini_DF,aes(best_of_bests,merged_model_score, color=convergence))+geom_point(shape=1, size=4)
b+theme_bw()
#b +scale_colour_gradient(low = "blue")

#***********************************************************************
# *************** for all patients together: plot linear model error vs total num networks
#***********************************************************************
c= ggplot(data=mini_DF,aes(best_of_bests,total_nws))+geom_point(size=4,color="grey")+geom_smooth(method=lm)
c+theme_bw()
# plot distributions 
d= ggplot(data=mini_DF,aes(reorder(name,total_nws),total_nws, fill=convergence))+geom_bar(stat="identity")
d  + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#***********************************************************************
# *************** for each patient: plot linear model error vs total num networks
#***********************************************************************
e= ggplot(data=total_patientDF,aes(best_score,num_nws,color=name))+geom_point(size=4)+guides(color=FALSE)#+geom_smooth(method=lm)
# e=e+theme_bw()+facet_grid(name~)
# e+theme_bw()
e+facet_wrap(~name,ncol=17)+theme_bw()
# plot distributions 
f= ggplot(data=total_patientDF,aes(reorder(run_id,num_nws),num_nws, fill=convergence))+geom_bar(stat="identity")
f=f  + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
f+facet_wrap(~name,ncol=17)

#***********************************************************************
# *************** plot num edges vs MSE
#***********************************************************************
g= ggplot(data=mini_DF,aes(edges,MSE))+geom_smooth(method=lm)+geom_point(size=4)+theme_bw()
g+theme_bw()

g= ggplot(data=mini_DF,aes(edges))#+geom_bar()+theme_bw()
g+geom_density(fill=NA,aes(y=..count..))+geom_histogram(aes(y=..count..))+theme_bw()

g=ggplot(data=mini_DF, aes(reorder(name, edges), edges,fill=MSE)) #+ scale_fill_manual(values=Greys)
g=g+geom_bar(stat="identity") + theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
g + scale_fill_continuous(low='grey',high='grey11')

#***********************************************************************
# *************** calculate similarity between patients with new median network
#***********************************************************************
distanceModels=jaccNello(allMedianNetworks[1,],allMedianNetworks[2,])

#***********************************************************************
# solve problem: why are there patients below the diagonal in plot "a" best_of_best vs merged_model_score?
# this are models better than found during optimization!
# answer found: this was a bug in computeScore, now fixed and benchmarked against own version of computeScoreT1.R
#***********************************************************************

###
#
#
#
# HERE, I GET DIFFERENT RESULTS (MAYBE BECAUSE I DON'T HAVE THE RIGHT CALCULATESCORE-FUNCTION?)
#
#
#
###
mini_DF[which(mini_DF$best_of_bests-mini_DF$merged_model_score>0.01),'name']
# "IB064" "IB069" "UZ010" "UZ011" "UZ017"
# from here, it runs again
networks_folder='../../files/all_solutions_models/'# "/Users/marti/Documents/R/combiMS/cluster/all/"
patient="IB064.RData"
load(paste0(networks_folder,patient)) # this loads PatientResults
which(PatientResults$AllScores==min(PatientResults$AllScores)) # there are several models with the exaxt same score!
#372303 372696 372706 372754 373726 373811 374029 
# get the first model out of the 7 with (same) best score

best_model=PatientResults$AllNwsPatient[which(PatientResults$AllScores==min(PatientResults$AllScores))[1],]
#plot simulation and data
fileName=strsplit(patient,"\\.")[[1]][1]
midas=CNOlist(paste0(data_folder,fileName,".csv"))
model=readSIF(model_path)  
sprintf("*********Before preprocessing:  %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))
#model=preprocessing(midas,model,expansion=FALSE)
#numInteractions=length(model$reacID)
#sprintf("*********After compressing: %s nodes and %s reactions",length(model$namesSpecies),length(model$reacID))
# Here it fails again
prova=calculateScore(model = model,CNOlist = midas,bString=best_model)
#plotModel(model,midas,bString=best_model)
prova=computeScoreT1(CNOlist=midas,model=model,bString=best_model)
computeScoreT1multiple(CNOlist=midas,model=model,bString=best_model)
