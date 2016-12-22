# ###############################################################################################################################################
# Based on script: analyzePatientSimilarity.R by Marti July 2015
#
# Minor modification for nicer plotting by Jakob in December 2016
# #################################################################################################################################################

library(apcluster)
library(gclus)
library(ggplot2)
library(reshape2)
library(NMF)
library(corrplot)
library(cowplot)
library(gdata)
  
# Use relative paths instead of absolut paths for the files
# in Rstudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ##########################################################################################################
# PART A OF FIGURE
# ##########################################################################################################

# ##########################################################################################################
# PART B OF FIGURE
# ##########################################################################################################

# ##########################################################################################################
# load annotation file and similarity BETWEEN matrix
# ####################################################################################################################
load("../../files/similarityMatrix.RData")
annot=read.csv("../../files/annot169pat_v2.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)

# ************************************ 
# add some fields to annotation for ggplot
annot2=annot
# correct that PPMS patients appear in group instead of by untreated and treatment
# Rename groups for consistency with other figures in the paper
annot2[which(annot2$Group=='PPMS' & annot2$Treatment=='no'),'Group']='Untreated'
annot2[which(annot2$Group=='PPMS' & annot2$Treatment!='no'),'Group']=annot2[which(annot2$Group=='PPMS' & annot2$Treatment!='no'),'Treatment']
annot2[which(annot2$Group=='Gilenya'),'Group']='FTY'
annot2[which(annot2$Group=='Fingolimod'),'Group']='FTY'
annot2[which(annot2$Group=='Glatiramer'),'Group']='GAL'
annot2[which(annot2$Group=='Interferon B'),'Group']='IFNb'
annot2[which(annot2$Group=='Natalizumab'),'Group']='NTZ'
annot=annot2
rm(annot2)

# ************************************ 
# include similarities WITHIN
load('../../files/similarityWITHINMeans.RData')
for (i in 1:dim(similarityMatrix)[1]){
  similarityMatrix[i,i]=means_within[i]
}
 

# ##########################################################################################################
# create dataframe of similarities to plot
# ##########################################################################################################
sim_df=melt(similarityMatrix)
sim_df$Group=rep(F,dim(similarityMatrix)[1])
sim_df$GroupSim=rep(F,dim(similarityMatrix)[1])
sim_df$Condition=rep(F,dim(similarityMatrix)[1])
sim_df$ConditionSim=rep(F,dim(similarityMatrix)[1])
sim_df$Self=rep(F,dim(similarityMatrix)[1])
sim_df$RR_Untreated=rep(F,dim(similarityMatrix)[1])

colnames(sim_df)=c('SimWith','Patient','Similarity','Group','GroupSim','Condition','ConditionSim','Self','RR_Untreated')  
sim_df$GroupSim=rep(annot$Group,length(annot$ID))
sim_df$ConditionSim=rep(annot$condition,length(annot$ID))
sim_df$all="All"

for (i in 1:length(annot$ID)){
  sim_df$Group[sim_df$Patient==annot[i,'ID']] = annot[i,'Group']
  sim_df$Condition[sim_df$Patient==annot[i,'ID']] = annot[i,'condition'] 
}
pb = txtProgressBar(min=0, max=length(sim_df$Patient))
for (i in 1:length(sim_df$Patient)){
  if (as.character(sim_df$Patient[i])==as.character(sim_df$SimWith[i])){
    sim_df$Self[i]='Self'
  } else {
     sim_df$Self[i]='Another'
  }
  
  if (as.character(sim_df$Patient[i]) %in% annot[which(annot$Disease.Subtype=='RR' & annot$Treatment=='no'),'ID']){
    sim_df$RR_Untreated[i]='RR_Untreated'
  } else {
    sim_df$RR_Untreated[i]='N'
  }
  setTxtProgressBar(pb, (pb$getVal()+1))
}

# ##########################################################################################################
# plot distribution of subsets where similarities are subsetting only between patients with same conditions
# and against all
# ############################################################################################################

sub_same_group_similarity=subset(sim_df,sim_df$Group==sim_df$GroupSim)
sub_self_similarity=subset(sub_same_group_similarity,sub_same_group_similarity$Self=='Self')
sub_IFN=subset(sub_same_group_similarity,sub_same_group_similarity$RR_Untreated=='RR_Untreated')

# delete Tysabri or Placebo, add All similarities, self similarities, and RR_untreated as reference
sub_same_group_similarity = sub_same_group_similarity[-which(sub_same_group_similarity$Group == 'Tysabri or Placebo'),]
sub_self_similarity$Group = 'Self'
sub_IFN$Group = 'RR_Untreated'
sub_all = sim_df
sub_all$Group = 'All'

# boxplot by Group
# In a notched box plot, the notches extend 1.58 * IQR / sqrt(n). This gives a roughly 95 interval for comparing medians
# So non-overlapping notches strongly suggest that medians are significantly different
# This is the case for RR untreated vis RR IFNb

plot_df = rbind(sub_same_group_similarity, sub_self_similarity, sub_IFN, sub_all)

## Reorder Groups
sorted_groups = c('All', 'Self', 'Healthy', 'EGCG', 'FTY', 'IFNb', 'GAL', 'NTZ', 'RR_Untreated', 'Untreated')
plot_df$Group = factor(plot_df$Group, sorted_groups)

boxplot_similarity = ggplot(plot_df) + 
  geom_boxplot(notch=T,aes(x=Group, y=Similarity, fill=Group), outlier.alpha = 0.2, outlier.stroke = 0) +
  theme_classic() + 
  scale_fill_manual(values=c('#ECEDED', '#646567', '#407FB7', '#2D7F83', '#00B1B7', '#8DC060', '#D0D95C', '#FABE50', '#A8859E', '#834E75')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  theme(legend.position = 'none') + 
  ylim(c(0.5, 1))

# ###################################################################################################################################################
# Divert output into pdf
# ###################################################################################################################################################
pdf("../../figures/figure_similarity.pdf", width=7, height=3.5, onefile = FALSE)

plot_grid(boxplot_similarity, boxplot_similarity, labels=c('A', 'B'))

dev.off()
