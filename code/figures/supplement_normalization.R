####################################################################################
## PLOT SUPPLEMENTARY FIGURE ABOUT NORMALISATION PIPELINE FOR  THECOMBIMS-PAPER
## Plot densities for the different proteins and how they change 
## after each normalization step
## 
## Only the structure of the figure, annotations will be added
## with an image manipulation program
## 
## Panel A: Densities of log2 fold changes for all phosphos for one patient
## Panel B-D: Densities of up/down/non-regulated phosphos
## Panel E-G: Densities of up/down/non-regulated phosphos after normalization
## Panel H: Densities of all normalized phosphos
##
####################################################################################

# Load libraries
library(ggplot2)
library(reshape2)
library(CellNOptR)
library(cowplot)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Define colours for all phosphos
colour_phosphos = c("#00549F",
                    "dodgerblue",
                    "turquoise",
                    "saddlebrown",
                    "orangered",
                    "lightsalmon",
                    "#7A6FAC",
                    "#612158",
                    "#A11035",
                    "#CC071E",
                    "#F6A800",
                    "#BDCD00",
                    "#57AB27",
                    "#0098A1",
                    "#006165",
                    "#E30066",
                    'navyblue')

####################################################################################
# READ DATA

# Read MIDAS file, convert into CNOlist, and extract data columns
M = readMIDAS(paste0("../../data/phosphos_merged/IB015_phosphos_midas.csv"), verbose=FALSE)
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

####################################################################################
# PANEL A
plot_df_a = melt(data_fc)
colnames(plot_df_a) = c('stimulus', 'phospho', 'logFC')
levels(plot_df_a$phospho) = sort(levels(plot_df_a$phospho))
p_a = ggplot(plot_df_a, aes(x=logFC, color=phospho)) + 
  geom_density() + theme_cowplot() + 
  scale_color_manual(values=colour_phosphos) + 
  guides(color=guide_legend(ncol=2))

####################################################################################
# ************ for each protein, find out which distribuitions are skewed, and to which side
# calculate t-test and keep pvalue. for each protein, a patient has a distribution
# prepare a matrix to hold the pvalue for each patient, that this protein is skewed to the right
numProt=length(colnames(data_fc))
pvalueGreater=vector()
pvalueLess=vector()
# to the right?  (i.e. for that prot, if significant the patient has a signfivantly positive fold change)
for (i in 1:numProt){
  pvalueGreater[i]=t.test(data_fc[,i],mu=0,alternative="greater")$p.value
}
positiveProt=which(pvalueGreater<0.05)
# to the left?  (i.e. for that prot, if significant the patient has a signfivantly negative fold change)
for (i in 1:numProt){
  pvalueLess[i]=t.test(data_fc[,i],mu=0,alternative="less")$p.value
}
negativeProt=which(pvalueLess<0.05)
# proteins that fail both tests are not significantly pos or neg phosphorylated
for (i in 1:numProt){
  noiseProt=intersect(which(pvalueGreater>0.05),which(pvalueLess>0.05))
}

####################################################################################
# PANEL B-D
# noise proteins
nameProts=colnames(data_fc)[noiseProt]
plot_df_b=plot_df_a[which(plot_df_a$phospho %in% nameProts),]
p_b = ggplot(plot_df_b, aes(logFC, colour=phospho)) + 
  geom_density() +
  geom_vline(xintercept=0, colour='grey') +
  theme(legend.position = 'none') + scale_color_manual(values=colour_phosphos[noiseProt])

# negatively phosphorylated
nameProts=colnames(data_fc)[negativeProt]
plot_df_c=plot_df_a[which(plot_df_a$phospho %in% nameProts),]
p_c = ggplot(plot_df_c, aes(logFC, colour=phospho)) + 
  geom_density() +
  geom_vline(xintercept=0, colour='grey') +
  theme(legend.position = "none") + 
  scale_color_manual(values=colour_phosphos[negativeProt])

# positively phosphorylated
nameProts=colnames(data_fc)[positiveProt]
plot_df_d=plot_df_a[which(plot_df_a$phospho %in% nameProts),]
p_d = ggplot(plot_df_d, aes(logFC, colour=phospho)) + 
  geom_density() +
  geom_vline(xintercept=0, colour='grey') +
  theme(legend.position = 'none')+
  scale_color_manual(values=colour_phosphos[positiveProt])

#************** step3: correct pos and neg distributions (non-significant still need to be corrected)
# in positive distributions: leave positives unchanged. between 0 and -1, replace by 0. below -1 replace by NA
data_fc_c = data_fc
numExp=dim(data_fc)[1]
if (length(positiveProt)>0){
  for(i in 1:numExp){
    for (j in 1:length(positiveProt)){
      if (data_fc_c[i,positiveProt[j]]>=(-1) && data_fc_c[i,positiveProt[j]]<0 && !is.na(data_fc_c[i,positiveProt[j]])){
        data_fc_c[i,positiveProt[j]]=0
      }
      
      if(data_fc_c[i,positiveProt[j]]<(-1) && !is.na(data_fc_c[i,positiveProt[j]])){
        data_fc_c[i,positiveProt[j]]=NA
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
      if (data_fc_c[i,negativeProt[j]]>0 && data_fc_c[i,negativeProt[j]]<=1 && !is.na(data_fc_c[i,negativeProt[j]])){
        data_fc_c[i,negativeProt[j]]=0
      } else if(data_fc_c[i,negativeProt[j]]>=1 && !is.na(data_fc_c[i,negativeProt[j]])){
        data_fc_c[i,negativeProt[j]]=NA
      } else if (data_fc_c[i,negativeProt[j]]!=0  && !is.na(data_fc_c[i,negativeProt[j]])){
        data_fc_c[i,negativeProt[j]]=abs(data_fc_c[i,negativeProt[j]]) 
      }
    }
  }
}

# ************* step 4: correct non-significant distributions (not done before because it crashes normalization by 0)  
if(length(noiseProt)>0){
  for(i in 1:numExp){
    for (j in 1:length(noiseProt)){
      if(data_fc_c[i,noiseProt[j]]>=(-1) && data_fc_c[i,noiseProt[j]]<1 && !is.na(data_fc_c[i,noiseProt[j]])){
        data_fc_c[,noiseProt[j]]=0} else if ( (data_fc_c[i,noiseProt[j]]<(-1) || data_fc_c[i,noiseProt[j]]>1) && !is.na(data_fc_c[i,positiveProt[j]])){
          data_fc_c[,noiseProt[j]]=NA
      }
    }
  }
}

# ************* step 5: Hill function on FCz, we are using as EC50 the median of each signal
# taula2 = 1 / (1 + (EC50 / taula2)^HillCoef)
data_h=data_fc_c
indexAllButNoise=setdiff(seq(1:ncol(data_fc_c)), noiseProt)
data_h[,indexAllButNoise]=apply(data_fc_c[,indexAllButNoise], 2, function(column){
  1 / (1 + (median(column,na.rm=TRUE) / column) ^ HillCoef)
})

#   #************* step6: shift negatives
if (length(negativeProt)>0){
  for (i in 1:length(negativeProt)){
    data_h[,negativeProt[i]] = 1- abs(data_h[,negativeProt[i]])
  }
}

####################################################################################
# PANEL E-G
# noise proteins
plot_df_c = melt(data_fc_c)
colnames(plot_df_c) = c('stimulus', 'phospho', 'logFC')
nameProts=colnames(data_fc_c)[noiseProt]
plot_df_e=plot_df_c[which(plot_df_c$phospho %in% nameProts),]
p_e = ggplot(plot_df_e, aes(logFC, colour=phospho)) + 
  geom_density(adjust=.01) +
  geom_vline(xintercept=0, colour='grey') +
  theme(legend.position = 'none') + scale_color_manual(values=colour_phosphos[noiseProt])

# negatively phosphorylated
nameProts=colnames(data_h)[negativeProt]
plot_df_f=plot_df_c[which(plot_df_c$phospho %in% nameProts),]
p_f = ggplot(plot_df_f, aes(logFC, colour=phospho)) + 
  geom_density() +
  geom_vline(xintercept=0, colour='grey') +
  theme(legend.position = "none") + 
  scale_color_manual(values=colour_phosphos[negativeProt])

# positively phosphorylated
nameProts=colnames(data_h)[positiveProt]
plot_df_g=plot_df_c[which(plot_df_c$phospho %in% nameProts),]
p_g = ggplot(plot_df_g, aes(logFC, colour=phospho)) + 
  geom_density() +
  geom_vline(xintercept=0, colour='grey') +
  theme(legend.position = 'none')+
  scale_color_manual(values=colour_phosphos[positiveProt])


####################################################################################
# PANEL H

####################################################################################
# PLOT EVERYTHING TOGETHER
# pdf('../figures/supplementary_normalization.pdf', size=c(8.3, 11.7))
top_row = plot_grid(p_a, NULL, ncol=2)
upper_middle_row = plot_grid(p_b, p_c, p_d, ncol=3)
lower_middle_row = plot_grid(p_e, p_f, p_g, ncol=3)
bottom_row = plot_grid(NULL, NULL, p_h, ncol=3)

plot_grid(top_row, upper_middle_row, lower_middle_row, top_row, nrow=4)
# dev.off()