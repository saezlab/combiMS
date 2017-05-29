####################################################################################
## PLOT FIGURE ABOUT NORMALISATION FOR COMBIMS-PAPER
## Plot boxplots for different steps in data pre-processing
## And normalisation statistics
##
## Panel A: Boxplot with the log2-transformed Phospho-Data for all Phosphos for all Donors throughout all Conditions
## Panel B: Boxplot with the log2 fold change in respect to the medium condition
## Panel C: Boxplot with Hill-transformed data
## 
## In these three Panels, highlight the data points by patient with identifier KI044
##
## Panel D: Barplot depicting, in how many patients the respective Phospho was found to
##          be upregulated, downregulated or not significantly regulated at all
##          plus a last bar with overall-values
####################################################################################


## Load Packages
library(CellNOptR) # Version 1.16
library(reshape2) # Version 1.4.1
library(ggplot2) # Version 2.1.0
library(cowplot) # Version 0.6.2


####################################################################################################################
## READ DATA #######################################################################################################

## Set working directory - make sure to update this line for use on your own computer
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Vector with patients to include into analysis
Patients = as.character(read.csv("../../files/annotations/annot169pat_v2.csv")$ID)

###################################################################################################################
## merged raw data

## Function to read the raw merged data, yield matrix as output
read_phospho_data_raw_merged <- function(Patient){
  
  ## Read MIDAS file and turn it into CNO list
  M = readMIDAS(paste("../../data/phosphos_merged/", Patient, "_phosphos_midas.csv", sep=""), verbose=F)
  CNO_List = makeCNOlist(M, subfield = T, verbose=F)
  
  ## Extract Data from CNOlist (wihtout BSA and PE, which were used as internal technical control)
  Matrix = CNO_List$valueSignals[[2]][,-c(18,19)]
  
  ## Name rows and columns with the phospho proteins and the stimuli, without BSA and PE
  rownames(Matrix) <- CNO_List$namesStimuli
  colnames(Matrix) <- CNO_List$namesSignals[-c(18,19)]
  
  return(Matrix)
}


## Prepare Matrix for data
Temp = read_phospho_data_raw_merged(Patients[1])
M = melt(Temp)
Stimuli = M$Var1
Phosphos = M$Var2
# Hack: rename MKO3 to MK03
levels(Phosphos)[13] = 'MK03'

Matrix_merged = matrix(data = NA, nrow = length(Patients), ncol = length(Temp))
rownames(Matrix_merged) = Patients

pb = txtProgressBar(min = 0, max =length(Patients), initial = 0) 
## Fill Matrices with data
for (p in 1:length(Patients)){
  Temp = read_phospho_data_raw_merged(Patients[p])
  Matrix_merged[p,] = as.vector(Temp)
  setTxtProgressBar(pb, (pb$getVal()+1))
}

####################################################################################################################
## merged log2


Matrix_log = log(Matrix_merged, base=2)

####################################################################################################################
## log2FC


Matrix_FC = Matrix_log

## Calculate Fold Change compared to Medium
medium = 1
for (i in 1:dim(Matrix_FC)[2]){
  if (Stimuli[i] == 'medium'){
    medium = i
  } else {
    Matrix_FC[,i] = Matrix_FC[,i] - Matrix_FC[,medium]
  } 
}
Matrix_FC[, which(Stimuli == 'medium')] = rep(0, dim(Matrix_FC)[1])

####################################################################################################################
## Hill Transformed


## Function to read the normalised, Hill transformed Data
read_phospho_data_normalised <- function(Patient){
  
  ## Read MIDAS file and turn it into CNO list
  M = readMIDAS(paste("../../data/phosphos_normalised/", Patient, ".csv", sep=""), verbose=F)
  CNO_List = makeCNOlist(M, subfield = T, verbose=F)
  
  ## Extract Data from CNOlist
  Matrix = rbind(CNO_List$valueSignals[[1]][1,], CNO_List$valueSignals[[2]])
  
  ## Name rows and columns with the phospho proteins and the stimuli, without BSA and PE
  rownames(Matrix) <- c('medium', CNO_List$namesStimuli)
  colnames(Matrix) <- CNO_List$namesSignals
  
  return(Matrix)
}

## Prepare Matrix for data
Matrix_Hill = Matrix_FC
pb = txtProgressBar(min = 0, max = length(Patients), initial = 0) 
## Fill Matrices with data
for (p in 1:length(Patients)){
  Temp = read_phospho_data_normalised(Patients[p])
  Matrix_Hill[p,] = as.vector(Temp)
  setTxtProgressBar(pb, (pb$getVal()+1))
}
## For Barplots later we need an unchanged matrix with the Hill transformed data
Matrix_Hill_2 = Matrix_Hill

## Invert negative regulation again, to show up and down-regulation in boxplots
medium = 1
for (i in 1:dim(Matrix_Hill)[2]){
  if (Stimuli[i] == 'medium'){
    medium = i
  } else {
    Matrix_Hill[,i] = Matrix_Hill[,i] - Matrix_Hill[,medium]
  } 
}
Matrix_Hill[, which(Stimuli == 'medium')] = rep(0, dim(Matrix_Hill)[1])


####################################################################################################################
## Dalek says: "ORDER!!!!" 
## For the boxplots
Melted = melt(Matrix_log)
names(Melted) = c('Stimulus', 'Phospho', 'Value')
Melted$Stimulus = rep(Stimuli, each = length(Patients))
Melted$Phospho = rep(Phosphos, each = length(Patients))
Melted$Patients = rep(Patients, dim(Matrix_log)[2])  

## Reorder according to mean
m = sapply(levels(Melted$Phospho), function(x){mean(subset(Melted$Value, Melted$Phospho == x), na.rm=T)})
sorted_by_means = sort(m)



####################################################################################################################
## PREPARE BARPLOTS ################################################################################################


## New Matrix to save the information about up-regulation or down-regulation
Regulation = matrix(data=0, ncol = length(levels(Phosphos)), nrow = 3)
colnames(Regulation) = levels(Phosphos)
rownames(Regulation) = c('phosphorylated', 'dephosphorylated', 'None')

Melted = melt(Matrix_Hill_2)
names(Melted) = c('Stimulus', 'Phospho', 'Value')
Melted$Stimulus = rep(Stimuli, each = length(Patients))
Melted$Phospho = rep(Phosphos, each = length(Patients))
Melted$Patients = rep(Patients, dim(Matrix_Hill)[2])  

pb = txtProgressBar(min = 0, max = length(levels(Phosphos))*length(Patients), initial = 0) 
for (p in levels(Phosphos)){
  for (p2 in Patients){
    Temp = subset(Melted, Melted$Phospho == p & Melted$Patients == p2)
    
    if (!is.na(Temp$Value[1]) & Temp$Value[1] == 1){
      Regulation[2,which(colnames(Regulation)==p)] = Regulation[2,which(colnames(Regulation)==p)] + 1
    } else {
      if (sd(Temp$Value, na.rm=T) == 0){
        Regulation[3,which(colnames(Regulation)==p)] = Regulation[3,which(colnames(Regulation)==p)] + 1
      } else {
        Regulation[1,which(colnames(Regulation)==p)] = Regulation[1,which(colnames(Regulation)==p)] + 1
      }
    }
    setTxtProgressBar(pb, (pb$getVal()+1))
  }
}

## Get Percentages
Regulation = Regulation/169


## Add last column with Total
Regulation = cbind(Regulation, rowSums(Regulation))
colnames(Regulation)[18] = 'Total'
Regulation[,18] = Regulation[,18]/sum(Regulation[,18])

####################################################################################################################
## Function to plot Boxplots

plotBoxplot = function(Matrix, ylim, ylabel){
  
  ## Melt Matrix and change Information in melted Dataframe
  Melted = melt(Matrix)
  names(Melted) = c('Stimulus', 'Phospho', 'Value')
  Melted$Stimulus = rep(Stimuli, each = length(Patients))
  Melted$Phospho = rep(Phosphos, each = length(Patients))
  Melted$Patients = rep(Patients, dim(Matrix)[2])
  Melted$KI044 = rep(FALSE, length.out = length(Melted$Patients))
  Melted$KI044[which(Melted$Patients == 'KI044')] = TRUE
  
  ## Relevel to order it according to mean
  Melted$Phospho = factor(Melted$Phospho, names(sorted_by_means))
  
  ## Plot in ggplot Object
  p = ggplot(Melted, aes(x = Phospho, y = Value)) +
    theme_classic() + xlab("") + ylab(ylabel) + ylim(ylim) + theme(panel.grid.major=element_blank(), 
                                                                   panel.grid.minor=element_blank(), 
                                                                   panel.border = element_blank(), 
                                                                   axis.line.x = element_line(), axis.line.y = element_line()) + 
    ## Tilted x Axis ticks
    theme(axis.text.x = element_text(family = 'sans', size=8, angle = 45, hjust = 1), 
          axis.title.y=element_text(family='sans', size=8), axis.text.y=element_text(family='sans', size=8)) +
    ## Remove Outliers from the boxplots
    geom_boxplot(outlier.shape=NA) + 
    ## Jitter data points on boxplots
    geom_point(position=position_jitter(width=0.2), colour='black', alpha=0.2, size=0.01, shape=16)
  
  ## Plot Patient KI044 on top
  Sub_Melted = subset(Melted, Melted$KI044 == TRUE)
  p = p + geom_point(data= Sub_Melted, colour = 'orange', position = position_jitter(width=0.15), size= .25)
  
  ## Output
  return(p)
}



####################################################################################################################
## PLOT BOXPLOTS ###################################################################################################

####################################################################################################################
## Panel A: log2 Data
####################################################################################################################
p1 = plotBoxplot(Matrix_log, ylim=c(6, 16), ylabel=expression('log'['2']*'(MFI)'))
####################################################################################################################
## Panel B: log2 FC
####################################################################################################################
p2 = plotBoxplot(Matrix_FC, ylim=c(-7,7), ylabel=expression('log'['2']*'(Fold Change)'))
####################################################################################################################
## Panel C: Hill transformed
####################################################################################################################
p3 = plotBoxplot(Matrix_Hill, ylim=c(-1,1), ylabel='Normalized Value')
####################################################################################################################
## PLOT BARPLOT ####################################################################################################

## Transform into Dataframe
Melted = melt(Regulation*100)
names(Melted) = c('Significance', 'Phospho', 'Value')

## Relevel according to mean
Melted$Phospho = factor(Melted$Phospho, c(names(sorted_by_means), 'Total'))
Melted$Significance = factor(Melted$Significance, c('phosphorylated', 'dephosphorylated', 'None'))
Melted$cut = c(rep('Phospho-Proteins', 51), rep('Total',3))

p4 = ggplot(Melted, aes(x = Phospho, y = Value, fill = Significance)) + 
  geom_bar(stat='identity', alpha=0.6) + 
  xlab("") + ylab("Percentage [%]") + 
  scale_fill_manual(values = c('steelblue1', 'gold', 'grey75')) +
  facet_grid(~cut, scales = "free_x", space='free', margins=F) +
  theme_classic() + theme(strip.background=element_blank(), strip.text.x = element_blank()) + 
  theme(axis.text.x = element_text(family = 'sans', size=8, angle = 45, hjust = 1), 
        legend.position = 'top', legend.title=element_blank(), 
        legend.text=element_text(family='sans', size=8), legend.margin=unit(.1, 'mm'),
        legend.key.size=unit(10, 'pt'), axis.title.y=element_text(family='sans', size=8), axis.text.y=element_text(family='sans', size=8))

####################################################################################################################
## PLOT ALL WITH COWPLOT ###########################################################################################
####################################################################################################################

## Divert output into pdf
pdf("../../figures/figure_data_normalisation.pdf", width=7, height=5.6, onefile = FALSE)

plot_grid(p1, p2, p3, p4, labels=c('A', 'B', 'C', 'D'))

dev.off()
