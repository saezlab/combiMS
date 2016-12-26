# #######################################################################################################################################
# Volcano plots - phosphos and cytos
# #######################################################################################################################################

# ################################
# CellNOptR to read the MIDAS files
library(CellNOptR)
## limma for regression model
library(limma)
## for plotting
library(ggplot2)
library(reshape2)
library(gridBase)
library(igraph)

# #######################################################################################################################################
# # Logic Model Interactions
# #######################################################################################################################################

model_interactions = list(FTY = data.frame(t(matrix(data=c(c('ANTICD3', 'STAT5'), c('EGCG', 'STAT5'), c('NACL', 'MK12'), c('NACL', 'HSPB1'), 
                                                           c('TNFA', 'MK12'), c('TNFA', 'HSPB1'), c('TNFA', 'IKBA'), c('TNFA', 'TF65'), c('TNFA', 'MKO3'),  
                                                           c('IL1A', 'MK12'), c('IL1A', 'HSPB1'), c('IL1A', 'IKBA'), c('IL1A', 'TF65'), c('IL1A', 'MKO3'), 
                                                           c('POLYIC', 'MK12'), c('POLYIC', 'HSPB1'), c('POLYIC', 'IKBA'), c('POLYIC', 'TF65'), c('POLYIC', 'MKO3'), 
                                                           c('LPS', 'MK12'), c('LPS', 'HSPB1'), c('LPS', 'IKBA'), c('LPS', 'TF65'), c('LPS', 'MKO3'),
                                                           c('INS', 'AKT1'), c('INS', 'MP2K1'),
                                                           c('REBIF', 'STAT6')), nrow=2))), 
                          GAL = data.frame(t(matrix(data=c(c('TNFA', 'MK12'), c('TNFA', 'MKO3'), 
                                                           c('LPS', 'MKO3'), c('POLYIC', 'MKO3'), c('IL1A', 'MKO3'),
                                                           c('INS', 'AKT1'), c('BDNF', 'AKT1'),
                                                           c('CONA', 'AKT1'), c('CONA', 'STAT1'),
                                                           c('INFG', 'STAT1'), c('INFG', 'STAT3'),
                                                           c('REBIF', 'STAT1'), c('REBIF', 'STAT3'),
                                                           c('IL6', 'STAT1'), c('IL6', 'STAT3')), nrow=2))), 
                          IFNb = data.frame(t(matrix(data=c(c('NACL', 'MP2KA'), c('NACL', 'MK12'), c('NACL', 'HSPB1'),
                                                           c('TNFA', 'MK12'), c('TNFA', 'HSPB1'), c('TNFA', 'MKO3'),
                                                           c('IL1A', 'MKO3'), c('LPS', 'MKO3'), c('POLYIC', 'MKO3'),
                                                           c('BDNF', 'AKT1'), c('INS', 'AKT1'),
                                                           c('ANTICD3', 'STAT1'), c('REBIF', 'STAT1'), c('INFG', 'STAT1'), c('IL6', 'STAT1')), nrow=2))),
                          EGCG = data.frame(t(matrix(data=c(c('NACL', 'MK12'), c('NACL', 'MKO3'),
                                                            c('TNFA', 'IKBA'), c('TNFA', 'TF65'), c('TNFA', 'MK12'), c('TNFA', 'MKO3'),
                                                            c('POLYIC', 'IKBA'), c('POLYIC', 'TF65'), c('POLYIC', 'MK12'), c('POLYIC', 'MKO3'),
                                                            c('LPS', 'IKBA'), c('LPS', 'TF65'), c('LPS', 'MK12'), c('LPS', 'MKO3'),
                                                            c('IL1A', 'IKBA'), c('IL1A', 'TF65'), c('IL1A', 'MK12'), c('IL1A', 'MKO3'),
                                                            c('CONA', 'MK12'), c('CONA', 'MKO3'),
                                                            c('BDNF', 'MK12'), c('BDNF', 'MKO3'),
                                                            c('IL6', 'MK12'), c('IL6', 'MKO3'), c('IL6', 'STAT3'),
                                                            c('EGCG', 'STAT1'), c('EGCG', 'STAT5'), c('EGCG', 'STAT3'),
                                                            c('ANTICD3', 'STAT1'), c('ANTICD3', 'STAT5'), c('ANTICD3', 'STAT3')), nrow=2))),
                          NTZ = data.frame(t(matrix(data=c(c('NACL', 'MK12'), c('NACL', 'HSPB1'),
                                                           c('TNFA', 'MK12'), c('TNFA', 'HSPB1'), c('TNFA', 'MKO3'),
                                                           c('IL1A', 'MK12'), c('IL1A', 'HSPB1'), c('IL1A', 'MKO3'),
                                                           c('POLYIC', 'MK12'), c('POLYIC', 'HSPB1'), c('POLYIC', 'MKO3'),
                                                           c('LPS', 'MK12'), c('LPS', 'HSPB1'), c('LPS', 'MKO3'),
                                                           c('INS', 'AKT1'),
                                                           c('CONA', 'STAT1'), c('ANTICD3', 'STAT1'), c('REBIF', 'STAT1'), c('INFG', 'STAT1'), c('IL6', 'STAT1')), nrow=2))))

# #######################################################################################################################################
# # Prepare Stuff
# #######################################################################################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Not all donors have phospho-Data available
DonorsWithPhospho = list.files('../../../data/phosphos_normalised/')
DonorsWithPhospho = gsub(".csv", '', DonorsWithPhospho)

## Annotations
annotations = read.table('../../../files/annotations_169_patients.csv', header = T, row.names=1, sep=',')

# Hack for nicer annotation
annotations$Group = as.character(annotations$Group)
group_temp = list('GAL', 'H', 'U', 'U', 'IFNb', 'FTY', 'EGCG', 'NTZ')
names(group_temp) = c("Glatiramer", "Healthy", "Untreated", "PPMS", "Interferon B", "Fingolimod", "EGCG", "Natalizumab")

for (i in 1:dim(annotations)[1]){
  annotations[i, 'Group'] = group_temp[[annotations[i, 'Group']]]
}

# #####################################################################################################################################
# Convert pkn into igraph object and find all possible Stimulus-Readout pairs
# #####################################################################################################################################

# pknmodel = readSIF('../../../files/combiMS_pkn_cut.sif')
# cnolist = readMIDAS('../../../data/phosphos_normalised/CH003.csv')
# cnolist = makeCNOlist(cnolist, subfield = FALSE)
# 
# # convert to igraph object
# edges <- pknmodel$interMat
# adjacency = matrix(0, nrow(edges), nrow(edges))
# for(i in 1:ncol(edges)){
#   startIdx = which(edges[, i]==-1)
#   endIdx = which(edges[, i]==1)
#   
#   adjacency[startIdx, endIdx] = 1
# }
# rownames(adjacency) <- colnames(adjacency) <- rownames(edges)
# 
# netGraph <- graph.adjacency(adjacency, mode = "directed")
# 
# connectivity = matrix(data=NA, ncol=length(cnolist$namesSignals), nrow=length(cnolist$namesStimuli))
# row.names(connectivity) = toupper(cnolist$namesStimuli)
# colnames(connectivity) = cnolist$namesSignals
# 
# for (stim in toupper(cnolist$namesStimuli)){
#   for (sig in cnolist$namesSignals){
#     if (length(all_shortest_paths(netGraph, from=stim, to=sig)$res) == 0){
#       connectivity[which(row.names(connectivity) == stim), which(colnames(connectivity) == sig)] = 0
#     } else {
#       connectivity[which(row.names(connectivity) == stim), which(colnames(connectivity) == sig)] = length(all_shortest_paths(netGraph, from=stim, to=sig)$res[[1]])
#     }
#   }
# }
# colnames(connectivity)[13] = 'MK03'
# possible_interactions = melt(connectivity)
# possible_interactions[['interaction']] = paste(possible_interactions$Var1, possible_interactions$Var2, sep='_')

# #######################################################################################################################################
# # Prepare Metadata, such as Stimuli and Phospho-information
# #######################################################################################################################################

Temp = readMIDAS('../../../data/phosphos_merged/CH003_phosphos_midas.csv', verbose=F)$dataMatrix
Phosphos = as.vector(sapply(colnames(Temp)[which(grepl('DV', colnames(Temp)))], FUN = function(x){strsplit(x, split=':')[[1]][2]}))[1:17]
Stimuli = toupper(as.vector(sapply(colnames(Temp)[which(grepl('TR', colnames(Temp)))], FUN = function(x){strsplit(x, split=':')[[1]][2]}))[2:21])
rm(Temp)

# ######################################################################################################################################
# # Load Phospho data, log2(FC)
# ######################################################################################################################################
ListPhosphos = lapply(DonorsWithPhospho, FUN = function(donor){
  ## Load Data
  Temp = readMIDAS(paste0("../../../data/phosphos_merged/", donor, '_phosphos_midas.csv'), verbose=F)$dataMatrix;
  Temp = Temp[-1,]
  Temp = as.matrix(Temp[,which(grepl('DV:', colnames(Temp)))]);
  Temp = Temp[,-c(18,19)]
  ## Log-transform data
  Temp = log(Temp, base=2)
  ## Compute fold change compared to the medium condition (first row)
  Temp = sapply(1:dim(Temp)[2], FUN = function(column){Temp[,column] = Temp[,column] - Temp[1,column]})
  ## Delete medium condition
  Temp = Temp[-1,]
  ## Compute Z-scores
  # Temp = sapply(1:dim(Temp)[2], FUN = function(column){Temp[,column] = (Temp[,column] - mean(Temp[,column], na.rm=T))/sd(Temp[,column], na.rm=T)})
  ## name rows and columns
  rownames(Temp) = Stimuli
  colnames(Temp) = Phosphos
  return(Temp)
})
# ################# 
# ListPhospho_zs = lapply(DonorsWithPhospho, FUN = function(donor){
#   ## Load Data
#   Temp = readMIDAS(paste0("../../../data/phosphos_merged/", donor, '_phosphos_midas.csv'), verbose=F)$dataMatrix;
#   Temp = Temp[-1,]
#   Temp = as.matrix(Temp[,which(grepl('DV:', colnames(Temp)))]);
#   Temp = Temp[,-c(18,19)]
#   ## Log-transform data
#   Temp = log(Temp, base=2)
#   ## Compute fold change compared to the medium condition (first row)
#   Temp = sapply(1:dim(Temp)[2], FUN = function(column){Temp[,column] = Temp[,column] - Temp[1,column]})
#   ## Delete medium condition
#   Temp = Temp[-1,]
#   ## Compute Z-scores
#   Temp = sapply(1:dim(Temp)[2], FUN = function(column){Temp[,column] = (Temp[,column] - mean(Temp[,column], na.rm=T))/sd(Temp[,column], na.rm=T)})
#   ## name rows and columns
#   rownames(Temp) = Stimuli
#   colnames(Temp) = Phosphos
#   return(Temp)
# })

# # Load normalised phospho-data
# ListPhosphos_n = lapply(DonorsWithPhospho, FUN = function(donor){
#   ## Load Data
#   Temp = readMIDAS(paste("~/Documents/combiMS/data/phosphos_normalised/", donor, '.csv', sep=""), verbose=F)$dataMatrix;
#   Temp = as.matrix(Temp[,which(grepl('DV:', colnames(Temp)))]);
#   ## Delete medium condition
#   Temp = Temp[-1,]
#   ## name rows and columns
#   rownames(Temp) = Stimuli
#   colnames(Temp) = Phosphos
#   return(Temp)
# })
# 
names(ListPhosphos)  <- DonorsWithPhospho
# names(ListPhospho_zs)  <- DonorsWithPhospho
# names(ListPhosphos_n) <- DonorsWithPhospho

# #######################################################################################################################################
# # Function for wilcoxon-test
# #######################################################################################################################################
one_sample = function(DataList, annot, wilcoxon = TRUE){
  res = data.frame(row.names = c('a', 'b', 'group', 'logFC', 'p.value'))
  Temp = DataList[[1]]
  iter_1 = colnames(Temp)
  iter_2 = row.names(Temp)
  pb = txtProgressBar(min = 0, max = length(iter_1)*length(iter_2), initial = 0)
  for (i1 in iter_1){
    for (i2 in iter_2){
      if (i1 == i2){
        next()
      }
      Y = as.matrix(sapply(names(DataList), FUN=function(d, a, b){DataList[[d]][a, b]}, a=i2, b=i1))
      
      for (g in unique(annot$Group)){
        test = Y[row.names(annot)[which(annot$Group == g)],]
        test = test[!is.na(test)]
        if (length(test) == 0){next()}
        
        w = wilcox.test(test)
        res[[(length(res)+1)]] = c(i1, i2, g, median(test), w$p.value)
        
      }
      setTxtProgressBar(pb, (pb$getVal()+1))
    }
  }
  res = as.data.frame(t(res))
  res$logFC = as.numeric(as.character(res$logFC))
  res$p.value = as.numeric(as.character(res$p.value))
  res$adj.p.value = p.adjust(res$p.value, method='fdr')
  res$log.p.value = -log10(res$adj.p.value)
  res[['threshold']] = res$log.p.value > 1.3
  g = ggplot(res, aes(x=logFC, y=log.p.value, colour=threshold)) + 
    geom_point() + 
    theme_bw() + 
    facet_grid(~group) + 
    scale_colour_manual(values=c('#95a3a6', '#fd3c06')) +
    theme(legend.position='none')
  return(list(plot=g, data=res))
}

# #######################################################################################################################################
# # Run one-sample tests
# #######################################################################################################################################
g1 = one_sample(ListPhosphos, annotations)
# g1 = one_sample(ListPhospho_zs, annotations)
# g2 = one_sample(ListPhosphos_n, annotations)

# Compare with model interactions 
# So far its still manual and crude, but maybe worth looking into...?
temp = g1$data
temp$interaction_in_model = FALSE
for (g in names(model_interactions)){
  temp_model = model_interactions[[g]]
  for (i in 1:dim(temp_model)[1]){
    temp[which(temp$b == as.character(temp_model[i, 1]) & temp$a == as.character(temp_model[i, 2]) & temp$group == g),'interaction_in_model'] = TRUE
  }
}

temp$group = factor(temp$group, c('H', 'EGCG', 'FTY', 'IFNb', 'GAL', 'NTZ', 'U'))

g = ggplot(temp, aes(x=logFC, y=log.p.value, colour=interaction_in_model)) + 
  facet_grid(~group) + 
  geom_point(alpha=0.8) + 
  theme(legend.position='none') +
  scale_colour_manual(values=c('#d8dcd6', '#fd3c06')) +
  theme_bw()
print(g)

# ############################################################################################################################################
# Compute Confusion matrices for the different groups
# ############################################################################################################################################
results_all = data.frame(row.names=c('Sensitivity', 'Specificity', 'Group'))

for (cut_off in seq(.05, .95, 0.01)){
  for (g in unique(g1$data$group)){
    if (g == 'H'){next()}
    if (g == 'U'){next()}
    temp = g1$data[which(g1$data$group == g),]
    temp$interaction = FALSE
    temp_model = model_interactions[[g]]
    for (i in 1:dim(temp_model)[1]){
      temp[which(temp$b == as.character(temp_model[i, 1]) & temp$a == as.character(temp_model[i, 2]) & temp$group == g),'interaction'] = TRUE
    }
    temp$p.value.rank = rank(temp$adj.p.value, ties.method = 'average')
    temp$tp = temp$p.value.rank < dim(temp)[1] * cut_off
  
    true_positives = length(which(temp$interaction & temp$tp))
    false_positives = length(which(!temp$interaction & temp$tp))
    false_negative = length(which(temp$interaction & !temp$tp))
    true_negatives = length(which(!temp$interaction & !temp$tp))
  
  
    Sensitivity = true_positives/(true_positives+false_negative)
    Specificity = true_negatives/(true_negatives+false_positives)
    
    results_all[[(length(results_all)+1)]] = c(Sensitivity, 1-Specificity, g)
  }
}

res = as.data.frame(t(results_all))
res$Sensitivity = as.numeric(as.character(res$Sensitivity))
res$Specificity = as.numeric(as.character((res$Specificity)))
g = ggplot(res, aes(y=Sensitivity, x=Specificity, color=Group)) + geom_line() + xlim(0, 1) + ylim(0, 1) + geom_abline(slope = 1, intercept = 0, linetype=2) + theme_bw()
print(g)

cut_off=.2
for (g in unique(g1$data$group)){
  if (g == 'H'){next()}
  if (g == 'U'){next()}
  temp = g1$data[which(g1$data$group == g),]
  temp$interaction = FALSE
  temp_model = model_interactions[[g]]
  for (i in 1:dim(temp_model)[1]){
    temp[which(temp$b == as.character(temp_model[i, 1]) & temp$a == as.character(temp_model[i, 2]) & temp$group == g),'interaction'] = TRUE
  }
  temp$p.value.rank = rank(temp$adj.p.value, ties.method = 'average')
  temp$tp = temp$p.value.rank < dim(temp)[1] * cut_off
  
  true_positives = length(which(temp$interaction & temp$tp))
  false_positives = length(which(!temp$interaction & temp$tp))
  false_negative = length(which(temp$interaction & !temp$tp))
  true_negatives = length(which(!temp$interaction & !temp$tp))
  
  
  Sensitivity = true_positives/(true_positives+false_negative)
  Specificity = true_negatives/(true_negatives+false_positives)
  print(g)
  print(Sensitivity)
  print(Specificity)

}
