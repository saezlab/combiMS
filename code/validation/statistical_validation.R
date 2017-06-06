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
# # Logic Model Interactions, manually extracted from the optimized networks in Figure 5 of the main manuscript
# #######################################################################################################################################

# model_interactions_m = list(FTY = data.frame(t(matrix(data=c(c('ANTICD3', 'STAT5'), c('EGCG', 'STAT5'), c('NACL', 'MK12'), c('NACL', 'HSPB1'),
#                                                            c('TNFA', 'MK12'), c('TNFA', 'HSPB1'), c('TNFA', 'IKBA'), c('TNFA', 'TF65'),
#                                                            c('TNFA', 'MKO3'),  c('IL1A', 'MK12'), c('IL1A', 'HSPB1'), c('IL1A', 'IKBA'),
#                                                            c('IL1A', 'TF65'), c('IL1A', 'MKO3'), c('POLYIC', 'MK12'), c('POLYIC', 'HSPB1'),
#                                                            c('POLYIC', 'IKBA'), c('POLYIC', 'TF65'), c('POLYIC', 'MKO3'), c('LPS', 'MK12'),
#                                                            c('LPS', 'HSPB1'), c('LPS', 'IKBA'), c('LPS', 'TF65'), c('LPS', 'MKO3'),
#                                                            c('INS', 'AKT1'), c('INS', 'MP2K1'), c('REBIF', 'STAT6')), nrow=2))),
#                           GA = data.frame(t(matrix(data=c(c('TNFA', 'MK12'), c('TNFA', 'MKO3'), c('LPS', 'MKO3'), c('POLYIC', 'MKO3'),
#                                                           c('IL1A', 'MKO3'), c('INS', 'AKT1'), c('BDNF', 'AKT1'), c('CONA', 'AKT1'),
#                                                           c('CONA', 'STAT1'), c('INFG', 'STAT1'), c('INFG', 'STAT3'), c('REBIF', 'STAT1'),
#                                                           c('REBIF', 'STAT3'), c('IL6', 'STAT1'), c('IL6', 'STAT3')), nrow=2))),
#                           IFNb = data.frame(t(matrix(data=c(c('NACL', 'MP2K1'), c('NACL', 'MK12'), c('NACL', 'HSPB1'), c('TNFA', 'MK12'),
#                                                             c('TNFA', 'HSPB1'), c('TNFA', 'MKO3'), c('IL1A', 'MKO3'), c('LPS', 'MKO3'),
#                                                             c('POLYIC', 'MKO3'), c('BDNF', 'AKT1'), c('INS', 'AKT1'), c('ANTICD3', 'STAT1'),
#                                                             c('REBIF', 'STAT1'), c('INFG', 'STAT1'), c('IL6', 'STAT1')), nrow=2))),
#                           EGCG = data.frame(t(matrix(data=c(c('NACL', 'MK12'), c('NACL', 'MKO3'),
#                                                             c('TNFA', 'IKBA'), c('TNFA', 'TF65'), c('TNFA', 'MK12'), c('TNFA', 'MKO3'),
#                                                             c('POLYIC', 'IKBA'), c('POLYIC', 'TF65'), c('POLYIC', 'MK12'), c('POLYIC', 'MKO3'),
#                                                             c('LPS', 'IKBA'), c('LPS', 'TF65'), c('LPS', 'MK12'), c('LPS', 'MKO3'),
#                                                             c('IL1A', 'IKBA'), c('IL1A', 'TF65'), c('IL1A', 'MK12'), c('IL1A', 'MKO3'),
#                                                             c('CONA', 'MK12'), c('CONA', 'MKO3'),
#                                                             c('BDNF', 'MK12'), c('BDNF', 'MKO3'),
#                                                             c('IL6', 'MK12'), c('IL6', 'MKO3'), c('IL6', 'STAT3'),
#                                                             c('EGCG', 'STAT1'), c('EGCG', 'STAT5'), c('EGCG', 'STAT3'),
#                                                             c('ANTICD3', 'STAT1'), c('ANTICD3', 'STAT5'), c('ANTICD3', 'STAT3')), nrow=2))),
#                           NTZ = data.frame(t(matrix(data=c(c('NACL', 'MK12'), c('NACL', 'HSPB1'),
#                                                            c('TNFA', 'MK12'), c('TNFA', 'HSPB1'), c('TNFA', 'MKO3'),
#                                                            c('IL1A', 'MK12'), c('IL1A', 'HSPB1'), c('IL1A', 'MKO3'),
#                                                            c('POLYIC', 'MK12'), c('POLYIC', 'HSPB1'), c('POLYIC', 'MKO3'),
#                                                            c('LPS', 'MK12'), c('LPS', 'HSPB1'), c('LPS', 'MKO3'),
#                                                            c('INS', 'AKT1'),
#                                                            c('CONA', 'STAT1'), c('ANTICD3', 'STAT1'), c('REBIF', 'STAT1'), c('INFG', 'STAT1'), c('IL6', 'STAT1')), nrow=2))))

# #######################################################################################################################################
# # Prepare Stuff
# #######################################################################################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Not all donors have phospho-Data available
DonorsWithPhospho = list.files('../../data/phosphos_normalised/')
DonorsWithPhospho = gsub(".csv", '', DonorsWithPhospho)

## Annotations
annotations = read.table('../../files/annotations/annotations_169_patients.csv', header = T, row.names=1, sep=',')

# Hack for nicer annotation
annotations$Group = as.character(annotations$Group)
group_temp = list('GA', 'H', 'U', 'U', 'IFNb', 'FTY', 'EGCG', 'NTZ')
names(group_temp) = c("Glatiramer", "Healthy", "Untreated", "PPMS", "Interferon B", "Fingolimod", "EGCG", "Natalizumab")

for (i in 1:dim(annotations)[1]){
  annotations[i, 'Group'] = group_temp[[annotations[i, 'Group']]]
}

# #######################################################################################################################################
# # Prepare Metadata, such as Stimuli and Phospho-information
# #######################################################################################################################################

Temp = readMIDAS('../../data/phosphos_merged/CH003_phosphos_midas.csv', verbose=F)$dataMatrix
Phosphos = as.vector(sapply(colnames(Temp)[which(grepl('DV', colnames(Temp)))], FUN = function(x){strsplit(x, split=':')[[1]][2]}))[1:17]
Stimuli = toupper(as.vector(sapply(colnames(Temp)[which(grepl('TR', colnames(Temp)))], FUN = function(x){strsplit(x, split=':')[[1]][2]}))[2:21])
rm(Temp)

# ######################################################################################################################################
# # Load Phospho data, log2(FC)
# ######################################################################################################################################
ListPhosphos = lapply(DonorsWithPhospho, FUN = function(donor){
  ## Load Data
  Temp = readMIDAS(paste0("../../data/phosphos_merged/", donor, '_phosphos_midas.csv'), verbose=F)$dataMatrix;
  Temp = Temp[-1,]
  Temp = as.matrix(Temp[,which(grepl('DV:', colnames(Temp)))]);
  Temp = Temp[,-c(18,19)]
  ## Log-transform data
  Temp = log(Temp, base=2)
  ## Compute fold change compared to the medium condition (first row)
  Temp = sapply(1:dim(Temp)[2], FUN = function(column){Temp[,column] = Temp[,column] - Temp[1,column]})
  ## Delete medium condition
  Temp = Temp[-1,]
  ## name rows and columns
  rownames(Temp) = Stimuli
  colnames(Temp) = Phosphos
  return(Temp)
})

names(ListPhosphos)  <- DonorsWithPhospho

# #####################################################################################################################################
# Convert pkn into igraph object and find all possible Stimulus-Readout pairs
# #####################################################################################################################################

pknmodel = readSIF('../../files/model/combiMSplaneCUT.sif')
cnolist = readMIDAS('../../data/phosphos_normalised/CH003.csv')
cnolist = makeCNOlist(cnolist, subfield = FALSE)

# #######################################################################################################################################
# # Logic Model Interactions, automatically extracted from the optimized networks
# #######################################################################################################################################

find_model_interactions = function(group){
  name_conversion = list('FTY'='Gilenya',
                         "GA"='Copaxone',
                         "IFNb"='IFNb',
                         "EGCG"='EGCG',
                         "NTZ"='Tysabri')
  group = name_conversion[[group]]
  model_group = read.csv(paste0('../../files/group_models/', group, 'median.csv'))
  edges = pknmodel$interMat[,row.names(model_group)[which(model_group$x == 1)]]
  adjacency = matrix(0, nrow(edges), nrow(edges))
  for(i in 1:ncol(edges)){
    startIdx = which(edges[, i]==-1)
    endIdx = which(edges[, i]==1)

    adjacency[startIdx, endIdx] = 1
  }
  rownames(adjacency) <- colnames(adjacency) <- rownames(edges)

  netGraph <- graph.adjacency(adjacency, mode = "directed")

  connectivity = matrix(data=NA, ncol=length(cnolist$namesSignals), nrow=length(cnolist$namesStimuli))
  row.names(connectivity) = toupper(cnolist$namesStimuli)
  colnames(connectivity) = cnolist$namesSignals

  for (stim in toupper(cnolist$namesStimuli)){
    for (sig in cnolist$namesSignals){
      if (length(all_shortest_paths(netGraph, from=stim, to=sig)$res) == 0){
        connectivity[which(row.names(connectivity) == stim), which(colnames(connectivity) == sig)] = 0
      } else {
        connectivity[which(row.names(connectivity) == stim), which(colnames(connectivity) == sig)] = length(all_shortest_paths(netGraph, from=stim, to=sig)$res[[1]])
      }
    }
  }
    
  possible_interactions = melt(connectivity)
  possible_interactions[['interaction']] = paste(possible_interactions$Var1, possible_interactions$Var2, sep='_')
  model_int = possible_interactions[which(possible_interactions$value != 0),c(1,2)]
  return(model_int)
}

model_interactions = list()
for (group in c('FTY', 'NTZ', 'GA', 'IFNb', 'EGCG')){
  model_interactions[[group]] = find_model_interactions(group)
}

# convert to igraph object
edges <- pknmodel$interMat
adjacency = matrix(0, nrow(edges), nrow(edges))
for(i in 1:ncol(edges)){   
  startIdx = which(edges[, i]==-1)
  endIdx = which(edges[, i]==1)
  
  adjacency[startIdx, endIdx] = 1
}
rownames(adjacency) <- colnames(adjacency) <- rownames(edges)

netGraph <- graph.adjacency(adjacency, mode = "directed")

connectivity = matrix(data=NA, ncol=length(cnolist$namesSignals), nrow=length(cnolist$namesStimuli))
row.names(connectivity) = toupper(cnolist$namesStimuli)
colnames(connectivity) = cnolist$namesSignals

for (stim in toupper(cnolist$namesStimuli)){
  for (sig in cnolist$namesSignals){
    if (length(all_shortest_paths(netGraph, from=stim, to=sig)$res) == 0){
      connectivity[which(row.names(connectivity) == stim), which(colnames(connectivity) == sig)] = 0
    } else {
      connectivity[which(row.names(connectivity) == stim), which(colnames(connectivity) == sig)] = length(all_shortest_paths(netGraph, from=stim, to=sig)$res[[1]])
    }
  }
}
colnames(connectivity)[13] = 'MK03'
possible_interactions = melt(connectivity)
possible_interactions[['interaction']] = paste(possible_interactions$Var1, possible_interactions$Var2, sep='_')
impossible_interactions = possible_interactions[which(possible_interactions$value == 0),c(1,2)]

# #######################################################################################################################################
# # Function for wilcoxon-test
# #######################################################################################################################################

res = data.frame(row.names = c('a', 'b', 'group', 'logFC', 'p.value'))
pb = txtProgressBar(min = 0, max = length(Stimuli)*length(Phosphos), initial = 0)
for (stim in Stimuli){
  for (phos in Phosphos){
    if (stim == phos){
      next()
    }
    if (length(which(impossible_interactions$Var1 == stim & impossible_interactions$Var2 == phos)) != 0){
      next()
    }
    Y = as.matrix(sapply(names(ListPhosphos), FUN=function(d, a, b){ListPhosphos[[d]][a, b]}, b=phos, a=stim))
      
    for (g in unique(annotations$Group)){
      test = Y[row.names(annotations)[which(annotations$Group == g)],]
      test = test[!is.na(test)]
      
      if (length(test) == 0){next()}
      
      w = wilcox.test(test)
      res[[(length(res)+1)]] = c(stim, phos, g, median(test), w$p.value)
    }
  setTxtProgressBar(pb, (pb$getVal()+1))
  }
}
res = as.data.frame(t(res))

res$logFC = as.numeric(as.character(res$logFC))
res$p.value = as.numeric(as.character(res$p.value))

res$adj.p.value = res$p.value
for (g in unique(annotations$Group)){
  res[which(res$group == g),'adj.p.value'] = p.adjust(res[which(res$group == g),'adj.p.value'], method='fdr')
}

res$log.p.value = -log10(res$adj.p.value)
res[['threshold']] = res$log.p.value > 1.3

res$interaction_in_model = FALSE
res$players_in_model = FALSE
for (g in names(model_interactions)){
  temp_model = model_interactions[[g]]
  # record which reactions contain the players of the model
  for (stim in unique(temp_model[,1])){
    for (phos in unique(temp_model[,2])){
      res[which((res$a == as.character(stim)) & ((res$b == as.character(phos)) | (res$group == g))), 'players_in_model'] = TRUE
    }
  }
  # record which reactions could be in the model
  for (i in 1:dim(temp_model)[1]){
    res[which(res$a == as.character(temp_model[i, 1]) & res$b == as.character(temp_model[i, 2]) & res$group == g),'interaction_in_model'] = TRUE
  }
}

res$group = factor(res$group, c('H', 'EGCG', 'FTY', 'IFNb', 'GA', 'NTZ', 'U'))


# #######################################################################################################################################
# Fisher's exact test to test the ratio between statistically significant yes/no and interaction in model yes/no
# #######################################################################################################################################

# split by group
# IFN = res[which((res$group == 'IFNb')),]
# EGCG = res[which((res$group == 'EGCG')),]
# FTY = res[which((res$group == 'FTY')),]
# GAL = res[which((res$group == 'GA')),]
# NTZ = res[which((res$group == 'NTZ')),]
IFN = res[which((res$group == 'IFNb') & res$players_in_model == TRUE),]
EGCG = res[which((res$group == 'EGCG') & res$players_in_model == TRUE),]
FTY = res[which((res$group == 'FTY') & res$players_in_model == TRUE),]
GAL = res[which((res$group == 'GA') & res$players_in_model == TRUE),]
NTZ = res[which((res$group == 'NTZ') & res$players_in_model == TRUE),]


test_with_fisher = function(df){
  matrix_test = matrix(c(length(which((df$threshold == TRUE) & df$interaction_in_model == TRUE)),
                         length(which((df$threshold == FALSE) & df$interaction_in_model == TRUE)),
                         length(which((df$threshold == TRUE) & df$interaction_in_model == FALSE)),
                         length(which((df$threshold == FALSE) & df$interaction_in_model == FALSE))), nrow=2)
  fisher_res = fisher.test(matrix_test)
  print(matrix_test)
  return(fisher_res)
}

test_with_fisher(IFN)
test_with_fisher(GAL)
test_with_fisher(NTZ)
test_with_fisher(FTY)
test_with_fisher(EGCG)


all_groups = rbind(IFN, EGCG, FTY, GAL, NTZ)

# #############################################################################################################################################
# nice Volcano plots
# ############################################################################################################################################

all_groups$group = factor(all_groups$group, c('IFNb', 'NTZ', 'FTY', 'GA', 'EGCG'))  

g = ggplot(all_groups, aes(x=logFC, y=log.p.value, color=interaction_in_model, fill=interaction_in_model, shape=factor(group))) + 
  geom_point(alpha=0.8, stroke=0.5) + 
  geom_point(data=all_groups[which(all_groups$interaction_in_model),], size=2.5, show.legend = TRUE, stroke=0.5) + 
  xlab(expression('log'['2']*'(FC)')) + ylab(expression('-log'['10']*'(Adj. p-value)')) +
  xlim(-3,3) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype='dotted', alpha=.8) +
  geom_hline(yintercept = 2, linetype='dotted', alpha=.8) +
  geom_hline(yintercept = 1.3, linetype='dotted', alpha=.5) +
  scale_shape_manual('Subgroup', values=c(21, 22, 23, 24, 25)) +
  scale_fill_manual(values=c('#8EBAE5', '#D85C41'), guide=FALSE) +
  scale_color_manual('Interaction', values=c('#8EBAE5', '#D85C41'), labels=c('not in\nmodel', 'in model'))

plot(g)

ann_text <- data.frame(log.p.value = c(2.95, 2.25, 1.15, 1, .5), y = c(.3, .3, .5, 2, 2), lab = c("INFb", 'NTZ','FTY', 'GA', 'EGCG'), 
                       interaction_in_model=c('none', 'none', 'none', 'none', 'none'),
                       group = factor(c('IFNb', 'NTZ', 'FTY', 'GA', 'EGCG'), levels = c('IFNb', 'NTZ', 'FTY', 'GA', 'EGCG')))


g1 = ggplot(all_groups, aes(x=log.p.value, colour=interaction_in_model, fill=interaction_in_model)) + 
  geom_density() + 
  facet_wrap(~group, ncol=1, scales = "free") +
  theme_classic() + theme(legend.position='none', axis.text=element_text(size=6)) + 
  xlab(expression('-log'['10']*'(Adj. p-value)')) + ylab('Density') + 
  scale_color_manual('Interaction', values=c('#8EBAE5','black', '#D85C41'), labels=c('not in\nmodel', 'in model', 'none')) +
  theme(strip.background=element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual('Interaction', values=alpha(c( '#8EBAE5', 'black', '#D85C41'),0.3), labels=c('not in\nmodel', 'in model', 'none')) +
  geom_text(data=ann_text, y=ann_text$y, label=ann_text$lab)

plot(g1)

save(all_groups, file='../../files/statistics/all_groups_wilcoxon_results.RData')
