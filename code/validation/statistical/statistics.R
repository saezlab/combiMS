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

# #######################################################################################################################################
# # Logic Model Interactions
# #######################################################################################################################################

model_interactions = list(FTY = data.frame(t(matrix(data=c(c('antiCD3', 'STAT5'), c('EGCG', 'STAT5'), c('NaCl', 'MK12'), c('NaCl', 'HSPB1'), 
                                                           c('TNFA', 'MK12'), c('TNFA', 'HSPB1'), c('TNFA', 'IKBA'), c('TNFA', 'TF65'), c('TNFA', 'MKO3'),  
                                                           c('IL1A', 'MK12'), c('IL1A', 'HSPB1'), c('IL1A', 'IKBA'), c('IL1A', 'TF65'), c('IL1A', 'MKO3'), 
                                                           c('PolyIC', 'MK12'), c('PolyIC', 'HSPB1'), c('PolyIC', 'IKBA'), c('PolyIC', 'TF65'), c('PolyIC', 'MKO3'), 
                                                           c('LPS', 'MK12'), c('LPS', 'HSPB1'), c('LPS', 'IKBA'), c('LPS', 'TF65'), c('LPS', 'MKO3'),
                                                           c('INS', 'AKT1'), c('INS', 'MP2K1'),
                                                           c('rebif', 'STAT6')), nrow=2))), 
                          GAL = data.frame(t(matrix(data=c(c('TNFA', 'MK12'), c('TNFA', 'MKO3'), 
                                                           c('LPS', 'MKO3'), c('PolyIC', 'MKO3'), c('IL1A', 'MKO3'),
                                                           c('INS', 'AKT1'), c('BDNF', 'AKT1'),
                                                           c('conA', 'AKT1'), c('conA', 'STAT1'),
                                                           c('INFG', 'STAT1'), c('INFG', 'STAT3'),
                                                           c('rebif', 'STAT1'), c('rebif', 'STAT3'),
                                                           c('IL6', 'STAT1'), c('IL6', 'STAT3')), nrow=2))), 
                          IFNb = data.frame(t(matrix(data=c(c('NaCl', 'MP2KA'), c('NaCl', 'MK12'), c('NaCl', 'HSPB1'),
                                                           c('TNFA', 'MK12'), c('TNFA', 'HSPB1'), c('TNFA', 'MKO3'),
                                                           c('IL1A', 'MKO3'), c('LPS', 'MKO3'), c('PolyIC', 'MKO3'),
                                                           c('BDNF', 'AKT1'), c('INS', 'AKT1'),
                                                           c('antiCD3', 'STAT1'), c('rebif', 'STAT1'), c('INFG', 'STAT1'), c('IL6', 'STAT1')), nrow=2))),
                          EGCG = data.frame(t(matrix(data=c(c('NaCl', 'MK12'), c('NaCl', 'MKO3'),
                                                            c('TNFA', 'IKBA'), c('TNFA', 'TF65'), c('TNFA', 'MK12'), c('TNFA', 'MKO3'),
                                                            c('PolyIC', 'IKBA'), c('PolyIC', 'TF65'), c('PolyIC', 'MK12'), c('PolyIC', 'MKO3'),
                                                            c('LPS', 'IKBA'), c('LPS', 'TF65'), c('LPS', 'MK12'), c('LPS', 'MKO3'),
                                                            c('IL1A', 'IKBA'), c('IL1A', 'TF65'), c('IL1A', 'MK12'), c('IL1A', 'MKO3'),
                                                            c('conA', 'MK12'), c('conA', 'MKO3'),
                                                            c('BDNF', 'MK12'), c('BDNF', 'MKO3'),
                                                            c('IL6', 'MK12'), c('IL6', 'MKO3'), c('IL6', 'STAT3'),
                                                            c('EGCG', 'STAT1'), c('EGCG', 'STAT5'), c('EGCG', 'STAT3'),
                                                            c('antiCD3', 'STAT1'), c('antiCD3', 'STAT5'), c('antiCD3', 'STAT3')), nrow=2))),
                          NTZ = data.frame(t(matrix(data=c(c('NaCl', 'MK12'), c('NaCl', 'HSPB1'),
                                                           c('TNFA', 'MK12'), c('TNFA', 'HSPB1'), c('TNFA', 'MKO3'),
                                                           c('IL1A', 'MK12'), c('IL1A', 'HSPB1'), c('IL1A', 'MKO3'),
                                                           c('PolyIC', 'MK12'), c('PolyIC', 'HSPB1'), c('PolyIC', 'MKO3'),
                                                             c('LPS', 'MK12'), c('LPS', 'HSPB1'), c('LPS', 'MKO3'),
                                                           c('INS', 'AKT1'),
                                                           c('conA', 'STAT1'), c('antiCD3', 'STAT1'), c('rebif', 'STAT1'), c('INFG', 'STAT1'), c('IL6', 'STAT1')), nrow=2))))

# #######################################################################################################################################
# # Prepare Stuff
# #######################################################################################################################################

Path = '~/Dropbox/combiMS/'

## Which Donors are to be included? See Donor_Selection.R for more details, that must be run to generate the file Donor_IDs.txt
DonorsCytos = read.table(paste(Path, "Donor_IDs.txt", sep=""))
DonorsCytos = as.vector(unlist(DonorsCytos))

## Not all donors have phospho-Data available
DonorsWithPhospho = list.files(paste(Path, 'phosphos/', sep=''))
DonorsWithPhospho = gsub(".csv", '', DonorsWithPhospho)

## Annotations
annotations = read.table(paste0(Path, 'annotation.csv'), sep=',', header = T, row.names=1)
annot_c = annotations[DonorsCytos, ]
annot_p = annotations[DonorsWithPhospho, ]
rm(annotations)

# Hack for nicer annotation
annot_p$Group = as.character(annot_p$Group)
group_temp = list('GAL', 'H', 'U', 'U', 'IFNb', 'FTY', 'EGCG', 'NTZ')
names(group_temp) = c("Glatiramer", "Healthy", "Untreated", "PPMS", "Interferon B", "Fingolimod", "EGCG", "Natalizumab")

for (i in 1:dim(annot_p)[1]){
  annot_p[i, 'Group'] = group_temp[[annot_p[i, 'Group']]]
}

# #######################################################################################################################################
# # Prepare Metadata, such as Stimuli and Cytokine/Phospho-information
# #######################################################################################################################################
Temp = readMIDAS(paste(Path, 'midas_raw_log/CH003_cytos_midas.csv', sep=""), verbose=F)$dataMatrix
Cytokines = as.vector(sapply(colnames(Temp)[which(grepl('DV', colnames(Temp)))], FUN = function(x){strsplit(x, split=':')[[1]][2]}))
Stimuli = as.vector(sapply(colnames(Temp)[which(grepl('TR', colnames(Temp)))], FUN = function(x){strsplit(x, split=':')[[1]][2]}))[2:21]
Temp = readMIDAS(paste(Path, 'phosphos/CH003.csv', sep=""), verbose=F)$dataMatrix
Phosphos = as.vector(sapply(colnames(Temp)[which(grepl('DV', colnames(Temp)))], FUN = function(x){strsplit(x, split=':')[[1]][2]}))
rm(Temp)

# ######################################################################################################################################
# # Load Phospho data, log2(FC)
# ######################################################################################################################################
ListPhosphos = lapply(DonorsWithPhospho, FUN = function(donor){
  ## Load Data
  Temp = readMIDAS(paste(Path, "phosphos/", list.files(paste(Path, "phosphos/", sep=""))[grep(donor, list.files(paste(Path, "phosphos/", sep="")))], sep=""), verbose=F)$dataMatrix;
  Temp = as.matrix(Temp[,which(grepl('DV:', colnames(Temp)))]);
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

# Load normalised phospho-data
ListPhosphos_n = lapply(DonorsWithPhospho, FUN = function(donor){
  ## Load Data
  Temp = readMIDAS(paste("~/Documents/combiMS/data/phosphos_normalised/", donor, '.csv', sep=""), verbose=F)$dataMatrix;
  Temp = as.matrix(Temp[,which(grepl('DV:', colnames(Temp)))]);
  ## Delete medium condition
  Temp = Temp[-1,]
  ## name rows and columns
  rownames(Temp) = Stimuli
  colnames(Temp) = Phosphos
  return(Temp)
})

names(ListPhosphos) <- names(ListPhosphos_n) <- DonorsWithPhospho

# #######################################################################################################################################
# # Function for wilcoxon-test
# #######################################################################################################################################
one_sample = function(DataList, annot, wilcoxon = TRUE){
  res = data.frame(row.names = c('a', 'b', 'group', 'logFC', 'p.value'))
  Temp = DataList[[1]]
  iter_1 = colnames(Temp)
  iter_2 = row.names(Temp)
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
        if (wilcoxon){
          w = wilcox.test(test)
          res[[(length(res)+1)]] = c(i1, i2, g, median(test), w$p.value)
        } else {
          w = t.test(test)
          res[[(length(res)+1)]] = c(i1, i2, g, w$estimate, w$p.value)
        }
        
      }
    }
  }
  res = as.data.frame(t(res))
  res$logFC = as.numeric(as.character(res$logFC))
  res$p.value = as.numeric(as.character(res$p.value))
  res$adj.p.value = res$p.value
  for (g in unique(annot$Group)){
    res$adj.p.value[which(res$group == g)] = p.adjust(res$adj.p.value[which(res$group == g)], method='BH')
  }
  res$log.p.value = -log10(res$adj.p.value)
  res[['threshold']] = res$log.p.value > 2
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
g1 = one_sample(ListPhosphos, annot_p)
g2 = one_sample(ListPhosphos, annot_p, wilcoxon = FALSE)
g3 = one_sample(ListPhosphos_n, annot_p)
g4 = one_sample(ListPhosphos_n, annot_p, wilcoxon = FALSE)


# Compare with model interactions 
# So far its still manual and crude, but maybe worth looking into...?
temp = g2$data
temp$interaction = FALSE
for (g in names(model_interactions)){
  temp_model = model_interactions[[g]]
  for (i in 1:dim(temp_model)[1]){
    temp[which(temp$b == as.character(temp_model[i, 1]) & temp$a == as.character(temp_model[i, 2]) & temp$group == g),'interaction'] = TRUE
  }
}
g = ggplot(temp, aes(x=logFC, y=log.p.value, colour=interaction)) + 
  facet_grid(~group) + 
  geom_point() + 
  theme(legend.position='none') +
  scale_colour_manual(values=c('#d8dcd6', '#fd3c06')) +
  theme_bw()
print(g)

# #######################################################################################################################################
# # Function for Marti's approach (email)
# #######################################################################################################################################
#
#
#
# About the statistical test to validate group models - not prediction of combination therapies, importantly 
# The question is to find which phosphos are active in certain treatment groups, and then look if they fall in active pathways of the merged group model. 
# Therefore, this poses a second question: active with respect to what? I think one option would be against both the control and the ms untreated state. 
# So this could be a simple way:
#   • First test: which prots are diffferentially phosphorylated between Ms untreated and control? From the result, take those that fail the test, 
#     i.e. those that are in a similar state (rather, not significantly different).
#   • For the subset of the prots above, merge the control and Ms untreated values, or take the lowest. Then second-test this against MS treated 
#     for each in vivo group sepparately (5 tests). 
#   • For those that are found, which ones are active in MS treated? 
#   • Do they fall on an active pathway according to the models, thereby validating our pathways predicted to be active upon each treatment?
#
#
#
#######################################################################################################################################
stat_comparison = function(DataList, annot, wilcoxon = TRUE){
  res = data.frame(row.names = c('a', 'b', 'logFC', 'p.value'))
  Temp = DataList[[1]]
  iter_1 = colnames(Temp)
  iter_2 = row.names(Temp)
  for (i1 in iter_1){
    for (i2 in iter_2){
      if (i1 == i2){
        next()
      }
      Y = as.matrix(sapply(names(DataList), FUN=function(d, a, b){DataList[[d]][a, b]}, a=i2, b=i1))
      
      healthy = Y[row.names(annot)[which(annot$Group == 'H')],]
      healthy = healthy[!is.na(healthy)]
      
      untreated = Y[row.names(annot)[which(annot$Group == 'U')],]
      untreated = untreated[!is.na(untreated)]
      
      if (wilcoxon){
        w = wilcox.test(healthy, untreated, conf.int = TRUE)
        res[[(length(res)+1)]] = c(i1, i2, w$estimate, w$p.value)
      } else {
        w = t.test(healthy, untreated)
        res[[(length(res)+1)]] = c(i1, i2, w$estimate[1] - w$estimate[2], w$p.value)
      }
      
    }
  }
  res = as.data.frame(t(res))
  res$logFC = as.numeric(as.character(res$logFC))
  res$p.value = as.numeric(as.character(res$p.value))
  res$adj.p.value = res$p.value
  res$log.p.value = -log10(res$adj.p.value)
  res[['threshold']] = res$log.p.value > 2
  g = ggplot(res, aes(x=logFC, y=log.p.value, colour=threshold)) + 
    geom_point() + 
    theme_bw() + 
    # facet_grid(~group) + 
    scale_colour_manual(values=c('#95a3a6', '#fd3c06')) +
    theme(legend.position='none')
  return(list(plot=g, data=res))
}

######################################################################################################################################
stats_1 = stat_comparison(ListPhosphos, annot_p)
stats_2 = stat_comparison(ListPhosphos, annot_p, wilcoxon = FALSE)
stats_3 = stat_comparison(ListPhosphos_n, annot_p)
stats_4 = stat_comparison(ListPhosphos_n, annot_p, wilcoxon = FALSE)

# which are differentially phosphorylated between MS untreated and control
print(stats_1$data[which(stats_1$data[['threshold']] == TRUE),]) # AKT1-vitD3, CREB-IFNG, CREB1-EGCG, STAT6-vitD3
print(stats_2$data[which(stats_2$data[['threshold']] == TRUE),]) # CREB1-EGCG, STAT6-vitD3
print(stats_3$data[which(stats_3$data[['threshold']] == TRUE),]) # WNK1-BDNF, STAT3-INFG
print(stats_4$data[which(stats_4$data[['threshold']] == TRUE),]) # WNK1-BDNF, STAT3-INFG, STAT3-INS

# These are not so many --> worth investigating in more detail?
######################################################################################################################################
stat_comparison_second = function(DataList, annot, wilcoxon = TRUE){
  res = data.frame(row.names = c('a', 'b', 'group','logFC', 'p.value'))
  Temp = DataList[[1]]
  iter_1 = colnames(Temp)
  iter_2 = row.names(Temp)
  for (i1 in iter_1){
    for (i2 in iter_2){
      if (i1 == i2){
        next()
      }
      Y = as.matrix(sapply(names(DataList), FUN=function(d, a, b){DataList[[d]][a, b]}, a=i2, b=i1))
      
      control = Y[row.names(annot)[which(annot$Group == 'H' | annot$Group == 'U')],]
      control = control[!is.na(control)]
      
      for (g in names(model_interactions)){
        test = Y[row.names(annot)[which(annot$Group == g)],]
        test = test[!is.na(test)]
        if (wilcoxon){
          w = wilcox.test(control, test, conf.int = TRUE)
          res[[(length(res)+1)]] = c(i1, i2, g, w$estimate, w$p.value)
        } else {
          w = t.test(control, test)
          res[[(length(res)+1)]] = c(i1, i2, g, w$estimate[1] - w$estimate[2], w$p.value)
        }
      }
    }
  }
  res = as.data.frame(t(res))
  res$logFC = as.numeric(as.character(res$logFC))
  res$p.value = as.numeric(as.character(res$p.value))
  res$adj.p.value = res$p.value
  res$log.p.value = -log10(res$adj.p.value)
  res[['threshold']] = res$log.p.value > 2
  g = ggplot(res, aes(x=logFC, y=log.p.value, colour=threshold)) + 
    geom_point() + 
    theme_bw() + 
    facet_grid(~group) + 
    scale_colour_manual(values=c('#95a3a6', '#fd3c06')) +
    theme(legend.position='none')
  return(list(plot=g, data=res))
}

stats_1_2 = stat_comparison_second(ListPhosphos, annot_p)
stats_2_2 = stat_comparison_second(ListPhosphos, annot_p, wilcoxon = FALSE)
stats_3_2 = stat_comparison_second(ListPhosphos_n, annot_p)
stats_4_2 = stat_comparison_second(ListPhosphos_n, annot_p, wilcoxon = FALSE)

# Compare with model interactions 
# So far its still manual and crude, but maybe worth looking into...?
temp = stats_4_2$data
temp$interaction = FALSE
for (g in names(model_interactions)){
  temp_model = model_interactions[[g]]
  for (i in 1:dim(temp_model)[1]){
    temp[which(temp$b == as.character(temp_model[i, 1]) & temp$a == as.character(temp_model[i, 2]) & temp$group == g),'interaction'] = TRUE
  }
}
g = ggplot(temp, aes(x=logFC, y=log.p.value, colour=interaction)) + 
  facet_grid(~group) + 
  geom_point() + 
  theme(legend.position='none') +
  scale_colour_manual(values=c('#d8dcd6', '#fd3c06')) +
  theme_bw()
print(g)
