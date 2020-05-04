# ###############################################################################################################################################
# Based on script: analyzePatientSimilarity.R by Marti July 2015
#
# Minor modification for nicer plotting by Jakob in December 2016
# New figure 3B by Jakob and Marti March 2020
# #################################################################################################################################################

library(apcluster)
library(gclus)
library(ggplot2)
library(reshape2)
library(NMF)
library(corrplot)
library(cowplot)
library(gdata)
library(ggsignif)
  
# Use relative paths instead of absolut paths for the files
# in Rstudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("~/GITtrunk/combiMS/code/figures/")
# ###########################################################################################################
# Load and clean annotation data
# ###########################################################################################################
annot=read.csv("../../files/annotations/annot169pat_v2.csv",header=TRUE,dec=".",check.names=TRUE, stringsAsFactors=FALSE)

# ************************************ 
# add some fields to annotation for ggplot
annot2=annot
# correct that PPMS patients appear in group instead of by untreated and treatment
# Rename groups for consistency with other figures in the paper
annot2[which(annot2$Group=='PPMS' & annot2$Treatment=='no'),'Group']='Untreated'
annot2[which(annot2$Group=='PPMS' & annot2$Treatment!='no'),'Group']=annot2[which(annot2$Group=='PPMS' & annot2$Treatment!='no'),'Treatment']
annot2[which(annot2$Group=='Gilenya'),'Group']='FTY'
annot2[which(annot2$Group=='Fingolimod'),'Group']='FTY'
annot2[which(annot2$Group=='Glatiramer'),'Group']='GA'
annot2[which(annot2$Group=='Interferon B'),'Group']='IFNb'
annot2[which(annot2$Group=='Natalizumab'),'Group']='NTZ'
annot2[which(annot2$Disease.Subtype == ''), 'Disease.Subtype'] = 'healthy'
annot2[which(annot2$Disease.Subtype == 'SP'), 'Disease.Subtype'] = substring(annot2[which(annot2$Disease.Subtype == 'SP'), 'Category'], 1, 2)
annot2[which(annot2$Disease.Subtype == 'PR'), 'Disease.Subtype'] = substring(annot2[which(annot2$Disease.Subtype == 'PR'), 'Category'], 1, 2)
annot=annot2
rm(annot2)


###################################
### FIGURE 3B New Approach!
###################################

load('../../files/median_models/allMeanModels.RData')
load('../../files/median_models/allMedianModels.RData')

# We calculate the group model as the mean of the median models of each patient
# Then, we calculate the non-binary Jaccard distance between each
# patient model and the mean group model


library(tidyverse)
library(vegan)
library(matrixStats)


#Calculate similarity between all donors within a given group and the 
#MEAN group network. This includes the untreated ms patients, whose distance
#is calculataed against the given subgroup, not the untreated mean network


df.annot <- annot %>% 
  as_tibble() %>% 
  mutate(Group2=case_when(Group=='Untreated'~paste0(Group, Disease.Subtype),
                          TRUE~Group)) %>% 
  filter(Group2!='Tysabri or Placebo')

df.plot <- tibble(group=character(0), jacc=double(0), type=character(0),
                  id.group=character(0), id=character(0))
for (x in unique(df.annot$Group)){
  group.model <- colMeans(allMedianNetworks[df.annot %>% 
                                              filter(Group==x) %>% 
                                              pull(ID),])
  temp <- rbind(allMeanNetworks, group.model)
  dis <- vegdist(temp, method='jaccard', binary = FALSE)
  
  df.plot <- bind_rows(df.plot,
                       enframe(as.matrix(dis)['group.model',], 
                               value='jacc', name='ID') %>% 
                         right_join(df.annot %>% select(ID, Group)) %>% 
                         mutate(type=case_when(Group==x~'Same', 
                                               TRUE~'Different')) %>% 
                         transmute(group=x, jacc=jacc, type=type,
                                   id=ID, id.group=Group))
}

### Discarded attempt for figure 3B where we compared subsets of the data to untreated 

g <- df.plot %>% 
  ggplot(aes(x=group, y=jacc, fill=type)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(col=type), 
              position = position_jitterdodge(jitter.width = 0.07)) + 
  scale_fill_manual(values=c('#009F4D75', '#FFA30075')) + 
  scale_colour_manual(values=c('#009F4D', '#FFA300')) + 
  theme_bw()
ggsave(g, filename = '~/Desktop/sim_test_1.pdf',
       width = 6, height = 4, useDingbats=FALSE)

# new Figure 3

library(ggpubr) #necessary for p-value

# calculate within distances, i.e. all patients against mean group models
x <- df.plot %>% 
  filter(group%in%c('EGCG', 'FTY', "GA", 'IFNb', 'NTZ', "Healthy")) %>%
  filter(group==id.group | id.group=='Untreated')

# calculate all distances between all patients. This will be the grey "all" bar
dist.mat <- vegdist(allMeanNetworks, method='jaccard', binary=FALSE)
temp <- as.matrix(dist.mat)
diag(temp) <- NA
temp[upper.tri(temp)] <- NA

# transform "all" distances into a tibble so that it can be included to the 
# the tible we already have
#after creating y, confirm that x and y look the same so that they can be joined togeher

y <- temp %>% 
  as_tibble(rownames = 'id') %>% 
  pivot_longer(-id, names_to = 'id2', values_to = 'jacc') %>% 
  filter(!is.na(jacc)) %>% 
  transmute(id=id, jacc=jacc, type='All', group='All', id.group='All')


symnum.args <-list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))
g <- bind_rows(x,y) %>% 
  mutate(type2=case_when(group=='Healthy'&type=='Same'~'Healthy',
                         type=='Same'~'Treatment', 
                         type=="Different"~'Untreated MS', TRUE~'All')) %>% 
  mutate(type2=factor(type2, levels = c('All', "Healthy",
                                        "Treatment", "Untreated MS"))) %>% 
  mutate(group=factor(group, levels = c('Healthy', "All", 
                                        'EGCG', 'FTY', 'GA', 
                                        'NTZ', 'IFNb'))) %>%
  ggplot(aes(x=group, y=jacc, fill=type2)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes(col=type2), 
              position = position_jitterdodge(jitter.width = 0.07)) + 
    scale_fill_manual(values=c('#A9A9A9BF', "#307FE275", 
                               '#009F4D75', '#FFA30075'), 
                      name='Donor') + 
    scale_colour_manual(values=c('darkgrey', "#307FE2", 
                                 '#009F4D', '#FFA300'), 
                        name='Donor') + 
    theme_bw() + ylab('Jaccard distance') + 
    xlab('') + 
    theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line()) + 
    stat_compare_means(method='wilcox.test', 
                       aes(label = ..p.signif..))
ggsave(g, filename = '~/Desktop/new_figure_3b.pdf',
       width = 6, height = 4, useDingbats=FALSE)

#######################
### p-values in final figure 3B. Drug donors compared to 
#drug network vs untreated ms patients compared to drug network
#######################
x %>% 
  group_by(group) %>% 
  summarise(p=wilcox.test(jacc~type, exact=FALSE)$p.val)

# # A tibble: 6 x 2
# group          p
# <chr>      <dbl>
#   1 EGCG    0.00676 
# 2 FTY     0.000418
# 3 GA      0.107   
# 4 Healthy 0.195   
# 5 IFNb    0.0462  
# 6 NTZ     0.0730 


#######################
### p-values in final figure 3B. healthy donors compared to 
#healthy network vs untreated ms patients compared to healthy network
#######################

wilcox.test(df.plot %>% filter(group=='Untreated', type=='Same') %>% pull(jacc),
            y$jacc)
# Wilcoxon rank sum test with continuity correction
# 
# data:  df.plot %>% filter(group == "Untreated", type == "Same") %>%  and y$jacc    pull(jacc) and y$jacc
# W = 346260, p-value = 0.0215
# alternative hypothesis: true location shift is not equal to 0
wilcox.test(df.plot %>% filter(group=='Healthy', type=='Same') %>% pull(jacc),
            y$jacc)
# Wilcoxon rank sum test with continuity correction
# 
# data:  df.plot %>% filter(group == "Healthy", type == "Same") %>% pull(jacc) and y$jacc
# W = 209038, p-value = 0.0318
# alternative hypothesis: true location shift is not equal to 0  

#########################################
### Discarded attempt for figure 3B where we compared subsets of the data to untreated 
### main difference: here the comparsisons were pairwise, i.e. healthy to untreated and untreated to healthy
#########################################

g <- df.plot %>% 
  filter(group%in%c('Healthy', 'Untreated')) %>% 
  filter(id.group%in%c('Healthy', 'Untreated')) %>% 
  ggplot(aes(x=group, y=jacc, fill=id.group)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.07),
                aes(colour=id.group)) + 
    scale_fill_manual(values=c('#307FE275', '#FFA30075'), guide=FALSE) + 
    scale_colour_manual(values=c('#307FE2', '#FFA300'), name='Group') + 
    xlab('') + ylab('Jaccard similarity') + 
    theme_bw()
ggsave(g, filename = '~/Desktop/sim_healthy_untreated.pdf',
       width = 4, height = 4, useDingbats=FALSE)

