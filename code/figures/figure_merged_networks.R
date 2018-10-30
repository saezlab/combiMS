# ###############################################################################################################################################
# SCRIPT FOR FIGURE 5
#
# 
# PANEL A-E: Networks made in Cytoscape, to be inserted by Inkscape
# PANEL F: Density-Distribution of P-values for the statistical validation
#
# 
# #################################################################################################################################################

# Libraries
library(reshape2)
library(ggplot2)
library(cowplot)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load results from statistical test
load('../../files/statistics/all_groups_wilcoxon_results.RData') # c(2.95, 2.25, 1.15, 1, .5), y = c(.3, .3, .5, 2, 2)

ann_text <- data.frame(log.p.value = c(.5, 1.15, 2.95, 1, 2.25), y = c(5, .5, .25, 3, .3), lab = c('EGCG', 'FTY', 'IFNb', 'GA', 'NTZ'), 
                       interaction_in_model=c('none', 'none', 'none', 'none', 'none'),
                       group = factor(c('EGCG', 'FTY', 'IFNb', 'GA', 'NTZ'), levels = c('EGCG', 'FTY', 'IFNb', 'GA', 'NTZ')))
all_groups$group = factor(all_groups$group, levels=c('EGCG', 'FTY', 'IFNb', 'GA', 'NTZ'))

g = ggplot(all_groups, aes(x=log.p.value, colour=interaction_in_model, fill=interaction_in_model)) + 
  geom_density() + 
  facet_wrap(~group, ncol=1, scales = "free") +
  theme_classic() + theme(legend.position='none', axis.text=element_text(size=6)) + 
  xlab(expression('-log'['10']*'(Adj. p-value)')) + ylab('Density') + 
  scale_color_manual('Interaction', values=c('#8EBAE5','black', '#D85C41'), labels=c('not in\nmodel', 'in model', 'none')) +
  theme(strip.background=element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual('Interaction', values=alpha(c( '#8EBAE5', 'black', '#D85C41'),0.3), labels=c('not in\nmodel', 'in model', 'none')) +
  geom_text(data=ann_text, y=ann_text$y, label=ann_text$lab)

pdf('../../figures/figure_merged_networks.pdf', width = 7, height = 10)
plot_grid(NULL, NULL, NULL, NULL, NULL, g, ncol=2, labels='AUTO')
dev.off()
