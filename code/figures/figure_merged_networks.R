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
load('../../files/statistics/all_groups_wilcoxon_results.RData')

ann_text <- data.frame(log.p.value = c(2.95, 2.25, 1.15, 1, .5), y = c(.3, .3, .5, 2, 2), lab = c("INFb", 'NTZ','FTY', 'GA', 'EGCG'), 
                       interaction_in_model=c('none', 'none', 'none', 'none', 'none'),
                       group = factor(c('IFNb', 'NTZ', 'FTY', 'GA', 'EGCG'), levels = c('IFNb', 'NTZ', 'FTY', 'GA', 'EGCG')))

g = ggplot(all_groups, aes(x=log.p.value, colour=interaction_in_model, fill=interaction_in_model)) + 
  geom_density() + 
  facet_wrap(~group, ncol=1, scales = "free") +
  theme_classic() + theme(legend.position='none', axis.text=element_text(size=6)) + 
  xlab(expression('-log'['10']*'(Adj. p-value)')) + ylab('Density') + 
  scale_color_manual('Interaction', values=c('#8EBAE5','black', '#D85C41'), labels=c('not in\nmodel', 'in model', 'none')) +
  theme(strip.background=element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual('Interaction', values=alpha(c( '#8EBAE5', 'black', '#D85C41'),0.3), labels=c('not in\nmodel', 'in model', 'none')) +
  geom_text(data=ann_text, y=ann_text$y, label=ann_text$lab)

plot(g)

pdf('../../figures/figure_merged_networks.pdf', width = 7, height = 10)
plot_grid(NULL, NULL, NULL, NULL, NULL, g, ncol=2, labels='AUTO')
dev.off()
