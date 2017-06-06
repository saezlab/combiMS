####################################################################################
## PLOT FIGURE PREDICTION AND VALIDATION OF COMBINATION THERAPIES
##
## Panel A: Drug combination scores
## Panel B: FTY network with co-druggable reactions, to be inserted later with Inkscape
## Panel C: Validation FTY-EGCG
## Panel D: Validation FTY-TAK1i
## 
####################################################################################

## Load Packages
library(CellNOptR) # Version 1.16
library(reshape2) # Version 1.4.1
library(ggplot2) # Version 2.1.0
# library(cowplot) # Version 0.6.2 # loaded later, otherwise it messes up the facet grid plot
library(gdata) # Version 2.17.0

## Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ####################################################################################
# Panel A
# ####################################################################################

# ####################################################################################
# Panel C and D
# ####################################################################################

## get data
plot_df_c = read.xls('../../data/validation_experiments/EAE_FTY_EGCG.xls', sheet=1, header=TRUE)
colnames(plot_df_c) = c('days', 'FTY', 'Placebo', 'FTY+\nEGCG', 'EGCG')
plot_df_c = melt(plot_df_c, id='days')
colnames(plot_df_c) = c('days', 'treatment', 'clinical score')

## reorder treatment
plot_df_c$treatment = factor(plot_df_c$treatment, levels=c('Placebo', 'FTY', 'EGCG', 'FTY+\nEGCG'))
plot_df_c$comb = 'FTY + EGCG'

## get data
plot_df_d = read.xls('../../data/validation_experiments/EAE_FTY_TAK1i.xls', sheet=1, header=TRUE)[,c(1,2,3,4,5)]
colnames(plot_df_d) = c('days', 'FTY+\nTAKi', 'FTY', 'Placebo', 'TAKi')
plot_df_d = melt(plot_df_d, id='days')
colnames(plot_df_d) = c('days', 'treatment', 'clinical score')

## reorder treatment
plot_df_d$treatment = factor(plot_df_d$treatment, levels=c('Placebo', 'FTY', 'TAKi', 'FTY+\nTAKi'))
plot_df_d$comb = 'FTY + TAKi'

# Combine and plot
temp = rbind(plot_df_c, plot_df_d)
temp$treatment = factor(temp$treatment,  levels=c('Placebo', 'FTY', 'EGCG', 'TAKi', 'FTY+\nEGCG', 'FTY+\nTAKi'))

g_v = ggplot(temp, aes(x=days, y=`clinical score`, color=treatment)) + 
  geom_line() + facet_grid(~comb, scales='free_x') + theme_classic() + 
  geom_point(shape=15) + theme(strip.background=element_blank(), legend.key.size = unit(1.7, 'lines')) + 
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  scale_color_manual(values=c('#9C9E9F', '#57AB27', '#0098A1',  '#006165', '#CC071E', '#A11035'))

# ####################################################################################
# Combine everything
# ####################################################################################

library(cowplot) # Version 0.6.2 # loaded later, otherwise it messes up the facet grid plot

pdf('../../figures/figure_combination_therapy.pdf', width=7, height=7)
top_row = plot_grid(NULL, NULL, ncol=2, labels=c('A', 'B'))
plot_grid(top_row, g_v, nrow=2, labels=c('', 'C'))
dev.off()
