####################################################################################
## PLOT SUPPLEMENTARY FIGURE ABOUT MODEL MERGING QUALITY CONTROL
##
## Panel A: Clustering of the similarity matrix for median models of all patients
## Panel B: Scatter plot of the scores for the best model vs the merged model error
## Panel C: Scatter plot of number of edges against the MSE
## 
## Based on the script plotMedianModels.R
####################################################################################

# *************** libraries
library(ggplot2)
library(ggrepel)
library(corrplot)
library(cowplot)
library(reshape2)

# *************** set working directory for relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# *************** get data saved before
load('../../files/model_merging/mini_DF.RData')
load("../../files/similarity/similarityMatrix.RData")

#***********************************************************************
# *************** Clustering of similarity Matrix
#***********************************************************************
ordering = hclust(dist(similarityMatrix))
dendro = as.dendrogram(ordering)
plot_dendro = dendro_data(dendro, type='rectangle')
p = ggplot(segment(plot_dendro)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  theme_nothing()
upper_dendrogramm = p
lower_dendrogramm = p + coord_flip()

plot_df = melt(similarityMatrix)

plot_df$Var2 <- with(plot_df, factor(Var2, levels = levels(Var2)[ordering$order]))
plot_df$Var1 <- with(plot_df, factor(Var1, levels = levels(Var1)[ordering$order]))

clustermap = ggplot(plot_df, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color='white', size=0.0001) +
  scale_fill_gradient2(low='blue', 
                       high = 'red', 
                       mid = 'whitesmoke', midpoint = 0.78, guide = FALSE) +
  coord_equal() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank())


test = plot_grid(upper_dendrogramm, NULL, NULL, clustermap, lower_dendrogramm, legend,
                 nrow=2, rel_heights = c(0.1, 1), rel_widths = c(1, .1, .3))

#***********************************************************************
# *************** plot linear model best error vs merged model error
#***********************************************************************
best_vs_merged = ggplot(data=mini_DF, aes(best_of_bests,merged_model_score)) + 
  geom_point(shape=1) + 
  geom_smooth(method=lm) + 
  theme_bw() + 
  ylab('Score of Merged Model') + 
  xlab('Score of Best Model') + 
  geom_text_repel(aes(label=ifelse(abs(best_of_bests-merged_model_score)>0.01,as.character(name),'')))

#***********************************************************************
# *************** plot num edges vs MSE
#***********************************************************************
mse = ggplot(data=mini_DF,aes(edges,MSE)) + 
  geom_smooth(method=lm) + 
  geom_point(size=2) + 
  xlab('Number of Edges in Model') + 
  ylab('Mean Squared Error') + 
  theme_bw()

#***********************************************************************
# *************** combine all plots
#***********************************************************************
upper_panel = plot_grid(NULL, clustermap, NULL, ncol=3, labels=c('', 'A', ''), rel_widths = c(0.5, 1, 0.5))
lower_panel = plot_grid(best_vs_merged, mse, ncol=2, labels=c('B', 'C'))
plot_grid(upper_panel, lower_panel, nrow=2)
