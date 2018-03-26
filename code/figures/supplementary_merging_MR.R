
# supplementary_merging_MR.R

#supplementary_merging_MR__using_Cluster_MR_Results_based_on_using_MIDAS_phosphos_processed__original_MB_JW__AND__combiMSplaneCUT__gaBinaryT1_postCellNOptRupdate__OKpostCellNOptRupdate



# Some notes added
# 
# AND
# 
# paths adapted for CombiMS rerun
# 
# AND
# 
# some modifications with no influence on the results
# 
# AND
# 



####################################################################################
# ERROR correction concerning :
# 
# The script supplementary_merging.R has some ERROR PROBLEMS:

# Error: could not find function "dendro_data"
# Error in ggplot(segment(plot_dendro)) : could not find function "segment"
# ...
####################################################################################

# 
# by 
# Melanie Rinas
# December 2017







Results_storage_name = "OptCombiMSclusterFork10_MR__using_MIDAS_phosphos_processed__original_MB_JW__AND__combiMSplaneCUT__gaBinaryT1_postCellNOptRupdate"        # MR inserted




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
# load('../../files/model_merging/mini_DF.RData')               # MR modified
# load("../../files/similarity/similarityMatrix.RData")      # MR modified





files_Similarity_MR_folder = "../../files/Similarity_MR"                                                             # MR inserted
sub_dir__Results_Similarity_MR = file.path(files_Similarity_MR_folder,Results_storage_name)                          # MR inserted

load(paste(sub_dir__Results_Similarity_MR,"/similarityMatrix.RData",sep="")) # called similarityMatrix               # MR inserted

# wrong_similarityMatrix
load(paste(sub_dir__Results_Similarity_MR,"/wrong_similarityMatrix.RData",sep="")) # called wrong_similarityMatrix               # MR inserted






files_model_merging_MR_folder = "../../files/model_merging_MR"   
sub_dir__Results__model_merging_MR = file.path(files_model_merging_MR_folder,Results_storage_name)              # MR inserted

load(paste(sub_dir__Results__model_merging_MR,"/mini_DF.RData",sep="")) # called similarityMatrix               # MR inserted

#mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR
load(paste(sub_dir__Results__model_merging_MR,"/mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR.RData",sep="")) # called mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR               # MR inserted




figures_MR_folder = "../../figures_MR"   


ifelse(!dir.exists(file.path(figures_MR_folder,"supplementary_merging_MR")), dir.create(file.path(figures_MR_folder,"supplementary_merging_MR")), FALSE)              # MR inserted
sub_dir__Figures_folder_of_current_script_class = file.path(figures_MR_folder,"supplementary_merging_MR")                                                                              # MR inserted

ifelse(!dir.exists(file.path(sub_dir__Figures_folder_of_current_script_class,Results_storage_name)), dir.create(file.path(sub_dir__Figures_folder_of_current_script_class,Results_storage_name)), FALSE)              # MR inserted
sub_dir__Figures_folder_of_current_script = file.path(sub_dir__Figures_folder_of_current_script_class,Results_storage_name)                                                                              # MR inserted





#***********************************************************************
# *************** Clustering of similarity Matrix
#***********************************************************************

dist_similarityMatrix = dist(similarityMatrix)   # MR inserted

ordering = hclust(dist(similarityMatrix))       # Hierarchical Clustering
dendro = as.dendrogram(ordering)




# START MR commented ................................................................................................................

# plot_dendro = dendro_data(dendro, type='rectangle')
# 
# 
# p = ggplot(segment(plot_dendro)) + 
#    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
#    theme_nothing()
# upper_dendrogramm = p
# lower_dendrogramm = p + coord_flip()


# START MR commented ................................................................................................................


plot_df = melt(similarityMatrix)

plot_df_original = plot_df   # MR inserted


plot_df$Var2 <- with(plot_df, factor(Var2, levels = levels(Var2)[ordering$order]))
plot_df$Var1 <- with(plot_df, factor(Var1, levels = levels(Var1)[ordering$order]))





FigureWidth__map =5   # MR inserted
FigureHeight__map=8   # MR inserted

clustermap = ggplot(plot_df, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color='white', size=0.0001) +
  scale_fill_gradient2(low='blue', 
                       high = 'red', 
                       mid = 'whitesmoke', midpoint = 0.78, guide = FALSE) +
  coord_equal() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank())
ggsave(file.path(sub_dir__Figures_folder_of_current_script,paste("clustermap__similarityMatrix__based_on_hclust.pdf",sep = "")), width = FigureWidth__map, height = FigureHeight__map)             # MR inserted


# test = plot_grid(upper_dendrogramm, NULL, NULL, clustermap, lower_dendrogramm, legend,               # MR commented
#                 nrow=2, rel_heights = c(0.1, 1), rel_widths = c(1, .1, .3))












#***********************************************************************
# *************** plot linear model best error vs merged model error
#***********************************************************************

FigureWidth__geom_point = 3.5   # MR inserted
FigureHeight__geom_point = 3.5   # MR inserted

best_vs_merged = ggplot(data=mini_DF, aes(best_of_bests,merged_model_score)) + 
  geom_point(shape=1) + 
  geom_smooth(method=lm) + 
  theme_bw() + 
  ylab('Score of Merged Model') + 
  xlab('Score of \nBest Model') + 
  geom_text_repel(aes(label=ifelse(abs(best_of_bests-merged_model_score)>0.01,as.character(name),'')))
ggsave(file.path(sub_dir__Figures_folder_of_current_script,paste("geom_point__x__best_of_bests__y__merged_model_score.pdf",sep = "")), width = FigureWidth__geom_point, height = FigureHeight__geom_point)             # MR inserted


#***********************************************************************
# *************** plot num edges vs MSE
#***********************************************************************
mse = ggplot(data=mini_DF,aes(edges,MSE_MR)) +           # MR modified
  geom_smooth(method=lm) + 
  geom_point(size=2) + 
  #   xlab('Number of Edges in Model') + 
  xlab('Number of Edges \nin Merged Model') +             # MR modified
  #   ylab('Mean Squared Error') + 
  ylab('Mean Squared Error \nof Merged Model') +          # MR modified
  theme_bw()
ggsave(file.path(sub_dir__Figures_folder_of_current_script,paste("geom_point__x__edges__y__MSE_MR.pdf",sep = "")), width = FigureWidth__geom_point, height = FigureHeight__geom_point)             # MR inserted

#***********************************************************************
# *************** combine all plots
#***********************************************************************

pdf(file.path(sub_dir__Figures_folder_of_current_script,paste("supplementary_merging_MR__A_clustermap_similarityMatrix_B_best_of_bests_VS_merged_model_score_C_edges_VS_MSE_MR.pdf",sep = "")),width = FigureWidth__map, height = FigureHeight__map,onefile=FALSE) # MR inserted

upper_panel = plot_grid(NULL, clustermap, NULL, ncol=3, labels=c('', 'A', ''), rel_widths = c(0.5, 1, 0.5))
lower_panel = plot_grid(best_vs_merged, mse, ncol=2, labels=c('B', 'C'))
plot_grid(upper_panel, lower_panel, nrow=2)

dev.off()
















# START MR inserted ................................................................................................................


####################################################################################################################################################

#***********************************************************************
# *************** Repeat everything with wrong results...
#***********************************************************************

####################################################################################################################################################






#***********************************************************************
# *************** Clustering of wrong similarity Matrix
#***********************************************************************

dist_wrong_similarityMatrix = dist(wrong_similarityMatrix)   # MR inserted

wrong_ordering = hclust(dist(wrong_similarityMatrix))       # Hierarchical Clustering
wrong_dendro = as.dendrogram(wrong_ordering)




# START MR commented ................................................................................................................

# plot_dendro = dendro_data(dendro, type='rectangle')
# 
# 
# p = ggplot(segment(plot_dendro)) + 
#    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
#    theme_nothing()
# upper_dendrogramm = p
# lower_dendrogramm = p + coord_flip()


# START MR commented ................................................................................................................


wrong_plot_df = melt(wrong_similarityMatrix)

wrong_plot_df_original = wrong_plot_df   # MR inserted


wrong_plot_df$Var2 <- with(wrong_plot_df, factor(Var2, levels = levels(Var2)[wrong_ordering$order]))
wrong_plot_df$Var1 <- with(wrong_plot_df, factor(Var1, levels = levels(Var1)[wrong_ordering$order]))





FigureWidth__map =5   # MR inserted
FigureHeight__map=8   # MR inserted

wrong_clustermap = ggplot(wrong_plot_df, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color='white', size=0.0001) +
  scale_fill_gradient2(low='blue', 
                       high = 'red', 
                       mid = 'whitesmoke', midpoint = 0.78, guide = FALSE) +
  coord_equal() +
  labs(title = "Clustermap of wrong_similarityMatrix") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank())
ggsave(file.path(sub_dir__Figures_folder_of_current_script,paste("clustermap__wrong_similarityMatrix__based_on_hclust.pdf",sep = "")), width = FigureWidth__map, height = FigureHeight__map)             # MR inserted


# test = plot_grid(upper_dendrogramm, NULL, NULL, clustermap, lower_dendrogramm, legend,
#                  nrow=2, rel_heights = c(0.1, 1), rel_widths = c(1, .1, .3))












#***********************************************************************
# *************** plot linear model wrong best error vs merged model error
#***********************************************************************

FigureWidth__geom_point = 3.5   # MR inserted
FigureHeight__geom_point = 3.5   # MR inserted

wrong_best_vs_merged = ggplot(data=mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR, aes(wrong_best_of_bests,merged_model_score)) + 
  geom_point(shape=1) + 
  geom_smooth(method=lm) + 
  theme_bw() + 
  ylab('Score of Merged Model') + 
  xlab('Wrong Score \nof Best Model') + 
  geom_text_repel(aes(label=ifelse(abs(wrong_best_of_bests-merged_model_score)>0.01,as.character(name),'')))
ggsave(file.path(sub_dir__Figures_folder_of_current_script,paste("geom_point__x__wrong_best_of_bests__y__merged_model_score.pdf",sep = "")), width = FigureWidth__geom_point, height = FigureHeight__geom_point)             # MR inserted


#***********************************************************************
# *************** plot num edges vs MSEwrong
#***********************************************************************
wrong_mse = ggplot(data=mini_DF__with_1wrong_best_of_bests_2wrongMSE_3MSEs_MR,aes(edges,MSEwrong)) +           # MR modified
  geom_smooth(method=lm) + 
  geom_point(size=2) + 
  #   xlab('Number of Edges in Model') + 
  xlab('Number of Edges \nin Merged Model') +             # MR modified
  #   ylab('Mean Squared Error') + 
  ylab('Wrong Mean Squared Error \nof Merged Model') +          # MR modified
  theme_bw()
ggsave(file.path(sub_dir__Figures_folder_of_current_script,paste("geom_point__x__edges__y__MSEwrong.pdf",sep = "")), width = FigureWidth__geom_point, height = FigureHeight__geom_point)             # MR inserted

#***********************************************************************
# *************** combine all plots
#***********************************************************************

pdf(file.path(sub_dir__Figures_folder_of_current_script,paste("supplementary_merging_MR__A_clustermap_wrong_similarityMatrix_B_wrong_best_of_bests_VS_merged_model_score_C_edges_VS_MSEwrong.pdf",sep = "")),width = FigureWidth__map, height = FigureHeight__map,onefile=FALSE) # MR inserted

upper_panel = plot_grid(NULL, wrong_clustermap, NULL, ncol=3, labels=c('', 'A', ''), rel_widths = c(0.5, 1, 0.5))
lower_panel = plot_grid(wrong_best_vs_merged, wrong_mse, ncol=2, labels=c('B', 'C'))
plot_grid(upper_panel, lower_panel, nrow=2)

dev.off()

# END MR inserted ................................................................................................................


print("Script finished!")   # MR inserted

