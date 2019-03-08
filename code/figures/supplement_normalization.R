####################################################################################
## PLOT SUPPLEMENTARY FIGURE ABOUT NORMALISATION PIPELINE FOR  THECOMBIMS-PAPER
## Plot densities for the different proteins and how they change 
## after each normalization step
## 
## Only the structure of the figure, annotations will be added
## with Inkscape
## 
## Panel A: Densities of log2 fold changes for all phosphos for one patient
## Panel B-D: Densities of up/down/non-regulated phosphos
## Panel E-G: Densities of up/down/non-regulated phosphos after normalization
## Panel H: Densities of all normalized phosphos
##
####################################################################################

# Load libraries
library(ggplot2)
library(reshape2)
library(CellNOptR)
library(cowplot)

library(latex2exp)
library(ggridges)
library(gridExtra)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Define colours for all phosphos
colour_phosphos = c("#00549F",
                    "dodgerblue",
                    "turquoise",
                    "saddlebrown",
                    "orangered",
                    "lightsalmon",
                    "#7A6FAC",
                    "#612158",
                    "#A11035",
                    "#CC071E",
                    "#F6A800",
                    "#BDCD00",
                    "#57AB27",
                    "#0098A1",
                    "#006165",
                    "#E30066",
                    'navyblue')

####################################################################################
# READ DATA

# Read MIDAS file, convert into CNOlist, and extract data columns
M = readMIDAS(paste0("../../data/phosphos_merged/KI044_phosphos_midas.csv"), verbose=FALSE)
M2 = makeCNOlist(M, subfield = TRUE, verbose=FALSE)
data = M2$valueSignals[[2]]
row.names(data) = M2$namesStimuli
colnames(data) = M2$namesSignals


# remove BSA and PE (internal controls)
data = data[,-c(which(colnames(data) == "BSA"),which(colnames(data) == "PE"))]


names(colour_phosphos) = sort(colnames(data))


#************* step1: X= log2 (rawx)
data_log_2=log2(data)

#************* step2: compute Fold Change, since we have log2 we just substract
data_fc = t(apply(data_log_2, 1, FUN=function(x){return(x-data_log_2[1,])}))
data_fc = data_fc[-1,]
####################################################################################
####################################################################################
ylim_max__data_density = 3.5
# PANEL A
plot_df_a = melt(data_fc)
colnames(plot_df_a) = c('stimulus', 'phospho', 'logFC')
plot_df_a$phospho = factor(plot_df_a$phospho, levels=sort(levels(plot_df_a$phospho)))
p_a = ggplot(plot_df_a, aes(x=logFC, color=phospho)) + 
   geom_density() + theme_cowplot() + 
   scale_color_manual(values=colour_phosphos) + 
   xlab(TeX("$ \\log_{2}$(FC)"))+
   xlim(round(min(plot_df_a$logFC),2),round(max(plot_df_a$logFC),2))+
   ylim(0,ylim_max__data_density)+
   ylab('Density') +
   guides(color=guide_legend(ncol=2))+
   geom_vline(xintercept=0, colour='grey') +
   theme(legend.title = element_blank(),legend.key.size = unit(.375, "cm"),legend.text = element_text(size = 11)) +
   #theme(legend.justification = c(1, 1), legend.position = c(1.1, 1.11),legend.title = element_blank(),legend.key.size = unit(.375, "cm"),legend.text = element_text(size = 11)) +
   #theme(legend.justification = c(1, 1), legend.position = c(1.05, 1.11),legend.title = element_blank(),legend.key.size = unit(.4, "cm"))+
   #theme_bw()+
   #ggtitle("Donor KI044 data, 20 experiments (stimuli), 17 measured phosphoproteins")
   ggtitle("Measured fold change phosphorylation intensity \n")

####################################################################################
# ************ for each protein, find out which distribuitions are skewed, and to which side
# calculate t-test and keep pvalue. for each protein, a patient has a distribution
# prepare a matrix to hold the pvalue for each patient, that this protein is skewed to the right
numProt=length(colnames(data_fc))
pvalueGreater=vector()
pvalueLess=vector()
# to the right?  (i.e. for that prot, if significant the patient has a signfivantly positive fold change)
for (i in 1:numProt){
   pvalueGreater[i]=t.test(data_fc[,i],mu=0,alternative="greater")$p.value
}
positiveProt=which(pvalueGreater<0.05)
# to the left?  (i.e. for that prot, if significant the patient has a signfivantly negative fold change)
for (i in 1:numProt){
   pvalueLess[i]=t.test(data_fc[,i],mu=0,alternative="less")$p.value
}
negativeProt=which(pvalueLess<0.05)
# proteins that fail both tests are not significantly pos or neg phosphorylated
for (i in 1:numProt){
   noiseProt=intersect(which(pvalueGreater>=0.05),which(pvalueLess>=0.05))
}

####################################################################################
# PANEL B-D
# noise proteins
nameProts=colnames(data_fc)[noiseProt]
plot_df_b=plot_df_a[which(plot_df_a$phospho %in% nameProts),]

p_c = ggplot(plot_df_b, aes(logFC)) + 
   geom_density(aes(color=phospho)) +
   geom_vline(xintercept=0, colour='grey') +
   theme(legend.justification = c(1, 1), legend.position = c(1.05, 1.11),legend.title = element_blank(),legend.key.size = unit(.375, "cm"),legend.text = element_text(size = 11)) +
   #theme(legend.justification = c(1, 1), legend.position = c(1.05, 1.11),legend.title = element_blank(),legend.key.size = unit(.4, "cm"))+
   #theme(legend.position = 'none') + 
   scale_color_manual(values=colour_phosphos) +
   #scale_color_manual(values=colour_phosphos[noiseProt]) +
   xlab(TeX("$ \\log_{2}$(FC)"))+
   xlim(round(min(plot_df_a$logFC),2),round(max(plot_df_a$logFC),2))+
   ylim(0,ylim_max__data_density)+
   ylab('Density') +
   ggtitle('Not phosphorylated')+
   theme(plot.title=element_text(face='plain',hjust = 0.5))

# negatively phosphorylated
nameProts=colnames(data_fc)[negativeProt]
plot_df_c=plot_df_a[which(plot_df_a$phospho %in% nameProts),]
p_b = ggplot(plot_df_c, aes(logFC, colour=phospho)) + 
   geom_density() +
   geom_vline(xintercept=0, colour='grey') +
   theme(legend.justification = c(1, 1), legend.position = c(1.05, 1.11),legend.title = element_blank(),legend.key.size = unit(.375, "cm"),legend.text = element_text(size = 11)) +
   #theme(legend.justification = c(1, 1), legend.position = c(1.05, 1.11),legend.title = element_blank(),legend.key.size = unit(.4, "cm"))+
   #theme(legend.justification = c(1, 1), legend.position = c(1.11, 1.11),legend.title = element_blank())+
   #theme(legend.position = "none") + 
   scale_color_manual(values=colour_phosphos) + 
   #scale_color_manual(values=colour_phosphos[negativeProt]) + 
   xlab(TeX("$ \\log_{2}$(FC)"))+
   xlim(round(min(plot_df_a$logFC),2),round(max(plot_df_a$logFC),2))+
   ylim(0,ylim_max__data_density)+
   ylab('Density') +
   ggtitle('Negatively phosphorylated')+
   theme(plot.title=element_text(face='plain',hjust = 0.5))


# positively phosphorylated
nameProts=colnames(data_fc)[positiveProt]
plot_df_d=plot_df_a[which(plot_df_a$phospho %in% nameProts),]
p_d = ggplot(plot_df_d, aes(logFC, colour=phospho)) + 
   geom_density() +
   geom_vline(xintercept=0, colour='grey') +
   guides(color=guide_legend(ncol=2))+
   theme(legend.justification = c(1, 1), legend.position = c(1.1, 1.11),legend.title = element_blank(),legend.key.size = unit(.375, "cm"),legend.text = element_text(size = 11)) +
   #theme(legend.position = 'none')+  # KI044_setting
   #theme(legend.justification = c(1, 1), legend.position = c(1.11, 1.11),legend.title = element_blank())+
   scale_color_manual(values=colour_phosphos)+
   #scale_color_manual(values=colour_phosphos[positiveProt])+
   xlab(TeX("$ \\log_{2}$(FC)"))+
   xlim(round(min(plot_df_a$logFC),2),round(max(plot_df_a$logFC),2))+
   ylim(0,ylim_max__data_density)+
   ylab('Density') +
   ggtitle('Positively phosphorylated')+
   theme(plot.title=element_text(face='plain',hjust = 0.5))

#************** step3: correct pos and neg distributions (non-significant still need to be corrected)
# in positive distributions: leave positives unchanged. between 0 and -1, replace by 0. below -1 replace by NA
data_fc_c = data_fc
numExp=dim(data_fc)[1]
if (length(positiveProt)>0){
   for(i in 1:numExp){
      for (j in 1:length(positiveProt)){
         if (data_fc_c[i,positiveProt[j]]>=(-1) && data_fc_c[i,positiveProt[j]]<0 && !is.na(data_fc_c[i,positiveProt[j]])){
            data_fc_c[i,positiveProt[j]]=0
         }
         
         if(data_fc_c[i,positiveProt[j]]<(-1) && !is.na(data_fc_c[i,positiveProt[j]])){
            data_fc_c[i,positiveProt[j]]=NA
         }
      }
   }
}
# in negative distributions: between 0 and 1, replace by 0. above 1 replace by NA. 
# Last, turn remaining negatives into positives.
if (length(negativeProt)>0){
   print("there's negative distributions")
   for(i in 1:numExp){
      for (j in 1:length(negativeProt)){
         if (data_fc_c[i,negativeProt[j]]>0 && data_fc_c[i,negativeProt[j]]<=1 && !is.na(data_fc_c[i,negativeProt[j]])){
            data_fc_c[i,negativeProt[j]]=0
         } else if(data_fc_c[i,negativeProt[j]]>1 && !is.na(data_fc_c[i,negativeProt[j]])){
            data_fc_c[i,negativeProt[j]]=NA
         } else if (data_fc_c[i,negativeProt[j]]!=0  && !is.na(data_fc_c[i,negativeProt[j]])){
            data_fc_c[i,negativeProt[j]]= abs(data_fc_c[i,negativeProt[j]])
         }
      }
   }
}

# ************* step 4: correct non-significant distributions (not done before because it crashes normalization by 0)  
if(length(noiseProt)>0){
   for(i in 1:numExp){
      for (j in 1:length(noiseProt)){
         if(data_fc_c[i,noiseProt[j]]>=(-1) && data_fc_c[i,noiseProt[j]]<=1 && !is.na(data_fc_c[i,noiseProt[j]])){
            data_fc_c[,noiseProt[j]]=0} else if ( (data_fc_c[i,noiseProt[j]]<(-1) || data_fc_c[i,noiseProt[j]]>1 || is.na(data_fc_c[i,noiseProt[j]])) && !is.na(data_fc_c[i,positiveProt[j]])){
               data_fc_c[,noiseProt[j]]=NA
            }
      }
   }
}


data_fc_c_sign = data_fc_c


# ************* step 5: Hill function on FCz, we are using as EC50 the median of each signal
# taula2 = 1 / (1 + (EC50 / taula2)^HillCoef)
HillCoef=2
data_h=data_fc_c
indexAllButNoise=setdiff(seq(1:ncol(data_fc_c)), noiseProt)
data_h[,indexAllButNoise]=apply(data_fc_c[,indexAllButNoise], 2, function(column){
   1 / (1 + (median(column,na.rm=TRUE) / column) ^ HillCoef)
})


data_h_HILL_of_logFC = data_h


#   #************* step6: shift negatives
if (length(negativeProt)>0){
   for (i in 1:length(negativeProt)){
      data_h[,negativeProt[i]] = 1- abs(data_h[,negativeProt[i]])
   }
}

####################################################################################
# PANEL E-G
# noise proteins
plot_df_c = melt(data_h)
colnames(plot_df_c) = c('stimulus', 'phospho', 'logFC')
nameProts=colnames(data_fc_c)[noiseProt]
plot_df_f=plot_df_c[which(plot_df_c$phospho %in% nameProts),]
p_f = ggplot(plot_df_f, aes(logFC, colour=phospho)) + 
   geom_density(adjust=.000000001) +
   #geom_vline(xintercept=0, colour='grey') +
   theme(legend.position = 'none')+   # KI044_setting
   #theme(legend.justification = c(1, 1), legend.position = c(1, 1.11),legend.title = element_blank())+
   scale_color_manual(values=colour_phosphos)+
   #scale_color_manual(values=colour_phosphos[noiseProt])+
   xlab(TeX("x_{norm}"))+
   #xlab("Normalized data")+
   xlim(0,1)+
   #xlim(round(min(plot_df_a$logFC),2),round(max(plot_df_a$logFC),2))+
   ylim(0,ylim_max__data_density)+
   ylab('Density') +
   #ggtitle('Not phosphorylated')
   ggtitle(' \n ')

# negatively phosphorylated
nameProts=colnames(data_h)[negativeProt]
plot_df_e=plot_df_c[which(plot_df_c$phospho %in% nameProts),]
p_e = ggplot(plot_df_e, aes(logFC, colour=phospho)) + 
   geom_density(na.rm = TRUE,trim=TRUE) +
   #geom_vline(xintercept=0, colour='grey') +
   theme(legend.position = 'none')+   # KI044_setting
   #theme(legend.justification = c(1, 1), legend.position = c(1, 1.11),legend.title = element_blank())+
   scale_color_manual(values=colour_phosphos) +
   #scale_color_manual(values=colour_phosphos[negativeProt]) + 
   xlab(TeX("x_{norm}"))+
   #xlab("Normalized data")+
   xlim(0,1)+
   #xlim(round(min(plot_df_a$logFC),2),round(max(plot_df_a$logFC),2))+
   ylim(0,ylim_max__data_density)+
   ylab('Density') +
   #+ggtitle('Negatively phosphorylated')
   ggtitle("Normalized phosphorylation intensity \n")

# positively phosphorylated
nameProts=colnames(data_h)[positiveProt]
plot_df_g=plot_df_c[which(plot_df_c$phospho %in% nameProts),]
p_g = ggplot(plot_df_g, aes(logFC, colour=phospho)) + 
   geom_density(na.rm = TRUE,trim=TRUE) +
   #geom_vline(xintercept=0, colour='grey') +
   #theme(legend.justification = c(0, 1), legend.position = c(0, 1.11),legend.title = element_blank())+
   scale_color_manual(values=colour_phosphos) + 
   theme(legend.position = 'none')+   # KI044_setting
   #scale_color_manual(values=colour_phosphos[positiveProt]) + 
   xlab(TeX("x_{norm}"))+
   #xlab("Normalized data")+
   xlim(0,1)+
   #xlim(round(min(plot_df_a$logFC),2),round(max(plot_df_a$logFC),2))+
   ylim(0,ylim_max__data_density)+
   ylab('Density') +
#   ggtitle('Positively phosphorylated')
   ggtitle(' \n')


####################################################################################
# PANEL H


# plot all normalized proteins

nameProts=colnames(data_h)
plot_df_h=plot_df_c[which(plot_df_c$phospho %in% nameProts),]


p_h = ggplot(plot_df_h, aes(logFC, colour=phospho)) + 
   geom_density(trim=TRUE,na.rm = TRUE) +
   #geom_bar(position='dodge') +
   # facet_grid(protein~., scales='free') +
   #geom_vline(xintercept=0, colour='grey') +
   theme(legend.position = 'none')+
   scale_color_manual(values=colour_phosphos) + 
   xlab(TeX("x_{norm}"))+
   xlim(0,1)+
   #xlim(round(min(plot_df_a$logFC),2),round(max(plot_df_a$logFC),2))+
   ylim(0,ylim_max__data_density)+
   ylab('Density')
   #theme_bw()+
   # ggtitle("Donor KI044 normalized data, 20 experiments (stimuli), 17 measured phosphoproteins")
  # ggtitle("Normalized data")



####################################################################################
# ANNOTATION PANEL
annotation = data.frame(matrix(c(2, 1, 1, 3), nrow=2))




annotation = data.frame(matrix(c(2, 1, .7, 3), nrow=2))

p_annotation_pos = 
   ggplot(annotation, aes(X1, X2)) + theme_nothing() + geom_point(color='white') +
   annotate('text', x=1.1, y=2.5, label =TeX("(1) -1 $\\leq$ $ \\log_{2}$(FC) < 0 $\\rightarrow$ 0", output = 'character'), parse=TRUE,hjust ='left') + 
   annotate('text', x=1.1, y=2, label =TeX("(2) $ \\log_{2}$(FC) < -1 $\\rightarrow$ NA", output = 'character'), parse=TRUE,hjust ='left') + 
   annotate('text', x=1.1, y=1.5, label =TeX("(3) Non-linear normalization", output = 'character'), parse=TRUE,hjust ='left')


p_annotation_not = 
   ggplot(annotation, aes(X1, X2)) + theme_nothing() + geom_point(color='white') +
   annotate('text', x=1.1, y=2.5, label =TeX("(1) -1 $\\leq$ $ \\log_{2}$(FC) $\\leq$ 1 $\\rightarrow$ 0", output = 'character'), parse=TRUE,hjust ='left') + 
   annotate('text', x=1.1, y=2, label =TeX("(2) $ \\log_{2}$(FC) < -1 $\\rightarrow$ NA", output = 'character'), parse=TRUE,hjust ='left') + 
   annotate('text', x=1.1, y=1.5, label =TeX("(3) $ \\log_{2}$(FC) > 1 $\\rightarrow$ NA", output = 'character'), parse=TRUE,hjust ='left') +
   annotate('text', x=1.1, y=1, label =TeX("(4) $ \\log_{2}$(FC) $\\rightarrow$ x_{norm}", output = 'character'), parse=TRUE,hjust ='left')

p_annotation_neg = 
   ggplot(annotation, aes(X1, X2)) + theme_nothing() + geom_point(color='white') +
   annotate('text', x=1.1, y=2.5, label =TeX("(1) 0 < $ \\log_{2}$(FC) $\\leq$ 1 $\\rightarrow$ 0", output = 'character'), parse=TRUE,hjust ='left') + 
   annotate('text', x=1.1, y=2, label =TeX("(2) $ \\log_{2}$(FC) > 1 $\\rightarrow$ NA", output = 'character'), parse=TRUE,hjust ='left') + 
   annotate('text', x=1.1, y=1.5, label =TeX("(3) Non-linear normalization", output = 'character'), parse=TRUE,hjust ='left')+
   annotate('text', x=1.1, y=1, label =TeX("(4) 1 - abs(x) $\\rightarrow$ x_{norm}", output = 'character'), parse=TRUE,hjust ='left')

#  annotate('text', x=1.1, y=1, label =TeX("(4) $ \\log_{2}$(FC) $\\rightarrow$ 1 - abs($ \\log_{2}$(FC))", output = 'character'), parse=TRUE,hjust ='left')

####################################################################################
# PLOT EVERYTHING TOGETHER
# pdf('../figures/supplementary_normalization.pdf', size=c(8.3, 11.7))
top_row = plot_grid(NULL, p_a, NULL, ncol=3, rel_widths=c(0.3, 1, 0.2))
upper_middle_row = plot_grid(p_d, p_b, p_c, ncol=3)                                    #upper_middle_row = plot_grid(p_b, p_c, p_d, ncol=3)
middle_row = plot_grid(p_annotation_pos, p_annotation_neg, p_annotation_not, ncol=3)   #middle_row = plot_grid(p_annotation_neg, p_annotation_not, p_annotation_pos, ncol=3)
lower_middle_row = plot_grid(p_g, p_e, p_f, ncol=3)                                    #lower_middle_row = plot_grid(p_e, p_f, p_g, ncol=3)
bottom_row = plot_grid(NULL,p_h, NULL, ncol=3)

pdf('../../figures/supplement_normalization_KI044.pdf', height = 11, width = 8.5)
plot_grid(top_row, upper_middle_row, middle_row, lower_middle_row, bottom_row, nrow=5, 
          rel_heights = c(1, 1, .4, 1, 1))
dev.off()


# 
# # PLOT EVERYTHING TOGETHER
# # pdf('../figures/supplementary_normalization.pdf', size=c(8.3, 11.7))
# top_row = plot_grid(NULL, p_a, NULL, ncol=3, rel_widths=c(0.3, 1, 0.2))
# upper_middle_row = plot_grid(p_d, p_b, p_c, ncol=3)                                    #upper_middle_row = plot_grid(p_b, p_c, p_d, ncol=3)
# middle_row = plot_grid(p_annotation_pos, p_annotation_neg, p_annotation_not, ncol=3)   #middle_row = plot_grid(p_annotation_neg, p_annotation_not, p_annotation_pos, ncol=3)
# lower_middle_row = plot_grid(p_g, p_e, p_f, ncol=3)                                    #lower_middle_row = plot_grid(p_e, p_f, p_g, ncol=3)
# bottom_row = plot_grid(NULL,p_h, NULL, ncol=3)
# 
# pdf('../../figures/supplement_normalization_KI044v2labels.pdf', height = 11, width = 8.5)
# plot_grid(top_row, upper_middle_row, middle_row, lower_middle_row, bottom_row, nrow=5, 
#           rel_heights = c(1, 1, .4, 1, 1), labels=c('A', 'B', 'C', 'D','E'))
# 
# dev.off()


