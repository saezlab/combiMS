####################################################################################
## PLOT FIGURE PREDICTION AND VALIDATION OF COMBINATION THERAPIES
##
## Panel A: Drug combination scores,
##          to be added later as already generated with 
##          the script calculateDefective.R during the run of pathDrugTargetsFinalv3.R
## Panel B: FTY network with co-druggable reactions, 
##          to be inserted later with Inkscape as generated using Cytoscape
## Panel C: Validation FTY-TAK1i
## 
####################################################################################


# By the contribution of
# Marti Bernardo-Faura, 2016
# Jakob Wirbel, July 2017
# Melanie Rinas, 2019


## Load Packages
library(CellNOptR) # Version 1.16
library(reshape2) # Version 1.4.1
library(ggplot2) # Version 2.1.0
# library(cowplot) # Version 0.6.2 # loaded later, otherwise it messes up the facet grid plot
library(gdata) # Version 2.17.0
library(MESS)




## Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

working_directory = getwd()
figure_folder_of_current_script_ = '../../figures/figure_combination_therapy_validation'



# ####################################################################################
### Panel A: Drug combination scores,
##          to be added later as already generated with 
##          the script calculateDefective.R during the run of pathDrugTargetsFinalv3.R
# ####################################################################################

# ####################################################################################
## Panel B: FTY network with co-druggable reactions, 
##          to be inserted later with Inkscape as generated using Cytoscape
####################################################################################






# ####################################################################################
#  Panel C: Validation FTY-TAK1i
# ####################################################################################




# The standard error of the mean (SEM) see e.g. https://en.wikipedia.org/wiki/Standard_error
SE_fct <- function(x) (sd(x,na.rm = T)) / (sqrt(length(x[!is.na(x)])))




Placebo__Exp1_FTY_TAKi = read.xls('../../data/validation_experiments/EAE_FTY_TAK1i__raw_data__clinical_score_per_mouse_per_day.xls', sheet=1, header=TRUE)
FTY__Exp1_FTY_TAKi = read.xls('../../data/validation_experiments/EAE_FTY_TAK1i__raw_data__clinical_score_per_mouse_per_day.xls', sheet=2, header=TRUE)
TAKi__Exp1_FTY_TAKi = read.xls('../../data/validation_experiments/EAE_FTY_TAK1i__raw_data__clinical_score_per_mouse_per_day.xls', sheet=3, header=TRUE)
FTY_TAKi__Exp1_FTY_TAKi = read.xls('../../data/validation_experiments/EAE_FTY_TAK1i__raw_data__clinical_score_per_mouse_per_day.xls', sheet=4, header=TRUE)



mean_Placebo__Exp1_FTY_TAKi = rowMeans(as.matrix(Placebo__Exp1_FTY_TAKi),na.rm = T)
#std_Placebo__Exp1_FTY_TAKi = apply(Placebo__Exp1_FTY_TAKi, 1, function(x) sd(x,na.rm = T))
SE_Placebo__Exp1_FTY_TAKi = apply(Placebo__Exp1_FTY_TAKi, 1,  SE_fct)



mean_FTY__Exp1_FTY_TAKi = rowMeans(as.matrix(FTY__Exp1_FTY_TAKi),na.rm = T)
#std_FTY__Exp1_FTY_TAKi = apply(FTY__Exp1_FTY_TAKi, 1, function(x) sd(x,na.rm = T))  #apply(FTY__Exp1_FTY_TAKi,1,sd)
SE_FTY__Exp1_FTY_TAKi = apply(FTY__Exp1_FTY_TAKi, 1, SE_fct)  #apply(FTY__Exp1_FTY_TAKi,1,sd)

mean_TAKi__Exp1_FTY_TAKi = rowMeans(as.matrix(TAKi__Exp1_FTY_TAKi),na.rm = T)
#std_TAKi__Exp1_FTY_TAKi =apply(TAKi__Exp1_FTY_TAKi, 1, function(x) sd(x,na.rm = T)) # apply(TAKi__Exp1_FTY_TAKi,1,sd)
SE_TAKi__Exp1_FTY_TAKi =apply(TAKi__Exp1_FTY_TAKi, 1, SE_fct) # apply(TAKi__Exp1_FTY_TAKi,1,sd)


mean_FTY_TAKi__Exp1_FTY_TAKi = rowMeans(as.matrix(FTY_TAKi__Exp1_FTY_TAKi),na.rm = T)
#std_FTY_TAKi__Exp1_FTY_TAKi = apply(FTY_TAKi__Exp1_FTY_TAKi, 1, function(x) sd(x,na.rm = T)) # apply(FTY_TAKi__Exp1_FTY_TAKi,1,sd)
SE_FTY_TAKi__Exp1_FTY_TAKi = apply(FTY_TAKi__Exp1_FTY_TAKi, 1, SE_fct) # apply(FTY_TAKi__Exp1_FTY_TAKi,1,sd)








day_vec = 1:length(Placebo__Exp1_FTY_TAKi[,1])

FigureHeight_replicates = 5
FigureWidth_replicates= 7

line_strength=2

x_min = 0
x_max = length(Placebo__Exp1_FTY_TAKi[,1])


y_min=0
y_max=6.5



#   scale_color_manual(values=c('#9C9E9F', '#57AB27', '#0098A1',  '#006165', '#CC071E', '#A11035'), 

col_placebo = "#9C9E9F"
col_FTY = "#57AB27"
col_TAKi = "#0098A1"
col_FTY_TAKi = "#CC071E"



col_TAKi = '#006165'
col_FTY_TAKi = '#A11035'










# PART  FTY- TAKi: Combine ALL  -----------------------------------------------------------------------------------------------------------------------------------------------------------------

FigureHeight_replicates = 5 #5
FigureWidth_replicates= 7 # 7

pch_point_type = 15
y_max = 4.25


pdf(file.path(figure_folder_of_current_script_,paste("figure_combination_therapy_validation_SE_lines.pdf",sep = "")),width = FigureWidth_replicates, height = FigureHeight_replicates)


par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(c(x_min,x_max),c(y_min,y_max), type = "n",
     #main ="Treated mice (mean \u00B1 SE)",
     xlab = "Days",
     ylab = paste("Clinical score "),
     font=1.8,cex.lab=1.2,
     cex.axis=1.2,
     #bty='L',   # make space for legend at the right
     frame.plot=FALSE)

lines(day_vec,mean_Placebo__Exp1_FTY_TAKi, type="o", pch=pch_point_type,col = col_placebo,bg = col_placebo)


arrows(day_vec,mean_Placebo__Exp1_FTY_TAKi, 
       day_vec,mean_Placebo__Exp1_FTY_TAKi+SE_Placebo__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_placebo) 
arrows(day_vec,mean_Placebo__Exp1_FTY_TAKi, 
       day_vec,mean_Placebo__Exp1_FTY_TAKi-SE_Placebo__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_placebo) 




lines(day_vec,mean_FTY__Exp1_FTY_TAKi, type="o", pch=pch_point_type,col = col_FTY,bg = col_FTY)


arrows(day_vec,mean_FTY__Exp1_FTY_TAKi, 
       day_vec,mean_FTY__Exp1_FTY_TAKi+SE_FTY__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_FTY) 
arrows(day_vec,mean_FTY__Exp1_FTY_TAKi, 
       day_vec,mean_FTY__Exp1_FTY_TAKi-SE_FTY__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_FTY) 



lines(day_vec,mean_TAKi__Exp1_FTY_TAKi, type="o", pch=pch_point_type,col = col_TAKi,bg = col_TAKi)


arrows(day_vec,mean_TAKi__Exp1_FTY_TAKi, 
       day_vec,mean_TAKi__Exp1_FTY_TAKi+SE_TAKi__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_TAKi) 
arrows(day_vec,mean_TAKi__Exp1_FTY_TAKi, 
       day_vec,mean_TAKi__Exp1_FTY_TAKi-SE_TAKi__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_TAKi) 

lines(day_vec,mean_FTY_TAKi__Exp1_FTY_TAKi, type="o", pch=pch_point_type,col = col_FTY_TAKi,bg = col_FTY_TAKi)


arrows(day_vec,mean_FTY_TAKi__Exp1_FTY_TAKi, 
       day_vec,mean_FTY_TAKi__Exp1_FTY_TAKi+SE_FTY_TAKi__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_FTY_TAKi) 
arrows(day_vec,mean_FTY_TAKi__Exp1_FTY_TAKi, 
       day_vec,mean_FTY_TAKi__Exp1_FTY_TAKi-SE_FTY_TAKi__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_FTY_TAKi) 



# abline(h = 1,
#        col = 'darkred',
#        lwd = line_strength)
par(xpd=TRUE)
legend("topright", inset=c(-0.3,0),
       legend=c("Placebo", "FTY","TAKi","FTY+TAKi"), 
       pch=pch_point_type,
       #lty=2, 
       col=c(col_placebo,col_FTY,col_TAKi,col_FTY_TAKi),lwd=1, bty='n', cex=1)

dev.off()











pdf(file.path(figure_folder_of_current_script_,paste("figure_combination_therapy_validation_SE_lines_and_area.pdf",sep = "")),width = FigureWidth_replicates, height = FigureHeight_replicates)


par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(c(x_min,x_max),c(y_min,y_max), type = "n",
     #main ="Treated mice (mean \u00B1 SE)",
     xlab = "Days",
     ylab = paste("Clinical score "),
     font=1.8,cex.lab=1.2,
     cex.axis=1.2,
     #bty='L',   # make space for legend at the right
     frame.plot=FALSE)

lines(day_vec,mean_Placebo__Exp1_FTY_TAKi, type="o", pch=pch_point_type,col = col_placebo,bg = col_placebo)


#make polygon where coordinates start with lower limit and 
# then upper limit in reverse order
polygon(c(day_vec,rev(day_vec)),c(mean_Placebo__Exp1_FTY_TAKi-SE_Placebo__Exp1_FTY_TAKi,rev(mean_Placebo__Exp1_FTY_TAKi+SE_Placebo__Exp1_FTY_TAKi)),col = adjustcolor( col_placebo, alpha.f = 0.2), border = FALSE)
#lines(day_vec, df$F, lwd = 2)
#add red lines on borders of polygon
# lines(day_vec, mean_Placebo__Exp1_FTY_TAKi+SE_Placebo__Exp1_FTY_TAKi, col="red",lty=2)
# lines(day_vec, mean_Placebo__Exp1_FTY_TAKi-SE_Placebo__Exp1_FTY_TAKi, col="red",lty=2)

arrows(day_vec,mean_Placebo__Exp1_FTY_TAKi, 
       day_vec,mean_Placebo__Exp1_FTY_TAKi+SE_Placebo__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_placebo) 
arrows(day_vec,mean_Placebo__Exp1_FTY_TAKi, 
       day_vec,mean_Placebo__Exp1_FTY_TAKi-SE_Placebo__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_placebo) 




lines(day_vec,mean_FTY__Exp1_FTY_TAKi, type="o", pch=pch_point_type,col = col_FTY,bg = col_FTY)
polygon(c(day_vec,rev(day_vec)),c(mean_FTY__Exp1_FTY_TAKi-SE_FTY__Exp1_FTY_TAKi,rev(mean_FTY__Exp1_FTY_TAKi+SE_FTY__Exp1_FTY_TAKi)),col = adjustcolor( col_FTY, alpha.f = 0.2), border = FALSE)


arrows(day_vec,mean_FTY__Exp1_FTY_TAKi, 
       day_vec,mean_FTY__Exp1_FTY_TAKi+SE_FTY__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_FTY) 
arrows(day_vec,mean_FTY__Exp1_FTY_TAKi, 
       day_vec,mean_FTY__Exp1_FTY_TAKi-SE_FTY__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_FTY) 



lines(day_vec,mean_TAKi__Exp1_FTY_TAKi, type="o", pch=pch_point_type,col = col_TAKi,bg = col_TAKi)
polygon(c(day_vec,rev(day_vec)),c(mean_TAKi__Exp1_FTY_TAKi-SE_TAKi__Exp1_FTY_TAKi,rev(mean_TAKi__Exp1_FTY_TAKi+SE_TAKi__Exp1_FTY_TAKi)),col = adjustcolor( col_TAKi, alpha.f = 0.2), border = FALSE)


arrows(day_vec,mean_TAKi__Exp1_FTY_TAKi, 
       day_vec,mean_TAKi__Exp1_FTY_TAKi+SE_TAKi__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_TAKi) 
arrows(day_vec,mean_TAKi__Exp1_FTY_TAKi, 
       day_vec,mean_TAKi__Exp1_FTY_TAKi-SE_TAKi__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_TAKi) 

lines(day_vec,mean_FTY_TAKi__Exp1_FTY_TAKi, type="o", pch=pch_point_type,col = col_FTY_TAKi,bg = col_FTY_TAKi)
polygon(c(day_vec,rev(day_vec)),c(mean_FTY_TAKi__Exp1_FTY_TAKi-SE_FTY_TAKi__Exp1_FTY_TAKi,rev(mean_FTY_TAKi__Exp1_FTY_TAKi+SE_FTY_TAKi__Exp1_FTY_TAKi)),col = adjustcolor( col_FTY_TAKi, alpha.f = 0.2), border = FALSE)


arrows(day_vec,mean_FTY_TAKi__Exp1_FTY_TAKi, 
       day_vec,mean_FTY_TAKi__Exp1_FTY_TAKi+SE_FTY_TAKi__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_FTY_TAKi) 
arrows(day_vec,mean_FTY_TAKi__Exp1_FTY_TAKi, 
       day_vec,mean_FTY_TAKi__Exp1_FTY_TAKi-SE_FTY_TAKi__Exp1_FTY_TAKi, length=0.05, angle=90,
       col = col_FTY_TAKi) 



# abline(h = 1,
#        col = 'darkred',
#        lwd = line_strength)
par(xpd=TRUE)
legend("topright", inset=c(-0.3,0),
       legend=c("Placebo", "FTY","TAK1i","FTY+TAK1i"), 
       pch=pch_point_type,
       #lty=2, 
       col=c(col_placebo,col_FTY,col_TAKi,col_FTY_TAKi),lwd=1, bty='n', cex=1)

dev.off()





pdf(file.path(figure_folder_of_current_script_,paste("figure_combination_therapy_validation_SE_area.pdf",sep = "")),width = FigureWidth_replicates, height = FigureHeight_replicates)


par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(c(x_min,x_max),c(y_min,y_max), type = "n",
     #main ="Treated mice (mean \u00B1 SE)",
     xlab = "Days",
     ylab = paste("Clinical score "),
     font=1.8,cex.lab=1.2,
     cex.axis=1.2,
     #bty='L',   # make space for legend at the right
     frame.plot=FALSE)

lines(day_vec,mean_Placebo__Exp1_FTY_TAKi, type="o", pch=pch_point_type,col = col_placebo,bg = col_placebo)


#make polygon where coordinates start with lower limit and 
# then upper limit in reverse order
polygon(c(day_vec,rev(day_vec)),c(mean_Placebo__Exp1_FTY_TAKi-SE_Placebo__Exp1_FTY_TAKi,rev(mean_Placebo__Exp1_FTY_TAKi+SE_Placebo__Exp1_FTY_TAKi)),col = adjustcolor( col_placebo, alpha.f = 0.2), border = FALSE)
#lines(day_vec, df$F, lwd = 2)
#add red lines on borders of polygon
# lines(day_vec, mean_Placebo__Exp1_FTY_TAKi+SE_Placebo__Exp1_FTY_TAKi, col="red",lty=2)
# lines(day_vec, mean_Placebo__Exp1_FTY_TAKi-SE_Placebo__Exp1_FTY_TAKi, col="red",lty=2)

# arrows(day_vec,mean_Placebo__Exp1_FTY_TAKi, 
#        day_vec,mean_Placebo__Exp1_FTY_TAKi+SE_Placebo__Exp1_FTY_TAKi, length=0.05, angle=90,
#        col = col_placebo) 
# arrows(day_vec,mean_Placebo__Exp1_FTY_TAKi, 
#        day_vec,mean_Placebo__Exp1_FTY_TAKi-SE_Placebo__Exp1_FTY_TAKi, length=0.05, angle=90,
#        col = col_placebo) 




lines(day_vec,mean_FTY__Exp1_FTY_TAKi, type="o", pch=pch_point_type,col = col_FTY,bg = col_FTY)
polygon(c(day_vec,rev(day_vec)),c(mean_FTY__Exp1_FTY_TAKi-SE_FTY__Exp1_FTY_TAKi,rev(mean_FTY__Exp1_FTY_TAKi+SE_FTY__Exp1_FTY_TAKi)),col = adjustcolor( col_FTY, alpha.f = 0.2), border = FALSE)


# arrows(day_vec,mean_FTY__Exp1_FTY_TAKi, 
#        day_vec,mean_FTY__Exp1_FTY_TAKi+SE_FTY__Exp1_FTY_TAKi, length=0.05, angle=90,
#        col = col_FTY) 
# arrows(day_vec,mean_FTY__Exp1_FTY_TAKi, 
#        day_vec,mean_FTY__Exp1_FTY_TAKi-SE_FTY__Exp1_FTY_TAKi, length=0.05, angle=90,
#        col = col_FTY) 



lines(day_vec,mean_TAKi__Exp1_FTY_TAKi, type="o", pch=pch_point_type,col = col_TAKi,bg = col_TAKi)
polygon(c(day_vec,rev(day_vec)),c(mean_TAKi__Exp1_FTY_TAKi-SE_TAKi__Exp1_FTY_TAKi,rev(mean_TAKi__Exp1_FTY_TAKi+SE_TAKi__Exp1_FTY_TAKi)),col = adjustcolor( col_TAKi, alpha.f = 0.2), border = FALSE)


# arrows(day_vec,mean_TAKi__Exp1_FTY_TAKi, 
#        day_vec,mean_TAKi__Exp1_FTY_TAKi+SE_TAKi__Exp1_FTY_TAKi, length=0.05, angle=90,
#        col = col_TAKi) 
# arrows(day_vec,mean_TAKi__Exp1_FTY_TAKi, 
#        day_vec,mean_TAKi__Exp1_FTY_TAKi-SE_TAKi__Exp1_FTY_TAKi, length=0.05, angle=90,
#        col = col_TAKi) 

lines(day_vec,mean_FTY_TAKi__Exp1_FTY_TAKi, type="o", pch=pch_point_type,col = col_FTY_TAKi,bg = col_FTY_TAKi)
polygon(c(day_vec,rev(day_vec)),c(mean_FTY_TAKi__Exp1_FTY_TAKi-SE_FTY_TAKi__Exp1_FTY_TAKi,rev(mean_FTY_TAKi__Exp1_FTY_TAKi+SE_FTY_TAKi__Exp1_FTY_TAKi)),col = adjustcolor( col_FTY_TAKi, alpha.f = 0.2), border = FALSE)


# arrows(day_vec,mean_FTY_TAKi__Exp1_FTY_TAKi, 
#        day_vec,mean_FTY_TAKi__Exp1_FTY_TAKi+SE_FTY_TAKi__Exp1_FTY_TAKi, length=0.05, angle=90,
#        col = col_FTY_TAKi) 
# arrows(day_vec,mean_FTY_TAKi__Exp1_FTY_TAKi, 
#        day_vec,mean_FTY_TAKi__Exp1_FTY_TAKi-SE_FTY_TAKi__Exp1_FTY_TAKi, length=0.05, angle=90,
#        col = col_FTY_TAKi) 



# abline(h = 1,
#        col = 'darkred',
#        lwd = line_strength)
par(xpd=TRUE)
legend("topright", inset=c(-0.3,0),
       legend=c("Placebo", "FTY","TAK1i","FTY+TAK1i"), 
       pch=pch_point_type,
       #lty=2, 
       col=c(col_placebo,col_FTY,col_TAKi,col_FTY_TAKi),lwd=1, bty='n', cex=1)

dev.off()





print("Script finished!")

