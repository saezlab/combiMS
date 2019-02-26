
# SMC_of_2networks.R

# Script to calculate the similarity of 2 Boolean networks represented as bStrings 
# by using the SMC


# SMC = simple matching coefficient
# (e.g. see https://en.wikipedia.org/wiki/Simple_matching_coefficient)
# 
# is a measure to calculate similarity
# 

# by 
# Melanie Rinas
# December 2017






SMC_of_2networks <- function(bString_network1, bString_network2) {     
   
   sum(bString_network1 == bString_network2) / length(bString_network1)
   
}






# 
# # Testing SMC_fct:
# 
# n1=c(0,0,1,1)
# n2=c(0,1,0,1)
# 
# number_of_matching_attributes__M00_plus_M11 = sum(n1 == n2)
# number_of_attributes__M00_plus_M01_plus_M10_plus_M11 = length(n1)
# 
# SMC_fct(n1,n2)
