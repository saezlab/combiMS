
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
