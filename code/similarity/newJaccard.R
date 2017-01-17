newJaccard <- function(d1, d2) {
  apply(d1, 1, function(row1) {
    apply(d2, 1, function(row2) {
      prova=jaccNello_2(row1, row2)
      # cat(prova,", ", sep="")
      return(prova)
      })
  })

}

jaccNello <- function(x, y) {
  sum(x == y) / length(x)
} 

jaccNello_2 = function(x, y){
  temp = rbind(x, y)
  temp = temp[,-which(colSums(temp == 0) == 2)]
  return(sum(temp[1,] == temp[2,])/dim(temp)[2])
}
