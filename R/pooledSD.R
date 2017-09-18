pooledSD=function(sd,n){
  #function that pools SDs into a single value
  SD=sqrt(sum((n-1)*sd*sd)/(sum(n)-length(n)))
  SD
}



