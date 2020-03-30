example.NMAJags=function (){
  cat("\n","\n")
  cat("##########################################################################","\n")
cat(" E X A M P L E using JAGS to fit NMA for a continuous outcome","\n")
cat("##########################################################################","\n")
cat("\n","\n")
cat("First, have a look at the data typing head(CPPS)" ,"\n" )
head(CPPS)
  cat(
"\n","\n", "Then, transform the data into a list suitable for JAGS analysis typing:"
,"\n", "NMAdataContinuous=make.jagsNMA.data(studyid=id,t=t,y=y,sd=sd,n=n,data=CPPS,type=\"cont\")" ,"\n")
  cat("\n","This might also re-name the drugs - if so, note them down!","\n")
  NMAdataContinuous=make.jagsNMA.data(studyid=id,t=t,y=y,sd=sd,n=n,data=CPPS,type="cont")
  cat("\n","Look at the object you created by typing NMAdataContinuous","\n")


cat("\n","\n", "Then run Jags and create a jags object typing"
,"\n","NMAinJAGS<- jags(data = NMAdataContinuous, inits = NULL,
                 parameters.to.save = c(\"SMD\",\"tau\"), n.chains = 2, n.iter = 10000,
                 n.burnin = 1000,DIC=F,
                 model.file = modelNMAContinuous)"
,"\n","\n","and print the results asking for  print(NMAinJAGS)","\n","\n")

  #run Jags and create a jags object
  NMAinJAGS<- jags(data = NMAdataContinuous, inits = NULL,
                   parameters.to.save = c("SMD","tau"), n.chains = 2, n.iter = 10000,
                   n.burnin = 1000,DIC=F,
                   model.file = modelNMAContinuous)

  cat("\n","\n", "Then print the object by typing print(NMAinJAGS)","\n","\n")
  print(NMAinJAGS)
  cat("\n","\n", "Check chain mixing by typing traceplot(NMAinJAGS)","\n","\n")

 traceplot(NMAinJAGS)
 cat("\n","\n", "Then, get the league table and a forest plot by typing","\n",
     "leaguetable=out.jagsNMA.results(NMAinJAGS,\"SMD\")","\n")
  leaguetable=out.jagsNMA.results(NMAinJAGS,"SMD")

  cat("\n","\n", "The name of the object (leaguetable) to see the leaguetable","\n","\n")
      leaguetable
      cat("##########################################################################","\n")
      cat("END of  E X A M P L E using JAGS to fit NMA for a continuous outcome","\n")
      cat("##########################################################################","\n")
}
