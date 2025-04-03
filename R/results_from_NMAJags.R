
results_from_NMAJAGS = function(JAGSobject, reference.trt.ind, trt.names, name.comp="OR", PI=TRUE, plots = TRUE, higher.desirable=TRUE, 
                                rounding=2){
  # parameters: 
  # JAGSobject = what running, for example, jags.parallel returns
  # reference.trt.ind = index (= number) of the reference treatment in the JAGSobject
  # trt.names = names of treatments in the order that they appear in the JAGSobject
  # name.comp = name of the comparison between treatments (for example OR, SMD. inspect JAGSobject to see what you have)
  # PI = should prediction intervals be returned? Only works, if JAGSobject has paste("pred", name.comp, sep="") parameters saved
  # plots = should forest plots be made for CI and PI against reference
  # higher.desirable = are higher numbers in the comparison desirable (for example is an odds ratio over 1 desirable?) PLots will be ranked this way
  # returns: 
  # mean = leaguetable (= every treatment compared to every other treatment) with mean of comparison per cell
  # median = leaguetable with median of comparison per cell
  # upper.CI = leaguetable with upper bound of CI of comparison per cell
  # lower.CI = leaguetable with lower bound of CI of comparison per cell
  # upper.PI = leaguetable with upper bound of PI of comparison per cell
  # lower.PI = leaguetable with lower bound of PI of comparison per cell
  # print.CI = leaguetable with median and CI of comparison per cell
  # print.PI = leaguetable with median and PI of comparison per cell
  # SUCRA = SUCRA values for the treatmetns
  # tau = heterogeneity variance tau, if tau in JAGSobject summary output  
  # if plots == "yes" plots forest plots with confidence intervals and prediction intervals
  
  
  sm = JAGSobject$BUGSoutput$summary
  inds = which(grepl(paste("\\b", name.comp, "\\[", sep=""), rownames(sm)))
  sm = sm[inds,]
  # dimension of resulting matrix = nt x nt (nt = number of treatments)
  nt = (1+sqrt(1+8*length(inds)))/2 # from midnight formula
  
  # Credible intervals
  ltable.mean = matrix(nrow=nt, ncol=nt)
  ltable.median = matrix(nrow=nt, ncol=nt)
  ltable.ub = matrix(nrow=nt, ncol=nt) # upper bound CI
  ltable.lb = matrix(nrow=nt, ncol=nt) # lower bound CI
  
  for (j in 1:(nt-1)){
    for (i in (j+1):nt){
      ind = which(grepl(paste(name.comp, "\\[", i, ",", j, "\\]", sep=""), rownames(sm)))
      ltable.mean[i,j] = sm[ind,"mean"]
      ltable.mean[j,i] = 1/(sm[ind,"mean"])
      ltable.median[i,j] = sm[ind,"50%"]
      ltable.median[j,i] = 1/(sm[ind,"50%"])
      ltable.lb[i,j] = sm[ind,"2.5%"]
      ltable.lb[j,i] = 1/sm[ind,"97.5%"] 
      ltable.ub[i,j] = sm[ind,"97.5%"]
      ltable.ub[j,i] = 1/sm[ind,"2.5%"]
      
    }
  }
  ltable.mean = round(ltable.mean, rounding)
  ltable.median = round(ltable.median, rounding)
  ltable.ub = round(ltable.ub, rounding)
  ltable.lb = round(ltable.lb, rounding)
  
  ltable.print = matrix(nrow=nt, ncol=nt) # gives for each treatment against each other treatment median and CI in a cell
  for (i in 1:nt){
    for (j in 1:nt){
      ltable.print[i,j] = paste(ltable.median[i,j], " (", ltable.lb[i,j], ", ", ltable.ub[i,j], ")", sep="")
    }
  }
  diag(ltable.print) = trt.names
  
  # prediction intervals
  sm = JAGSobject$BUGSoutput$summary
  if (PI){ # only make PI, if we have tau
    ltable.ub.PI = matrix(nrow=nt, ncol=nt) # upper bound PI
    ltable.lb.PI = matrix(nrow=nt, ncol=nt) # lower bound PI
    for (j in 1:(nt-1)){
      for (i in (j+1):nt){
        ind = which(grepl(paste("pred", name.comp, "\\[", i, ",", j, "\\]", sep=""), rownames(sm)))
        ltable.lb.PI[i,j] = sm[ind,"2.5%"]
        ltable.lb.PI[j,i] = 1/sm[ind,"97.5%"] 
        ltable.ub.PI[i,j] = sm[ind,"97.5%"]
        ltable.ub.PI[j,i] = 1/sm[ind,"2.5%"]
      }
    }
    
    ltable.lb.PI=round(ltable.lb.PI, rounding)
    ltable.ub.PI=round(ltable.ub.PI, rounding)
    ltable.print.PI = matrix(nrow=nt, ncol=nt) # gives for each treatment against each other treatment median and PI in a cell
    for (i in 1:nt){
      for (j in 1:nt){
        ltable.print.PI[i,j] = paste(ltable.median[i,j], " (", ltable.lb.PI[i,j], ", ", ltable.ub.PI[i,j], ")", sep="")
      }
    }
    diag(ltable.print.PI) = trt.names
  }else{
    ltable.ub.PI = "was not requested"
    ltable.lb.PI = "was not requested"
    ltable.print.PI = "was not requested"
  }
  
  # SUCRA
  sm = JAGSobject$BUGSoutput$summary
  inds = which(grepl("SUCRA", rownames(sm)))
  median.SUCRA = sm[inds, "50%"]
  names(median.SUCRA) = trt.names
  median.SUCRA = sort(median.SUCRA, decreasing=TRUE)
  median.SUCRA = round(median.SUCRA, rounding)
  
  # plots
  if (plots){
    library(metafor)
    forest(ltable.median[-reference.trt.ind,reference.trt.ind], ci.lb=ltable.lb[-reference.trt.ind,reference.trt.ind], 
           ci.ub=ltable.ub[-reference.trt.ind,reference.trt.ind],
           slab=trt.names[-reference.trt.ind], xlab=paste("Network", name.comp, "vs", trt.names[reference.trt.ind]),
           header=c("Treatment", "Median [95% CI]"), refline=ifelse(name.comp=="OR", 1, 0), 
           order=order(median.SUCRA[trt.names[-reference.trt.ind]], decreasing=higher.desirable), 
           psize=0.75)
    if (PI){
      forest(ltable.median[-reference.trt.ind,reference.trt.ind], ci.lb=ltable.lb.PI[-reference.trt.ind,reference.trt.ind], 
             ci.ub=ltable.ub.PI[-reference.trt.ind,reference.trt.ind], slab=trt.names[-reference.trt.ind], 
             xlab=paste("Network", name.comp, "vs", trt.names[reference.trt.ind]), header=c("Treatment", "Median [95% PI]"), 
             refline=ifelse(name.comp=="OR", 1, 0), order=order(median.SUCRA[trt.names[-reference.trt.ind]], decreasing=higher.desirable),
             psize=0.75)
    }
    
  }
  
  if ((length(which(grepl("tau", rownames(sm))))>0)){
    tau = round(sm["tau","mean"], rounding)
  }else{
    tau = "tau not in summary table of JAGSobject"
  }
  
  toreturn = list("mean"=ltable.mean, "median"=ltable.median, "median.ci"=ltable.print, "upper.ci"=ltable.ub, "lower.ci"=ltable.lb,
                  "upper.PI"=ltable.ub.PI, "lower.PI"=ltable.lb.PI, "print.CI" = ltable.print, "print.PI" = ltable.print.PI, 
                  "SUCRA"=median.SUCRA, "tau"=tau)
  return(toreturn)
  
}
