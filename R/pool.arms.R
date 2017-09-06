pool.arms=function(data, studlab,treat, r, n){
  
  # This is a function that looks for arms within a study that have the same treatment code and pools data 
  #on a binary outcome
  #The required data format has to be one study arm per row
  #studyid: the study id
  #t the treatment
  #r the events
  #n the numbers randomised
  
  #Get the data
  initialid=eval(substitute(studlab), data)
  data$idd = eval(substitute(studlab), data)
  tt=data$tt = eval(substitute(treat), data)
  n = data$n = eval(substitute(n), data)
  r = data$n = eval(substitute(r), data)
  idd = data$idd
  idd = as.numeric(as.factor(idd))

  #Check for and identify duplicates
  treatsbyarmmat=table(idd,tt)
  studies.with.problem=unique(idd)[apply(treatsbyarmmat>1,1,sum)>0]
  
  #identify the studies that do not need pooling
  clean=is.na(match(idd,studies.with.problem))

  #Pool duplicates
  allr=alln=allt=allid=c()
  for(i in studies.with.problem){
    nonproblematicarm=colnames(treatsbyarmmat)[treatsbyarmmat[i,]==1]
    if(length(nonproblematicarm)==0)
      {cat("\n", paste("the study with id",  unique(initialid)[i], "is deleted because it compares the same treatements", "\n"))}
    else
      {
    problematicarm=colnames(treatsbyarmmat)[treatsbyarmmat[i,]>1]
    rnew=sum(r[idd==i & tt==problematicarm])
    nnew=sum(n[idd==i & tt==problematicarm])
    rold=r[idd==i & tt==nonproblematicarm]
    nold=n[idd==i & tt==nonproblematicarm]
    allr=c(allr,rnew,rold)
    alln=c(alln,nnew,nold)
    allid=c(allid,i,i)
    allt=c(allt,problematicarm,nonproblematicarm)
      }
  }
  
  #put all data together
  data=cbind.data.frame(id=c(initialid[clean],unique(initialid)[allid] ),
  treat=c(as.character(tt)[clean],allt),
  r=c(r[clean], allr),
  n=c(n[clean], alln))
  
  #push them out
  return(data)
  
}

