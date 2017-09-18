pool.arms=function(data, studlab,treat, r, n, y,sd, type = "binary"){
  
  # This is a function that looks for arms within a study that have the same treatment code and pools data 
  #on a binary outcome
  #The required data format has to be one study arm per row
  #studyid: the study id
  #t the treatment
  #r the events
  #n the numbers randomised
  # y the means
  #sd the SDs
  #type should be "binary"or "cont"

  #Get the data
  initialid=eval(substitute(studlab), data)
  data$idd = eval(substitute(studlab), data)
  tt=data$tt = eval(substitute(treat), data)
  
  if(type=="binary"){
     n = data$n = eval(substitute(n), data)
     r = data$n = eval(substitute(r), data)}
  else{
     y=  data$y = eval(substitute(y), data)
     n = data$n = eval(substitute(n), data)
     sd = data$sd = eval(substitute(sd), data)
        }
  idd = data$idd
  idd = as.numeric(as.factor(idd))

  #Check for and identify duplicates
  treatsbyarmmat=table(idd,tt)
  studies.with.problem=unique(idd)[apply(treatsbyarmmat>1,1,sum)>0]
  
  #identify the studies that do not need pooling
  clean=is.na(match(idd,studies.with.problem))

  #Pool duplicates
  if(type=="binary"){allr=alln=allt=allid=c()}
  if(type=="cont"){ally=alln=allt=allid=allsd=c()}

  for(i in studies.with.problem)
    {
      nonproblematicarm=colnames(treatsbyarmmat)[treatsbyarmmat[i,]==1]
    
      if(length(nonproblematicarm)==0){cat("\n", paste("the study with id",  unique(initialid)[i], "is deleted because it compares the same treatements", "\n"))}
      else{
          problematicarm=colnames(treatsbyarmmat)[treatsbyarmmat[i,]>1]
            if(type=="binary") { 
              rnew=sum(r[idd==i & tt==problematicarm])
              rold=r[idd==i & tt==nonproblematicarm]
              allr=c(allr,rnew,rold)}
            else{
              select=c(idd==i & tt==problematicarm)
              noselect=c(idd==i & tt==nonproblematicarm)
              sdnew=pooledSD(sd[select],n[select])
              sdold=sd[noselect]
              ynew=sum(y[select]*n[select])/sum(n[select])
              yold=y[noselect]
              ally=c(ally,ynew,yold)
              allsd=c(allsd,sdnew,sdold)
               }
        nnew=sum(n[idd==i & tt==problematicarm])
        nold=n[idd==i & tt==nonproblematicarm]
        alln=c(alln,nnew,nold)
        allid=c(allid,i,i)
        allt=c(allt,problematicarm,nonproblematicarm)
        }
    }
  
  #put all data together
  if(type=="binary") 
    { 
      data=cbind.data.frame(
        id=c(initialid[clean],unique(initialid)[allid] ),
        treat=c(as.character(tt)[clean],allt),
        r=c(r[clean], allr),
        n=c(n[clean], alln))
      }
  else
    {
      data=cbind.data.frame(
          id=c(initialid[clean],unique(initialid)[allid] ),
          treat=c(as.character(tt)[clean],allt),
          y=c(y[clean], ally),
          sd=c(sd[clean], allsd),
          n=c(n[clean], alln))
      }
  
  #push them out
  return(data)
  
}

