pool.arms=function (data, studlab, treat, r, n, y, sd, type = "binary")
{


  initialid = eval(substitute(studlab), data)
  data$idd = eval(substitute(studlab), data)
  tt = data$tt = eval(substitute(treat), data)
  if(type == "binary"){
    n = data$n = eval(substitute(n), data)
    r = data$n = eval(substitute(r), data)
  }
  if(type == "cont"){
    y = data$y = eval(substitute(y), data)
    n = data$n = eval(substitute(n), data)
    sd = data$sd = eval(substitute(sd), data)
  }
  idd = data$idd
  idd = as.numeric(as.factor(idd))
  treatsbyarmmat = table(idd, tt)
  studies.with.problem = unique(idd)[apply(treatsbyarmmat >
                                             1, 1, sum) > 0]
  clean = is.na(match(idd, studies.with.problem))
  if (type == "binary") {
    allr = alln = allt = allid = c()
  }
  if (type == "cont") {
    ally = alln = allt = allid = allsd = c()
  }

  for (i in studies.with.problem) {
    nonproblematicarm = colnames(treatsbyarmmat)[treatsbyarmmat[i,] == 1]
    problematicarm = colnames(treatsbyarmmat)[treatsbyarmmat[i,] > 1]

      if (type == "binary") {
        rnew=c()
        nnew=c()
        for(j in problematicarm){
          select=c(idd == i & tt == j)
          rnew = c(rnew,sum(r[select]))
          nnew = c(nnew,sum(n[select]))}

        rold = r[idd == i & tt == nonproblematicarm]
        nold = n[idd == i & tt == nonproblematicarm]
        allr = c(allr, rnew, rold)
        alln = c(alln, nnew, nold)
      }

      if (type == "cont"){
        sdnew=c()
        ynew=c()
        nnew=c()
        for(j in problematicarm){
        select = c(idd == i & tt == j)
        sdnew = c(sdnew,pooledSD(sd[select], n[select]))
        ynew = c(ynew,sum(y[select] * n[select])/sum(n[select]))
        nnew = c(nnew,sum(n[select]))
                                }

        noselect = c(idd == i & tt == nonproblematicarm)
        sdold = sd[noselect]
        yold = y[noselect]
        nold=n[noselect]
        ally = c(ally, ynew, yold)
        allsd = c(allsd, sdnew, sdold)
        alln = c(alln, nnew, nold)

      }


      allt = c(allt, problematicarm, nonproblematicarm)
      allid = c(allid, rep(i,times=length(problematicarm)), rep(i,length(nonproblematicarm)))
  }

  if (type == "binary") {
    data = cbind.data.frame(id = c(initialid[clean], unique(initialid)[allid]),
                            treat = c(as.character(tt)[clean], allt), r = c(r[clean],
                                                                            allr), n = c(n[clean], alln))
                        }
  if (type == "cont") {
    data = cbind.data.frame(id = c(initialid[clean], unique(initialid)[allid]),
                            treat = c(as.character(tt)[clean], allt), y = c(y[clean],
                                                                            ally), sd = c(sd[clean], allsd), n = c(n[clean],
                                                                                                                   alln))
                      }
  nr.of.arms=table(data$id)
  oxo=names(table(data$id))[table(data$id)<2]
  cat(paste("\n","The studies with the following IDs have been excluded because they comapre the same treatments:","\n"))
  cat(paste(oxo,sep=", ","\n"))

  delete.arms=match(oxo,data$id)
  data=data[-delete.arms,]

  return(data)
}
