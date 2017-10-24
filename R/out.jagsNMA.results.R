out.jagsNMA.results=function(JAGSobject,parameter=parameter,forestplot=F,treatnames=NA, rounding=3){
  resultstable=JAGSobject$BUGSoutput$summary
  allvariablenames=rownames(resultstable)
  rowsmatching=substr(allvariablenames,1,nchar(parameter))
  rowstokeep=startsWith(rowsmatching,parameter)

  resultstabletokeep=resultstable[rowstokeep,c(1,3,7,2)]
  
  library(stringr)
  tosplit=unlist(strsplit(rownames(resultstabletokeep),","))
  tosplit2 <- as.numeric(str_extract(tosplit, "[0-9]+"))
  nroftreatments=max(tosplit2)
  location=matrix(tosplit2,ncol=2,byrow=T)
  meanmat=CImat=sdmat=matrix(NA,nrow=nroftreatments,ncol=nroftreatments)

for(i in 1:nrow(location)){
  meanmat[location[i,1],location[i,2]]=resultstabletokeep[i,1]
  sdmat[location[i,1],location[i,2]]=resultstabletokeep[i,4]
  CImat[location[i,1],location[i,2]]=resultstabletokeep[i,3]#high CI
  CImat[location[i,2],location[i,1]]=resultstabletokeep[i,2]#low CI
  
}
  if(forestplot){
    library(metafor)
    slab1=rep(1:(nroftreatments-1),(nroftreatments-1):1)
    a=t(sapply(1:nroftreatments,rep,nroftreatments))
    slab2=a[lower.tri(a,F)]
    slab=paste(slab1,"vs",slab2,sep="")
  forest(x=meanmat[upper.tri(meanmat)], ci.lb=CImat[lower.tri(CImat)],ci.ub=CImat[upper.tri(CImat)], slab=slab,xlab="Network meta-analysis results")
  }
  meanmat=round(meanmat,rounding)
 CImat=round(CImat,rounding)

 #create to print
 Ttreat=dim(meanmat)[1]
 toprintmat=matrix(nrow=Ttreat,ncol=Ttreat)

 for(i in c(1:(Ttreat-1)))
 {for (j in c((i+1):Ttreat))
       {
   toprintmat[i,j]=paste(meanmat[i,j],"(",CImat[j,i],",",CImat[i,j],")",sep="")
   toprintmat[j,i]=paste(c(-meanmat[i,j]),"(",c(-CImat[i,j]),",",c(-CImat[j,i]),")",sep="")
 }}
 if(!missing(treatnames)){
   diag(meanmat)=treatnames
   diag(CImat)=treatnames
   diag(toprintmat)=treatnames
 }

list(Means=meanmat,CI=CImat, leaguetable=toprintmat)
}
