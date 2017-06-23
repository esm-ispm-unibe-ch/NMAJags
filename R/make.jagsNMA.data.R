make.jagsNMA.data=
function (studyid, r, n, y, sd, t, type = "cont", data, reference = 1, 
          othervar = NA) 
{
  data$idd = eval(substitute(studyid), data)
  data$tt = eval(substitute(t), data)
  n = data$n = eval(substitute(n), data)
  idd = data$idd
  idd = as.numeric(as.factor(idd))
  print(idd)
  tt = data$tt
  ns = length(unique(idd))
  nt = length(unique(tt))
  na = table(idd)
  idd = rep(1:ns, table(idd))
  t = as.numeric(as.factor(tt))
  if (!identical(t, tt)) {
    print("Note: the treatments have been renamed as follows")
    out <- cbind.data.frame(`old names` = sort(unique(tt)), 
                            `new names` = sort(unique(t)))
    refer = sort(unique(t))[sort(unique(tt)) == reference]
    print(out)
  }
  maxnrofarms = max(table(idd))
  nofarms <- length(idd)
  armsenumerate = unlist(sapply(na, seq))
  tmat <- matrix(999, nrow = ns, ncol = maxnrofarms)
  for (i in 1:nofarms) {
    tmat[idd[i], armsenumerate[i]] = t[i]
  }
  tmat2 = t(apply(tmat, 1, sort))
  tmat2[tmat2 == 999] <- NA
  tmat[tmat == 999] <- NA
  nmat <- matrix(-99, nrow = ns, ncol = nt)
  for (i in 1:nofarms) {
    nmat[idd[i], t[i]] <- n[i]
  }
  nmat[nmat == -99] <- NA
  if (!missing(othervar)) {
    variable = data$variable = eval(substitute(othervar), 
                                    data)
    varmat <- matrix(NA, nrow = ns, ncol = nt)
    for (i in 1:nofarms) {
      varmat[idd[i], t[i]] <- variable[i]
    }
    varmat[is.na(varmat)] <- 0
  }
  if (type == "cont") {
    y = data$y = eval(substitute(y), data)
    sd = data$sd = eval(substitute(sd), data)
    se = sd/sqrt(n)
    prec = 1/(se * se)
    ymat <- matrix(9999, nrow = ns, ncol = nt)
    precmat <- matrix(-99, nrow = ns, ncol = nt)
    for (i in 1:nofarms) {
      ymat[idd[i], t[i]] = y[i]
      precmat[idd[i], t[i]] = prec[i]
    }
    ymat[ymat == 9999] <- NA
    precmat[precmat == -99] <- NA
    nominator = sqrt(tapply(n * sd * sd, idd, sum))
    denominator = sqrt(tapply(n, idd, sum) - na)
    pooled.sd = nominator/denominator
    if (!missing(othervar)) {
      toreturn = list(ns = ns, nt = nt, na = na, t = tmat2, 
                      y = ymat, prec = precmat, pooled.sd = pooled.sd, 
                      ref = refer, variab = apply(varmat,1,max))
    }
    else {
      toreturn = list(ns = ns, nt = nt, na = na, t = tmat2, 
                      y = ymat, prec = precmat, pooled.sd = pooled.sd, 
                      ref = refer)
    }
  }
  if (type == "binary") {
    r = data$r = eval(substitute(r), data)
    rmat <- matrix(-99, nrow = ns, ncol = nt)
    for (i in 1:nofarms) {
      rmat[idd[i], t[i]] = r[i]
    }
    rmat[rmat == -99] <- NA
    if (!missing(othervar)) {
      toreturn = list(ns = ns, nt = nt, na = na, t = tmat2, 
                      r = rmat, n = nmat, ref = refer, variab =apply(varmat,1,max))
    }
    else {
      toreturn = list(ns = ns, nt = nt, na = na, t = tmat2, 
                      r = rmat, n = nmat, ref = refer)
    }
  }
  return(toreturn)
}