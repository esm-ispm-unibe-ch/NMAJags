modelIPDNMA_randomEffects <- function(){
  
  # 1) Likelihood for IPD
  for (i in 1:Np){
    outcome[i] ~ dbern(p[i])
    logit(p[i]) <- u[studyid[i]] + d[studyid[i], arm[i]]
  }
  
  # 2) Multi-arm specification with random effects
  for (i in 1:Nstudies) {
    d[i,1] <- 0
    w[i,1] <- 0
    for (k in 2:na[i]) {
      # Random effects: each study-arm difference
      d[i,k] ~ dnorm(md[i,k], precD[i,k])
      
      # The "fixed" mean difference + multi-arm shift
      md[i,k] <- (delta[treat[i,k]] - delta[treat[i,1]]) + sw[i,k]
      
      # Multi-arm correction
      w[i,k]  <- d[i,k] - (delta[treat[i,k]] - delta[treat[i,1]])
      sw[i,k] <- sum(w[i,1:(k-1)]) / (k - 1)
      precD[i,k] <- prec * 2*(k-1)/k  # arm-specific precision adjustment
      
    }
  }
  
  # 3) Study intercepts
  for (i in 1:Nstudies){
    u[i] ~ dnorm(0,0.001)
  }
  
  # 4) Priors for treatment effects
  delta[ref] <- 0
  for(k in 1:(ref-1)) {
    delta[k] ~ dnorm(0, 0.001)
  }
  for(k in (ref+1):nt) {
    delta[k] ~ dnorm(0, 0.001)
  }
  
  # 5) Random-effects variance
  tau ~  dnorm(0, 1)%_%T(0,)    # for example
  prec <- 1 / pow(tau, 2)
  
  # OR for each comparison
  for(i in 1:(nt-1)) {
    for (j in (i+1):nt) {
      OR[j,i]<- exp(delta[j] - delta[i])
      LOR[j,i]<- delta[j] - delta[i]}}
  
  # 6) Derived ORs vs reference
  for(j in 1:(ref-1)){ORref[j]<- exp(delta[j] - delta[ref])}
  for(j in (ref+1):nt) {ORref[j]<- exp(delta[j] - delta[ref])}
  
  #Predictions
  for (i in 1:(nt - 1)) {
    for (j in (i + 1):nt) {
      predLOR[j, i] ~ dnorm(LOR[j, i], prec)
      predOR[j, i] = exp(predLOR[j, i])
    }
  }
  
  #Ranking of treatments
  order[1:nt] <- rank(-1*delta[1:nt]) # if OR > 1 desirable, else rank(d[1:nt])
  for (i in 1:nt) {
    most.effective[i] <- equals(order[i], 1)
    for (j in 1:nt) {
      effectiveness[i, j] <- equals(order[i], j)
    }
  }
  for (i in 1:nt) {
    for (j in 1:nt) {
      cumeffectiveness[i, j] <- sum(effectiveness[i, 1:j])
    }
  }
  for (i in 1:nt) {
    SUCRA[i] <- sum(cumeffectiveness[i, 1:(nt - 1)])/(nt-1)
  }
  
}
