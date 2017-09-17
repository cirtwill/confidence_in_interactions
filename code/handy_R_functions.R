calculate_parameters<-function(priordata,n,k){
  pars=fitdistr(x=priordata,"beta",start=list(shape1=1,shape2=1),lower=c(0,0))$estimate
  alpha=pars[[1]]
  beta=pars[[2]]
  alpha_prime=alpha+k
  beta_prime=beta+n-k
  pars2=c(alpha_prime,beta_prime)
  return(pars2)
}


calculate_distribution<-function(pars){
  alpha=pars[[1]]
  beta=pars[[2]]

  mu_num=alpha
  mu_den=alpha+beta
  mu=mu_num/mu_den

  sig_num=alpha*beta
  den1=alpha+beta
  den2=den1**2
  sig_den=den1*den2
  sigma2=sig_num/sig_den

  return(c(mu,sigma2))
}


credible_interval<-function(pars,p_lower,p_upper){
  alpha=pars[[1]]
  beta=pars[[2]]
  lowCI=qbeta(p=p_lower,shape1=alpha,shape2=beta)
  highCI=qbeta(p=p_upper,shape1=alpha,shape2=beta)
  return(c(lowCI,highCI))
}

