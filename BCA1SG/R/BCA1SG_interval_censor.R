BCA1SG_interval_censor<-function(input_data,initial_beta,initial_Lambda=function(x){x},threshold=1e-5,max_iter=5000,max_stepsize=1e4,xi=0.3,contraction=0.5){
  if(sum(is.na(input_data))>0){
    stop("The input data frame contains NA.")
  }
  eps<-10^(-8)
  eps1<-10^(-300)
  tmp_index<-which(input_data[,1]==input_data[,2])
  if(length(tmp_index)>0){
    for(i in 1:length(tmp_index)){
      input_data[tmp_index[i],1]<-input_data[tmp_index[i],1]-0.01
      input_data[tmp_index[i],2]<-input_data[tmp_index[i],2]+0.01
    }
  }
  tmp_index<-which((input_data[,1]==0)&(input_data[,2]==Inf))
  if(length(tmp_index)>0){
    input_data<-input_data[-tmp_index,]
  }
  Lowerbound<-input_data[,1]
  Upperbound<-input_data[,2]
  if(max(Upperbound[Upperbound!=Inf])>max(Lowerbound)){
    stop("The MLE of the nonparametric cumulative hazard contain Infinity.")
  }
  samplesize<-nrow(input_data)
  cov_dim<-ncol(input_data)-2
  if(cov_dim!=length(initial_beta)){
    stop("The dimension of initial_beta does not comply with the dimension of the covariates.")
  }
  covariate<-as.matrix(input_data[,-c(1,2)])
  uniquepoints<-unique(sort(c(Lowerbound,Upperbound)))
  if(uniquepoints[1]==0){
    uniquepoints<-uniquepoints[-1]
  }
  if(max(uniquepoints)==Inf){
    uniquepoints<-uniquepoints[-length(uniquepoints)]
  }
  n<-length(uniquepoints)
  transL<-rep(0,samplesize)
  transU<-rep(-1,samplesize)
  for(i in 1:samplesize){
    if(Lowerbound[i]>0){
      transL[i]<-which(abs(Lowerbound[i]-uniquepoints)<eps)
    }
    if(Upperbound[i]<Inf){
      transU[i]<-which(abs(Upperbound[i]-uniquepoints)<eps)
    }
  }
  backupU<-transU
  backupL<-transL
  backupL[transL==0]<-1
  backupU[transU==-1]<-n
  indicator<-rep(2,samplesize)
  indicator[Lowerbound==0]<-1
  indicator[Upperbound==Inf]<-3
  loglikelihood<-function(Lambda,beta=numeric(0)){
    for(i in 1:n){
      Lambda[n-i+1]<-sum(Lambda[1:(n-i+1)])
    }
    term1<-as.numeric(exp(beta%*%t(covariate)))
    term2<-exp(-Lambda[backupL]*term1)
    term3<-exp(-Lambda[backupU]*term1)
    lk<-sum((indicator==1)*pmax(log(1-term3),-1e200)+(indicator==2)*pmax(log(term2-term3),-1e200)+(indicator==3)*pmax(log(term2),-1e200))
    return(lk)
  }
  loglikelihood0<-function(beta,Lambda=numeric(0)){
    for(i in 1:n){
      Lambda[n-i+1]<-sum(Lambda[1:(n-i+1)])
    }
    term1<-as.numeric(exp(beta%*%t(covariate)))
    term2<-exp(-Lambda[backupL]*term1)
    term3<-exp(-Lambda[backupU]*term1)
    lk<-sum((indicator==1)*pmax(log(1-term3),-1e200)+(indicator==2)*pmax(log(term2-term3),-1e200)+(indicator==3)*pmax(log(term2),-1e200))
    return(-lk)
  }
  explc_grad<-function(Lambda,beta=numeric(0)){
    for(i in 1:n){
      Lambda[n-i+1]<-sum(Lambda[1:(n-i+1)])
    }
    mygrad<-rep(0,n)
    term1<-as.numeric(exp(beta%*%t(covariate)))
    term2<-exp(-Lambda[backupL]*term1)
    term3<-exp(-Lambda[backupU]*term1)
    cond1<-term1*term3/pmax((1-term3),1e-8)
    cond2<-term1*term3/pmax((term2-term3),1e-8)
    for(i in 1:n){
      mygrad[i]<-sum((indicator==1)*(i<=transU)*cond1+(indicator==2)*((i<=transL)*(-term1)+(i>transL)*(i<=transU)*pmin(cond2,1e100))+(indicator==3)*(i<=transL)*(-term1))
    }
    return(mygrad)
  }
  explc_hess<-function(Lambda,beta=numeric(0)){
    for(i in 1:n){
      Lambda[n-i+1]<-sum(Lambda[1:(n-i+1)])
    }
    mygrad<-rep(0,n)
    term1<-as.numeric(exp(beta%*%t(covariate)))
    term2<-exp(-Lambda[backupL]*term1)
    term3<-exp(-Lambda[backupU]*term1)
    cond3<--term1^2*term3/pmax((1-term3)^2,1e-8)
    cond4<--term1^2*term2*term3/pmax((term2-term3)^2,1e-8)
    for(i in 1:n){
      mygrad[i]<-sum((indicator==1)*(i<=transU)*cond3+(indicator==2)*(i>transL)*(i<=transU)*pmax(cond4,-1e100))
    }
    return(mygrad)
  }

  Int_Lambda<-initial_Lambda(uniquepoints)
  ###################################################BCA1SG Algorithm
  old_Lambda<-Int_Lambda
  old_Lambda<-c(old_Lambda[1],diff(old_Lambda))
  old_Lambda[old_Lambda<=0]<-eps
  old_beta<-initial_beta
  active_set<-numeric(0)
  ptm<-proc.time()
  flag<-TRUE
  count<-0
  while((flag)&(count<max_iter)){
    movement<-Inf
    while((movement>threshold)&(count<max_iter)){
      ##########first update the parametric part
      res<-stats::optim(par=old_beta,fn=loglikelihood0,Lambda=old_Lambda)$par
      movement1<-max(abs(old_beta-res))
      last_beta<-old_beta
      old_beta<-res
      ###########then update the nonparametric part
      hess<-explc_hess(old_Lambda,beta = old_beta)
      mygrad<-explc_grad(old_Lambda,beta = old_beta)
      hess<-hess-rep(eps,n)
      ###############compute the movement direction
      direct<--1/hess*mygrad
      count<-count+1
      ###############project the direction onto the region defined by the active constraints
      direct[active_set]<-0
      ###############compute the largest step length
      d<-direct[direct<0]
      old<-old_Lambda[direct<0]
      if(length(d)>0){
        alpha<-min((eps-old)/d)
      }else{
        alpha<-1
      }
      alpha<-min(alpha,max_stepsize)
      direct<-direct*alpha
      ###############compute the new_Theta candidate
      new_Lambda<-old_Lambda+direct
      ########line search
      line_count<-0
      while((loglikelihood(new_Lambda,beta = old_beta)<(loglikelihood(old_Lambda,beta = old_beta)+xi*mygrad%*%direct))&&(line_count<50)){
        alpha<-alpha*contraction
        direct<-direct*contraction
        new_Lambda<-old_Lambda+direct
        line_count<-line_count+1
      }
      ###############update the active set
      active_set<-which(((new_Lambda-eps)<=10^(-6)))
      increment<-loglikelihood(new_Lambda,beta = old_beta)-loglikelihood(old_Lambda,beta = last_beta)
      movement<-max(abs(c(new_Lambda-old_Lambda,movement1)))
      old_Lambda<-new_Lambda
      # print(c(count,movement))
    }
    ###############decide whether to remove some point from the boundary
    if(length(active_set)>0){
      lagrg<-rep(0,length(active_set))
      for(i in 1:length(active_set)){
        lagrg[i]<-mygrad[active_set[i]]
      }
      if(max(lagrg)>10^(-6)){
        # print("remove")
        # print(max(lagrg))
        active_set<-active_set[-which.max(lagrg)]
      }else{
        flag<-FALSE
      }
    }else{
      flag<-FALSE
    }
  }
  for(i in 1:n){
    old_Lambda[n-i+1]<-sum(old_Lambda[1:(n-i+1)])
  }
  timecost<-as.numeric((proc.time()-ptm)[3])
  if(count>=max_iter){
    warning("The algorithm fails to converge. Please try to increase the threshold or the max_iter.")
  }
  output_MLE<-list(distinct_time=uniquepoints,est_Lambda=old_Lambda,est_beta=old_beta,iteration=count,timecost=timecost)
  return(output_MLE)
}
