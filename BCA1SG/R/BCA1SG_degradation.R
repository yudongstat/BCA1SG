BCA1SG_degradation<-function(input_data,initial_delta,initial_r,initial_Lambda=function(x){x},threshold=1e-5,max_iter=5000,max_stepsize=1e5,xi=0.3,contraction=0.5){
  # requireNamespace()
  # requireNamespace("Matrix")
  # requireNamespace("logOfGamma")
  if(sum(is.na(input_data))>0){
    stop("The input data frame contains NA.")
  }
  eps<-10^(-8)
  samplesize<-length(unique(input_data[,1]))
  measure_time<-rep(0,samplesize)
  for(i in 1:samplesize){
    tmp<-input_data[input_data[,1]==i,]
    measure_time[i]<-length(tmp[,1])
  }
  max_measure<-max(measure_time)
  timepoints<-matrix(rep(-1,samplesize*max_measure),samplesize,max_measure)
  transpoints<-matrix(rep(-1,samplesize*max_measure),samplesize,max_measure)
  measure_results<-matrix(rep(-1,samplesize*max_measure),samplesize,max_measure)
  for(i in 1:samplesize){
    tmp<-input_data[input_data[,1]==i,]
    timepoints[i,1:measure_time[i]]<-tmp[,2]
    measure_results[i,1:measure_time[i]]<-tmp[,3]
  }
  uniquepoints<-sort(unique(as.vector(timepoints)))[-1]
  for(i in 1:samplesize){
    for(j in 1:measure_time[i]){
      transpoints[i,j]<-which(timepoints[i,j]==uniquepoints)
    }
  }
  n<-length(uniquepoints)
  ###################################################################################################################
  explc_grad<-function(Theta,delta=numeric(0),r=numeric(0)){
    for(i in 1:n){
      Theta[n-i+1]<-sum(Theta[1:(n-i+1)])
    }
    mygrad<-rep(0,n)
    nonzero<-diff(AA@p)
    for(q in 1:n){
      if(q>1){
        for(i in 1:(q-1)){
          if(nonzero[i]>0){
            for(j in (AA@p[i]+1):AA@p[i+1]){
              if((AA@i[j]+1)>=q){
                mygrad[q]<-mygrad[q]+AA@x[j]/(Theta[AA@i[j]+1]-Theta[i])
              }
            }
          }
        }
      }
      for(i in q:n){
        mygrad[q]<-mygrad[q]+A[i,i]/Theta[i]
      }
    }
    for(i in 1:samplesize){
      deltai<-delta+measure_time[i]/2
      ri<-r+(measure_results[i,1]-Theta[transpoints[i,1]])^2/(2*measure_results[i,1])
      if(measure_time[i]>1){
        for(j in 2:measure_time[i]){
          ri<-ri+(measure_results[i,j]-measure_results[i,j-1]-(Theta[transpoints[i,j]]-Theta[transpoints[i,j-1]]))^2/(2*(measure_results[i,j]-measure_results[i,j-1]))
        }
      }

      for(j in 1:transpoints[i,1]){
        mygrad[j]<-mygrad[j]-(deltai*(Theta[transpoints[i,1]])/(measure_results[i,1]))/ri
      }
      if(measure_time[i]>1){
        for(j in 2:measure_time[i]){
          for(k in (transpoints[i,j-1]+1):transpoints[i,j]){
            mygrad[k]<-mygrad[k]-(deltai*(Theta[transpoints[i,j]]-Theta[transpoints[i,j-1]])/(measure_results[i,j]-measure_results[i,j-1]))/ri
          }
        }
      }

      for(j in 1:transpoints[i,measure_time[i]]){
        mygrad[j]<-mygrad[j]+deltai/(ri)
      }
    }
    return(mygrad)
  }
  explc_hess<-function(Theta,delta=numeric(0),r=numeric(0)){
    for(i in 1:n){
      Theta[n-i+1]<-sum(Theta[1:(n-i+1)])
    }
    mygrad<-rep(0,n)
    nonzero<-diff(AA@p)
    for(q in 1:n){
      if(q>1){
        for(i in 1:(q-1)){
          if(nonzero[i]>0){
            for(j in (AA@p[i]+1):AA@p[i+1]){
              if((AA@i[j]+1)>=q){
                mygrad[q]<-mygrad[q]-AA@x[j]/(Theta[AA@i[j]+1]-Theta[i])^2
              }
            }
          }

        }
      }

      for(i in q:n){
        mygrad[q]<-mygrad[q]-A[i,i]/Theta[i]^2
      }
    }
    for(i in 1:samplesize){
      deltai<-delta+measure_time[i]/2
      ri<-r+(measure_results[i,1]-Theta[transpoints[i,1]])^2/(2*measure_results[i,1])
      if(measure_time[i]>1){
        for(j in 2:measure_time[i]){
          ri<-ri+(measure_results[i,j]-measure_results[i,j-1]-(Theta[transpoints[i,j]]-Theta[transpoints[i,j-1]]))^2/(2*(measure_results[i,j]-measure_results[i,j-1]))
        }
      }

      for(j in 1:transpoints[i,1]){
        mygrad[j]<-mygrad[j]-(deltai*ri/measure_results[i,1]-deltai*(Theta[transpoints[i,1]]/measure_results[i,1]-1)^2)/ri^2
      }
      if(measure_time[i]>1){
        for(j in 2:measure_time[i]){
          for(k in (transpoints[i,j-1]+1):transpoints[i,j]){
            mygrad[k]<-mygrad[k]-(deltai*ri/(measure_results[i,j]-measure_results[i,j-1])-deltai*((Theta[transpoints[i,j]]-Theta[transpoints[i,j-1]])/(measure_results[i,j]-measure_results[i,j-1])-1)^2)/ri^2
          }
        }
      }

    }
    return(mygrad)
  }
  loglikelihood1<-function(par,Theta=numeric(0)){
    if(min(par)<0){
      return(Inf)
    }else{
      delta<-par[1]
      r<-par[2]
      for(i in 1:n){
        Theta[n-i+1]<-sum(Theta[1:(n-i+1)])
      }
      lk<-0
      for(i in 1:samplesize){
        deltai<-delta+measure_time[i]/2
        ri<-r+(measure_results[i,1]-Theta[transpoints[i,1]])^2/(2*measure_results[i,1])
        if(measure_time[i]>1){
          for(j in 2:measure_time[i]){
            ri<-ri+(measure_results[i,j]-measure_results[i,j-1]-(Theta[transpoints[i,j]]-Theta[transpoints[i,j-1]]))^2/(2*(measure_results[i,j]-measure_results[i,j-1]))
          }
        }

        lk<-lk+logOfGamma::gammaln(deltai)-logOfGamma::gammaln(delta)+delta*log(r)-deltai*log(ri)

      }
      return(-lk)
    }

  }
  loglikelihood<-function(Theta,delta=numeric(0),r=numeric(0)){
    for(i in 1:n){
      Theta[n-i+1]<-sum(Theta[1:(n-i+1)])
    }
    lk<-0
    for(i in 1:samplesize){
      deltai<-delta+measure_time[i]/2
      ri<-r+(measure_results[i,1]-Theta[transpoints[i,1]])^2/(2*measure_results[i,1])
      if(measure_time[i]>1){
        for(j in 2:measure_time[i]){
          ri<-ri+(measure_results[i,j]-measure_results[i,j-1]-(Theta[transpoints[i,j]]-Theta[transpoints[i,j-1]]))^2/(2*(measure_results[i,j]-measure_results[i,j-1]))
        }
      }

      lk<-lk+logOfGamma::gammaln(deltai)-logOfGamma::gammaln(delta)+delta*log(r)-deltai*log(ri)+log(Theta[transpoints[i,1]])
      if(measure_time[i]>1){
        for(j in 2:measure_time[i]){
          lk<-lk+log(Theta[transpoints[i,j]]-Theta[transpoints[i,j-1]])
        }
      }
    }
    return(lk)
  }
  A<-matrix(rep(0,n^2),n,n)
  for(i in 1:samplesize){
    A[transpoints[i,1],transpoints[i,1]]<-A[transpoints[i,1],transpoints[i,1]]+1
    if(measure_time[i]>1){
      for(j in 2:measure_time[i]){
        A[transpoints[i,j-1],transpoints[i,j]]<-A[transpoints[i,j-1],transpoints[i,j]]+1
      }
    }

  }
  AA<-Matrix::Matrix(t(A))
  Int_Lambda<-initial_Lambda(uniquepoints)
  ###################################################BCA1SG Algorithm
  old_Theta<-Int_Lambda
  old_Theta<-c(old_Theta[1],diff(old_Theta))
  old_Theta[old_Theta<=0]<-eps
  old_delta<-initial_delta
  old_r<-initial_r
  increment<-Inf
  active_set<-numeric(0)
  count<-0
  ptm<-proc.time()
  flag<-TRUE
  while((flag)&(count<max_iter)){
    line_count<-0
    movement<-Inf
    while((movement>threshold)&(count<max_iter)){
      res<-stats::optim(par=c(old_delta,old_r),fn=loglikelihood1,Theta=old_Theta)$par
      movement1<-old_delta-res[1]
      movement2<-old_r-res[2]
      last_delta<-old_delta
      last_r<-old_r
      old_delta<-res[1]
      old_r<-res[2]
      hess<-explc_hess(old_Theta,delta=old_delta,r=old_r)
      mygrad<-explc_grad(old_Theta,delta=old_delta,r=old_r)
      ###############compute the movement direction
      direct<--1/hess*mygrad
      ###############project the direction
      if(length(active_set)>0){
        for(i in 1:length(active_set)){
          direct[active_set[i]]<-0
        }
      }
      #direct_norm1<-max(abs(c(direct,movement1,movement2)))
      count<-count+1
      ###############compute the largest step length
      d<-direct[direct<0]
      old<-old_Theta[direct<0]
      if(length(d)>0){
        beta<-min((eps-old)/d)
      }else{
        beta<-1
      }
      beta<-min(beta,max_stepsize)
      direct<-direct*beta
      ###############compute the new_Theta candidate
      new_Theta<-old_Theta+direct
      ########line search
      line_count<-0
      while((loglikelihood(new_Theta,delta=old_delta,r=old_r)<(loglikelihood(old_Theta,delta=old_delta,r=old_r)+xi*mygrad%*%direct))&&(line_count<100)){
        beta<-beta*contraction
        direct<-direct*contraction
        new_Theta<-old_Theta+direct
        line_count<-line_count+1
      }
      ###############update the active set
      if(sum(((new_Theta-eps)<=10^(-6)))>0){
        active_set<-which(((new_Theta-eps)<=10^(-6)))
      }
      movement<-max(abs(c(new_Theta-old_Theta,old_delta-last_delta,old_r-last_r)))
      # print(c(count,movement))
      old_Theta<-new_Theta
    }


    ###############decide whether to remove some point from the boundary
    if(length(active_set)>0){
      lagrg<-rep(0,length(active_set))
      for(i in 1:length(active_set)){
        lagrg[i]<-mygrad[active_set[i]]
      }
      if(max(lagrg)>0){
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
  for(i in 1:length(new_Theta)){
    new_Theta[length(new_Theta)-i+1]<-sum(new_Theta[1:(length(new_Theta)-i+1)])
  }
  timecost<-as.numeric((proc.time()-ptm)[3])
  if(count>=max_iter){
    warning("The algorithm fails to converge. Please try to increase the threshold or the max_iter.")
  }
  output_MLE<-list(distinct_time=uniquepoints,est_Lambda=new_Theta,est_delta=old_delta,est_r=old_r,iteration=count,timecost=timecost)
  return(output_MLE)
}
