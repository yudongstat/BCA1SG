#' @title BCA1SG algorithm for panel count data which are modeled by the nonhomogeneous Poisson process.
#'
#' @description This function implements the BCA1SG algorithm on the semiparametric NHPP model for panel count data. Details about the algorithm and the semiparametric NHPP model can be found in Wang et al.(2020, Biometrics).
#'
#' @param input_data An object of class data.frame.
#                    The structure of the data.frame must be
#                    \{patient ID, time of measurement, measurement(cumulative counts),covariate_1,...,covariate_p\}.
#                    This data frame cannot contain missing values. See the dataset skiTum for an example.
#' @param initial_beta The initial value of the regression coefficients.
#                      The dimension of this input should comply with the dimension of the covariates.
#' @param initial_Lambda An R function which serves as the initial value of the baseline mean cumulative function.
#' @param threshold Convergence threshold. The algorithm is terminated when the infinity norm of the difference between successive iterates is less than the convergence threshold.
#' @param max_iter Maximum number of iterations allowed.
#' @param max_stepsize Maximum stepsize allowed.
#' @param xi The xi parameter in the inexact backtracking line search algorithm. See Wang et al.(2020) for details.
#' @param contraction The contraction parameter in the inexact backtracking line search algorithm. See Wang et al.(2020) for details.
#' @return A list of length 5：
#          distinct_time is the set of distinct observation time points;
#          est_Lambda is the estimated baseline mean cumulative function at the set of distinct observation time points;
#          est_beta is the estimated regression coefficients;
#          iteration is the number of iterations;
#          timecost is the computational time in seconds.
#'
#' @export stock_predict
#'
#' @examples
#'
#' # Example 1: Application on the skin tumor data set
#' data(adapt_skiTum)
#' res<-BCA1SG_NHPP(adapt_skiTum,initial_beta = rep(0,4))
#' res$est_beta
#' res$iteration
#' res$timecost
#' plot(res$distinct_time,res$est_Lambda,type="s",lwd=3)
#'




BCA1SG_NHPP<-function(input_data,initial_beta,initial_Lambda=function(x){x},threshold=1e-5,max_iter=5000,max_stepsize=1e4,xi=0.3,contraction=0.5){
  if(sum(is.na(input_data))>0){
    stop("The input data frame contains NA.")
  }
  eps<-10^(-8)
  samplesize<-length(unique(input_data[,1]))
  cov_dim<-ncol(input_data)-3
  if(cov_dim!=length(initial_beta)){
    stop("The dimension of initial_beta does not comply with the dimension of the covariates.")
  }
  measure_time<-rep(0,samplesize)
  for(i in 1:samplesize){
    tmp<-input_data[input_data[,1]==i,]
    measure_time[i]<-length(tmp[,1])
  }
  max_measure<-max(measure_time)
  timepoints<-matrix(rep(-1,samplesize*max_measure),samplesize,max_measure)
  transpoints<-matrix(rep(-1,samplesize*max_measure),samplesize,max_measure)
  measure_results<-matrix(rep(-1,samplesize*max_measure),samplesize,max_measure)
  covariate<-matrix(rep(0,cov_dim*samplesize),samplesize,cov_dim)
  for(i in 1:samplesize){
    tmp<-input_data[input_data[,1]==i,]
    timepoints[i,1:measure_time[i]]<-tmp[,2]
    measure_results[i,1:measure_time[i]]<-tmp[,3]
    covariate[i,]<-as.numeric(tmp[1,-c(1,2,3)])
  }
  uniquepoints<-sort(unique(as.vector(timepoints)))[-1]
  for(i in 1:samplesize){
    for(j in 1:measure_time[i]){
      transpoints[i,j]<-which(timepoints[i,j]==uniquepoints)
    }
  }
  n<-length(uniquepoints)


  loglikelihood<-function(Lambda,beta){
    for(i in 1:n){
      Lambda[n-i+1]<-sum(Lambda[1:(n-i+1)])
    }
    term1<-covariate%*%beta
    lk<-0
    for(i in 1:samplesize){
      lk<-lk+log(Lambda[transpoints[i,1]])*measure_results[i,1]
      if(measure_time[i]>1){
        for(j in 2:measure_time[i]){
          lk<-lk+log(Lambda[transpoints[i,j]]-Lambda[transpoints[i,j-1]])*(measure_results[i,j]-measure_results[i,j-1])
        }
      }
      lk<-lk+term1[i]*measure_results[i,measure_time[i]]-exp(term1[i])*Lambda[transpoints[i,measure_time[i]]]
    }
    return(lk)
  }
  explc_hess<-function(Lambda,beta=numeric(0)){
    mygrad<-rep(0,n)
    for(i in 1:n){
      Lambda[n-i+1]<-sum(Lambda[1:(n-i+1)])
    }
    for(i in 1:samplesize){
      for(k in 1:transpoints[i,1]){
        mygrad[k]<-mygrad[k]-measure_results[i,1]/(Lambda[transpoints[i,1]])^2
      }
      if(measure_time[i]>1){
        for(j in 2:measure_time[i]){
          for(k in (transpoints[i,j-1]+1):transpoints[i,j])
            mygrad[k]<-mygrad[k]-(measure_results[i,j]-measure_results[i,j-1])/(Lambda[transpoints[i,j]]-Lambda[transpoints[i,j-1]])^2
        }
      }

    }
    return(mygrad)
  }
  explc_grad<-function(Lambda,beta=numeric(0)){
    mygrad<-rep(0,n)
    term1<-covariate%*%beta
    for(i in 1:n){
      Lambda[n-i+1]<-sum(Lambda[1:(n-i+1)])
    }
    for(i in 1:samplesize){
      for(k in 1:transpoints[i,1]){
        mygrad[k]<-mygrad[k]+measure_results[i,1]/(Lambda[transpoints[i,1]])-exp(term1[i])
      }
      if(measure_time[i]>1){
        for(j in 2:measure_time[i]){
          for(k in (transpoints[i,j-1]+1):transpoints[i,j])
            mygrad[k]<-mygrad[k]+(measure_results[i,j]-measure_results[i,j-1])/(Lambda[transpoints[i,j]]-Lambda[transpoints[i,j-1]])-exp(term1[i])
        }
      }
    }
    return(mygrad)
  }
  loglikelihood0<-function(beta,LL=numeric(0),NN=numeric(0)){
    term1<-covariate%*%beta
    lk<-sum(NN*term1-exp(term1)*LL)
    return(-lk)
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
  while((flag)&&(count<max_iter)){
    movement<-Inf
    while((movement>threshold)&&(count<max_iter)){
      ##########first update the parametric part
      LL<-sum(old_Lambda[1:transpoints[1,measure_time[1]]])
      NN<-measure_results[1,measure_time[1]]
      for(i in 2:samplesize){
        LL<-c(LL,sum(old_Lambda[1:transpoints[i,measure_time[i]]]))
        NN<-c(NN,measure_results[i,measure_time[i]])
      }
      res<-stats::optim(par=old_beta,fn=loglikelihood0,LL=LL,NN=NN)$par
      movement1<-max(abs(old_beta-res))
      old_beta<-res
      ###########then update the nonparametric part
      hess<-explc_hess(old_Lambda,beta=old_beta)
      hess<-hess-rep(eps,n)
      mygrad<-explc_grad(old_Lambda,beta=old_beta)
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
      if(sum(((new_Lambda-eps)<=10^(-6)))>0){
        active_set<-which(((new_Lambda-eps)<=10^(-6)))
      }
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
      if(max(lagrg)>10^(-5)){
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
