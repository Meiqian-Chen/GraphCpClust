#' Basic description
#'
#' @description Given the starting point, ending point, and the shortest Hamiltonian path, this function returns the change points estimation introduced in Shi, Chen, Dong and Rao (2020)
#' @usage CpCOVID(start,end,data)
#' @param start,end starting point, ending point
#' @param data input data 
#' @return the change points estimation
#' @examples
#' 
#' # The data example of change points estimation in Shi, Chen, Dong and Rao (2020).
#' data(origin_data,package="GraphCpClust")
#' CpCOVID(1,142,Data.new[,"HB"])
#' @seealso
#' @export

# ICSS algorithm
CpCOVID <- function(start,end,data){
  data = log(1+data)
  ini.par=function(beta00,beta01,beta10,beta11,x,y){
    seq1=seq(beta00,beta01,0.1)
    seq2=seq(beta10,beta11,0.1)
    M=NULL
    for(i in 1:length(seq1))
      for(j in 1:length(seq2))
      {
        fit=lm(y~pnorm(seq1[i]+x*seq2[j]))
        if(length(summary(fit)$coefficients[,1])>1)
          #the residuals, that is response minus fitted values,coefficients a named vector of coefficients
          M=rbind(M,c(sum(fit$residuals^2),seq1[i],seq2[j],summary(fit)$coefficients[,1]))
      }
    Re=M[which.min(M[,1]),2:5]
    return(Re)
  }


  # BIC
  BIC_base <- function(origindata,fitdata,para_num){
    # calculate BIC with m number
    BIC_num=log(sum((origindata-fitdata)^2)/length(origindata))+para_num*log(length(origindata))/length(origindata)
    return(BIC_num)
  }


  # fitting data
  fun_fit <- function(start,end,data){
    # fitting the origin data and return fitting results
    para=GridSearch1(start:end,data,-20,20,-20,20,100)
    # para=ini.par(beta00=-10,beta01=10,beta10=-10,beta11=10,x=c(start:end),y=data)
    # cat('this is para:',c(as.numeric(para[1]),as.numeric(para[2]),as.numeric(para[3]),as.numeric(para[4]),as.numeric(para[5]),as.numeric(para[6])),'\n')
    y <- data[3:(end-start+1)]
    z1<-data[2:(end-start+1-1)]#t-1
    z2<-data[1:(end-start+1-2)]#t-2
    x <- seq(start+2,end,1)
    fit<-try(nlm1<- nls(y ~ alpha1*z1+alpha2*z2 + b*pnorm(beta1+beta2*x)+a, start=list(a=para$a ,b=para$b 
                                                                 ,beta1=para$beta1 ,beta2=para$beta2, alpha1=para$alpha1,alpha2=para$alpha2),
                        control=nls.control(maxiter = 2000, tol = 1e-09, minFactor = 1/20240000, printEval = FALSE, warnOnly = T))
    )
    return(list(fit=fit,residual=para$residual,para=para))
  }


  # find change point
  fun_changepoint <- function(start,end,data){
    all_var = c()
    num = 10
    cat('find changepoint\n')
    for(i in (start+num):end){
      if((end-i)<=10){
        break
      }
      cat(i,'\n')
      fit<-fun_fit(start,i,data[1:(i-start+1)])
      cat('get new fitdata 1\n')
      if("try-error" %in% class(fit$fit)){
        result_a <- fit$residual
        cat('a',result_a,'\n')
      }else{
        result_a <- sum((c(data[1:2],fitted(fit$fit))-data[1:(i-start+1)])^2)
      }

      fit1<-fun_fit(i+1,end,data[(i+1-start+1):(end-start+1)])
      cat('get new fitdata 2\n')
      if("try-error" %in% class(fit1$fit)){
        result_b <- fit1$residual
        cat('b',result_b,'\n')
      }else{
        result_b<-sum((c(data[(i+1-start+1):(i+1-start+1+1)],fitted(fit1$fit))-data[(i+1-start+1):(end-start+1)])^2)
      }
      cat(result_a,result_b,'\n')
      a<-result_a+result_b
      all_var <- c(all_var, a)
    }
    if(length(all_var)==0){
      return(1)
    }else{
      point=which(all_var==min(all_var),arr.ind=TRUE)
      changepoint=point+start+num-1
      cat('new point',changepoint,'\n')
      result_a_final<-fun_fit(start,changepoint,data[1:(point+num)])$fit
      result_b_final<-fun_fit(changepoint+1,end,data[(point+num+1):(end-start+1)])$fit
      if((is.atomic(result_a_final))||(is.atomic(result_b_final))){
        return(1)
      }else{
        return(list(changepoint=changepoint, fitdata=c(data[1:2],fitted(result_a_final), data[(point+num+1):(point+num+1+1)],fitted(result_b_final))))
      }
    }
  }

  # main ICSS
  k=1
  t1=start
  t2=end
  CP <- c()
  while (k>0) {
    #step2
    cat('t1=',t1,'t2=',t2,'\n')
    result_fit <- fun_fit(t1,t2,data[t1:t2])$fit
    result_CP <- fun_changepoint(t1,t2,data[t1:t2])
    # if t1,t2 can not find change point, the program doesn't have to run
    if(is.atomic(result_CP)){
      k=-1
    }else{
      BIC0<-BIC_base(origindata=data[t1:t2],fitdata=c(data[t1:(t1+1)],fitted(result_fit)),para_num=6)
      BIC1<-BIC_base(origindata=data[t1:t2],fitdata=result_CP$fitdata,para_num=12)
      cat(BIC0,BIC1,'\n')
      if(BIC0<=BIC1){
        k = -1
      }
      else{
        t_a_low = t1
        t_a_up = result_CP$changepoint
        k1 = 1
        while(k1>0){
          # Complete the first half of the iteration process
          # calculate process
          if((t_a_up-t_a_low)>10){
            cat(t_a_low,t_a_up,'\n')
            result_a_fit <- fun_fit(t_a_low,t_a_up,data[t_a_low:t_a_up])$fit
            result_a_CP <- fun_changepoint(t_a_low,t_a_up,data[t_a_low:t_a_up])
            # cat(result_a_CP)
            # if t_a_up t_a_low can not find change point, this part finish
            if(is.atomic(result_a_CP)){
              change_point_first = t_a_up
              k1 = -1
            }else{
              BICa0<-BIC_base(origindata=data[t_a_low:t_a_up],fitdata=c(data[t_a_low:(t_a_low+1)],fitted(result_a_fit)),para_num=6)
              BICa1<-BIC_base(origindata=data[t_a_low:t_a_up],fitdata=result_a_CP$fitdata,para_num=12)
              cat(BICa0,BICa1,'\n')
              # Determine process
              if(BICa0<=BICa1){
                change_point_first = t_a_up
                k1 = -1
              }else{
                t_a_up = result_a_CP$changepoint
              }
            }
          }else{
            cat('save left change point',t_a_up - 1)
            change_point_first = t_a_up
            k1 = -1
          }
        }
        t_b_low = result_CP$changepoint + 1
        t_b_up = t2
        k2=1
        while (k2>0) {
          if((t_b_up-t_b_low)>10){
            cat(t_b_low,t_b_up,'\n')
            result_b_fit <- fun_fit(t_b_low,t_b_up,data[t_b_low:t_b_up])$fit
            result_b_CP <- fun_changepoint(t_b_low,t_b_up,data[t_b_low:t_b_up])
            if(is.atomic(result_b_CP)){
              change_point_last = t_b_low - 1
              k2 = -1
            }else{
              BICb0<-BIC_base(origindata=data[t_b_low:t_b_up],fitdata=c(data[t_b_low:(t_b_low+1)],fitted(result_b_fit)),para_num=6)
              BICb1<-BIC_base(origindata=data[t_b_low:t_b_up],fitdata=result_b_CP$fitdata,para_num=12)
              cat(BICb0,BICb1,'\n')
              # Determine process
              if(BICb0<=BICb1){
                change_point_last = t_b_low - 1
                k2 = -1
              }else{
                t_b_low = result_b_CP$changepoint + 1 # +1???
              }
            }
          }else{
            cat('save right change point',t_b_low - 1)
            change_point_last = t_b_low - 1
            k2 = -1
          }
        }
        # step4
        if(change_point_first==change_point_last){
          # return current changepoint
          CP<-c(CP,change_point_first)
          k=-1
        }else{
          CP<-c(CP,change_point_first,change_point_last)
          t1 = change_point_first
          t2 = change_point_last
          cat('next repeat\n')
        }
      }
    }
  }
  
  # step5&6
  cat(CP,'\n')
  CP_new <- CP
  if(length(CP)>=2){
    CP<-sort(c(CP,start,end))
    k3=1
    while (k3>0) {
      CP_new <- c()
      for(j in 2:(length(CP)-1)){
        # NO SURE

        fit0 = fun_fit(CP[j-1]+1,CP[j+1],data[(CP[j-1]+1):CP[j+1]])
        if("try-error" %in% class(fit0$fit)){
          result_fit = fit0$pred
        }else{
          result_fit <- fitted(fit0$fit)
        }


        fit1 = fun_fit(CP[j-1]+1,CP[j],data[(CP[j-1]+1):CP[j]])
        if("try-error" %in% class(fit1$fit)){
          result_fit1 = fit1$pred
        }else{
          result_fit1 <- fitted(fit1$fit)
        }


        fit2 = fun_fit(CP[j]+1,CP[j+1],data[(CP[j]+1):CP[j+1]])
        if("try-error" %in% class(fit2$fit)){
          result_fit2 = fit2$pred
        }else{
          result_fit2 <- fitted(fit2$fit)
        }
        

        myfit <- c(c(data[(CP[j-1]+1):(CP[j-1]+1+1)],result_fit1),c(data[(CP[j]+1):(CP[j]+1+1)],result_fit2))
        cat(myfit)
        BICc0<-BIC_base(origindata=data[(CP[j-1]+1):CP[j+1]],fitdata=c(data[(CP[j-1]+1):(CP[j-1]+1+1)],result_fit),para_num=6)
        BICc1<-BIC_base(origindata=data[(CP[j-1]+1):CP[j+1]],fitdata=myfit,para_num=12)
        cat(BICc0,BICc1,'\n')
        if(BICc1<=BICc0){
          CP_new <- c(CP_new,CP[j])
          cat(CP_new,'\n')
        }
      }
      if(length(CP)==(length(CP_new)+2)||length(CP_new)==0){
        k3 = -1
      }else{
        CP<-sort(c(CP_new,start,end))
      }
    }
  }
  return(CP_new)
}


