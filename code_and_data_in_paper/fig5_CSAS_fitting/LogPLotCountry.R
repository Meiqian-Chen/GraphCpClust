library(readr)
library(ggplot2)
library(purrr)
library(stringr)
library(dplyr)
library(grid)
library(gridExtra)
library(tidyr)
# library(gridsearch)


setwd("")

# initial parameter satting
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("test_new_a.cpp")

load('cluster_data.RData')
load('origin_data.RData')

#' Basic description
#'
#' @description This function is used to fit the log-transformed infection counts of 33 countries in the rest of world that introduced in Shi, Chen, Dong and Rao (2020)
#' @usage LogPLotCountry(data_clust_C, final_C)
#' @param data_clust_C 
#' @param final_C 
#' @examples
#' 
#' 
#' 
#' LogPLotCountry(data_clust_C, final_C)
#' @seealso 
#' @export
plot_res_C=LogPLotCountry(data_clust_C, final_C)
residual_list_C_origin=plot_res_C$residual_list
pred_plot_C=plot_res_C$pred_plot
pred_plot_sim_C=plot_res_C$pred_plot_sim
# save(pred_plot_C,pred_plot_sim_C,data_clust_C,final_C,file='')
# save(residual_list_C_origin,file='')

LogPLotCountry <- function(data_clust_C, final_C){
  residual_list=list()
  pred_plot=list()
  pred_plot_sim=list()
  mytheme <- theme(plot.title = element_text(hjust = 0.5, size = 22),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.text.x = element_text(size = 18,color='black'),
                 axis.text.y = element_text(size = 18,color='black'),
                 axis.title.x=element_text(size=25),
                 axis.title.y=element_text(size=25),
                 plot.margin=unit(c(2,2,2,2), 'lines'))
  # fitting data
  fun_fit <- function(start,end,data){
    # fitting the origin data and return fitting results
    para=GridSearch1(start:end,data,-20,20,-20,20,100)
    cat(c(as.numeric(para[1]),as.numeric(para[2]),as.numeric(para[3]),as.numeric(para[4]),as.numeric(para[5]),as.numeric(para[6])),'\n')
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

  #-----------------------------------------------------------------------------------------------#
  # this section aim to predict the next part according to previous part
  predict_fun<-function(a, b, beta1, beta2, alpha1,alpha2,x.new,z1_num,z2_num){
    y_list = c()
    for (x in x.new) {
      cat('x=',x,'this is z2_num:',z2_num,'this is z1_num',z1_num,'\n')
      y=alpha1*z1_num+alpha2*z2_num+b*pnorm(beta1+beta2*x)+a
      cat('prediction result is:',y,'\n')
      y_list=c(y_list,y)
      z2_num = z1_num
      z1_num = y
    }
    cat('res is:',y_list,'\n')
    return(y_list)
  }
  #-----------------------------------------------------------------------------------------------#


  Date = rownames(data_clust_C)
  # Cluster 1(CN) AND Cluster 2(IR,IT)
  # cluster1data
  data=log(1+data_clust_C[1:dim(data_clust_C)[1],'cluster1'])
  data1<-data
  qq <- c(1,sort(final_C[[1]]),dim(data_clust_C)[1])
  # plot fitting data
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  phase<-c('Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1')
  phase_num<-c(13,34,10,16,16,53)
  fig12<-ggplot(dataplot,aes(Date))
  fig12 <- fig12 +
    geom_point(data = dataplot, aes(y=y)) +
    labs(title = expression(italic(paste('Cluster 1(CN) and Cluster 2(IR,IT)'))),
         x='', y=expression(paste("log(1+infection)")),fill="") +
    mytheme
  fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
  df.delta_data1 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                               lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
  
  #-----------------------------------------------------------------------------------------------#
  # this section aim to predict the next part according to previous part
  pred_x.new_1=list()
  pred_f.new_1=list()
  if(length(final)!=0){
    for (item in 1:(length(qq)-2)) {
      sat=item
      a<-as.numeric(qq[sat])+1
      b=as.numeric(qq[sat+1])
      x.new <- seq((b+1), qq[sat+2], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1_num<-data[b]#t-1
      z2_num<-data[b-1]#t-2
      pred_x.new_1[[item]] <- dataplot[(b+1):qq[sat+2],2]
      cat('part',item+1,'a=',beta2.est[1],'b=',beta2.est[2],'beta1=',beta2.est[3],'beta2=',beta2.est[4],'alpha1=', beta2.est[5],'alpha2=',beta2.est[6],'\n')
      pred_f.new_1[[item]] <- predict_fun(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1_num,z2_num)
    }
    
  }
  #-----------------------------------------------------------------------------------------------#
  
  
  residual_list_1=list()
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
      nlm <- fun_fit(a,b,data[a:b])$fit
      residual_list_1[[1]]=c(data[a:(a+1)],fitted(nlm))-data[a:b]
      data_one <- data.frame(x.new=dataplot[a:b,2],f.new=c(data[a:(a+1)],fitted(nlm)),lwr.conf=c(data[a:(a+1)],fitted(nlm))
                             ,upr.conf=c(data[a:(a+1)],fitted(nlm)),lwr.pred=c(data[a:(a+1)],fitted(nlm)),upr.pred=c(data[a:(a+1)],fitted(nlm)))
      data_one['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data1 <- rbind(df.delta_data1,data_one)
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
      x.new <- seq((a+2), qq[i], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1<-data[(a+1):(b-1)]#t-1
      z2<-data[(a):(b-2)]#t-2
      f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1,z2)
      residual_list_1[[i-1]]=c(data[a:(a+1)],f.new)-data[a:b]
      print(class(f.new))
      g.new <- attr(f.new,"gradient")
      # f.new <-as.numeric(c(f.new[1:2],f.new))
      # g.new<-rbind(g.new[1:2,],g.new)
      V.beta2 <- vcov(nlm2)
      GS=rowSums((g.new[,1:6]%*%V.beta2)%*%t(g.new[,1:6]))
      head(GS)
      alpha <- 0.05
      df <- length(data[a:b])-length(beta2.est)#freedom degree
      deltaf <- sqrt(GS)*qt(1-alpha/2,df)
      if(max(deltaf)>20){
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new), upr.conf=c(data[a:(a+1)],f.new))
      }else{
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new-deltaf), upr.conf=c(data[a:(a+1)],f.new+deltaf))
      }
      head(df.delta)
      sigma2.est <- summary(nlm2)$sigma
      deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
      line_dataframa<-data.frame(x.new=dataplot[a:b,2], f.new=as.numeric(c(data[a:(a+1)],f.new)))
      if(max(deltay)>20){
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new),c(data[a:(a+1)],f.new))
      }else{
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new - deltay),c(data[a:(a+1)],f.new + deltay))
      }
      df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data1<-rbind(df.delta_data1,df.delta)
    }
  }
  residual_list[[1]]=residual_list_1
  #cluster2 data
  data=log(1+data_clust_C[1:dim(data_clust_C)[1],'cluster2'])
  data2<-data
  qq <- c(1,sort(final_C[[2]]),dim(data_clust_C)[1])
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  phase<-c('Cluster 2','Cluster 2','Cluster 2')
  phase_num<-c(80,12,50)
  fig12 <- fig12 +
    geom_point(data = dataplot, aes(y=y)) +
    mytheme
  fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
  df.delta_data2 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                               lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
  #-----------------------------------------------------------------------------------------------#
  # this section aim to predict the next part according to previous part
  pred_x.new_2=list()
  pred_f.new_2=list()
  if(length(final)!=0){
    for (item in 1:(length(qq)-2)) {
      sat=item
      a<-as.numeric(qq[sat])+1
      b=as.numeric(qq[sat+1])
      x.new <- seq((b+1), qq[sat+2], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1_num<-data[b]#t-1
      z2_num<-data[b-1]#t-2
      pred_x.new_2[[item]] <- dataplot[(b+1):qq[sat+2],2]
      cat('part',item+1,'a=',beta2.est[1],'b=',beta2.est[2],'beta1=',beta2.est[3],'beta2=',beta2.est[4],'alpha1=', beta2.est[5],'alpha2=',beta2.est[6],'\n')
      pred_f.new_2[[item]] <- predict_fun(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1_num,z2_num)
    }
    
  }
  #-----------------------------------------------------------------------------------------------#
  
  residual_list_2=list()
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
    }
    x.new <- seq((a+2), qq[i], by=1)
    nlm2 <- fun_fit(a,b,data[a:b])$fit
    beta2.est <- coef(nlm2)
    z1<-data[(a+1):(b-1)]#t-1
    z2<-data[(a):(b-2)]#t-2
    f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1,z2)
    residual_list_2[[i-1]]=c(data[a:(a+1)],f.new)-data[a:b]
    print(class(f.new))
    g.new <- attr(f.new,"gradient")
    # f.new <-as.numeric(c(f.new[1:2],f.new))
    # g.new<-rbind(g.new[1:2,],g.new)
    V.beta2 <- vcov(nlm2)
    GS=rowSums((g.new[,1:6]%*%V.beta2)%*%t(g.new[,1:6]))
    head(GS)
    alpha <- 0.05
    df <- length(data[a:b])-length(beta2.est)#freedom degree
    deltaf <- sqrt(GS)*qt(1-alpha/2,df)
    if(max(deltaf)>20){
      df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new), upr.conf=c(data[a:(a+1)],f.new))
    }else{
      df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new-deltaf), upr.conf=c(data[a:(a+1)],f.new+deltaf))
    }
    head(df.delta)
    sigma2.est <- summary(nlm2)$sigma
    deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
    line_dataframa<-data.frame(x.new=dataplot[a:b,2], f.new=as.numeric(c(data[a:(a+1)],f.new)))
    if(max(deltay)>20){
      df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new),c(data[a:(a+1)],f.new))
    }else{
      df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new - deltay),c(data[a:(a+1)],f.new + deltay))
    }
    df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
    
    df.delta_data2<-rbind(df.delta_data2,df.delta)
  }
  pred_plot_sim[[1]]=df.delta_data2
  residual_list[[2]]=residual_list_2
  df.delta_data<-rbind(df.delta_data1,df.delta_data2)
  fig12 <- fig12 +
    # geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.pred, ymax=upr.pred, group=phase), alpha=0.2, fill='black') +
    geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.conf, ymax=upr.conf, group=phase,fill=phase), alpha=0.4) +
    geom_line(data=df.delta_data, aes(x=x.new, y=f.new,colour=phase, group=phase), size=1)+
    ylim(-0.1,14.5)+
    scale_fill_manual(values=c('#FF8C00','#006400'),name='Province')+
    scale_colour_manual(values=c('#FF8C00','#006400'),name='Province')+
    #cluster1
    geom_line(data=data.frame(x=c(date[13],date[13]),y=c(data1[13]-1,data1[13]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#FF8C00')+
    geom_line(data=data.frame(x=c(date[47],date[47]),y=c(data1[47]-1,data1[47]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#FF8C00')+
    geom_line(data=data.frame(x=c(date[57],date[57]),y=c(data1[57]-1,data1[57]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#FF8C00')+
    geom_line(data=data.frame(x=c(date[73],date[73]),y=c(data1[73]-1,data1[73]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#FF8C00')+
    geom_line(data=data.frame(x=c(date[89],date[89]),y=c(data1[89]-1,data1[89]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#FF8C00')+
    #cluster2
    geom_line(data=data.frame(x=c(date[80],date[80]),y=c(data2[80]-1,data2[80]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#006400')+
    geom_line(data=data.frame(x=c(date[92],date[92]),y=c(data2[92]-1,data2[92]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#006400')
  fig12 <- fig12 +
    scale_x_date(breaks =as.Date(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),date_labels="%y/%m/%d")+
    # cluster1
    annotate(geom='text',x=as.Date(c("2019-12-29")),y=1.7,label="2019-12-13",size=7)+
    annotate(geom='text',x=as.Date(c("2020-02-01")),y=3.5,label="2020-01-16",size=7)+
    annotate(geom='text',x=as.Date(c("2020-01-12")),y=8,label="2020-01-26",size=7)+
    annotate(geom='text',x=as.Date(c("2020-01-26")),y=10.9,label="2020-02-11",size=7)+
    annotate(geom='text',x=as.Date(c("2020-03-12")),y=12,label="2020-02-27",size=7)+
    #cluster2
    annotate(geom='text',x=as.Date(c("2020-03-05")),y=1,label="2020-02-18",size=7)+
    annotate(geom='text',x=as.Date(c("2020-03-17")),y=7.5,label="2020-03-01",size=7)+
    theme_bw() +
    mytheme+
    theme(legend.position=c(0.15,0.85))+theme(legend.title = element_blank())+
    theme(legend.background = element_blank())+
    theme(legend.text = element_text(size=20))+
    theme(legend.key.size = unit(40, "pt"))+
    labs(tag = "A") +
    theme(plot.tag.position = c(0.05, 1))+
    theme(plot.tag=element_text(size = 30))
  
  # if(length(pred_x.new_1)!=0){
  #   for (i in 1:length(pred_x.new_1)) {
  #     cat(pred_x.new_1[[i]],'\n')
  #     pred_data=data.frame(x.new=pred_x.new_1[[i]],f.new=pred_f.new_1[[i]])
  #     fig12 <- fig12 +
  #       geom_line(data=pred_data, aes(x=x.new, y=f.new),shape = 2, size=1,colour='#F8766D')
  #   }
  # }
  # 
  # if(length(pred_x.new_2)!=0){
  #   pred_plot_test=list()
  #   for (i in 1:length(pred_x.new_2)) {
  #     cat(pred_x.new_2[[i]],'\n')
  #     pred_data=data.frame(x.new=pred_x.new_2[[i]],f.new=pred_f.new_2[[i]])
  #     pred_plot_test[[i]]=pred_data
  #     fig12 <- fig12 +
  #       geom_line(data=pred_data, aes(x=x.new, y=f.new),shape = 2, size=1,colour='#B79F00')
  #   }
  #   pred_plot[[1]]=pred_plot_test
  # }

  fig12

  # cluster3data
  data=log(1+data_clust_C[1:dim(data_clust_C)[1],'cluster3'])
  data3<-data
  qq <- c(1,sort(final_C[[3]]),dim(data_clust_C)[1])
  # plot fitting data
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  phase<-c('Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3')
  phase_num<-c(54,11,19,11,47)
  fig34<-ggplot(dataplot,aes(Date))
  fig34 <- fig34 +
    geom_point(data = dataplot, aes(y=y)) +
    labs(title = expression(italic(paste('Cluster 3(KR) and Cluster 4(JP)'))),
         x='', y=expression(paste("log(1+infection)")),fill="") +
    mytheme
  fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
  df.delta_data3 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                               lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
  
  #-----------------------------------------------------------------------------------------------#
  # this section aim to predict the next part according to previous part
  pred_x.new_3=list()
  pred_f.new_3=list()
  if(length(final)!=0){
    for (item in 1:(length(qq)-2)) {
      sat=item
      a<-as.numeric(qq[sat])+1
      b=as.numeric(qq[sat+1])
      x.new <- seq((b+1), qq[sat+2], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1_num<-data[b]#t-1
      z2_num<-data[b-1]#t-2
      pred_x.new_3[[item]] <- dataplot[(b+1):qq[sat+2],2]
      cat('part',item+1,'a=',beta2.est[1],'b=',beta2.est[2],'beta1=',beta2.est[3],'beta2=',beta2.est[4],'alpha1=', beta2.est[5],'alpha2=',beta2.est[6],'\n')
      pred_f.new_3[[item]] <- predict_fun(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1_num,z2_num)
    }
    
  }
  #-----------------------------------------------------------------------------------------------#
  
  residual_list_3=list()
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
      nlm <- fun_fit(a,b,data[a:b])$fit
      residual_list_3[[1]]=c(data[a:(a+1)],fitted(nlm))-data[a:b]
      data_one <- data.frame(x.new=dataplot[a:b,2],f.new=c(data[a:(a+1)],fitted(nlm)),lwr.conf=c(data[a:(a+1)],fitted(nlm))
                             ,upr.conf=c(data[a:(a+1)],fitted(nlm)),lwr.pred=c(data[a:(a+1)],fitted(nlm)),upr.pred=c(data[a:(a+1)],fitted(nlm)))
      data_one['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data3 <- rbind(df.delta_data3,data_one)
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
      x.new <- seq((a+2), qq[i], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1<-data[(a+1):(b-1)]#t-1
      z2<-data[(a):(b-2)]#t-2
      f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1,z2)
      residual_list_3[[i-1]]=c(data[a:(a+1)],f.new)-data[a:b]
      print(class(f.new))
      g.new <- attr(f.new,"gradient")
      # f.new <-as.numeric(c(f.new[1:2],f.new))
      # g.new<-rbind(g.new[1:2,],g.new)
      V.beta2 <- vcov(nlm2)
      GS=rowSums((g.new[,1:6]%*%V.beta2)%*%t(g.new[,1:6]))
      head(GS)
      alpha <- 0.05
      df <- length(data[a:b])-length(beta2.est)#freedom degree
      deltaf <- sqrt(GS)*qt(1-alpha/2,df)
      if(max(deltaf)>20){
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new), upr.conf=c(data[a:(a+1)],f.new))
      }else{
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new-deltaf), upr.conf=c(data[a:(a+1)],f.new+deltaf))
      }
      head(df.delta)
      sigma2.est <- summary(nlm2)$sigma
      deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
      line_dataframa<-data.frame(x.new=dataplot[a:b,2], f.new=as.numeric(c(data[a:(a+1)],f.new)))
      if(max(deltay)>20){
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new),c(data[a:(a+1)],f.new))
      }else{
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new - deltay),c(data[a:(a+1)],f.new + deltay))
      }
      df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data3<-rbind(df.delta_data3,df.delta)
    }
  }
  residual_list[[3]]=residual_list_3
  #cluster4 data
  data=log(1+data_clust_C[1:dim(data_clust_C)[1],'cluster4'])
  data4<-data
  qq <- c(1,sort(final_C[[4]]),dim(data_clust_C)[1])
  # plot fitting data
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  phase<-c('Cluster 4','Cluster 4','Cluster 4','Cluster 4','Cluster 4')
  phase_num<-c(50,16,11,20,45)#   c(66,11,20,45)  c(50,16,11,20,45)
  fig34 <- fig34 +
    geom_point(data = dataplot, aes(y=y)) +
    mytheme
  fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
  df.delta_data4 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                               lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
  
  #-----------------------------------------------------------------------------------------------#
  # this section aim to predict the next part according to previous part
  pred_x.new_4=list()
  pred_f.new_4=list()
  if(length(final)!=0){
    for (item in 1:(length(qq)-2)) {
      sat=item
      a<-as.numeric(qq[sat])+1
      b=as.numeric(qq[sat+1])
      x.new <- seq((b+1), qq[sat+2], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1_num<-data[b]#t-1
      z2_num<-data[b-1]#t-2
      pred_x.new_4[[item]] <- dataplot[(b+1):qq[sat+2],2]
      cat('part',item+1,'a=',beta2.est[1],'b=',beta2.est[2],'beta1=',beta2.est[3],'beta2=',beta2.est[4],'alpha1=', beta2.est[5],'alpha2=',beta2.est[6],'\n')
      pred_f.new_4[[item]] <- predict_fun(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1_num,z2_num)
    }
    
  }
  #-----------------------------------------------------------------------------------------------#
  
  residual_list_4=list()
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
      nlm <- fun_fit(a,b,data[a:b])$fit
      residual_list_4[[1]]=c(data[a:(a+1)],fitted(nlm))-data[a:b]
      data_one <- data.frame(x.new=dataplot[a:b,2],f.new=c(data[a:(a+1)],fitted(nlm)),lwr.conf=c(data[a:(a+1)],fitted(nlm))
                             ,upr.conf=c(data[a:(a+1)],fitted(nlm)),lwr.pred=c(data[a:(a+1)],fitted(nlm)),upr.pred=c(data[a:(a+1)],fitted(nlm)))
      data_one['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data4 <- rbind(df.delta_data4,data_one)
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
      x.new <- seq((a+2), qq[i], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1<-data[(a+1):(b-1)]#t-1
      z2<-data[(a):(b-2)]#t-2
      f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1,z2)
      residual_list_4[[i-1]]=c(data[a:(a+1)],f.new)-data[a:b]
      print(class(f.new))
      g.new <- attr(f.new,"gradient")
      # f.new <-as.numeric(c(f.new[1:2],f.new))
      # g.new<-rbind(g.new[1:2,],g.new)
      V.beta2 <- vcov(nlm2)
      GS=rowSums((g.new[,1:6]%*%V.beta2)%*%t(g.new[,1:6]))
      head(GS)
      alpha <- 0.05
      df <- length(data[a:b])-length(beta2.est)#freedom degree
      deltaf <- sqrt(GS)*qt(1-alpha/2,df)
      if(max(deltaf)>20){
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new), upr.conf=c(data[a:(a+1)],f.new))
      }else{
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new-deltaf), upr.conf=c(data[a:(a+1)],f.new+deltaf))
      }
      head(df.delta)
      sigma2.est <- summary(nlm2)$sigma
      deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
      line_dataframa<-data.frame(x.new=dataplot[a:b,2], f.new=as.numeric(c(data[a:(a+1)],f.new)))
      if(max(deltay)>20){
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new),c(data[a:(a+1)],f.new))
      }else{
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new - deltay),c(data[a:(a+1)],f.new + deltay))
      }
      df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data4<-rbind(df.delta_data4,df.delta)
    }
  }
  residual_list[[4]]=residual_list_4
  df.delta_data<-rbind(df.delta_data3,df.delta_data4)
  fig34 <- fig34 +
    # geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.pred, ymax=upr.pred, group=phase), alpha=0.2, fill='black') +
    geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.conf, ymax=upr.conf, group=phase,fill=phase), alpha=0.4) +
    geom_line(data=df.delta_data, aes(x=x.new, y=f.new,colour=phase, group=phase), size=1)+
    ylim(-0.1,14.5)+
    scale_fill_manual(values=c('#FF8C00','#006400'),name='Province')+
    scale_colour_manual(values=c('#FF8C00','#006400'),name='Province')+
    #cluster3
    geom_line(data=data.frame(x=c(date[54],date[54]),y=c(data3[54]-0.5,data3[54]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#FF8C00')+
    geom_line(data=data.frame(x=c(date[65],date[65]),y=c(data3[65]-1,data3[65]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#FF8C00')+
    geom_line(data=data.frame(x=c(date[84],date[84]),y=c(data3[84]-1,data3[84]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#FF8C00')+
    geom_line(data=data.frame(x=c(date[95],date[95]),y=c(data3[95]-1,data3[95]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#FF8C00')+
    #cluster4
    geom_line(data=data.frame(x=c(date[50],date[50]),y=c(data4[50]-0.1,data4[50]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#006400')+
    geom_line(data=data.frame(x=c(date[66],date[66]),y=c(data4[66]-1,data4[66]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#006400')+
    geom_line(data=data.frame(x=c(date[77],date[77]),y=c(data4[77]-1,data4[77]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#006400')+
    geom_line(data=data.frame(x=c(date[97],date[97]),y=c(data4[97]-1,data4[97]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#006400')
  fig34 <- fig34 +
    scale_x_date(breaks =as.Date(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),date_labels="%y/%m/%d")+
    # cluster3KR
    annotate(geom='text',x=as.Date(c("2020-02-07")),y=0.8,label="2020-01-23",size=7)+
    annotate(geom='text',x=as.Date(c("2020-02-19")),y=2,label="2020-02-03",size=7)+
    annotate(geom='text',x=as.Date(c("2020-03-07")),y=5.5,label="2020-02-22",size=7)+
    annotate(geom='text',x=as.Date(c("2020-03-18")),y=10,label="2020-03-04",size=7)+
    #cluster4JP
    annotate(geom='text',x=as.Date(c("2020-01-05")),y=1,label="2020-01-19",size=7)+
    annotate(geom='text',x=as.Date(c("2020-01-19")),y=2.9,label="2020-02-04",size=7)+
    annotate(geom='text',x=as.Date(c("2020-01-31")),y=5.7,label="2020-02-15",size=7)+
    annotate(geom='text',x=as.Date(c("2020-03-12")),y=7.7,label="2020-03-06",size=7)+
    theme_bw() +
    mytheme+
    theme(legend.position=c(0.15,0.85))+theme(legend.title = element_blank())+
    theme(legend.background = element_blank())+
    theme(legend.text = element_text(size=20))+
    theme(legend.key.size = unit(40, "pt"))+
    labs(tag = "B") +
    theme(plot.tag.position = c(0.05, 1))+
    theme(plot.tag=element_text(size = 30))

  # if(length(pred_x.new_3)!=0){
  #   for (i in 1:length(pred_x.new_3)) {
  #     cat(pred_x.new_3[[i]],'\n')
  #     pred_data=data.frame(x.new=pred_x.new_3[[i]],f.new=pred_f.new_3[[i]])
  #     fig34 <- fig34 +
  #       geom_line(data=pred_data, aes(x=x.new, y=f.new),shape = 2, size=1,colour='#F8766D')
  #   }
  # }
  # 
  # if(length(pred_x.new_4)!=0){
  #   for (i in 1:length(pred_x.new_4)) {
  #     cat(pred_x.new_4[[i]],'\n')
  #     pred_data=data.frame(x.new=pred_x.new_4[[i]],f.new=pred_f.new_4[[i]])
  #     fig34 <- fig34 +
  #       geom_line(data=pred_data, aes(x=x.new, y=f.new),shape = 2, size=1,colour='#B79F00')
  #   }
  # }
  fig34

  # cluster5data
  data=log(1+data_clust_C[1:dim(data_clust_C)[1],'cluster5'])
  data5<-data
  qq <- c(1,sort(final_C[[5]]),dim(data_clust_C)[1])
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  phase<-c('Cluster 5','Cluster 5','Cluster 5','Cluster 5')
  phase_num<-c(65,18,11,48) # c(65,18,11,48) c(65,18,12,47)
  fig56<-ggplot(dataplot,aes(Date))
  fig56 <- fig56 +
    geom_point(data = dataplot, aes(y=y)) +
    labs(title = expression(italic(paste('Cluster 5(US,DE,FR,GB,CA) and Cluster 6(ES)'))),
         x='', y=expression(paste("log(1+infection)")),fill="") +
    mytheme
  fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
  df.delta_data5 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                               lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())

  #-----------------------------------------------------------------------------------------------#
  # this section aim to predict the next part according to previous part
  pred_x.new_5=list()
  pred_f.new_5=list()
  if(length(final)!=0){
    for (item in 1:(length(qq)-2)) {
      sat=item
      a<-as.numeric(qq[sat])+1
      b=as.numeric(qq[sat+1])
      x.new <- seq((b+1), qq[sat+2], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1_num<-data[b]#t-1
      z2_num<-data[b-1]#t-2
      pred_x.new_5[[item]] <- dataplot[(b+1):qq[sat+2],2]
      cat('part',item+1,'a=',beta2.est[1],'b=',beta2.est[2],'beta1=',beta2.est[3],'beta2=',beta2.est[4],'alpha1=', beta2.est[5],'alpha2=',beta2.est[6],'\n')
      pred_f.new_5[[item]] <- predict_fun(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1_num,z2_num)
    }
    
  }
  #-----------------------------------------------------------------------------------------------#
  
  residual_list_5=list()
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
    }
    x.new <- seq((a+2), qq[i], by=1)
    nlm2 <- fun_fit(a,b,data[a:b])$fit
    beta2.est <- coef(nlm2)
    z1<-data[(a+1):(b-1)]#t-1
    z2<-data[(a):(b-2)]#t-2
    f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1,z2)
    residual_list_5[[i-1]]=c(data[a:(a+1)],f.new)-data[a:b]
    print(class(f.new))
    g.new <- attr(f.new,"gradient")
    # f.new <-as.numeric(c(f.new[1:2],f.new))
    # g.new<-rbind(g.new[1:2,],g.new)
    V.beta2 <- vcov(nlm2)
    GS=rowSums((g.new[,1:6]%*%V.beta2)%*%t(g.new[,1:6]))
    head(GS)
    alpha <- 0.05
    df <- length(data[a:b])-length(beta2.est)#freedom degree
    deltaf <- sqrt(GS)*qt(1-alpha/2,df)
    if(max(deltaf)>20){
      df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new), upr.conf=c(data[a:(a+1)],f.new))
    }else{
      df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new-deltaf), upr.conf=c(data[a:(a+1)],f.new+deltaf))
    }
    head(df.delta)
    sigma2.est <- summary(nlm2)$sigma
    deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
    line_dataframa<-data.frame(x.new=dataplot[a:b,2], f.new=as.numeric(c(data[a:(a+1)],f.new)))
    if(max(deltay)>20){
      df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new),c(data[a:(a+1)],f.new))
    }else{
      df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new - deltay),c(data[a:(a+1)],f.new + deltay))
    }
    df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
    df.delta_data5<-rbind(df.delta_data5,df.delta)
  }
  residual_list[[5]]=residual_list_5
  #cluster6 data
  data=log(1+data_clust_C[1:dim(data_clust_C)[1],'cluster6'])
  data6<-data
  qq <- c(1,sort(final_C[[6]]),dim(data_clust_C)[1])
  # plot fitting data
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  phase<-c('Cluster 6','Cluster 6')
  phase_num<-c(86,56)
  # fig2<-ggplot(dataplot,aes(Date))
  fig56 <- fig56 +
    geom_point(data = dataplot, aes(y=y)) +
    mytheme
  fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
  df.delta_data6 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                               lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())

  #-----------------------------------------------------------------------------------------------#
  # this section aim to predict the next part according to previous part
  pred_x.new_6=list()
  pred_f.new_6=list()
  if(length(final)!=0){
    for (item in 1:(length(qq)-2)) {
      sat=item
      a<-as.numeric(qq[sat])+1
      b=as.numeric(qq[sat+1])
      x.new <- seq((b+1), qq[sat+2], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1_num<-data[b]#t-1
      z2_num<-data[b-1]#t-2
      pred_x.new_6[[item]] <- dataplot[(b+1):qq[sat+2],2]
      cat('part',item+1,'a=',beta2.est[1],'b=',beta2.est[2],'beta1=',beta2.est[3],'beta2=',beta2.est[4],'alpha1=', beta2.est[5],'alpha2=',beta2.est[6],'\n')
      pred_f.new_6[[item]] <- predict_fun(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1_num,z2_num)
    }
    
  }
  #-----------------------------------------------------------------------------------------------#
  
  residual_list_6=list()
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
    }
    x.new <- seq((a+2), qq[i], by=1)
    nlm2 <- fun_fit(a,b,data[a:b])$fit
    beta2.est <- coef(nlm2)
    z1<-data[(a+1):(b-1)]#t-1
    z2<-data[(a):(b-2)]#t-2
    f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1,z2)
    residual_list_6[[i-1]]=c(data[a:(a+1)],f.new)-data[a:b]
    print(class(f.new))
    g.new <- attr(f.new,"gradient")
    # f.new <-as.numeric(c(f.new[1:2],f.new))
    # g.new<-rbind(g.new[1:2,],g.new)
    V.beta2 <- vcov(nlm2)
    GS=rowSums((g.new[,1:6]%*%V.beta2)%*%t(g.new[,1:6]))
    head(GS)
    alpha <- 0.05
    df <- length(data[a:b])-length(beta2.est)#freedom degree
    deltaf <- sqrt(GS)*qt(1-alpha/2,df)
    if(max(deltaf)>20){
      df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new), upr.conf=c(data[a:(a+1)],f.new))
    }else{
      df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new-deltaf), upr.conf=c(data[a:(a+1)],f.new+deltaf))
    }
    head(df.delta)
    sigma2.est <- summary(nlm2)$sigma
    deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
    line_dataframa<-data.frame(x.new=dataplot[a:b,2], f.new=as.numeric(c(data[a:(a+1)],f.new)))
    if(max(deltay)>20){
      df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new),c(data[a:(a+1)],f.new))
    }else{
      df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new - deltay),c(data[a:(a+1)],f.new + deltay))
    }
    df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
    df.delta_data6<-rbind(df.delta_data6,df.delta)
  }
  residual_list[[6]]=residual_list_6
  df.delta_data<-rbind(df.delta_data5,df.delta_data6)
  fig56 <- fig56 +
    # geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.pred, ymax=upr.pred, group=phase), alpha=0.2, fill='black') +
    geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.conf, ymax=upr.conf, group=phase,fill=phase), alpha=0.4) +
    geom_line(data=df.delta_data, aes(x=x.new, y=f.new,colour=phase, group=phase), size=1)+
    ylim(-0.1,14.5)+
    scale_fill_manual(values=c('#FF8C00','#006400'),name='Province')+
    scale_colour_manual(values=c('#FF8C00','#006400'),name='Province')+
    #cluster5
    geom_line(data=data.frame(x=c(date[65],date[65]),y=c(data5[65]-1,data5[65]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#FF8C00')+
    geom_line(data=data.frame(x=c(date[83],date[83]),y=c(data5[83]-1,data5[83]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#FF8C00')+
    geom_line(data=data.frame(x=c(date[94],date[94]),y=c(data5[94]-1,data5[94]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#FF8C00')+
    #cluster6
    geom_line(data=data.frame(x=c(date[86],date[86]),y=c(data6[86]-1,data6[86]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#006400')
  fig56 <- fig56 +
    scale_x_date(breaks =as.Date(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),date_labels="%y/%m/%d")+
    # cluster5
    annotate(geom='text',x=as.Date(c("2020-01-20")),y=3.5,label="2020-02-03",size=7)+
    annotate(geom='text',x=as.Date(c("2020-02-07")),y=4.8,label="2020-02-21",size=7)+
    annotate(geom='text',x=as.Date(c("2020-02-16")),y=6.8,label="2020-03-03",size=7)+
    #cluster6
    annotate(geom='text',x=as.Date(c("2020-03-10")),y=0.8,label="2020-02-24",size=7)+
    theme_bw() +
    mytheme+
    theme(legend.position=c(0.15,0.85))+theme(legend.title = element_blank())+
    theme(legend.background = element_blank())+
    theme(legend.text = element_text(size=20))+
    theme(legend.key.size = unit(40, "pt"))+
    labs(tag = "C") +
    theme(plot.tag.position = c(0.05, 1))+
    theme(plot.tag=element_text(size = 30))
  
  # if(length(pred_x.new_5)!=0){
  #   for (i in 1:length(pred_x.new_5)) {
  #     cat(pred_x.new_5[[i]],'\n')
  #     pred_data=data.frame(x.new=pred_x.new_5[[i]],f.new=pred_f.new_5[[i]])
  #     fig56 <- fig56 +
  #       geom_line(data=pred_data, aes(x=x.new, y=f.new),shape = 2, size=1,colour='#F8766D')
  #   }
  # }
  # 
  # if(length(pred_x.new_6)!=0){
  #   for (i in 1:length(pred_x.new_6)) {
  #     cat(pred_x.new_6[[i]],'\n')
  #     pred_data=data.frame(x.new=pred_x.new_6[[i]],f.new=pred_f.new_6[[i]])
  #     fig56 <- fig56 +
  #       geom_line(data=pred_data, aes(x=x.new, y=f.new),shape = 2, size=1,colour='#B79F00')
  #   }
  # }

  fig56

  # cluster7data
  data=log(1+data_clust_C[1:dim(data_clust_C)[1],'cluster7'])
  data7<-data
  qq <- c(1,sort(final_C[[7]]),dim(data_clust_C)[1])
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  phase<-c('Cluster 7','Cluster 7','Cluster 7')
  phase_num<-c(66,30,46)
  fig78<-ggplot(dataplot,aes(Date))
  fig78 <- fig78 +
    geom_point(data = dataplot, aes(y=y)) +
    labs(title = expression(italic(paste('Cluster 7(CH,NL,BE,AT,SE,DK,PT,BR,IE,RO,ID,PL,PE,\n EC,RU,IN,PH,EG,DZ,MK,DO) and Cluster 8(TR)'))),
         x='', y=expression(paste("log(1+infection)")),fill="") +
    mytheme
  fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
  df.delta_data7 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                               lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())

  #-----------------------------------------------------------------------------------------------#
  # this section aim to predict the next part according to previous part
  pred_x.new_7=list()
  pred_f.new_7=list()
  if(length(final)!=0){
    for (item in 1:(length(qq)-2)) {
      sat=item
      a<-as.numeric(qq[sat])+1
      b=as.numeric(qq[sat+1])
      x.new <- seq((b+1), qq[sat+2], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1_num<-data[b]#t-1
      z2_num<-data[b-1]#t-2
      pred_x.new_7[[item]] <- dataplot[(b+1):qq[sat+2],2]
      cat('part',item+1,'a=',beta2.est[1],'b=',beta2.est[2],'beta1=',beta2.est[3],'beta2=',beta2.est[4],'alpha1=', beta2.est[5],'alpha2=',beta2.est[6],'\n')
      pred_f.new_7[[item]] <- predict_fun(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1_num,z2_num)
    }
    
  }
  #-----------------------------------------------------------------------------------------------#
  
  residual_list_7=list()
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
      nlm <- fun_fit(a,b,data[a:b])$fit
      residual_list_7[[1]]=c(data[a:(a+1)],fitted(nlm))-data[a:b]
      data_one <- data.frame(x.new=dataplot[a:b,2],f.new=c(data[a:(a+1)],fitted(nlm)),lwr.conf=c(data[a:(a+1)],fitted(nlm))
                             ,upr.conf=c(data[a:(a+1)],fitted(nlm)),lwr.pred=c(data[a:(a+1)],fitted(nlm)),upr.pred=c(data[a:(a+1)],fitted(nlm)))
      data_one['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data7 <- rbind(df.delta_data7,data_one)
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
      x.new <- seq((a+2), qq[i], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1<-data[(a+1):(b-1)]#t-1
      z2<-data[(a):(b-2)]#t-2
      f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1,z2)
      residual_list_7[[i-1]]=c(data[a:(a+1)],f.new)-data[a:b]
      print(class(f.new))
      g.new <- attr(f.new,"gradient")
      # f.new <-as.numeric(c(f.new[1:2],f.new))
      # g.new<-rbind(g.new[1:2,],g.new)
      V.beta2 <- vcov(nlm2)
      GS=rowSums((g.new[,1:6]%*%V.beta2)%*%t(g.new[,1:6]))
      head(GS)
      alpha <- 0.05
      df <- length(data[a:b])-length(beta2.est)#freedom degree
      deltaf <- sqrt(GS)*qt(1-alpha/2,df)
      if(max(deltaf)>20){
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new), upr.conf=c(data[a:(a+1)],f.new))
      }else{
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new-deltaf), upr.conf=c(data[a:(a+1)],f.new+deltaf))
      }
      head(df.delta)
      sigma2.est <- summary(nlm2)$sigma
      deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
      line_dataframa<-data.frame(x.new=dataplot[a:b,2], f.new=as.numeric(c(data[a:(a+1)],f.new)))
      if(max(deltay)>20){
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new),c(data[a:(a+1)],f.new))
      }else{
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new - deltay),c(data[a:(a+1)],f.new + deltay))
      }
      df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data7<-rbind(df.delta_data7,df.delta)
    }
  }
  pred_plot_sim[[2]]=df.delta_data7
  residual_list[[7]]=residual_list_7
  #cluster8 data
  data=log(1+data_clust_C[1:dim(data_clust_C)[1],'cluster8'])
  data8<-data
  qq <- c(1,sort(final_C[[8]]),dim(data_clust_C)[1])
  # plot fitting data
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  phase<-c('Cluster 8')
  phase_num<-c(142)
  # fig2<-ggplot(dataplot,aes(Date))
  fig78 <- fig78 +
    geom_point(data = dataplot, aes(y=y)) +
    mytheme
  fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
  df.delta_data8 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                               lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())

  #-----------------------------------------------------------------------------------------------#
  # this section aim to predict the next part according to previous part
  # pred_x.new_8=list()
  # pred_f.new_8=list()
  # if(length(final)!=0){
  #   for (item in 1:(length(qq)-2)) {
  #     sat=item
  #     a<-as.numeric(qq[sat])+1
  #     b=as.numeric(qq[sat+1])
  #     x.new <- seq((b+1), qq[sat+2], by=1)
  #     nlm2 <- fun_fit(a,b,data[a:b])$fit
  #     beta2.est <- coef(nlm2)
  #     z1_num<-data[b]#t-1
  #     z2_num<-data[b-1]#t-2
  #     pred_x.new_8[[item]] <- dataplot[(b+1):qq[sat+2],2]
  #     cat('part',item+1,'a=',beta2.est[1],'b=',beta2.est[2],'beta1=',beta2.est[3],'beta2=',beta2.est[4],'alpha1=', beta2.est[5],'alpha2=',beta2.est[6],'\n')
  #     pred_f.new_8[[item]] <- predict_fun(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1_num,z2_num)
  #   }
  #   
  # }
  #-----------------------------------------------------------------------------------------------#
   
  residual_list_8=list()
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
      nlm <- fun_fit(a,b,data[a:b])$fit
      residual_list_8[[1]]=c(data[a:(a+1)],fitted(nlm))-data[a:b]
      data_one <- data.frame(x.new=dataplot[a:b,2],f.new=c(data[a:(a+1)],fitted(nlm)),lwr.conf=c(data[a:(a+1)],fitted(nlm))
                             ,upr.conf=c(data[a:(a+1)],fitted(nlm)),lwr.pred=c(data[a:(a+1)],fitted(nlm)),upr.pred=c(data[a:(a+1)],fitted(nlm)))
      data_one['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data8 <- rbind(df.delta_data8,data_one)
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
      x.new <- seq((a+2), qq[i], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1<-data[(a+1):(b-1)]#t-1
      z2<-data[(a):(b-2)]#t-2
      f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1,z2)
      residual_list_8[[i-1]]=c(data[a:(a+1)],f.new)-data[a:b]
      print(class(f.new))
      g.new <- attr(f.new,"gradient")
      # f.new <-as.numeric(c(f.new[1:2],f.new))
      # g.new<-rbind(g.new[1:2,],g.new)
      V.beta2 <- vcov(nlm2)
      GS=rowSums((g.new[,1:6]%*%V.beta2)%*%t(g.new[,1:6]))
      head(GS)
      alpha <- 0.05
      df <- length(data[a:b])-length(beta2.est)#freedom degree
      deltaf <- sqrt(GS)*qt(1-alpha/2,df)
      if(max(deltaf)>20){
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new), upr.conf=c(data[a:(a+1)],f.new))
      }else{
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new-deltaf), upr.conf=c(data[a:(a+1)],f.new+deltaf))
      }
      head(df.delta)
      sigma2.est <- summary(nlm2)$sigma
      deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
      line_dataframa<-data.frame(x.new=dataplot[a:b,2], f.new=as.numeric(c(data[a:(a+1)],f.new)))
      if(max(deltay)>20){
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new),c(data[a:(a+1)],f.new))
      }else{
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new - deltay),c(data[a:(a+1)],f.new + deltay))
      }
      df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data8<-rbind(df.delta_data8,df.delta)
    }
  }
  residual_list[[8]]=residual_list_8
  df.delta_data<-rbind(df.delta_data7,df.delta_data8)
  fig78 <- fig78 +
    # geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.pred, ymax=upr.pred, group=phase), alpha=0.2, fill='black') +
    geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.conf, ymax=upr.conf, group=phase,fill=phase), alpha=0.4) +
    geom_line(data=df.delta_data, aes(x=x.new, y=f.new,colour=phase, group=phase), size=1)+
    ylim(-0.1,14.5)+
    scale_fill_manual(values=c('#FF8C00','#006400'),name='Province')+
    scale_colour_manual(values=c('#FF8C00','#006400'),name='Province')+
    #cluster7
    geom_line(data=data.frame(x=c(date[66],date[66]),y=c(data7[66]-0.7,data7[66]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#FF8C00')+
    geom_line(data=data.frame(x=c(date[96],date[96]),y=c(data7[96]-1,data7[96]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='#FF8C00')

  fig78 <- fig78 +
    scale_x_date(breaks =as.Date(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),date_labels="%y/%m/%d")+
    # cluster8
    annotate(geom='text',x=as.Date(c("2020-01-20")),y=1.6,label="2020-02-04",size=7)+
    annotate(geom='text',x=as.Date(c("2020-02-18")),y=5.7,label="2020-03-05",size=7)+
    theme_bw() +
    mytheme+
    theme(legend.position=c(0.15,0.85))+theme(legend.title = element_blank())+
    theme(legend.background = element_blank())+
    theme(legend.text = element_text(size=20))+
    theme(legend.key.size = unit(40, "pt"))+
    labs(tag = "D") +
    theme(plot.tag.position = c(0.05, 1))+
    theme(plot.tag=element_text(size = 30))
  
  # if(length(pred_x.new_7)!=0){
  #   pred_plot_test=list()
  #   for (i in 1:length(pred_x.new_7)) {
  #     cat(pred_x.new_7[[i]],'\n')
  #     pred_data=data.frame(x.new=pred_x.new_7[[i]],f.new=pred_f.new_7[[i]])
  #     pred_plot_test[[i]]=pred_data
  #     fig78 <- fig78 +
  #       geom_line(data=pred_data, aes(x=x.new, y=f.new),shape = 2, size=1,colour='#F8766D')
  #   }
  #   pred_plot[[2]]=pred_plot_test
  # }

  # if(length(pred_x.new_8)!=0){
  #   for (i in 1:length(pred_x.new_8)) {
  #     cat(pred_x.new_8[[i]],'\n')
  #     pred_data=data.frame(x.new=pred_x.new_8[[i]],f.new=pred_f.new_8[[i]])
  #     fig78 <- fig78 +
  #       geom_line(data=pred_data, aes(x=x.new, y=f.new),shape = 2, size=1,colour='#B79F00')
  #   }
  # }

  fig78
  tu <- grid.arrange(fig12,fig34,fig56,fig78,nrow=1,ncol=4)
  # ggsave(filename='country_1.pdf', plot = tu, width = 40, height = 8, dpi = 300,limitsize = FALSE)
  return(list(pred_plot=pred_plot,residual_list=residual_list,pred_plot_sim=pred_plot_sim))
}

