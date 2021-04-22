library(readr)
library(ggplot2)

library(dplyr)
library(grid)
library(gridExtra)
library(tidyr)

setwd("")

# initial parameter satting
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("test_new_a1.cpp")

load('cluster_data.RData')
load('origin_data.RData')
load('final_new_formula2.RData')

residual_list=LogPLotChina(data_clust, final_new_formula2)

LogPLotChina <- function(data_clust, final){
  residual_list=list()
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
    para=GridSearch2(start:end,data,-20,20,-20,20,100)
    y <- data[1:(end-start+1)]
    x <- seq(start,end,1)
    fit<-try(nlm1<- nls(y ~ b*pnorm(beta1+beta2*x)+a, start=list(a=para$a ,b=para$b,beta1=para$beta1 ,beta2=para$beta2),
                        control=nls.control(maxiter = 2000, tol = 1e-09, minFactor = 1/20240000, printEval = FALSE, warnOnly = T))
    )
    return(list(fit=fit,residual=para$residual,para=para))
  }


  Date = rownames(data_clust)
  # interval_estimated
  #cluster 1(HB)
  data=log(1+data_clust[1:dim(data_clust)[1],'cluster1'])
  qq <- c(1,sort(final[[1]]),dim(data_clust)[1])
  # plot fitting data
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  fig1<-ggplot(dataplot,aes(Date))
  phase<-c(1,2,3,4,5,6)
  phase_num<-c(final[[1]][1])
  if(length(final[[1]])>1){
    for(i in 2:length(final[[1]])){
      phase_num <- c(phase_num, final[[1]][i]-final[[1]][i-1])
    }
  }
  phase_num <- c(phase_num, dim(data_clust)[1]-final[[1]][length(final[[1]])])

  fig1 <- fig1 +
    geom_point(data = dataplot, aes(y=y)) +
    labs(title = expression(italic(paste('Cluster 1(HB)'))),x='', y=expression(paste("log(1+infection)")),fill="") +
    mytheme
  fgh2 <- deriv(y ~ b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2"), function(a, b, beta1, beta2,x){} ) 
  df.delta_data <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                              lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
  
  
  residual_list_1 = list()
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
      nlm <- fun_fit(a,b,data[a:b])$fit
      residual_list_1[[1]]=c(fitted(nlm))-data[a:b]
      data_one <- data.frame(x.new=dataplot[a:b,2],f.new=c(fitted(nlm)),lwr.conf=c(fitted(nlm))
                             ,upr.conf=c(fitted(nlm)),lwr.pred=c(fitted(nlm)),upr.pred=c(fitted(nlm)))
      cat('phase_num now is',phase_num[i-1],'\n')
      data_one['phase']<-rep(phase[i-1],phase_num[i-1])
      
      df.delta_data <- rbind(df.delta_data,data_one)
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
      x.new <- seq((a), qq[i], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4],x.new)
      residual_list_1[[i-1]]=f.new-data[a:b]
      print(class(f.new))
      g.new <- attr(f.new,"gradient")
      V.beta2 <- vcov(nlm2)
      GS=rowSums((g.new[,1:4]%*%V.beta2)%*%t(g.new[,1:4]))
      head(GS)
      alpha <- 0.05
      df <- length(data[a:b])-length(beta2.est)#freedom degree
      deltaf <- sqrt(GS)*qt(1-alpha/2,df)
      if(max(deltaf)>20){
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(f.new), lwr.conf=c(f.new), upr.conf=c(f.new))
      }else{
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(f.new), lwr.conf=c(f.new-deltaf), upr.conf=c(f.new+deltaf))
      }
      head(df.delta)
      sigma2.est <- summary(nlm2)$sigma
      deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
      line_dataframa<-data.frame(x.new=dataplot[a:b,2], f.new=as.numeric(c(f.new)))
      if(max(deltay)>20){
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(f.new),c(f.new))
      }else{
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(f.new - deltay),c(f.new + deltay))
      }
      cat('phase_num now is',phase_num[i-1],'\n')
      df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data<-rbind(df.delta_data,df.delta)
    }
  }
  residual_list[[1]]=residual_list_1
  color <- c('red','blue','green','pink','yellow','red','blue','pink','green','yellow')
  xintercept<-c(date[14],date[47],date[57],date[72],date[90])
  fig1 <- fig1 +
    geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.conf, ymax=upr.conf), alpha=0.4,fill=color[1]) +#, group=phase,fill=phase
    geom_line(data=df.delta_data, aes(x=x.new, y=f.new), size=1,colour=color[1])+#, colour=phase
    ylim(0,12.5)+
    geom_line(data=data.frame(x=c(date[14],date[14]),y=c(data[14]-0.7,data[14]+1)),aes(x=x,y=y),linetype='longdash',size=1,color=color[1])+
    geom_line(data=data.frame(x=c(date[47],date[47]),y=c(data[47]-1,data[47]+1)),aes(x=x,y=y),linetype='longdash',size=1,color=color[1])
  fig1 <- fig1 +
    scale_x_date(breaks =as.Date(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),date_labels="%y/%m/%d")+
    annotate(geom='text',x=as.Date(c("2019-12-29")),y=1.5,label="2019-12-14",size=7)+
    annotate(geom='text',x=as.Date(c("2020-01-31")),y=3.8,label="2020-01-16",size=7)+
    theme_bw() +
    mytheme+
    labs(tag = "A") +
    theme(plot.tag.position = c(0.05, 1))+
    theme(plot.tag=element_text(size = 30))
  fig1


  #cluster 2(....)
  data=log(1+data_clust[1:dim(data_clust)[1],'cluster2'])
  qq <- c(1,sort(final[[2]]),dim(data_clust)[1])
  # plot fitting data
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  phase<-c(11,12,13,14)
  phase_num<-c(final[[2]][1])
  if(length(final[[2]])>1){
    for(i in 2:length(final[[2]])){
      phase_num <- c(phase_num, final[[2]][i]-final[[2]][i-1])
    }
  }
  phase_num <- c(phase_num, dim(data_clust)[1]-final[[2]][length(final[[2]])])
  cat('phase_num now is',phase_num,'\n')
  fig2<-ggplot(dataplot,aes(Date))
  fig2 <- fig2 +
    geom_point(data = dataplot, aes(y=y)) +
    labs(title = expression(italic(paste('Cluster 2(GD,ZJ,HA,HN,AH,JX,SD,JS,SC,CQ,BJ,SH,HL,\n HE,FJ,SN,GX,YN,HI,LN,SX,TJ,GZ,GS,NM,JL,XJ,NX)'))),
         x='', y=expression(paste("log(1+infection)")),fill="") +
    mytheme
  fgh2 <- deriv(y ~ b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2"), function(a, b, beta1, beta2,x){} ) 
  df.delta_data <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                              lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
  
  residual_list_2 = list()
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
      nlm <- fun_fit(a,b,data[a:b])$fit
      residual_list_2[[1]]=c(fitted(nlm))-data[a:b]
      data_one <- data.frame(x.new=dataplot[a:b,2],f.new=c(fitted(nlm)),lwr.conf=c(fitted(nlm))
                             ,upr.conf=c(fitted(nlm)),lwr.pred=c(fitted(nlm)),upr.pred=c(fitted(nlm)))
      cat('phase_num now is',phase_num[i-1],'\n')
      data_one['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data <- rbind(df.delta_data,data_one)
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
      x.new <- seq((a), qq[i], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4],x.new)
      residual_list_2[[i-1]]=f.new-data[a:b]
      print(class(f.new))
      g.new <- attr(f.new,"gradient")
      V.beta2 <- vcov(nlm2)
      GS=rowSums((g.new[,1:4]%*%V.beta2)%*%t(g.new[,1:4]))
      head(GS)
      alpha <- 0.05
      df <- length(data[a:b])-length(beta2.est)#freedom degree
      deltaf <- sqrt(GS)*qt(1-alpha/2,df)
      if(max(deltaf)>20){
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(f.new), lwr.conf=c(f.new), upr.conf=c(f.new))
      }else{
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(f.new), lwr.conf=c(f.new-deltaf), upr.conf=c(f.new+deltaf))
      }
      head(df.delta)
      sigma2.est <- summary(nlm2)$sigma
      deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
      line_dataframa<-data.frame(x.new=dataplot[a:b,2], f.new=as.numeric(c(f.new)))
      if(max(deltay)>20){
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(f.new),c(f.new))
      }else{
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(f.new - deltay),c(f.new + deltay))
      }
      cat('phase_num now is',phase_num[i-1],'\n')
      df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data<-rbind(df.delta_data,df.delta)
    }
  }
  residual_list[[2]]=residual_list_2
  fig2 <- fig2 +
    geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.conf, ymax=upr.conf), alpha=0.4, fill=color[2]) +
    geom_line(data=df.delta_data, aes(x=x.new, y=f.new), size=1, colour=color[2])+
    ylim(-0.1,12.5)+
    geom_line(data=data.frame(x=c(date[53],date[53]),y=c(data[53]-1,data[53]+1)),aes(x=x,y=y),linetype='longdash',size=1,color=color[2])+
    geom_line(data=data.frame(x=c(date[83],date[83]),y=c(data[83]-1,data[83]+1)),aes(x=x,y=y),linetype='longdash',size=1,color=color[2])

  fig2 <- fig2 +
    scale_x_date(breaks =as.Date(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),date_labels="%y/%m/%d")+
    annotate(geom='text',x=as.Date(c("2020-02-06")),y=5.3,label="2020-1-22",size=7)+
    annotate(geom='text',x=as.Date(c("2020-03-07")),y=10,label="2020-02-21",size=7)+
    theme_bw() +
    mytheme+
    labs(tag = "B") +
    theme(plot.tag.position = c(0.05, 1))+
    theme(plot.tag=element_text(size = 30))
  fig2

  # Cluster 3(TW) AND Cluster 4(HK)
  # cluster3data
  data=log(1+data_clust[1:dim(data_clust)[1],'cluster3'])
  data3<-data
  qq <- c(1,sort(final[[3]]),dim(data_clust)[1])
  # plot fitting data
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  phase<-c('Cluster 3','Cluster 3','Cluster 3','Cluster 3')
  phase_num<-c(final[[3]][1])
  if(length(final[[3]])>1){
    for(i in 2:length(final[[3]])){
      phase_num <- c(phase_num, final[[3]][i]-final[[3]][i-1])
    }
  }
  phase_num <- c(phase_num, dim(data_clust)[1]-final[[3]][length(final[[3]])])

  fig34<-ggplot(dataplot,aes(Date))
  fig34 <- fig34 +
    geom_point(data = dataplot, aes(y=y)) +
    labs(title = expression(italic(paste('Cluster 3(TW) and Cluster 4(HK)'))),
         x='', y=expression(paste("log(1+infection)")),fill="") +
    mytheme
  fgh2 <- deriv(y ~ b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2"), function(a, b, beta1, beta2,x){} ) 
  df.delta_data1 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                              lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
  
  residual_list_3=list()
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
      nlm <- fun_fit(a,b,data[a:b])$fit
      residual_list_3[[1]]=c(fitted(nlm))-data[a:b]
      data_one <- data.frame(x.new=dataplot[a:b,2],f.new=c(fitted(nlm)),lwr.conf=c(fitted(nlm))
                             ,upr.conf=c(fitted(nlm)),lwr.pred=c(fitted(nlm)),upr.pred=c(fitted(nlm)))
      data_one['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data1 <- rbind(df.delta_data1,data_one)
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
      x.new <- seq((a), qq[i], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4],x.new)
      residual_list_3[[i-1]]=f.new - data[a:b]
      print(class(f.new))
      g.new <- attr(f.new,"gradient")
      V.beta2 <- vcov(nlm2)
      GS=rowSums((g.new[,1:4]%*%V.beta2)%*%t(g.new[,1:4]))
      head(GS)
      alpha <- 0.05
      df <- length(data[a:b])-length(beta2.est)#freedom degree
      deltaf <- sqrt(GS)*qt(1-alpha/2,df)
      if(max(deltaf)>20){
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(f.new), lwr.conf=c(f.new), upr.conf=c(f.new))
      }else{
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(f.new), lwr.conf=c(f.new-deltaf), upr.conf=c(f.new+deltaf))
      }
      head(df.delta)
      sigma2.est <- summary(nlm2)$sigma
      deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
      line_dataframa<-data.frame(x.new=dataplot[a:b,2], f.new=as.numeric(c(f.new)))
      if(max(deltay)>20){
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(f.new),c(f.new))
      }else{
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(f.new - deltay),c(f.new + deltay))
      }
      df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data1<-rbind(df.delta_data1,df.delta)
    }
  }
  residual_list[[3]]=residual_list_3
  #cluster4 data
  data=log(1+data_clust[1:dim(data_clust)[1],'cluster4'])
  data4<-data
  qq <- c(1,sort(final[[4]]),dim(data_clust)[1])
  # plot fitting data
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  phase<-c('Cluster 4','Cluster 4','Cluster 4','Cluster 4')
  phase_num<-c(final[[4]][1])
  if(length(final[[4]])>1){
    for(i in 2:length(final[[4]])){
      phase_num <- c(phase_num, final[[4]][i]-final[[4]][i-1])
    }
  }
  phase_num <- c(phase_num, dim(data_clust)[1]-final[[4]][length(final[[4]])])

  fig34 <- fig34 +
    geom_point(data = dataplot, aes(y=y)) +
    mytheme
  fgh2 <- deriv(y ~ b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2"), function(a, b, beta1, beta2,x){} ) 
  df.delta_data2 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                               lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
  
  residual_list_4=list()
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
    }
    x.new <- seq((a), qq[i], by=1)
    nlm2 <- fun_fit(a,b,data[a:b])$fit
    beta2.est <- coef(nlm2)
    f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4],x.new)
    residual_list_4[[i-1]]=f.new-data[a:b]
    print(class(f.new))
    g.new <- attr(f.new,"gradient")
    V.beta2 <- vcov(nlm2)
    GS=rowSums((g.new[,1:4]%*%V.beta2)%*%t(g.new[,1:4]))
    head(GS)
    alpha <- 0.05
    df <- length(data[a:b])-length(beta2.est)#freedom degree
    deltaf <- sqrt(GS)*qt(1-alpha/2,df)
    if(max(deltaf)>20){
      df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(f.new), lwr.conf=c(f.new), upr.conf=c(f.new))
    }else{
      df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(f.new), lwr.conf=c(f.new-deltaf), upr.conf=c(f.new+deltaf))
    }
    head(df.delta)
    sigma2.est <- summary(nlm2)$sigma
    deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
    line_dataframa<-data.frame(x.new=dataplot[a:b,2], f.new=as.numeric(c(f.new)))
    if(max(deltay)>20){
      df.delta[c("lwr.pred","upr.pred")] <- cbind(c(f.new),c(f.new))
    }else{
      df.delta[c("lwr.pred","upr.pred")] <- cbind(c(f.new - deltay),c(f.new + deltay))
    }
    df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
    df.delta_data2<-rbind(df.delta_data2,df.delta)
  }
  residual_list[[4]]=residual_list_4
  df.delta_data<-rbind(df.delta_data1,df.delta_data2)
  fig34 <- fig34 +
    geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.conf, ymax=upr.conf, group=phase,fill=phase), alpha=0.4) +
    geom_line(data=df.delta_data, aes(x=x.new, y=f.new,colour=phase, group=phase), size=1)+
    ylim(-0.1,12.5)+
    scale_fill_manual(values=c('red','blue'),name='Province')+
    scale_colour_manual(values=c('red','blue'),name='Province')+
    #cluster3
    geom_line(data=data.frame(x=c(date[58],date[58]),y=c(data3[58]-1,data3[58]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='red')+
    geom_line(data=data.frame(x=c(date[68],date[68]),y=c(data3[68]-1,data3[68]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='red')+
    geom_line(data=data.frame(x=c(date[98],date[98]),y=c(data3[98]-1,data3[98]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='red')+
    #cluster4
    geom_line(data=data.frame(x=c(date[54],date[54]),y=c(data4[54]-1,data4[54]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='blue')+
    geom_line(data=data.frame(x=c(date[98],date[98]),y=c(data4[98]-1,data4[98]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='blue')
  fig34 <- fig34 +
    scale_x_date(breaks =as.Date(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),date_labels="%y/%m/%d")+
    # cluster3
    annotate(geom='text',x=as.Date(c("2020-02-12")),y=1.5,label="2020-01-27",size=7)+
    annotate(geom='text',x=as.Date(c("2020-02-20")),y=2.3,label="2020-02-06",size=7)+
    annotate(geom='text',x=as.Date(c("2020-03-22")),y=3.2,label="2020-03-07",size=7)+
    #cluster4
    annotate(geom='text',x=as.Date(c("2020-01-08")),y=2,label="2020-01-23",size=7)+
    annotate(geom='text',x=as.Date(c("2020-02-21")),y=5.3,label="2020-03-07",size=7)+
    theme_bw() +
    mytheme+
    theme(legend.position=c(0.15,0.85))+theme(legend.title = element_blank())+
    theme(legend.background = element_blank())+
    theme(legend.text = element_text(size=20))+
    theme(legend.key.size = unit(40, "pt"))+
    labs(tag = "C") +
    theme(plot.tag.position = c(0.05, 1))+
    theme(plot.tag=element_text(size = 30))
  fig34


  #cluster 5(MO),cluster 6(QH),cluster 7(XZ)
  # cluster5
  data=log(1+data_clust[1:dim(data_clust)[1],'cluster5'])
  data5<-data
  qq <- c(1,sort(final[[5]]),dim(data_clust)[1])
  # plot fitting data
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  phase<-c('Cluster 5','Cluster 5','Cluster 5')
  phase_num<-c(final[[5]][1])
  if(length(final[[5]])>1){
    for(i in 2:length(final[[5]])){
      phase_num <- c(phase_num, final[[5]][i]-final[[5]][i-1])
    }
  }
  phase_num <- c(phase_num, dim(data_clust)[1]-final[[5]][length(final[[5]])])

  fig567<-ggplot(dataplot,aes(Date))
  fig567 <- fig567 +
    geom_point(data = dataplot, aes(y=y)) +
    labs(title = expression(italic(paste('Cluster 5(MO), Cluster 6(QH) and Cluster 7(XZ)'))),
         x='', y=expression(paste("log(1+infection)")),fill="") +
    mytheme
  fgh2 <- deriv(y ~ b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2"), function(a, b, beta1, beta2,x){} ) 
  df.delta_data5 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                              lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
  
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
    x.new <- seq((a), qq[i], by=1)
    nlm2 <- fun_fit(a,b,data[a:b])$fit
    beta2.est <- coef(nlm2)
    f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4],x.new)
    residual_list_5[[i-1]]=f.new-data[a:b]
    print(class(f.new))
    g.new <- attr(f.new,"gradient")
    V.beta2 <- vcov(nlm2)
    GS=rowSums((g.new[,1:4]%*%V.beta2)%*%t(g.new[,1:4]))
    head(GS)
    alpha <- 0.05
    df <- length(data[a:b])-length(beta2.est)#freedom degree
    deltaf <- sqrt(GS)*qt(1-alpha/2,df)
    if(max(deltaf)>20){
      df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(f.new), lwr.conf=c(f.new), upr.conf=c(f.new))
    }else{
      df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(f.new), lwr.conf=c(f.new-deltaf), upr.conf=c(f.new+deltaf))
    }
    head(df.delta)
    sigma2.est <- summary(nlm2)$sigma
    deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
    line_dataframa<-data.frame(x.new=dataplot[a:b,2], f.new=as.numeric(c(f.new)))
    if(max(deltay)>20){
      df.delta[c("lwr.pred","upr.pred")] <- cbind(c(f.new),c(f.new))
    }else{
      df.delta[c("lwr.pred","upr.pred")] <- cbind(c(f.new - deltay),c(f.new + deltay))
    }
    df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
    df.delta_data5<-rbind(df.delta_data5,df.delta)
  }
  residual_list[[5]]=residual_list_5
  # cluster6
  data=log(1+data_clust[1:dim(data_clust)[1],'cluster6'])
  data6<-data
  qq <- c(1,sort(final[[6]]),dim(data_clust)[1])
  # plot fitting data
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  phase<-c('Cluster 6','Cluster 6')
  phase_num<-c(final[[6]][1])
  if(length(final[[6]])>1){
    for(i in 2:length(final[[6]])){
      phase_num <- c(phase_num, final[[6]][i]-final[[6]][i-1])
    }
  }
  phase_num <- c(phase_num, dim(data_clust)[1]-final[[6]][length(final[[6]])])

  fig567 <- fig567 +
    geom_point(data = dataplot, aes(y=y)) +
    mytheme
  fgh2 <- deriv(y ~ b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2"), function(a, b, beta1, beta2,x){} ) 
  df.delta_data6 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                               lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
  residual_list_6=list
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
    }
    x.new <- seq((a), qq[i], by=1)
    nlm2 <- fun_fit(a,b,data[a:b])$fit
    beta2.est <- coef(nlm2)
    f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4],x.new)
    residual_list_6=f.new -data[a:b]
    print(class(f.new))
    g.new <- attr(f.new,"gradient")
    V.beta2 <- vcov(nlm2)
    GS=rowSums((g.new[,1:4]%*%V.beta2)%*%t(g.new[,1:4]))
    head(GS)
    alpha <- 0.05
    df <- length(data[a:b])-length(beta2.est)#freedom degree
    deltaf <- sqrt(GS)*qt(1-alpha/2,df)
    if(max(deltaf)>20){
      df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(f.new), lwr.conf=c(f.new), upr.conf=c(f.new))
    }else{
      df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(f.new), lwr.conf=c(f.new-deltaf), upr.conf=c(f.new+deltaf))
    }
    head(df.delta)
    sigma2.est <- summary(nlm2)$sigma
    deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
    line_dataframa<-data.frame(x.new=dataplot[a:b,2], f.new=as.numeric(c(f.new)))
    if(max(deltay)>20){
      df.delta[c("lwr.pred","upr.pred")] <- cbind(c(f.new),c(f.new))
    }else{
      df.delta[c("lwr.pred","upr.pred")] <- cbind(c(f.new - deltay),c(f.new + deltay))
    }
    df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
    df.delta_data6<-rbind(df.delta_data6,df.delta)
  }
  residual_list[[6]]=residual_list_6
  #cluster7
  data=log(1+data_clust[1:dim(data_clust)[1],'cluster7'])
  data7<-data
  qq <- c(1,dim(data_clust)[1])
  # plot fitting data
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  phase<-c('Cluster 7')
  phase_num<-c(142)
  # phase_num<-c(final[[4]][1])
  # if(length(final[[4]])>1){
  #   for(i in 2:length(final[[4]])){
  #     phase_num <- c(phase_num, final[[4]][i]-final[[4]][i-1])
  #   }
  # }
  # phase_num <- c(phase_num, dim(data_clust)[4]-final[[4]][length(final[[4]])])

  fig567 <- fig567 +
    geom_point(data = dataplot, aes(y=y)) +
    mytheme
  fgh2 <- deriv(y ~ b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2"), function(a, b, beta1, beta2,x){} ) 
  df.delta_data7 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                               lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())

  residual_list_7=list()
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
      nlm <- fun_fit(a,b,data[a:b])$fit
      residual_list_7[[1]]=c(fitted(nlm))-data[a:b]
      data_one <- data.frame(x.new=dataplot[a:b,2],f.new=c(fitted(nlm)),lwr.conf=c(fitted(nlm))
                             ,upr.conf=c(fitted(nlm)),lwr.pred=c(fitted(nlm)),upr.pred=c(fitted(nlm)))
      data_one['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data7 <- rbind(df.delta_data7,data_one)
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
      x.new <- seq((a), qq[i], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4],x.new)
      residual_list_7[[1]]=f.new-data[a:b]
      print(class(f.new))
      g.new <- attr(f.new,"gradient")
      V.beta2 <- vcov(nlm2)
      GS=rowSums((g.new[,1:4]%*%V.beta2)%*%t(g.new[,1:4]))
      head(GS)
      alpha <- 0.05
      df <- length(data[a:b])-length(beta2.est)#freedom degree
      deltaf <- sqrt(GS)*qt(1-alpha/2,df)
      if(max(deltaf)>20){
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(f.new), lwr.conf=c(f.new), upr.conf=c(f.new))
      }else{
        df.delta <- data.frame(x.new=dataplot[a:b,2], f.new=c(f.new), lwr.conf=c(f.new-deltaf), upr.conf=c(f.new+deltaf))
      }
      head(df.delta)
      sigma2.est <- summary(nlm2)$sigma
      deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
      line_dataframa<-data.frame(x.new=dataplot[a:b,2], f.new=as.numeric(c(f.new)))
      if(max(deltay)>20){
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(f.new),c(f.new))
      }else{
        df.delta[c("lwr.pred","upr.pred")] <- cbind(c(f.new - deltay),c(f.new + deltay))
      }
      df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data7<-rbind(df.delta_data7,df.delta)
    }
  }
  residual_list[[7]]=residual_list_7
  df.delta_data<-rbind(df.delta_data5,df.delta_data6,df.delta_data7)
  fig567 <- fig567 +
    geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.conf, ymax=upr.conf, group=phase,fill=phase), alpha=0.4) +
    geom_line(data=df.delta_data, aes(x=x.new, y=f.new,colour=phase, group=phase), size=1)+
    ylim(-0.1,12.5)+
    scale_fill_manual(values=c('red','blue','#FF8C00'),name='Province')+
    scale_colour_manual(values=c('red','blue','#FF8C00'),name='Province')+
    #cluster5
    geom_line(data=data.frame(x=c(date[65],date[65]),y=c(data5[65]-1,data5[65]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='red')+
    #cluster6
    geom_line(data=data.frame(x=c(date[60],date[60]),y=c(data6[60]-1,data6[60]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='blue')
  fig567 <- fig567 +
    scale_x_date(breaks =as.Date(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),date_labels="%y/%m/%d")+
    #cluster5
    annotate(geom='text',x=as.Date(c("2020-02-19")),y=2,label="2020-02-03",size=7)+
    # cluster6
    annotate(geom='text',x=as.Date(c("2020-01-14")),y=2.5,label="2020-01-29",size=7)+
    theme_bw() +
    mytheme+
    theme(legend.position=c(0.15,0.8))+theme(legend.title = element_blank())+
    theme(legend.background = element_blank())+
    theme(legend.text = element_text(size=20))+
    theme(legend.key.size = unit(40, "pt"))+
    labs(tag = "D") +
    theme(plot.tag.position = c(0.05, 1))+
    theme(plot.tag=element_text(size = 30))
  fig567

  empty <- ggplot() + geom_point(aes(1, 1), colour = "white") +   
    theme(axis.ticks = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_blank(), 
          axis.text.x = element_blank(), axis.text.y = element_blank(), 
          axis.title.x = element_blank(), axis.title.y = element_blank())
  tu <- grid.arrange(fig1,fig2,fig34,fig567,nrow=1,ncol=4)
  ggsave(filename='E:/Desktop/R_test/covid_new_demand/test_china2.pdf', plot = tu, width = 40, height = 8, dpi = 300,limitsize = FALSE)
  # ggsave(filename='test_china.png', plot = tu, width = 40, height = 8, dpi = 300,limitsize = FALSE)
  return(residual_list)
}


