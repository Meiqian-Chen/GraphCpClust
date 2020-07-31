#' Basic description
#'
#' @description This function is used to fit the log-transformed infection counts of provinces/regions in China that introduced in Shi, Chen, Dong and Rao (2020)
#' @usage LogPLotChina(data_clust, final)
#' @param data_clust data_clust
#' @param final final
#' @examples
#' 
#' data("cluster_data",package = "GraphCpClust")
#' 
#' LogPLotChina(data_clust, final)
#' @seealso 
#' @export


LogPLotChina <- function(data_clust, final){
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
    # para=ini.par(beta00=-10,beta01=10,beta10=-10,beta11=10,x=c(start:end),y=data)
    cat(c(as.numeric(para[1]),as.numeric(para[2]),as.numeric(para[3]),as.numeric(para[4]),as.numeric(para[5]),as.numeric(para[6])),'\n')
    y <- data[3:(end-start+1)]
    z1<-data[2:(end-start+1-1)]#t-1
    z2<-data[1:(end-start+1-2)]#t-2
    x <- seq(start+2,end,1)
    fit<-try(nlm1<- nls(y ~ alpha1*z1+alpha2*z2 + b*pnorm(beta1+beta2*x)+a, start=list(a=para$a ,b=para$b
                                                                                       ,beta1=para$beta1 ,beta2=para$beta2, alpha1=para$alpha1,alpha2=para$alpha2),
                        control=nls.control(maxiter = 2000, tol = 1e-09, minFactor = 1/20240000, printEval = FALSE, warnOnly = T))
    )
    # fit<-try(nlm1<- nls(y ~ b*pnorm(beta1+beta2*x)+a, start=list(a=as.numeric(para[3]) ,b=as.numeric(para[4])
    #                                                     ,beta1=as.numeric(para[1]) ,beta2=as.numeric(para[2])),
    #            control=nls.control(maxiter = 100000, tol = 1e-08, minFactor = 1/20240000, printEval = FALSE, warnOnly = T))
    # )
    return(list(fit=fit,residual=para$residual,para=para))
  }


  Date = rownames(data_clust)
  # interval_estimated
  #cluster 1(HB)
  data=log(1+data_clust[1:dim(data_clust)[1],'cluster1'])
  qq <- c(1,sort(final[[1]]),dim(data_clust)[1])
  # plot fitting data
  # dataplot <- data.frame(y=as.vector(data),phase=c(rep('a',14),rep('b',33),rep('c',10),rep('d',15),rep('e',18),rep('f',52)))
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  fig1<-ggplot(dataplot,aes(Date))
  # color<-colorRampPalette(c('red','orange'))(10)
  # color<-colorRampPalette(c('#4B0082','#DB7093'))(6)
  phase<-c(1,2,3,4,5,6)
  phase_num<-c(14,33,10,15,18,52)
  xintercept<-c(date[14],date[47],date[57],date[72],date[90])
  fig1 <- fig1 +
    geom_point(data = dataplot, aes(y=y)) +
    labs(title = expression(italic(paste('Cluster 1(HB)'))),x='', y=expression(paste("log(1+infection)")),fill="") +
    mytheme
  fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
  df.delta_data <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                              lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
      nlm <- fun_fit(a,b,data[a:b])$fit
      data_one <- data.frame(x.new=dataplot[a:b,2],f.new=c(data[a:(a+1)],fitted(nlm)),lwr.conf=c(data[a:(a+1)],fitted(nlm))
                             ,upr.conf=c(data[a:(a+1)],fitted(nlm)),lwr.pred=c(data[a:(a+1)],fitted(nlm)),upr.pred=c(data[a:(a+1)],fitted(nlm)))
      data_one['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data <- rbind(df.delta_data,data_one)
      # fig1 <- fig1 +
      #   geom_line(data=data_one, aes(x=x, y=y), colour=color[i-1], size=0.7)
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
      x.new <- seq((a+2), qq[i], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1<-data[(a+1):(b-1)]#t-1
      z2<-data[(a):(b-2)]#t-2
      f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1,z2)
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
      df.delta_data<-rbind(df.delta_data,df.delta)
    }
  }
  color <- c('red','blue','green','pink','yellow','red','blue','pink','green','yellow')
  xintercept<-c(date[14],date[47],date[57],date[72],date[90])
  fig1 <- fig1 +
    # geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.pred, ymax=upr.pred), alpha=0.2, fill='black') +
    geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.conf, ymax=upr.conf), alpha=0.4,fill=color[1]) +#, group=phase,fill=phase
    geom_line(data=df.delta_data, aes(x=x.new, y=f.new), size=1,colour=color[1])+#, colour=phase
    ylim(0,12.5)+
    # scale_fill_gradient(low = 'blue', high = 'red')+
    # scale_color_gradient(low = 'blue', high = 'red',guide=FALSE)+
    geom_line(data=data.frame(x=c(date[14],date[14]),y=c(data[14]-0.7,data[14]+1)),aes(x=x,y=y),linetype='longdash',size=1,color=color[1])+
    geom_line(data=data.frame(x=c(date[47],date[47]),y=c(data[47]-1,data[47]+1)),aes(x=x,y=y),linetype='longdash',size=1,color=color[1])+
    geom_line(data=data.frame(x=c(date[57],date[57]),y=c(data[57]-1,data[57]+1)),aes(x=x,y=y),linetype='longdash',size=1,color=color[1])+
    geom_line(data=data.frame(x=c(date[72],date[72]),y=c(data[72]-1,data[72]+1)),aes(x=x,y=y),linetype='longdash',size=1,color=color[1])+
    geom_line(data=data.frame(x=c(date[90],date[90]),y=c(data[90]-1,data[90]+1)),aes(x=x,y=y),linetype='longdash',size=1,color=color[1])
  fig1 <- fig1 +
    scale_x_date(breaks =as.Date(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),date_labels="%y/%m/%d")+
    annotate(geom='text',x=as.Date(c("2019-12-28")),y=1.5,label="2019-12-14",size=7)+
    annotate(geom='text',x=as.Date(c("2020-01-30")),y=3.8,label="2020-01-16",size=7)+
    annotate(geom='text',x=as.Date(c("2020-02-09")),y=6.8,label="2020-01-26",size=7)+
    annotate(geom='text',x=as.Date(c("2020-01-25")),y=10.5,label="2020-02-10",size=7)+
    annotate(geom='text',x=as.Date(c("2020-03-19")),y=10.5,label="2020-02-28",size=7)+
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
  # dataplot <- data.frame(y=as.vector(data),phase=c(rep('a',14),rep('b',33),rep('c',10),rep('d',15),rep('e',18),rep('f',52)))
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  # color<-colorRampPalette(c('red','orange'))(10)
  phase<-c(11,12,13,14)
  phase_num<-c(54,10,17,61)
  fig2<-ggplot(dataplot,aes(Date))
  fig2 <- fig2 +
    geom_point(data = dataplot, aes(y=y)) +
    labs(title = expression(italic(paste('Cluster 2(GD,ZJ,HA,HN,AH,JX,SD,JS,SC,CQ,BJ,SH,HL,\n,HE,FJ,SN,GX,YN,HI,LN,SX,TJ,GZ,GS,NM,JL,XJ,NX)'))),
         x='', y=expression(paste("log(1+infection)")),fill="") +
    mytheme
  fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
  df.delta_data <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                              lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
      nlm <- fun_fit(a,b,data[a:b])$fit
      data_one <- data.frame(x.new=dataplot[a:b,2],f.new=c(data[a:(a+1)],fitted(nlm)),lwr.conf=c(data[a:(a+1)],fitted(nlm))
                             ,upr.conf=c(data[a:(a+1)],fitted(nlm)),lwr.pred=c(data[a:(a+1)],fitted(nlm)),upr.pred=c(data[a:(a+1)],fitted(nlm)))
      data_one['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data <- rbind(df.delta_data,data_one)
      # fig1 <- fig1 +
      #   geom_line(data=data_one, aes(x=x, y=y), colour=color[i-1], size=0.7)
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
      x.new <- seq((a+2), qq[i], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1<-data[(a+1):(b-1)]#t-1
      z2<-data[(a):(b-2)]#t-2
      f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1,z2)
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
      df.delta_data<-rbind(df.delta_data,df.delta)
    }
  }
  fig2 <- fig2 +
    # geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.pred, ymax=upr.pred), alpha=0.2, fill='black') +
    geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.conf, ymax=upr.conf), alpha=0.4, fill=color[2]) +
    geom_line(data=df.delta_data, aes(x=x.new, y=f.new), size=1, colour=color[2])+
    ylim(-0.1,12.5)+
    geom_line(data=data.frame(x=c(date[54],date[54]),y=c(data[54]-1,data[54]+1)),aes(x=x,y=y),linetype='longdash',size=1,color=color[2])+
    geom_line(data=data.frame(x=c(date[64],date[64]),y=c(data[64]-1,data[64]+1)),aes(x=x,y=y),linetype='longdash',size=1,color=color[2])+
    geom_line(data=data.frame(x=c(date[81],date[81]),y=c(data[81]-1,data[81]+1)),aes(x=x,y=y),linetype='longdash',size=1,color=color[2])

  fig2 <- fig2 +
    scale_x_date(breaks =as.Date(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),date_labels="%y/%m/%d")+
    annotate(geom='text',x=as.Date(c("2020-02-05")),y=5.3,label="2020-1-23",size=7)+
    annotate(geom='text',x=as.Date(c("2020-01-15")),y=8.2,label="2020-02-02",size=7)+
    annotate(geom='text',x=as.Date(c("2020-03-04")),y=10,label="2020-02-19",size=7)+
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
  # dataplot <- data.frame(y=as.vector(data),phase=c(rep('a',14),rep('b',33),rep('c',10),rep('d',15),rep('e',18),rep('f',52)))
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  # color<-colorRampPalette(c('red','orange'))(10)
  phase<-c('Cluster 3','Cluster 3','Cluster 3','Cluster 3')
  phase_num<-c(56,20,15,51)
  fig34<-ggplot(dataplot,aes(Date))
  fig34 <- fig34 +
    geom_point(data = dataplot, aes(y=y)) +
    labs(title = expression(italic(paste('Cluster 3(TW) and Cluster 4(HK)'))),
         x='', y=expression(paste("log(1+infection)")),fill="") +
    mytheme
  fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
  df.delta_data1 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                              lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
      nlm <- fun_fit(a,b,data[a:b])$fit
      data_one <- data.frame(x.new=dataplot[a:b,2],f.new=c(data[a:(a+1)],fitted(nlm)),lwr.conf=c(data[a:(a+1)],fitted(nlm))
                             ,upr.conf=c(data[a:(a+1)],fitted(nlm)),lwr.pred=c(data[a:(a+1)],fitted(nlm)),upr.pred=c(data[a:(a+1)],fitted(nlm)))
      data_one['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data1 <- rbind(df.delta_data1,data_one)
      # fig1 <- fig1 +
      #   geom_line(data=data_one, aes(x=x, y=y), colour=color[i-1], size=0.7)
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
      x.new <- seq((a+2), qq[i], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1<-data[(a+1):(b-1)]#t-1
      z2<-data[(a):(b-2)]#t-2
      f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1,z2)
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
  #cluster4 data
  data=log(1+data_clust[1:dim(data_clust)[1],'cluster4'])
  data4<-data
  qq <- c(1,sort(final[[4]]),dim(data_clust)[1])
  # plot fitting data
  # dataplot <- data.frame(y=as.vector(data),phase=c(rep('a',14),rep('b',33),rep('c',10),rep('d',15),rep('e',18),rep('f',52)))
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  # color<-colorRampPalette(c('red','orange'))(10)
  phase<-c('Cluster 4','Cluster 4','Cluster 4','Cluster 4')
  phase_num<-c(57,11,15,59)
  # fig2<-ggplot(dataplot,aes(Date))
  fig34 <- fig34 +
    geom_point(data = dataplot, aes(y=y)) +
    # labs(title = expression(italic(paste('Cluster 4()'))),
         # x='', y=expression(paste("log(1+infection)")),fill="") +
    mytheme
  fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
  df.delta_data2 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                               lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
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
  df.delta_data<-rbind(df.delta_data1,df.delta_data2)
  fig34 <- fig34 +
    # geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.pred, ymax=upr.pred, group=phase), alpha=0.2, fill='black') +
    geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.conf, ymax=upr.conf, group=phase,fill=phase), alpha=0.4) +
    geom_line(data=df.delta_data, aes(x=x.new, y=f.new,colour=phase, group=phase), size=1)+
    ylim(-0.1,12.5)+
    scale_fill_manual(values=c('red','blue'),name='Province')+
    scale_colour_manual(values=c('red','blue'),name='Province')+
    #cluster3
    geom_line(data=data.frame(x=c(date[56],date[56]),y=c(data3[56]-1,data3[56]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='red')+
    geom_line(data=data.frame(x=c(date[76],date[76]),y=c(data3[76]-1,data3[76]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='red')+
    geom_line(data=data.frame(x=c(date[91],date[91]),y=c(data3[91]-1,data3[91]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='red')+
    #cluster4
    geom_line(data=data.frame(x=c(date[57],date[57]),y=c(data4[57]-1,data4[57]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='blue')+
    geom_line(data=data.frame(x=c(date[68],date[68]),y=c(data4[68]-1,data4[68]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='blue')+
    geom_line(data=data.frame(x=c(date[83],date[83]),y=c(data4[83]-1,data4[83]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='blue')
  fig34 <- fig34 +
    scale_x_date(breaks =as.Date(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),date_labels="%y/%m/%d")+
    # cluster4HK
    annotate(geom='text',x=as.Date(c("2020-01-12")),y=2,label="2020-01-26",size=7)+
    annotate(geom='text',x=as.Date(c("2020-01-23")),y=3.5,label="2020-02-06",size=7)+
    annotate(geom='text',x=as.Date(c("2020-02-07")),y=4.8,label="2020-02-21",size=7)+
    #cluster4TW
    annotate(geom='text',x=as.Date(c("2020-02-08")),y=1.5,label="2020-01-25",size=7)+
    annotate(geom='text',x=as.Date(c("2020-02-28")),y=2.5,label="2020-02-14",size=7)+
    annotate(geom='text',x=as.Date(c("2020-03-14")),y=3.5,label="2020-02-29",size=7)+
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
  # dataplot <- data.frame(y=as.vector(data),phase=c(rep('a',14),rep('b',33),rep('c',10),rep('d',15),rep('e',18),rep('f',52)))
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  # color<-colorRampPalette(c('red','orange'))(10)
  phase<-c('Cluster 5','Cluster 5','Cluster 5')
  phase_num<-c(56,42,44)
  fig567<-ggplot(dataplot,aes(Date))
  fig567 <- fig567 +
    geom_point(data = dataplot, aes(y=y)) +
    labs(title = expression(italic(paste('Cluster 5(MO), Cluster 6(QH) and Cluster 7(XZ)'))),
         x='', y=expression(paste("log(1+infection)")),fill="") +
    mytheme
  fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
  df.delta_data5 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                              lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
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

  # cluster6
  data=log(1+data_clust[1:dim(data_clust)[1],'cluster6'])
  data6<-data
  qq <- c(1,sort(final[[6]]),dim(data_clust)[1])
  # plot fitting data
  # dataplot <- data.frame(y=as.vector(data),phase=c(rep('a',14),rep('b',33),rep('c',10),rep('d',15),rep('e',18),rep('f',52)))
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  # color<-colorRampPalette(c('red','orange'))(10)
  phase<-c('Cluster 6','Cluster 6')
  phase_num<-c(61,81)
  fig567 <- fig567 +
    geom_point(data = dataplot, aes(y=y)) +
    # labs(title = expression(italic(paste('Cluster 5(MO),Cluster 6(QH)','Cluster 7(XZ)'))),
    #      x='', y=expression(paste("log(1+infection)")),fill="") +
    mytheme
  fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
  df.delta_data6 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                               lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())
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
  #cluster7
  data=log(1+data_clust[1:dim(data_clust)[1],'cluster7'])
  data7<-data
  qq <- c(1,sort(final[[7]]),dim(data_clust)[1])
  # plot fitting data
  # dataplot <- data.frame(y=as.vector(data),phase=c(rep('a',14),rep('b',33),rep('c',10),rep('d',15),rep('e',18),rep('f',52)))
  dataplot <- data.frame(y=as.vector(data))
  date <- as.Date(Date,"%m/%d/%Y")
  dataplot[['Date']] <- date
  # color<-colorRampPalette(c('red','orange'))(10)
  phase<-c('Cluster 7')
  phase_num<-c(142)
  fig567 <- fig567 +
    geom_point(data = dataplot, aes(y=y)) +
    # labs(title = expression(italic(paste('Cluster 5(MO),Cluster 6(QH) and Cluster 7(XZ)'))),
    #      x='', y=expression(paste("log(1+infection)")),fill="") +
    mytheme
  fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
  df.delta_data7 <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                               lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())

  for(i in 2:length(qq)){
    b<-as.numeric(qq[i])
    if((i-1)==1){
      cat('a1\n')
      a <- as.numeric(qq[i-1])
      nlm <- fun_fit(a,b,data[a:b])$fit
      data_one <- data.frame(x.new=dataplot[a:b,2],f.new=c(data[a:(a+1)],fitted(nlm)),lwr.conf=c(data[a:(a+1)],fitted(nlm))
                             ,upr.conf=c(data[a:(a+1)],fitted(nlm)),lwr.pred=c(data[a:(a+1)],fitted(nlm)),upr.pred=c(data[a:(a+1)],fitted(nlm)))
      data_one['phase']<-rep(phase[i-1],phase_num[i-1])
      df.delta_data7 <- rbind(df.delta_data7,data_one)
      # fig1 <- fig1 +
      #   geom_line(data=data_one, aes(x=x, y=y), colour=color[i-1], size=0.7)
    }else{
      cat('b1\n')
      a<-as.numeric(qq[i-1])+1
      x.new <- seq((a+2), qq[i], by=1)
      nlm2 <- fun_fit(a,b,data[a:b])$fit
      beta2.est <- coef(nlm2)
      z1<-data[(a+1):(b-1)]#t-1
      z2<-data[(a):(b-2)]#t-2
      f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1,z2)
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
  df.delta_data<-rbind(df.delta_data5,df.delta_data6,df.delta_data7)
  fig567 <- fig567 +
    # geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.pred, ymax=upr.pred, group=phase), alpha=0.2, fill='black') +
    geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.conf, ymax=upr.conf, group=phase,fill=phase), alpha=0.4) +
    geom_line(data=df.delta_data, aes(x=x.new, y=f.new,colour=phase, group=phase), size=1)+
    ylim(-0.1,12.5)+
    scale_fill_manual(values=c('red','blue','#FF8C00'),name='Province')+
    scale_colour_manual(values=c('red','blue','#FF8C00'),name='Province')+
    #cluster5
    geom_line(data=data.frame(x=c(date[56],date[56]),y=c(data5[56]-1,data5[56]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='red')+
    geom_line(data=data.frame(x=c(date[98],date[98]),y=c(data5[98]-1,data5[98]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='red')+
    #cluster6
    geom_line(data=data.frame(x=c(date[61],date[61]),y=c(data6[61]-1,data6[61]+1)),aes(x=x,y=y),linetype='longdash',size=1,color='blue')
  fig567 <- fig567 +
    scale_x_date(breaks =as.Date(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),date_labels="%y/%m/%d")+
    #cluster5
    annotate(geom='text',x=as.Date(c("2020-01-10")),y=1,label="2020-1-25",size=7)+
    annotate(geom='text',x=as.Date(c("2020-03-22")),y=2,label="2020-03-07",size=7)+
    # cluster6
    annotate(geom='text',x=as.Date(c("2020-01-15")),y=2.5,label="2020-01-30",size=7)+
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
  ggsave(filename='china.pdf', plot = tu, width = 40, height = 8, dpi = 300,limitsize = FALSE)
}


