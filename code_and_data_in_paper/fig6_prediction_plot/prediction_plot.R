setwd('')

load('pred_plot_data_C.RData')
load('pred_plot_data.RData')

# china cluster1 2/10-4/20
# country cluster2 2/18-4-20, cluster7 2/4-4/20


par(mfrow=c(1,3),mar=c(3,3,4,2),oma=c(0,0,2,0))
# china cluster1 2/10-4/20
title_cluster_1=expression(italic(paste('Cluster 1(HB)')))
data=log(1+data_clust[1:dim(data_clust)[1],'cluster1'])
qq <- c(1,sort(final[[1]]),dim(data_clust)[1])
a=qq[5]
b=142
x=seq(a,b,by=1)
x_pred=seq(qq[6],b,by=1)
orgin_data=data[a:b]
pred_data=pred_plot[[1]][[5]]
plot(x,orgin_data,type='l',lwd=3,ylim = c(10.2,11.5),main = list(title_cluster_1),xaxt="n")
lines(x_pred,pred_data$f.new,lty=2,xaxt="n")
axis(1,c(1,32,63,92,123,142),as.character(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),cex.axis=1)

# country cluster2 2/18-4-20
title_cluster_2=expression(italic(paste('Cluster 2(IR,IT)')))
data=log(1+data_clust_C[1:dim(data_clust_C)[1],'cluster2'])
qq <- c(1,sort(final_C[[2]]),dim(data_clust_C)[1])
a=qq[2]
b=142
x=seq(a,b,by=1)
x_pred=seq(qq[3],b,by=1)
orgin_data=data[a:b]
pred_data=pred_plot_C[[1]][[2]]
plot(x,orgin_data,type='l',lwd=3,ylim = c(0,15),main = list(title_cluster_2),xaxt="n")
lines(x_pred,pred_data$f.new,lty=2,xaxt="n")
axis(1,c(1,32,63,92,123,142),as.character(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),cex.axis=1)


# country cluster7 2/4-4/20
title_cluster_7=expression(italic(paste('Cluster 7(CH,NL,BE,AT,SE,DK,PT,BR,IE,\n RO,ID,PL,PE,EC,RU,IN,PH,EG,DZ,MK,DO)')))
data=log(1+data_clust_C[1:dim(data_clust_C)[1],'cluster7'])
date=names(data)
qq <- c(1,sort(final_C[[7]]),dim(data_clust_C)[1])
a=qq[2]
b=142
x=seq(a,b,by=1)
x_pred=seq(qq[3],b,by=1)
orgin_data=data[a:b]
pred_data=pred_plot_C[[2]][[2]]
plot(x,orgin_data,type='l',lwd=3,ylim = c(0,15),main = list(title_cluster_7),xaxt="n")
lines(x_pred,pred_data$f.new,lty=2,xaxt="n")
axis(1,c(1,32,63,92,123,142),as.character(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),cex.axis=1)














library(readr)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(tidyr)


# initial parameter satting
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("test_new_a.cpp")

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

title_cluster_1=expression(italic(paste('Cluster 1(CN)')))


# one interval_estimated
qq <- c(1,sort(final[[1]]),dim(data_clust)[1])[5:7]
data=log(1+data_clust[1:dim(data_clust)[1],'cluster1'])
Date = as.Date(rownames(data_clust),"%m/%d/%Y")
# plot fitting data
dataplot <- data.frame(y=as.vector(data[(qq[1]+1):qq[3]]))
date <- as.Date(Date,"%m/%d/%Y")[(qq[1]+1):qq[3]]
dataplot[['Date']] <- date
fig1<-ggplot(dataplot,aes(Date))
phase<-c(1,2,3,4,5,6)
phase_num<-c(final[[1]][1])
if(length(final[[1]])>1){
  for(i in 2:length(final[[1]])){
    phase_num <- c(phase_num, final[[1]][i]-final[[1]][i-1])
  }
}
phase_num <- c(phase_num, dim(data_clust)[1]-final[[1]][length(final[[1]])])[5:6]
phase_num[1]=phase_num[1]
fig1 <- fig1 +
geom_point(data = dataplot, aes(y=y)) +
labs(title = expression(italic(paste('Cluster 1(HB)'))),x='', y=expression(paste("log(1+infection)")),fill="") +
mytheme
fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
df.delta_data <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                            lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())

#-----------------------------------------------------------------------------------------------#
# this section aim to predict the next part according to previous part
pred_x.new_1=list()
pred_f.new_1=list()
if((length(qq)-2)!=0){
  for (item in 1:(length(qq)-2)) {
    sat=item
    if((i-1)==1){
      cat('a1\n')
      a<-as.numeric(qq[sat])
    }else{
      cat('b1\n')
      a<-as.numeric(qq[sat])+1
    }
    b=as.numeric(qq[sat+1])
    x.new <- seq((b+1), qq[sat+2], by=1)
    nlm2 <- fun_fit(a,b,data[a:b])$fit
    beta2.est <- coef(nlm2)
    z1_num<-data[b]#t-1
    z2_num<-data[b-1]#t-2
    pred_x.new_1[[item]] <- Date[(b+1):qq[sat+2]]
    cat('part',item+1,'a=',beta2.est[1],'b=',beta2.est[2],'beta1=',beta2.est[3],'beta2=',beta2.est[4],'alpha1=', beta2.est[5],'alpha2=',beta2.est[6],'\n')
    pred_f.new_1[[item]] <- predict_fun(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1_num,z2_num)
  }
}
#-----------------------------------------------------------------------------------------------#

for(i in 2:length(qq)){
  b<-as.numeric(qq[i])
  if((i-1)==1){
    cat('a1\n')
    a <- as.numeric(qq[i-1])+1
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
    df.delta <- data.frame(x.new=Date[a:b], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new), upr.conf=c(data[a:(a+1)],f.new))
  }else{
    df.delta <- data.frame(x.new=Date[a:b], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new-deltaf), upr.conf=c(data[a:(a+1)],f.new+deltaf))
  }
  head(df.delta)
  sigma2.est <- summary(nlm2)$sigma
  deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
  line_dataframa<-data.frame(x.new=Date[a:b], f.new=as.numeric(c(data[a:(a+1)],f.new)))
  if(max(deltay)>20){
    df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new),c(data[a:(a+1)],f.new))
  }else{
    df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new - deltay),c(data[a:(a+1)],f.new + deltay))
  }
  df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
  df.delta_data<-rbind(df.delta_data,df.delta)
}

color <- c('red','blue','green','pink','yellow','red','blue','pink','green','yellow')
fig1 <- fig1 +
  # geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.pred, ymax=upr.pred), alpha=0.2, fill='black') +
  geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.conf, ymax=upr.conf), alpha=0.4,fill="red") +#, group=phase,fill=phase
  geom_line(data=df.delta_data, aes(x=x.new, y=f.new), size=1,colour="red")+#, colour=phase
  ylim(10,11.5)+
  # scale_fill_gradient(low = 'blue', high = 'red')+
  # scale_color_gradient(low = 'blue', high = 'red',guide=FALSE)+
  geom_line(data=data.frame(x=c(Date[90],Date[90]),y=c(10,11.5)),aes(x=x,y=y),linetype='solid',size=1,color='black')
  # geom_line(data=data.frame(x=c(Date[90],Date[90]),y=c(data[90]-0.1,data[90]+0.1)),aes(x=x,y=y),linetype='longdash',size=1,color="red")
fig1 <- fig1 +
  scale_x_date(breaks =as.Date(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),date_labels="%y/%m/%d")+
  annotate(geom='text',x=as.Date(c("2020-03-08")),y=11,label="2020-02-28",size=7)+
  theme_bw() +
  mytheme+
  labs(tag = "A") +
  theme(plot.tag.position = c(0.05, 1))+
  theme(plot.tag=element_text(size = 30))

if(length(pred_x.new_1)!=0){
  for (i in 1:length(pred_x.new_1)) {
    cat(pred_x.new_1[[i]],'\n')
    pred_data=data.frame(x.new=pred_x.new_1[[i]],f.new=pred_f.new_1[[i]])
    fig1 <- fig1 +
    geom_line(data=pred_data, aes(x=x.new, y=f.new),shape = 2, size=1,colour='black',linetype='dashed')
  }
}
fig1




# two interval_estimated
title_cluster_2=expression(italic(paste('Cluster 2(IR,IT)')))
qq <- c(1,sort(final_C[[2]]),dim(data_clust_C)[1])[2:4]
data=log(1+data_clust_C[1:dim(data_clust_C)[1],'cluster2'])
Date = as.Date(rownames(data_clust_C),"%m/%d/%Y")
# plot fitting data
dataplot <- data.frame(y=as.vector(data[(qq[1]+1):qq[3]]))
date <- as.Date(Date,"%m/%d/%Y")[(qq[1]+1):qq[3]]
dataplot[['Date']] <- date
fig2<-ggplot(dataplot,aes(Date))
phase<-c(1,2,3,4,5,6)
phase_num<-c(final_C[[2]][1])
if(length(final_C[[2]])>1){
  for(i in 2:length(final_C[[2]])){
    phase_num <- c(phase_num, final_C[[2]][i]-final_C[[2]][i-1])
  }
}
phase_num <- c(phase_num, dim(data_clust_C)[1]-final_C[[2]][length(final_C[[2]])])[2:3]
phase_num[1]=phase_num[1]
fig2 <- fig2 +
  geom_point(data = dataplot, aes(y=y)) +
  labs(title = title_cluster_2,x='', y=expression(paste("log(1+infection)")),fill="") +
  mytheme
fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
df.delta_data <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                            lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())

#-----------------------------------------------------------------------------------------------#
# this section aim to predict the next part according to previous part
pred_x.new_1=list()
pred_f.new_1=list()
if((length(qq)-2)!=0){
  for (item in 1:(length(qq)-2)) {
    sat=item
    if((item)==1){
      cat('a1\n')
      a<-as.numeric(qq[sat])+1
    }else{
      cat('b1\n')
      a<-as.numeric(qq[sat])+1
    }
    b=as.numeric(qq[sat+1])
    x.new <- seq((b+1), qq[sat+2], by=1)
    nlm2 <- fun_fit(a,b,data[a:b])$fit
    beta2.est <- coef(nlm2)
    z1_num<-data[b]#t-1
    z2_num<-data[b-1]#t-2
    pred_x.new_1[[item]] <- Date[(b+1):qq[sat+2]]
    cat('part',item+1,'a=',beta2.est[1],'b=',beta2.est[2],'beta1=',beta2.est[3],'beta2=',beta2.est[4],'alpha1=', beta2.est[5],'alpha2=',beta2.est[6],'\n')
    pred_f.new_1[[item]] <- predict_fun(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1_num,z2_num)
  }
}
#-----------------------------------------------------------------------------------------------#

for(i in 2:length(qq)){
  b<-as.numeric(qq[i])
  if((i-1)==1){
    cat('a1\n')
    a <- as.numeric(qq[i-1])+1
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
  V.beta2 <- vcov(nlm2)
  GS=rowSums((g.new[,1:6]%*%V.beta2)%*%t(g.new[,1:6]))
  head(GS)
  alpha <- 0.05
  df <- length(data[a:b])-length(beta2.est)#freedom degree
  deltaf <- sqrt(GS)*qt(1-alpha/2,df)
  if(max(deltaf)>20){
    df.delta <- data.frame(x.new=Date[a:b], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new), upr.conf=c(data[a:(a+1)],f.new))
  }else{
    df.delta <- data.frame(x.new=Date[a:b], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new-deltaf), upr.conf=c(data[a:(a+1)],f.new+deltaf))
  }
  head(df.delta)
  sigma2.est <- summary(nlm2)$sigma
  deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
  line_dataframa<-data.frame(x.new=Date[a:b], f.new=as.numeric(c(data[a:(a+1)],f.new)))
  if(max(deltay)>20){
    df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new),c(data[a:(a+1)],f.new))
  }else{
    df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new - deltay),c(data[a:(a+1)],f.new + deltay))
  }
  df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
  df.delta_data<-rbind(df.delta_data,df.delta)
}

color <- c('red','blue','green','pink','yellow','red','blue','pink','green','yellow')
fig2 <- fig2 +
  geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.conf, ymax=upr.conf), alpha=0.4,fill="red") +#, group=phase,fill=phase
  geom_line(data=df.delta_data, aes(x=x.new, y=f.new), size=1,colour="red")+#, colour=phase
  ylim(0,15)+
  geom_line(data=data.frame(x=c(Date[92],Date[92]),y=c(0,15)),aes(x=x,y=y),linetype='solid',size=1,color='black')
fig2 <- fig2 +
  scale_x_date(breaks =as.Date(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),date_labels="%y/%m/%d")+
  annotate(geom='text',x=as.Date(c("2020-03-09")),y=7,label="2020-03-01",size=7)+
  theme_bw() +
  mytheme+
  labs(tag = "B") +
  theme(plot.tag.position = c(0.05, 1))+
  theme(plot.tag=element_text(size = 30))

if(length(pred_x.new_1)!=0){
  for (i in 1:length(pred_x.new_1)) {
    cat(pred_x.new_1[[i]],'\n')
    pred_data=data.frame(x.new=pred_x.new_1[[i]],f.new=pred_f.new_1[[i]])
    fig2 <- fig2 +
      geom_line(data=pred_data, aes(x=x.new, y=f.new),shape = 2, size=1,colour='black',linetype='dashed')
  }
}
fig2




# three interval_estimated
title_cluster_7=expression(italic(paste('Cluster 7(CH,NL,BE,AT,SE,DK,PT,BR,IE,\n,RO,ID,PL,PE,EC,RU,IN,PH,EG,DZ,MK,DO)')))
qq <- c(1,sort(final_C[[7]]),dim(data_clust_C)[1])[2:4]
data=log(1+data_clust_C[1:dim(data_clust_C)[1],'cluster7'])
Date = as.Date(rownames(data_clust_C),"%m/%d/%Y")
# plot fitting data
dataplot <- data.frame(y=as.vector(data[(qq[1]+1):qq[3]]))
date <- as.Date(Date,"%m/%d/%Y")[(qq[1]+1):qq[3]]
dataplot[['Date']] <- date
fig3<-ggplot(dataplot,aes(Date))
phase<-c(1,2,3,4,5,6)
phase_num<-c(final_C[[7]][1])
if(length(final_C[[7]])>1){
  for(i in 2:length(final_C[[7]])){
    phase_num <- c(phase_num, final_C[[7]][i]-final_C[[7]][i-1])
  }
}
phase_num <- c(phase_num, dim(data_clust_C)[1]-final_C[[7]][length(final_C[[7]])])[2:3]
phase_num[1]=phase_num[1]
fig3 <- fig3 +
  geom_point(data = dataplot, aes(y=y)) +
  labs(title = title_cluster_7,x='', y=expression(paste("log(1+infection)")),fill="") +
  mytheme
fgh2 <- deriv(y ~ alpha1*z1+alpha2*z2+b*pnorm(beta1+beta2*x)+a, c("a", "b", "beta1", "beta2",'alpha1','alpha2'), function(a, b, beta1, beta2, alpha1,alpha2,x,z1,z2){} ) 
df.delta_data <- data.frame(x.new= c(), f.new= numeric(), lwr.conf=numeric(),upr.conf=numeric(),
                            lwr.pred=numeric(),upr.pred=numeric(),phase=numeric())

#-----------------------------------------------------------------------------------------------#
# this section aim to predict the next part according to previous part
pred_x.new_1=list()
pred_f.new_1=list()
if((length(qq)-2)!=0){
  for (item in 1:(length(qq)-2)) {
    sat=item
    if(item==1){
      cat('a1\n')
      a<-as.numeric(qq[sat])+1
    }else{
      cat('b1\n')
      a<-as.numeric(qq[sat])+1
    }
    b=as.numeric(qq[sat+1])
    x.new <- seq((b+1), qq[sat+2], by=1)
    nlm2 <- fun_fit(a,b,data[a:b])$fit
    beta2.est <- coef(nlm2)
    z1_num<-data[b]#t-1
    z2_num<-data[b-1]#t-2
    pred_x.new_1[[item]] <- Date[(b+1):qq[sat+2]]
    cat('part',item+1,'a=',beta2.est[1],'b=',beta2.est[2],'beta1=',beta2.est[3],'beta2=',beta2.est[4],'alpha1=', beta2.est[5],'alpha2=',beta2.est[6],'\n')
    pred_f.new_1[[item]] <- predict_fun(beta2.est[1],beta2.est[2],beta2.est[3],beta2.est[4], beta2.est[5],beta2.est[6],x.new,z1_num,z2_num)
  }
}
#-----------------------------------------------------------------------------------------------#

for(i in 2:length(qq)){
  b<-as.numeric(qq[i])
  if((i-1)==1){
    cat('a1\n')
    a <- as.numeric(qq[i-1])+1
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
  V.beta2 <- vcov(nlm2)
  GS=rowSums((g.new[,1:6]%*%V.beta2)%*%t(g.new[,1:6]))
  head(GS)
  alpha <- 0.05
  df <- length(data[a:b])-length(beta2.est)#freedom degree
  deltaf <- sqrt(GS)*qt(1-alpha/2,df)
  if(max(deltaf)>20){
    df.delta <- data.frame(x.new=Date[a:b], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new), upr.conf=c(data[a:(a+1)],f.new))
  }else{
    df.delta <- data.frame(x.new=Date[a:b], f.new=c(data[a:(a+1)],f.new), lwr.conf=c(data[a:(a+1)],f.new-deltaf), upr.conf=c(data[a:(a+1)],f.new+deltaf))
  }
  head(df.delta)
  sigma2.est <- summary(nlm2)$sigma
  deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,df)
  line_dataframa<-data.frame(x.new=Date[a:b], f.new=as.numeric(c(data[a:(a+1)],f.new)))
  if(max(deltay)>20){
    df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new),c(data[a:(a+1)],f.new))
  }else{
    df.delta[c("lwr.pred","upr.pred")] <- cbind(c(data[a:(a+1)],f.new - deltay),c(data[a:(a+1)],f.new + deltay))
  }
  df.delta['phase']<-rep(phase[i-1],phase_num[i-1])
  df.delta_data<-rbind(df.delta_data,df.delta)
}

color <- c('red','blue','green','pink','yellow','red','blue','pink','green','yellow')
fig3 <- fig3 +
  geom_ribbon(data=df.delta_data, aes(x=x.new, ymin=lwr.conf, ymax=upr.conf), alpha=0.4,fill="red") +#, group=phase,fill=phase
  geom_line(data=df.delta_data, aes(x=x.new, y=f.new), size=1,colour="red",linetype='solid')+#, colour=phase
  ylim(0,15)+
  geom_line(data=data.frame(x=c(Date[96],Date[96]),y=c(0,15)),aes(x=x,y=y),linetype='solid',size=1,color='black')
fig3 <- fig3 +
  scale_x_date(breaks =as.Date(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),date_labels="%y/%m/%d")+
  annotate(geom='text',x=as.Date(c("2020-03-15")),y=6,label="2020-03-05",size=7)+
  theme_bw() +
  mytheme+
  labs(tag = "C") +
  theme(plot.tag.position = c(0.05, 1))+
  theme(plot.tag=element_text(size = 30))

if(length(pred_x.new_1)!=0){
  for (i in 1:length(pred_x.new_1)) {
    cat(pred_x.new_1[[i]],'\n')
    pred_data=data.frame(x.new=pred_x.new_1[[i]],f.new=pred_f.new_1[[i]])
    fig3 <- fig3 +
      geom_line(data=pred_data, aes(x=x.new, y=f.new),size=1,colour='black',linetype='dashed')
  }
}

fig3


tu <- grid.arrange(fig1,fig2,fig3,nrow=1,ncol=3)
# ggsave(filename='pred.pdf', plot = tu, width = 30, height = 8, dpi = 300,limitsize = FALSE)


