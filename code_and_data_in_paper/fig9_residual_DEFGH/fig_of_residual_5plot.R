library(readr)
library(ggplot2)
# library(purrr)
# library(stringr)
library(dplyr)
library(grid)
library(gridExtra)
library(tidyr)

mytheme <- theme(plot.title = element_text(hjust = 0.5, size = 22),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.text.x = element_text(size = 18,color='black'),
                 axis.text.y = element_text(size = 18,color='black'),
                 axis.title.x=element_text(size=25),
                 axis.title.y=element_text(size=25),
                 plot.margin=unit(c(2,2,2,2), 'lines'))
setwd("")
load("cluster_data.RData")
load("origin_data.RData")
# initial parameter satting
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("test_new_a.cpp")

# initial parameter satting
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("test_new_a1.cpp")

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
# fitting data
fun_fit2 <- function(start,end,data){
  # fitting the origin data and return fitting results
  para=GridSearch2(start:end,data,-20,20,-20,20,100)
  # para=ini.par(beta00=-10,beta01=10,beta10=-10,beta11=10,x=c(start:end),y=data)
  y <- data
  x <- seq(start,end,1)
  fit<-try(nlm1<- nls(y ~ b*pnorm(beta1+beta2*x)+a, start=list(a=para$a ,b=para$b 
                                                               ,beta1=para$beta1 ,beta2=para$beta2),
                      control=nls.control(maxiter = 2000, tol = 1e-09, minFactor = 1/20240000, printEval = FALSE, warnOnly = T))
  )
  return(list(fit=fit,residual=para$residual))
}

data=log(1+data_clust[1:n,'cluster1'])

#all_four_residual_of_420000
qq<-c(14,47,57,72,90)
qq <- c(1,sort(qq),n)
# plot fitting data
dataplot <- data.frame(y=as.vector(data))
date <- as.Date(Date,"%m/%d/%Y")
dataplot[['Date']] <- date
all_summary <- vector(mode = "list")
all_residuals <- vector(mode = "list")
colnum = length(qq)-1
bbb<-c('A','A','A','A','B')
par(mfrow=c(1,5),mar = c(4,4,4,4),mai=c(0.3,0.3,0.6,0.3))
for(i in (length(qq)-2):(length(qq))){
  if(qq[i-1]==1){
    a<-as.numeric(qq[i-1])
  }else{
    a<-as.numeric(qq[i-1])+1
  }
  b<-as.numeric(qq[i])
  nlm3 <- fun_fit(a, b, data[a:b])$fit
  all_summary[[i-1]]<-summary(nlm3)
  residuals <- resid(nlm3)
  all_residuals[[i-1]]<-residuals
  var(residuals)
  # hist(residuals,main = paste("Histogram of residuals"),xlab = 'Residuals', ylab='Frequency')
  # mtext(paste("A.",i-1),side=3,adj=0,line = 0.5,font=3)
  # qqnorm(residuals,main = "Normal Q-Q Plot",xlab = "Theoretical Quantiles", ylab = "Sample Quantiles");qqline(residuals)
  # mtext(paste("B.",i-1),side=3,adj=0,line = 0.5,font=3)
  if(qq[i]==142){
    plot(residuals,main = '',xlab='',ylab='residuals',ylim = c(-0.006,0.006))#,ylim = c(-0.15,0.15)
    mtext("C",side=3,adj=0,line = 1,font=3)  #,Las=3
    }else{
    plot(residuals,main = '',xlab='',ylab='residuals')#,ylim = c(-0.15,0.15)
    mtext(paste(bbb[i-1]),side=3,adj=0,line = 1,font=3)  #,Las=3
  }
  # mtext(paste("C.",i-1),side=3,adj=0,line = 0.5,font=3,Las=3)
  # plot(fitted(nlm3),residuals,type='b')
}


#440000,changepoint c(1,54,77,142)
#420000,changepoint c(1,14,47,57,72,90,142)
four_pro<-list('HB'='HB','2'='')
four_pro_point<-list('HB'=c(1,14,47,57,72,90,142))
four_pro_name<-list('HB'='Scatter plot of residuals')
abcd<-c('D')
arima_residual<-c()
for(j in 1:1){
  cat(j,'\n')
  data=log(1+Data.new[1:n,four_pro[[j]]])
  qq1<-four_pro_point[[j]]
  x<-c()
  y_origin<-c()
  y_new<-c()
  for(i in 2:length(qq1)){
    y=data[(qq1[i-1]+1):qq1[i]]
    y1<-diff(y,1)
    am<-arima(y,order=c(2,1,0))
    coef<-am$coef
    cat(length(am$residuals),'\n')
    arima_residual<-c(arima_residual,am$residuals)
    y2=y1[2:(length(y)-1)]*coef[1]+y1[1:(length(y)-2)]*coef[2]
    cat(seq(qq1[i-1],qq1[i],1),'\n')
    # cat(y1,'\n')
    # cat(c(y1[qq1[i-1]:(qq1[i-1]+1)],y2),'\n')
    # if(i==2){
    #   x<-c(x,seq(qq1[i-1],qq1[i],1))
    # }else{
    #   x<-c(x,seq(qq1[i-1]+1,qq1[i],1))
    # }
    x<-c(x,seq(qq1[i-1]+1,qq1[i],1))
    y_origin<-c(y_origin,y1)
    y_new<-c(y_new,c(y2[1:1],y2))
  }
  cat('\n',length(x),length(arima_residual))
  plot(x,arima_residual,type='p',main='',ylab = 'residual',xlab = '',cex =0.8)  # main=four_pro_name[[j]]
  mtext(text=abcd[j],side=3,adj=0,line = 1,font=3) # ,Las=3
}


qq <- c(1,n)
# plot fitting data
dataplot <- data.frame(y=as.vector(data))
date <- as.Date(Date,"%m/%d/%Y")
dataplot[['Date']] <- date
all_summary <- vector(mode = "list")
all_residuals <- vector(mode = "list")
for(i in (length(qq)):(length(qq))){
  if(qq[i-1]==1){
    a<-as.numeric(qq[i-1])
  }else{
    a<-as.numeric(qq[i-1])+1
  }
  b<-as.numeric(qq[i])
  nlm3 <- fun_fit2(a, b, data[a:b])$fit
  all_summary[[i-1]]<-summary(nlm3)
  residuals <- resid(nlm3)
  all_residuals[[i-1]]<-residuals
  var(residuals)
  plot(residuals,main = '',xlab='',ylab='residuals')#,ylim = c(-0.15,0.15)
  mtext(paste("E"),side=3,adj=0,line = 1,font=3)  # ,Las=3
}
