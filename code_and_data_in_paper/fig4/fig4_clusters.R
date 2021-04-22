library(readr)
library(ggplot2)
library(igraph)
setwd("")
data <- read_csv("data/4.20.csv")
data_country <- read_csv("data/world_0420_33country.csv")


Data=as.matrix(data)
PC=unique(Data[,5])
Date=unique(Data[,1])
n=length(Date)
# The confirmed matrix is constructed by province
Data.new=matrix(0,n,length(PC))
colnames(Data.new)=c(PC)
rownames(Data.new)=Date
Data.new[,1]=c(1:n)
for(i in 1:n)
  for(j in 1:length(PC))
    if(length(intersect(which(Data[,1]==Date[i]),which(Data[,5]==PC[j])))!=0)
      Data.new[i,j]=as.numeric(Data[intersect(which(Data[,1]==Date[i]),which(Data[,5]==PC[j])),6])

# The death matrix is constructed by province
Data.new_dead=matrix(0,n,length(PC))
colnames(Data.new_dead)=c(PC)
rownames(Data.new_dead)=Date
Data.new_dead[,1]=c(1:n)
for(i in 1:n)
  for(j in 1:length(PC))
    if(length(intersect(which(Data[,1]==Date[i]),which(Data[,5]==PC[j])))!=0)
      Data.new_dead[i,j]=as.numeric(Data[intersect(which(Data[,1]==Date[i]),which(Data[,5]==PC[j])),9])

Data_country=as.matrix(data_country)
PC_country=unique(Data_country[,3])
Date=unique(Data_country[,1])
n=length(Date)
# The confirmed number matrix was constructed by country
Data.new_country=matrix(0,n,length(PC_country))
colnames(Data.new_country)=c(PC_country)
rownames(Data.new_country)=Date
Data.new_country[,1]=c(1:n)
for(i in 1:n)
  for(j in 1:length(PC_country))
    if(length(intersect(which(Data_country[,1]==Date[i]),which(Data_country[,3]==PC_country[j])))!=0)
      Data.new_country[i,j]=as.numeric(Data_country[intersect(which(Data_country[,1]==Date[i]),which(Data_country[,3]==PC_country[j])),4])
# There was a default in the January 5th section of the country confirmed
for(j in 2:length(PC_country)){
  for(i in 1:n){
    if(i==66){
      Data.new_country[i,j]=Data.new_country[65,j]
    }
    if(i==67){
      Data.new_country[i,j]=Data.new_country[68,j]
    }
  }
}
# Construct the death matrix by country
Data.new_country_dead=matrix(0,n,length(PC_country))
colnames(Data.new_country_dead)=c(PC_country)
rownames(Data.new_country_dead)=Date
Data.new_country_dead[,1]=c(1:n)
for(i in 1:n)
  for(j in 1:length(PC_country))
    if(length(intersect(which(Data_country[,1]==Date[i]),which(Data_country[,3]==PC_country[j])))!=0)
      Data.new_country_dead[i,j]=as.numeric(Data_country[intersect(which(Data_country[,1]==Date[i]),which(Data_country[,3]==PC_country[j])),7])
for(j in 2:length(PC_country)){
  for(i in 1:n){
    if(i==66){
      Data.new_country_dead[i,j]=Data.new_country_dead[65,j]
    }
    if(i==67){
      Data.new_country_dead[i,j]=Data.new_country_dead[68,j]
    }
  }
}


#Structural classification matrix
# Data_confirmed=t(log(1+Data.new))#Data.new_country
# Data_confirmed=t(log(1+Data.new))
source("Distance.data.R")
source("Hpath.R")
source("path.kruskal.R")

# Data_confirmed=t(log(1+Data.new))#Data.new_country
# Lnorm = function(x) sqrt(sum(t(x)*x))# Euclidean distance, but reader can define different distance such as max(abs(x)) and sum(abs(x)) 
# data.trans=Distance.data(Data_confirmed,1,dim(Data_confirmed)[1],Lnorm)# 34 province,30 country
# HD=data.trans$hd




clusters=function(HD){
  n=dim(HD)[1]
  SHP0=Hpath(1,n,HD)
  karate <-graph.edgelist(SHP0, directed=F)
  ii=sort(which(degree(karate)==1))[1]
  ie=ini=sort(which(degree(karate)==1))[2]
  path=all_simple_paths(karate,from=ii,to=ie)
  path=as.vector(path[[1]])
  dist.path=NULL
  for(j in 1:(n-1))
    dist.path=c(dist.path,HD[sort(path[j:(j+1)])[1],sort(path[j:(j+1)])[2]])

  #boxplot(dist.path)
  #outliers <- boxplot(dist.path, plot=FALSE)$out
  #dist.path
  #path
  #outliers
  #Out=which(dist.path>=min(outliers))
  T=n-1
  Clusters=list()
  Clusters[[1]]=c(1:n)
  BIC.all=NULL
  BIC.all[1]=log(var(dist.path))+2*log(T)/T
  #theta=min(outliers)
  ALL=sort(dist.path,TRUE)
  for(ell in 1:length(ALL)){
    theta=ALL[ell]
    #for 2.5, see equation (6) in
    #https://wis.kuleuven.be/stat/robust/papers/2011/rousseeuwhubert-robuststatisticsforoutlierdetectio.pdf

    adjm=matrix(0,n,n)
    for(j in 1:(n-1))
      if(HD[SHP0[j,1],SHP0[j,2]]<theta)
        adjm[SHP0[j,1],SHP0[j,2]]=adjm[SHP0[j,2],SHP0[j,1]]=1

    g1 <- graph_from_adjacency_matrix( adjm )
    Clusters[[ell+1]]=groups(components(g1))
    N=dim(Clusters[[ell+1]])
    dist.path.c=dist.path[dist.path<theta]
    BIC.all[ell+1]=log(var(dist.path.c))+2*N*log(T)/T
    if(BIC.all[ell+1]>BIC.all[ell]) break;
  }

  RE=Clusters[[ell]]
  return(list(path=path,groups=RE))
}


PLOT=function(theta1,theta2,delta,COL,NAME,zoom.in,zoom.out,SHIFT){ 
  xy1=c(radius*cos(theta1),radius*sin(theta1))
  xy2=c(radius*cos(theta2),radius*sin(theta2))
  newcenter=c((xy1[1]+xy2[1])/2,(xy1[2]+xy2[2])/2)
  dist.xy=sqrt(sum((xy1-xy2)^2))
  x.center=seq(center[1]-dist.xy/2,center[1]+dist.xy/2,0.01)
  y.center=dnorm(x.center,mean=0,sd=delta)
  v.center=c(x.center[which.max(y.center)],max(y.center)+1)###
  shift=newcenter-center
  x1=0;y1=1;x2=newcenter[1]-center[1];y2=newcenter[2]-center[2]
  dot = x1*x2 + y1*y2      # dot product between [x1, y1] and [x2, y2]
  det = x1*y2 - y1*x2      # determinant
  angle = 2*pi-atan2(det, dot) # atan2(y, x) or atan2(sin, cos)
  #theta <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
  rotate.theta=angle
  #atan((xy2[2]-xy1[2])/(xy2[1]-xy1[1]))
  #theta2-theta1
  xy.center=matrix(c(cos(rotate.theta),sin(rotate.theta),
                     -sin((rotate.theta)),cos(rotate.theta)),2,2,byrow=T)
  xy.shift=NULL
  for(j in 1:length(x.center)){
    xy.rotate=xy.center%*%matrix(c(x.center[j],y.center[j]),
                                 2,1)
    xy.shift=cbind(xy.shift,xy.rotate+shift)
  }
  v.shift.center=xy.center%*%matrix(v.center,2,1)+shift
  points(xy.shift[1,]-SHIFT,xy.shift[2,],type = "l",col="black")
  #points(v.shift.center[1,],v.shift.center[2,],pch=16,col=COL,cex=1)
  rr=sqrt(sum((v.shift.center-center)^2))
  x1=1;y1=0;x2=v.shift.center[1,];y2=v.shift.center[2,]
  dot = x1*x2 + y1*y2      # dot product between [x1, y1] and [x2, y2]
  det = x1*y2 - y1*x2      # determinant
  angle2 = atan2(det, dot) # atan2(y, x) or atan2(sin, cos)
  points((rr-zoom.in)*cos(angle2)-SHIFT,(rr-zoom.in)*sin(angle2),pch=16,col=COL,cex=1.5)
  text((rr+zoom.out)*cos(angle2)-SHIFT,(rr+zoom.out)*sin(angle2),NAME)
  #lines(c(center[1],newcenter[1]),c(center[2],newcenter[2]),col="green")
}

par(mar=c(3,2,3,2)+0.1,xpd=TRUE)
#china
Data_confirmed=t(log(1+Data.new))#Data.new_country
Lnorm = function(x) sqrt(sum(t(x)*x))# Euclidean distance, but reader can define different distance such as max(abs(x)) and sum(abs(x)) 
data.trans=Distance.data(Data_confirmed,1,dim(Data_confirmed)[1],Lnorm)# 34 province,30 country
HD=data.trans$hd
re.c=clusters(HD)
re.c$groups
PC
PC_country

n.c=length(re.c$group)
n.p=dim(Data_confirmed)[1]

center=c(0,0);radius=10;

plot(0,xlim=c(-1.5*radius,2*radius),ylim=2*c(-radius,radius),
     xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
#plot of circle
theta=seq(0,2*pi,0.01)
for(i in 1:length(theta))
  points(radius*cos(theta[i]),radius*sin(theta[i]),type="l")

re.c$path
COL=rainbow(n.c)
NAME <- read.csv("data/legend.csv",header = FALSE)

L=NULL
for(t1 in 1:n.c)   L=c(L,length(re.c$groups[[t1]]))
L=unique(L)
L=sort(L)
PL=seq(0.05,0.15,length.out=length(L))
ALL=matrix(NA,n.p,4)
ALL[,1]=c(1:n.p)
for(t1 in 1:n.c)
  for(t2 in 1:length(re.c$groups[[t1]])){
    temp=which(ALL[,1]==re.c$groups[[t1]][t2])
    ALL[temp,2]=COL[t1]
    temp.n=which(NAME[,1]==PC[re.c$groups[[t1]][t2]])
    ALL[temp,3]=as.vector(NAME[temp.n,3])
    temp.L=which(L==length(re.c$groups[[t1]]))
    ALL[temp,4]=PL[temp.L]
  }

delta.theta=(length(theta)-n.p*10)/n.p

center=c(0,0);radius=10;
plot(0,xlim=c(-3.5*radius,4*radius),ylim=c(-1.8*radius,1.4*radius),
     xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
#plot of circle
theta=seq(0,2*pi,0.01)
for(i in 1:length(theta))
  points(radius*cos(theta[i])-radius,radius*sin(theta[i]),type="l")
text(0-2*radius,0,"COVID-19 in China")

s.theta=1
for(k in 1:n.p){
  temp1=which(ALL[,1]==re.c$path[k])
  PLOT(theta[s.theta],theta[s.theta+9],as.numeric(ALL[temp1,4]),
       ALL[temp1,2],ALL[temp1,3],0.6,0.9,2*radius)
  theta1=theta[s.theta+9]
  s.theta=s.theta+9+delta.theta
  theta2=theta[s.theta]
  seg.theta=seq(theta1,theta2,0.01)
  xy.seg=NULL
  for(ell in 1:length(seg.theta))
    xy.seg=cbind(xy.seg,matrix(c(radius*cos(seg.theta[ell]),radius*sin(seg.theta[ell])),2,1))
  points(xy.seg[1,]-2*radius,xy.seg[2,],type="l")
}
center=c(0,0);radius=10;


#world
# Data were taken for 41 days from the date of diagnosis for each country, and those that did not were excluded
Data.new_test_41=matrix(0,41,33)
PC_new <- c()
rownames(Data.new_test_41)=seq(1,41,1)
# Data.new_test_41[,1]=c(1:41)
for(j in 1:length(PC_country)){
  for(k in 1:(n-40)){
    if(Data.new_country[k,j]!=0){
      for(h in k:(k+40)){
        Data.new_test_41[h-k+1,j]=Data.new_country[h,j]
      }
      PC_new <- c(PC_new,PC_country[j])
      break
    }
  }
}
colnames(Data.new_test_41)=c(PC_new)
Data.new_country_log=t(log(1+Data.new_country))#Data.new_country
Lnorm = function(x) sqrt(sum(t(x)*x))# Euclidean distance, but reader can define different distance such as max(abs(x)) and sum(abs(x)) 
data.trans=Distance.data(Data.new_country_log,1,dim(Data.new_country_log)[1],Lnorm)# 34 province,30 country
HD=data.trans$hd
re.c=clusters(HD)
re.c$groups
PC_country


n.c=n.c=length(re.c$groups)
n.p=dim(Data.new_country_log)[1]

re.c$path
COL=rainbow(n.c)
# NAME <- read.csv("legend.csv",header = FALSE)

L=NULL
for(t1 in 1:n.c)   L=c(L,length(re.c$groups[[t1]]))
L=unique(L)
L=sort(L)
PL=seq(0.05,0.15,length.out=length(L))
ALL=matrix(NA,n.p,4)
ALL[,1]=c(1:n.p)
for(t1 in 1:n.c)
  for(t2 in 1:length(re.c$groups[[t1]])){
    temp=which(ALL[,1]==re.c$groups[[t1]][t2])
    ALL[temp,2]=COL[t1]
    #temp.n=which(NAME[,1]==PC[re.c$groups[[t1]][t2]])
    ALL[temp,3]=PC_country[temp]
    temp.L=which(L==length(re.c$groups[[t1]]))
    ALL[temp,4]=PL[temp.L]
  }

theta=seq(0,2*pi,0.01)
for(i in 1:length(theta))
  points(radius*cos(theta[i])-radius,radius*sin(theta[i]),type="l")
text(0+2*radius,0,"COVID-19 in world")

s.theta=1
for(k in 1:n.p){
  temp1=which(ALL[,1]==re.c$path[k])
  PLOT(theta[s.theta],theta[s.theta+9],as.numeric(ALL[temp1,4]),
       ALL[temp1,2],ALL[temp1,3],0.6,0.9,-2*radius)
  theta1=theta[s.theta+9]
  s.theta=s.theta+9+delta.theta
  theta2=theta[s.theta]
  seg.theta=seq(theta1,theta2,0.01)
  xy.seg=NULL
  for(ell in 1:length(seg.theta))
    xy.seg=cbind(xy.seg,matrix(c(radius*cos(seg.theta[ell]),radius*sin(seg.theta[ell])),2,1))
  points(xy.seg[1,]+2*radius,xy.seg[2,],type="l")
}










