#' Basic description
#'
#' @description  This function returns the clustering result and the plot of clusters based on the graph-based clustering method that introduced in Shi, Chen, Dong and Rao (2020).
#' @usage PathClust(data)
#' @param data The input data must be a matrix, it's each column represents a cluster variable and this matrix must have names for rows and columns
#' @return the clustering result
#' @import igraph
#' @examples
#' 
#'data(origin_data,package="GraphCpClust")
#' #Example 1
#'d=150
#'s.e=0.5
#'Data=NULL
#'for(i in 1:1) Data=rbind(Data,0+10*pnorm(-4+0.05*c(1:d))+rnorm(d,0,s.e))#
#'for(i in 1:10) Data=rbind(Data,c(rep(0,50),20*pnorm(-3+0.03*c(1:(d-50)))+rnorm(d-50,0,s.e)))
#'for(i in 1:15) Data=rbind(Data,c(rep(0,100),5*pnorm(-2+0.07*c(1:(d-100)))+rnorm(d-100,0,s.e)))
#'Data<-t(Data)
#'simulation_data = Data
#'colnames(simulation_data)<-c('v1','v2','v3','v4','v5','v6','v7','v8','v9','v10',
#'                             'v11','v12','v13','v14','v15','v16','v17','v18','v19','v20',
#'                            'v21','v22','v23','v24','v25','v26')

#'PathClust(simulation_data)
#'
#' #Example 2
#' # The clustering example in Shi, Chen, Dong and Rao (2020)
#'PathClust(Data.new)
#'PathClust(Data.new_country)
#' @seealso Hpath_covid, Distance.data_covid, path.kruskal_covid
#' @export
#' @import ggplot2 igraph

PathClust=function(data){
  
  clusters=function(HD){
    n=dim(HD)[1]
    SHP0=Hpath_covid(1,n,HD)
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
    k=0
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
      if(length(dist.path.c)==0){
        k = 1
        break
      }
      # BIC.all[ell+1]=log(var(dist.path.c))+2*N*log(T)/T
      if(length(dist.path.c)==1){
        BIC.all[ell+1]=2*N*log(T)/T
      }else{
        BIC.all[ell+1]=log(var(dist.path.c))+2*N*log(T)/T
      }
      # BIC.all[ell+1]=log(var(dist.path.c))+2*N*log(T)/T
      if(BIC.all[ell+1]>BIC.all[ell]) break;
    }
    if(k==0){
      RE=Clusters[[ell]]
    }else{
      cat('now is one cluster',is.null(dim(Clusters[[1]])), '\n')
      RE=Clusters[[1]]
    }
    # RE=Clusters[[ell]]
    return(list(path=path,groups=RE))
  }

  PLOT=function(theta1,theta2,delta,COL,NAME,zoom.in,zoom.out,SHIFT){ 
    xy1=c(radius*cos(theta1),radius*sin(theta1))
    xy2=c(radius*cos(theta2),radius*sin(theta2))
    newcenter=c((xy1[1]+xy2[1])/2,(xy1[2]+xy2[2])/2)
    dist.xy=sqrt(sum((xy1-xy2)^2))
    x.center=seq(center[1]-dist.xy/2,center[1]+dist.xy/2,0.01)
    cat(x.center,delta,'\n')
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

  IOS_name = colnames(data)
  Date = rownames(data)
  n=length(Date)
  Data_confirmed=t(data)
  Lnorm = function(x) sqrt(sum(t(x)*x))# Euclidean distance, but reader can define different distance such as max(abs(x)) and sum(abs(x)) 
  data.trans=Distance.data_covid(Data_confirmed,1,dim(Data_confirmed)[1],Lnorm)
  HD=data.trans$hd
  re.c=clusters(HD)
  cluster_res = re.c$groups
  
  # save cluster result data
  cluster_res_ios_name = list()
  if(is.null(dim(re.c$groups))){
    cluster = rep(0,n)
    for(item in re.c$groups){
      cluster = cluster + data[1:n,IOS_name[item]]
    }
    data_clust <- matrix(cluster,nrow=1,ncol=n,byrow=TRUE)
    data_clust <- t(data_clust)
    colnames(data_clust)=c('cluster1')
    rownames(data_clust)=Date
    cluster_res_ios_name[[1]] = IOS_name

  }else{
    data_clust_cell = c()
    data_clust_name = c()
    k = 1
    for(item in re.c$groups){
      cat(item,'\n')
      data_clust_name = c(data_clust_name,paste(c("cluster",as.character(k)),collapse=""))
      k = k+1
      cluster = rep(0,n)
      for(i in item){
        cluster = cluster + data[1:n,IOS_name[i]]
      }
      data_clust_cell <- c(data_clust_cell,cluster)
    }
    data_clust <- matrix(data_clust_cell,nrow=dim(re.c$groups),ncol=n,byrow=TRUE)
    data_clust <- t(data_clust)
    colnames(data_clust)=data_clust_name
    rownames(data_clust)=Date
    # save(data_clust,data_clust_name,cluster_res, file='data_clust.RData')
    for(i in 1:dim(cluster_res)){
      clusteri = c()
      for(item in cluster_res[[i]]){
        clusteri = c(clusteri, IOS_name[item])
      }
      cat(clusteri,'\n')
      cluster_res_ios_name[[i]] = clusteri
    }
  }

  
  par(mar=c(3,2,7,2)+0.1,xpd=TRUE)
  if(is.null(dim(re.c$groups))){
    n.c=1
  }else{
    n.c=length(re.c$groups)
  }
  n.p=dim(Data_confirmed)[1]
  n.c_1 = length(re.c$groups)
  
  center=c(0,0);radius=10;
  
  plot(0,xlim=c(-1.5*radius,2*radius),ylim=2*c(-radius,radius),
       xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  # plot of circle
  theta=seq(0,2*pi,0.01)
  for(i in 1:length(theta))
    points(radius*cos(theta[i]),radius*sin(theta[i]),type="l")
  
  re.c$path
  COL_origin = c('#FF0000', '#FF1493', '#800080', '#0000FF', '#00FFFF', '#00FF7F', '#FFA500', '#FFFF00', '#778899')
  if(n.c != 1){
    COL = c()
    for(i in 1:n.c){
      COL = c(COL, COL_origin[i])
    }
  }else{
    COL = c('#FF0000')
  }
  # COL=rainbow(length(re.c$group))
  # COL=rainbow(n.c)
  if(n.c == 1){
    L=NULL
    for(t1 in 1:n.c_1)   L=c(L,length(re.c$groups[[t1]]))
    L=unique(L)
    L=sort(L)
    PL=seq(0.05,0.15,length.out=length(L))
    ALL=matrix(NA,n.p,4)
    ALL[,1]=c(1:n.p)
    for(t1 in 1:n.c){
      groups_res = length(re.c$groups)
      for(t2 in 1:groups_res){  # groups_res
        temp=which(ALL[,1]==re.c$groups[[t2]])
        ALL[temp,2]=COL[t1]
        # temp.n=which(NAME[,1]==IOS_name[re.c$groups[[t1]][t2]])
        # ALL[temp,3]=as.vector(NAME[temp.n,3])
        ALL[temp,3]=IOS_name[temp]
        temp.L=which(L==n.c)
        ALL[temp,4]=PL[temp.L]
      }
    }  
  }else{
    L=NULL
    for(t1 in 1:n.c)   L=c(L,length(re.c$groups[[t1]]))
    L=unique(L)
    L=sort(L)
    PL=seq(0.05,0.15,length.out=length(L))
    ALL=matrix(NA,n.p,4)
    ALL[,1]=c(1:n.p)
    for(t1 in 1:n.c){
      groups_res = length(re.c$groups[[t1]])
      for(t2 in 1:groups_res){  # groups_res
        temp=which(ALL[,1]==re.c$groups[[t1]][t2])
        ALL[temp,2]=COL[t1]
        # temp.n=which(NAME[,1]==IOS_name[re.c$groups[[t1]][t2]])
        # ALL[temp,3]=as.vector(NAME[temp.n,3])
        ALL[temp,3]=IOS_name[temp]
        temp.L=which(L==length(re.c$groups[[t1]]))
        ALL[temp,4]=PL[temp.L]
      }
    }
  }
  
  
  delta.theta=(length(theta)-n.p*10)/n.p
  
  center=c(0,0);radius=10;
  plot(0,xlim=c(-3.5*radius,4*radius),ylim=c(-1.8*radius,1.4*radius),
       xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  #plot of circle
  theta=seq(0,2*pi,0.01)
  for(i in 1:length(theta))
    points(radius*cos(theta[i])-radius,radius*sin(theta[i]),type="l")
  text(0,0,"clustering")  # 0-2*radius
  
  s.theta=1
  for(k in 1:n.p){
    temp1=which(ALL[,1]==re.c$path[k])
    # cat(theta[s.theta], radius,'\n')
    PLOT(theta[s.theta],theta[s.theta+9],as.numeric(ALL[temp1,4]),
         ALL[temp1,2],ALL[temp1,3],0.6,0.9,0)  # the last para 2*radius
    theta1=theta[s.theta+9]
    s.theta=s.theta+9+delta.theta
    theta2=theta[s.theta]
    seg.theta=seq(theta1,theta2,0.01)
    xy.seg=NULL
    for(ell in 1:length(seg.theta))
      xy.seg=cbind(xy.seg,matrix(c(radius*cos(seg.theta[ell]),radius*sin(seg.theta[ell])),2,1))
    points(xy.seg[1,],xy.seg[2,],type="l")  # xy.seg[1,]-2*radius
  }
  # center=c(0,0);radius=10;
  return(list(data_clust=data_clust,cluster_res=cluster_res_ios_name)) 
}