
#' Basic description
#'
#' @description  This function is used to generate simulation data introduced in Shi, Chen, Dong and Rao (2020)
#' @usage SimulationData()
#' @seealso Hpath_covid, Distance.data_covid, path.kruskal_covid
#' @export

SimulationData<-function(){
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
    return(RE)
  }
  
  
  Purity=function(L,N,K,re.c,c.t){
    purity=0
    for(ell in 1:L){
      s.c=NULL
      for(i in 1:N){
        s.c=c(s.c,length(intersect(re.c[[ell]],c.t[[i]])))
      }
      purity=purity+max(s.c)
    }
    purity=purity/K-abs(L-N)/max(c(L,N))
    return(purity)
  }
  
  
  
  ######simulation1
  set.seed(1) 
  #simulations for comparisons
  D=c(150);sigma=c(1:10)*0.1;loop.time=100;
  RE.g=RE.m=matrix(NA,length(D),length(sigma))
  for(r1 in 1:length(D))
    for(r2 in 1:length(sigma)){
      d=D[r1]
      s.e=sigma[r2]
      C1=C2=0;
      for(loop in 1:loop.time){
        Data=NULL
        for(i in 1:1) Data=rbind(Data,0+10*pnorm(-4+0.05*c(1:d))+rnorm(d,0,s.e))#
        for(i in 1:10) Data=rbind(Data,c(rep(0,50),20*pnorm(-3+0.03*c(1:(d-50)))+rnorm(d-50,0,s.e)))
        for(i in 1:15) Data=rbind(Data,c(rep(0,100),5*pnorm(-2+0.07*c(1:(d-100)))+rnorm(d-100,0,s.e)))
        sample.loc=sample(dim(Data)[1])
        Data=Data[sample.loc,]
        Lnorm = function(x) sqrt(sum(t(x)*x))# Euclidean distance, but reader can define different distance such as max(abs(x)) and sum(abs(x)) 
        
        data.trans=Distance.data(Data,1,dim(Data)[1],Lnorm)
        HD=data.trans$hd
        re.c=clusters(HD)
        re.c
        c.t=list()
        c.t[[1]]=1;c.t[[2]]=c(2:11);c.t[[3]]=c(12:26)
        
        L=dim(re.c);N=length(c.t);K=dim(Data)[1];
        re.c2=re.c
        for(ell in 1:dim(re.c))
          re.c2[[ell]]=sort(sample.loc[re.c[[ell]]])
        re.c2
        
        
        C1=C1+Purity(dim(re.c2),length(c.t),dim(Data)[1],re.c2,c.t)
        
        #sort(sample.loc[re.c[[1]]])
        #sort(sample.loc[re.c[[2]]])
        #sort(sample.loc[re.c[[3]]])
        # EM when d decreases to 150, sd 0.2 increases to 1 and the no. of first cluster 20 decreases to 1, it failed
        #https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html
        em<-Mclust(Data, G = NULL, modelNames = NULL, prior = NULL, 
                   control = emControl(), initialization = NULL, 
                   warn = mclust.options("warn"))
        
        L=length(unique(em$classification))
        
        EM=list()
        for(ell in 1:L)
          EM[[ell]]=sort(sample.loc[which(em$classification==ell)])
        EM
        #sort(sample.loc[which(em$classification==1)])
        #sort(sample.loc[which(em$classification==2)])
        #sort(sample.loc[which(em$classification==3)])
        
        C2=C2+Purity(length(unique(em$classification)),length(c.t),dim(Data)[1],EM,c.t)
      }
      
      
      RE.g[r1,r2]=C1/loop.time;RE.m[r1,r2]=C2/loop.time
    }
  
  saveRDS(RE.g, "REg.rds")
  saveRDS(RE.m, "REm.rds") 
  ##########end of simulation1
  
  
  ######simulation2
  set.seed(2) 
  #simulations for comparisons
  D=c(210);sigma=c(1:10)*0.1;loop.time=100;
  RE.g=RE.m=matrix(NA,length(D),length(sigma))
  for(r1 in 1:length(D))
    for(r2 in 1:length(sigma)){
      d=D[r1]
      s.e=sigma[r2]
      C1=C2=0; 
      for(loop in 1:loop.time){
        Data=NULL
        for(i in 1:5) Data=rbind(Data,0+10*pnorm(-4+0.05*c(1:d))+rnorm(d,0,s.e))#
        for(i in 1:20) Data=rbind(Data,c(rep(0,70),20*pnorm(-3+0.03*c(1:(d-70)))+rnorm(d-70,0,s.e)))
        for(i in 1:30) Data=rbind(Data,c(rep(0,140),5*pnorm(-2+0.07*c(1:(d-140)))+rnorm(d-140,0,s.e)))
        sample.loc=sample(dim(Data)[1])
        Data=Data[sample.loc,]
        Lnorm = function(x) sqrt(sum(t(x)*x))# Euclidean distance, but reader can define different distance such as max(abs(x)) and sum(abs(x)) 
        
        data.trans=Distance.data(Data,1,dim(Data)[1],Lnorm)
        HD=data.trans$hd
        re.c=clusters(HD)
        re.c
        c.t=list()
        c.t[[1]]=c(1:5);c.t[[2]]=c(6:25);c.t[[3]]=c(26:55)
        
        L=dim(re.c);N=length(c.t);K=dim(Data)[1];
        re.c2=re.c
        for(ell in 1:dim(re.c))
          re.c2[[ell]]=sort(sample.loc[re.c[[ell]]])
        re.c2
        
        
        
        C1=C1+Purity(dim(re.c2),length(c.t),dim(Data)[1],re.c2,c.t)
        
        #sort(sample.loc[re.c[[1]]])
        #sort(sample.loc[re.c[[2]]])
        #sort(sample.loc[re.c[[3]]])
        # EM when d decreases to 150, sd 0.2 increases to 1 and the no. of first cluster 20 decreases to 1, it failed
        #https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html
        em<-Mclust(Data, G = NULL, modelNames = NULL, prior = NULL, 
                   control = emControl(), initialization = NULL, 
                   warn = mclust.options("warn"))
        
        L=length(unique(em$classification))
        
        EM=list()
        for(ell in 1:L)
          EM[[ell]]=sort(sample.loc[which(em$classification==ell)])
        EM
        #sort(sample.loc[which(em$classification==1)])
        #sort(sample.loc[which(em$classification==2)])
        #sort(sample.loc[which(em$classification==3)])
        
        C2=C2+Purity(length(unique(em$classification)),length(c.t),dim(Data)[1],EM,c.t)
      }
      
      
      RE.g[r1,r2]=C1/loop.time;RE.m[r1,r2]=C2/loop.time
    }

  
  saveRDS(RE.g, "REg2.rds")
  saveRDS(RE.m, "REm2.rds") 
  ##########end of simulation2
  
  ######simulation3
  set.seed(3) 
  #simulations for comparisons
  D=c(300);sigma=c(1:10)*0.1;loop.time=100;
  RE.g=RE.m=matrix(NA,length(D),length(sigma))
  for(r1 in 1:length(D))
    for(r2 in 1:length(sigma)){
      d=D[r1]
      s.e=sigma[r2]
      C1=C2=0;
      for(loop in 1:loop.time){
        Data=NULL
        for(i in 1:20) Data=rbind(Data,0+10*pnorm(-4+0.05*c(1:d))+rnorm(d,0,s.e))#
        for(i in 1:100) Data=rbind(Data,c(rep(0,100),20*pnorm(-3+0.03*c(1:(d-100)))+rnorm(d-100,0,s.e)))
        for(i in 1:200) Data=rbind(Data,c(rep(0,200),5*pnorm(-2+0.07*c(1:(d-200)))+rnorm(d-200,0,s.e)))
        sample.loc=sample(dim(Data)[1])
        Data=Data[sample.loc,]
        Lnorm = function(x) sqrt(sum(t(x)*x))# Euclidean distance, but reader can define different distance such as max(abs(x)) and sum(abs(x)) 
        
        data.trans=Distance.data(Data,1,dim(Data)[1],Lnorm)
        HD=data.trans$hd
        re.c=clusters(HD)
        re.c
        c.t=list()
        c.t[[1]]=c(1:20);c.t[[2]]=c(21:120);c.t[[3]]=c(121:320)
        
        L=dim(re.c);N=length(c.t);K=dim(Data)[1];
        re.c2=re.c
        for(ell in 1:dim(re.c))
          re.c2[[ell]]=sort(sample.loc[re.c[[ell]]])
        re.c2
        
        
        C1=C1+Purity(dim(re.c2),length(c.t),dim(Data)[1],re.c2,c.t)
        
        #sort(sample.loc[re.c[[1]]])
        #sort(sample.loc[re.c[[2]]])
        #sort(sample.loc[re.c[[3]]])
        # EM when d decreases to 150, sd 0.2 increases to 1 and the no. of first cluster 20 decreases to 1, it failed
        #https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html
        em<-Mclust(Data, G = NULL, modelNames = NULL, prior = NULL, 
                   control = emControl(), initialization = NULL, 
                   warn = mclust.options("warn"))
        
        L=length(unique(em$classification))
        
        EM=list()
        for(ell in 1:L)
          EM[[ell]]=sort(sample.loc[which(em$classification==ell)])
        EM
        #sort(sample.loc[which(em$classification==1)])
        #sort(sample.loc[which(em$classification==2)])
        #sort(sample.loc[which(em$classification==3)])
        
        C2=C2+Purity(length(unique(em$classification)),length(c.t),dim(Data)[1],EM,c.t)
      }
      
      
      RE.g[r1,r2]=C1/loop.time;RE.m[r1,r2]=C2/loop.time
    }
  
  
  saveRDS(RE.g, "REg3.rds")
  saveRDS(RE.m, "REm3.rds") 
  ##########end of simulation3
}



