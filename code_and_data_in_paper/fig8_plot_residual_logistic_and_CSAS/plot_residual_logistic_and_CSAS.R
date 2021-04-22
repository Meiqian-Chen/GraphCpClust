setwd("")
load("residual_list2.RData")
load("residual_list_C2.RData")
load("residual_list_origin.RData")
load("residual_list_C_origin.RData")
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n' )
  legend(...)
}

country_clust_test=c(expression(italic(paste('Cluster 2(IR,IT)'))),
                     expression(italic(paste('Cluster 4(JP)'))),
                     expression(italic(paste('Cluster 6(ES)'))),
                     expression(italic(paste("Cluster 8(TR)"))))
country_clust_select=c(2,4,6,8)
china_clust_test=c(expression(italic(paste('Cluster 1(HB)'))),
                   expression(italic(paste('Cluster 3(TW)'))),
                   expression(italic(paste('Cluster 4(HK)'))),
                   expression(italic(paste('Cluster 5(MO)'))))
china_clust_select=c(1,3,4,5)
plot_target=c('A','B','C','D','E','F','G','H')


cluster_C_origin=list()
cluster_C=list()
for(i in 1:length(country_clust_select)){
  cluster_C_test=c()
  for(item in residual_list_C[[country_clust_select[i]]]){
    cluster_C_test=c(cluster_C_test,item)
  }
  cluster_C[[i]]=cluster_C_test
  cluster_C_test_origin=c()
    for(item in residual_list_C_origin[[country_clust_select[i]]]){
    cluster_C_test_origin=c(cluster_C_test_origin,item)
  }
  cluster_C_origin[[i]]=cluster_C_test_origin
}

cluster_origin=list()
cluster=list()
for(i in 1:length(country_clust_select)){
  cluster_test=c()
  for(item in residual_list[[china_clust_select[i]]]){
    cluster_test=c(cluster_test,item)
  }
  cluster[[i]]=cluster_test
  cluster_test_origin=c()
  for(item in residual_list_origin[[china_clust_select[i]]]){
    cluster_test_origin=c(cluster_test_origin,item)
  }
  cluster_origin[[i]]=cluster_test_origin
}

x=1:142
col_wheel<-c("blue","red")
par(mfrow=c(2,4),mar=c(3,3,4,2),oma=c(0,0,2,0))
for(i in 1:length(cluster_origin)){
  plot(x,cluster_origin[[i]],type='b',pch=2,ylim=c(-0.3,0.4),main=list(china_clust_test[i],cex=1.4),xaxt="n",ylab='',xlab='',col=col_wheel[1],lwd=1.8)
  points(x,cluster[[i]],type="b",pch=3,xaxt="n",ylab='',xlab='',col=col_wheel[2],lwd=1.8)
  axis(1,c(1,32,63,92,123,142),as.character(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),cex.axis=1)
  mtext(paste(plot_target[i]),side=3,adj=0,line = 1,font=3)
}

for(i in 1:length(cluster_C_origin)){
  plot(x,cluster_C_origin[[i]],type='b',pch=2,ylim=c(-0.9,0.9),main=list(country_clust_test[i],cex=1.4),xaxt="n",ylab='',xlab='',col=col_wheel[1],lwd=1.8)
  points(x,cluster_C[[i]],type="b",pch=3,xaxt="n",ylab='',xlab='',col=col_wheel[2],lwd=1.8)
  axis(1,c(1,32,63,92,123,142),as.character(c("2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-04-20")),cex.axis=1)
  mtext(paste(plot_target[i+4]),side=3,adj=0,line = 1,font=3)
}
add_legend("top", legend=c('CSAS','Logistic model'),
           cex=1.5, bty='n', horiz=TRUE,lwd=2, pch=2:3, col=col_wheel[1:2],text.width=0.1)




