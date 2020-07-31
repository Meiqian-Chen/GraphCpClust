#' Basic description
#'
#' @description Given the true starting point, number of points, and the shortest Hamiltonian path, this function returns the change-point estimation based on the ratio cut.
#' @usage CpBeeCell(ni,n1,n2,re.path)
#' @param ni true starting point
#' @param n1 always 1
#' @param n2 number of points
#' @param re.path the shortest Hamiltonian path returned from the function Hpath
#' @return the change-point estimate based on the ratio cut
#' @examples
#' data(Data_cell,Data_bee_A1,Data_bee_B1,package="GraphCpClust")
#' # Example in Shi, Wu and Rao (2017)
#' Lnorm_cell = function(x) sqrt(sum(t(x)*x))# Euclidean distance, but reader can define different distance such as max(abs(x)) and sum(abs(x)) 
#' Transf_cell = function(x) sqrt(x) #sqrt(x) log(1+x) # for cell images, squareroot transformation or log transformation
#'
#' # Note: Have the data in an image folder in your working directory and in the .png package subfolder img. 
#' # for whole sequence of data 
#' data.trans=Distance.data_cell_bee(Data_cell,1,285,Lnorm_cell)
#' HD=data.trans$hd
#' SHP=CpBeeCell(1,1,285,Hpath_cell_bee(1,285,HD))
#'
#' # cp1 = first change-point estimate 
#' # cp1 should be equal to SHP$Ratio.cut which is ratio cut change-point estiamte
#' cp1=SHP$Ratio.cp
#'
#' # subsequence from from first node to cp1 (the first change-point estiamte)
#' # repeat procdure to the change-point estimate
#' data.trans=Distance.data_cell_bee(Data_cell,1,cp1,Lnorm_cell)
#' HD=data.trans$hd
#' SHP=CpBeeCell(1,1,cp1,Hpath_cell_bee(1,cp1,HD))
#'
#' SHP$Sn # statistic - value shows whether our change-point test is signficant or not based on a comparison to the critical value from a table.
#' cp2=SHP$Ratio.cp
#' 
#' # subsequence fro cp1+1 to last node
#' # repeat procdure to the change-point estimate
#' data.trans=Distance.data_cell_bee(Data_cell,cp1+1,285,Lnorm_cell)
#' HD=data.trans$hd
#' SHP=CpBeeCell(cp1+1,1,285-cp1,Hpath_cell_bee(1,285-cp1,HD)) # first number is the true initial starting point location after the first change-point estimate
#'
#' cp3=SHP$Ratio.cp
#' # input data must be a data list
#' data.trans=Distance.data_cell_bee(Data_cell,1,285,Lnorm_cell)
#' HD=data.trans$hd
#' SHP=CpBeeCell(1,1,285,Hpath_cell_bee(1,285,HD))
#' print(cp1)
#' print(cp2)
#' print(cp3)
#' # Example in Shi, Wu and Rao (2018)
#' #####################################################start
#' # fig4A
#' Lnorm_bee = function(x) abs(sum(x))  #sqrt(sum(t(x)*x))# sqrt(sum(t(x)*x))#sqrt(sum(t(x)*x)) max(abs(x)) #sum(abs(x))
#' Transf_bee = function(x) sqrt(x) #sqrt(x) log(1+x)
#'
#' nb=1;ne=49;
#' nl=ne-nb+1
#' theta=0#thresheld
#' data.trans=Distance.data_cell_bee(Data_bee_A1,1,ne,Lnorm_bee)
#' HD=data.trans$hd
#' SHP1=CpBeeCell(1,1,ne,Hpath_cell_bee(1,ne,HD))
#'
#' par(mar=c(4,5,1,1)+0.1,fig=c(0,1,0,1))
#' plot(x=c(1:(ne-1)),y=SHP1$Ratio.cut,ylab=expression(italic(C[t]^{'SHP(w*)'}/(t(N-t)))),xlab=expression(italic(t)), type="b", lwd=1.5,lty=3, pch=20) 
#' abline(v=which.min(SHP1$Ratio.cut[1:20]),lty=5)
#' abline(v=which.min(SHP1$Ratio.cut[20:ne])+19,lty=5)
#' #####################################################end
#'
#' #####################################################start
#' # fig4B
#' Lnorm_bee = function(x) abs(sum(x))  #sqrt(sum(t(x)*x))# sqrt(sum(t(x)*x))#sqrt(sum(t(x)*x)) max(abs(x)) #sum(abs(x))
#' Transf_bee = function(x) sqrt(x) #sqrt(x) log(1+x)
#' nb=5;ne=49;
#' nl=ne-nb+1
#' theta=0#thresheld
#'
#' data.trans=Distance.data_cell_bee(Data_bee_B1,1,45,Lnorm_bee)
#' HD=data.trans$hd
#' SHP2=CpBeeCell(1,1,45,Hpath_cell_bee(1,45,HD)) 
#'
#' par(mar=c(4,5,1,1)+0.1,fig=c(0,1,0,1) )
#' plot(x=c(1:(45-1)),y=SHP2$Ratio.cut,ylab=expression(italic(C[t]^{'SHP(w*)'}/(t(N-t)))),xlab=expression(italic(t)), type="b", lwd=1.5,lty=3, pch=20) 
#' abline(v=which.min(SHP2$Ratio.cut[1:45]),lty=5)
#' @seealso 
#' @export

CpBeeCell=function(ni,n1,n2,re.path){
  n0=n2-n1+1
  Run.t=rep(0,(n0-1))
  for(t in 1:((n0-1)))
   for(i in 1:dim(re.path)[1])
    if(re.path[i,1]<=t &&re.path[i,2]>t) Run.t[t]=Run.t[t]+1;
  Z.t=Z2.t=NULL
  for(i in 1:(length(Run.t))){
    Z.t[i]=-(Run.t[i]-2*i*(n0-i)/n0)/sqrt(2*i*(n0-i)*(2*i*(n0-i)-n0)/(n0^3-n0^2))
    Z2.t[i]=(Z.t[i])^2*n0^2/(i*(n0-i))
    }
  D1=D2=NULL
  for(t in 1:(n0-1)) {
    D1[t]=min(c(t, n0-t));
    D2[t]=t*(n0-t);
  }
  return(list(Sn=sum(Z.t^2)/(n0-1),Zsum2=sum(abs(Z.t))/(n0-1),Zmax=max(Z.t),Run=Run.t,Zt=Z.t,Cheeger.cut=Run.t/D1,Ratio.cut=Run.t/D2,Cheeger.cp=which.min(Run.t/D1)+ni-1,Ratio.cp=which.min(Run.t/D2)+ni-1))
}
