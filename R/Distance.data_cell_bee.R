
#' Basic description
#'
#' @description Given the starting point, ending point and weight function, and this function returns the distance matrix
#' @usage Distance.data_cell_bee(Data1,n1,n2,fun.dist)
#' @param Data1 input data
#' @param n1 starting point/node
#' @param n2 ending point/node
#' @param fun.dist weight function (Euclidean, lnorm, or may be defined by the reader)
#' @return distance matrix - depends upon weight function what kind of distance it returns
#' @examples # See Hpath_cell_bee example.
#' @seealso Hpath_cell_bee
#' @export

Distance.data_cell_bee=function(Data1,n1,n2,fun.dist){
n0=n2-n1+1
hd=matrix(NA,nrow=n0,ncol=n0)
for(i in n1:(n2-1))
  for(j in (i+1):(n2)) 
    #hd[i-n1+1,j-n1+1]= fun.dist( Data1[[i]]-Data1[[j]])+fun.dist( Data2[[i]]-Data2[[j]])+fun.dist( Data3[[i]]-Data3[[j]])+fun.dist( Data4[[i]]-Data4[[j]])+fun.dist( Data5[[i]]-Data5[[j]])+fun.dist( Data6[[i]]-Data6[[j]])
 	hd[i-n1+1,j-n1+1]= fun.dist(Data1[[i]]-Data1[[j]])  

return(list(hd=hd))
}