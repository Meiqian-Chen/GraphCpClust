#' Basic description
#'
#' @description Given the starting point, ending point and weight function, and this function returns the distance matrix
#' @usage Distance.data_covid(Data,n1,n2,fun.dist)
#' @param Data input data
#' @param n1 starting point/node
#' @param n2 ending point/node
#' @param fun.dist weight function (Euclidean, lnorm, or may be defined by the reader)
#' @return distance matrix - depends upon weight function what kind of distance it returns
#' @examples # See Hpath_covid example.
#' @seealso Hpath_covid
#' @export
Distance.data_covid=function(Data,n1,n2,fun.dist){
  n0=n2-n1+1
  hd=matrix(NA,nrow=n0,ncol=n0)
  for(i in n1:(n2-1))
    for(j in (i+1):(n2))
      hd[i-n1+1,j-n1+1]= fun.dist( Data[i,]-Data[j,])

  return(list(hd=hd))
}

