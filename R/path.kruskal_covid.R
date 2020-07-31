#' Basic description
#'
#' @description Calculates the shortest Hamiltonian path based on the sorted edge weights and the nodes
#' @usage path.kruskal_covid(nodes,edge_cost)
#' @param nodes sequence of nodes 1,...,n from the graph which is based on the high-dimensional data that is provided by the reader
#' @param edge_cost sorted edge weights
#' @return the shortest Hamiltonian path
#' @examples
#' # path.kruskal is used in the Hpath function as follows:
#' Hpath_covid=function(n1,n2,mat){
#' n0=n2-n1+1
#' edge.cost=matrix(NA,nrow=n0*(n0-1)/2,ncol=3)
#'
#' temp=1;
#'
#' for(i in n1:(n2-1))
#' for(j in (i+1):(n2))
#' {
#' edge.cost[temp,3]=mat[i,j];edge.cost[temp,1]=i-n1+1;edge.cost[temp,2]=j-n1+1;temp=temp+1;}
#'
#' edge.cost=edge.cost[sort.list(edge.cost[,3]), ]
#'
#' #based on path
#'
#' opt.path=path.kruskal_covid(c(1:n0),edge.cost)
#'
#' return(opt.path$edge.cost[,1:2])
#' }
#' @seealso Hpath_covid
#' @export
path.kruskal_covid=function (nodes, edge_cost)
{   n0=length(nodes)
cost=0
components <- matrix(c(nodes, nodes), ncol = 2)
shp <- matrix(ncol = 3)[-1, ]
edge_cost=rbind(matrix(ncol = 3)[-1, ],edge_cost)
i <- 1
degrees=rep(0,length(nodes))
while (nrow(shp) < n0 - 1) {
  min.shp <- edge_cost[i, ]
  iComp <- components[components[, 1] == min.shp[1], 2]
  jComp <- components[components[, 1] == min.shp[2], 2]
  if (iComp != jComp && max(c(degrees[min.shp[1]],degrees[min.shp[2]]))<2) {
    shp <- rbind(shp, min.shp)
    cost=cost+edge_cost[i, 3]
    components[components[, 2] == jComp, 2] <- iComp
    degrees[min.shp[1]]=degrees[min.shp[1]]+1;
    degrees[min.shp[2]]=degrees[min.shp[2]]+1;
  }
  i <- i + 1
}

return(list(edge.cost=shp))
}
