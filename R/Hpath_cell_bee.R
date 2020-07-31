#' Basic description
#'
#' @description Applies the path.kruskal function based on the nodes and edge.cost (sorts the weights from minimum to maximum). Given the starting node, ending node, and the distance matrix, this function returns the list of nodes of each edge from the shortest Hamiltonian path. We have the Hamiltonian path from path.kruskal
#' @usage Hpath_cell_bee(n1,n2,mat)
#' @param n1 starting node
#' @param n2 ending node
#' @param mat distance matrix (distance type is determined by the reader)
#' @return list of nodes of each edge from the shortest Hamiltonian path
#' @examples
#' data(Data_cell,package="GraphCpClust")
#' Lnorm = function(x) sqrt(sum(t(x)*x))
#' data.trans=Distance.data_cell_bee(Data_cell,1,285,Lnorm)
#' HD=data.trans$hd # distance matrix
#' Hpath_cell_bee(1,285,HD)
#' @seealso Distance.data_cell_bee, path.kruskal_cell_bee
#' @export
Hpath_cell_bee=function(n1,n2,mat){
	n0=n2-n1+1
	edge.cost=matrix(NA,nrow=n0*(n0-1)/2,ncol=3)

	temp=1;
	for(i in n1:(n2-1))
	  for(j in (i+1):(n2))
	    {
	edge.cost[temp,3]=mat[i,j];edge.cost[temp,1]=i-n1+1;edge.cost[temp,2]=j-n1+1;temp=temp+1;}

	edge.cost=edge.cost[sort.list(edge.cost[,3]), ]
	#based on path
	 
	opt.path=path.kruskal_cell_bee(c(1:n0),edge.cost)
	return(opt.path$edge.cost[,1:2])
}