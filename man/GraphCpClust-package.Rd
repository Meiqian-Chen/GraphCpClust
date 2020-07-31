\name{GraphCpClust-package}
\alias{GraphCpClust-package}
\alias{GraphCpClust}
\docType{package}
\title{
\packageTitle{GraphCpClust}
}
\description{
\packageDescription{GraphCpClust}
}
\details{
}
\author{
\strong{Maintainer}: Meiqian Chen \email{mqchen@stu.scu.edu.cn}

\itemize{
  \item Meiqian Chen
  \item Xiaoping Shi
  \item Hao Li
}}


\references{

SHI, X.P., WU, Y.H. & RAO, C.R. (2017). Consistent and powerful graph-based change-point test for high-dimensional data. Proc Natl Acad Sci 114, 3873-8.

SHI, X.P., WU, Y.H. & RAO, C.R. (2018). Consistent and powerful non-Euclidean graph-based change-point test with applications to segmenting random interfered video data. Proc Natl Acad Sci 115, 5914-5919.

SHI, X.P., CHEN, M.Q., DONG, Y.C. & RAO, C.R. (2020). Exploring the space-time pattern of log-transformed infectious count of COVID-19: a clustering-segmented autoregressive sigmoid model.

\strong{Data support:}

Lihoreau M, Chittka L, Raine NE (2016) Monitoring flower visitation networks and interactions between pairs of bumble bees in a large outdoor flight cage. PLoS One 11:e0150844.

Ciccone M, Gutwein M, Gunsalus KC, Geiger D (2014) Label free cell-tracking and division detection based on 2D time-lapse images for lineage analysis of early embryo development. Comput Biol Med 51:24-34. http://celltracking.bio.nyu.edu/
}


\keyword{cluster  changepoint}
\section{Explanation}{
library(GraphCpClust)\cr\cr
help(CpBeeCell) #View the change points estimation method in Shi, Wu and Rao (2017) and Shi, Wu and Rao (2018).\cr\cr
example("CpBeeCell") # Run the example of data analysis in Shi, Wu and Rao (2017) and Shi, Wu and Rao (2018).\cr\cr
help(PathClust) # The clustering method of PathClust in Shi, Chen, Dong and Rao (2020)\cr\cr
example("PathClust") #Some examples of the PathClust clustering method that include the clustering results in Shi, Chen, Dong and Rao (2020)\cr\cr
SimulationPlot() #the Comparisons of PathClust clustering method and model-based clustering method.\cr\cr
example("CpCOVID") #The data example of change points estimation in Shi, Chen, Dong and Rao (2020).\cr\cr
}