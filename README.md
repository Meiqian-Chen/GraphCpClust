
# GraphCpClust

<!-- badges: start -->
<!-- badges: end -->

The goal of GraphCpClust is to implement the change point estimation method and graph-based clustering method introduced in the three papers, and provide related data analysis.

    SHI, X.P., WU, Y.H. & RAO, C.R. (2017). Consistent and powerful graph-based change-point test for high-dimensional data. Proc Natl Acad Sci 114, 3873-8.

    SHI, X.P., WU, Y.H. & RAO, C.R. (2018). Consistent and powerful non-Euclidean graph-based change-point test with applications to segmenting random interfered video data. Proc Natl Acad Sci 115, 5914-5919.

    SHI, X.P., CHEN, M.Q., DONG, Y.C. & RAO, C.R. (2020). Exploring the space-time pattern of log-transformed infectious count of COVID-19: a clustering-segmented autoregressive sigmoid model.

# Package authors

Maintainer: Meiqian Chen (mqchen@stu.scu.edu.cn)

    Meiqian Chen
    
    Xiaoping shi
    
    Hao Li

# Data support

For the cell image data used in Shi, Wu and Rao (2017), Xiaoping Shi was authorized by Dr. M. Cicconet (Ciccone et al., 2014). And for the video data used in Shi, Wu and Rao (2018), Xiaoping Shi was authorized by Dr. Mathieu Lihoreau (Lihoreau et al., 2014).

Cicconet M, Gutwein M, Gunsalus KC, Geiger D (2014) Label free cell-tracking and division detection based on 2D time-lapse images for lineage analysis of early embryo development. Comput Biol Med 51:24-34. http://celltracking.bio.nyu.edu/

Lihoreau M, Chittka L, Raine NE (2016) Monitoring flower visitation networks and interactions between pairs of bumble bees in a large outdoor flight cage. PLoS One 11:e0150844.

The data used in Shi, Chen, Dong, Wu and Rao (2020): Wuhan-2019-nCoV (https://github.com/canghailan/Wuhan-2019-nCoV/)

## Installation

You can use one of the following two methods to install the released version of GraphCpClust:

```r
devtools::install_github("Meiqian-Chen/GraphCpClust")
```
```r
install.packages("https://github.com/Meiqian-Chen/GraphCpClust/releases/download/v1.0/GraphCpClust_1.0.tar.gz",repos = NULL,type="source")
```

