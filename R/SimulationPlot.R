#' Basic description
#'
#' @description  This fuction is the comparisons of graph-based clustering method and model-based clustering method introduced in Shi, Chen, Dong and Rao (2020)
#' @usage SimulationPlot()
#' @examples
#' SimulationPlot()
#' @seealso 
#' @export
SimulationPlot<-function(){
	add_legend <- function(...) {
	  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
	              mar=c(0, 0, 0, 0), new=TRUE)
	  on.exit(par(opar))
	  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n' )
	  legend(...)
	}
	# layout(matrix(c(1,1,1,2,3,4),2,3,byrow=T))
	fig2 <- function(...){
		D=c(150);sigma=c(1:10)*0.1;loop.time=100;c(1,10,15)
		set.seed(1)
		s.e=1;d=D
		y1=10*pnorm(-4+0.05*c(1:d))+rnorm(d,0,s.e)#
		y2=c(rep(0,50),20*pnorm(-3+0.03*c(1:(d-50)))+rnorm(d-50,0,s.e))
		y3=c(rep(0,100),5*pnorm(-2+0.07*c(1:(d-100)))+rnorm(d-100,0,s.e))

		#figure2
		plotchar=c(0,3,1)
		linetype=c(1,3,5)
		COL=c("blue","red","green")

		plot(c(1:d),y1,ylab="",xlab=expression(t),type="b", lwd=1.5,lty=linetype[1],
		 pch=plotchar[1], col=COL[1]) 
		lines(c(1:d),y2,ylab="",xlab=expression(t),type="b", lwd=1.5,lty=linetype[2],
		 pch=plotchar[2], col=COL[2]) 
		lines(c(1:d),y3,ylab="",xlab=expression(t),type="b", lwd=1.5,lty=linetype[3],
		 pch=plotchar[3], col=COL[3]) 
		    
		add_legend("top", legend=c("Class 1","Class 2", "Class 3"), 
		cex=1, bty='n', horiz=TRUE, lty=linetype,lwd=1.5, pch=plotchar, col=COL  ) 
	}


  fig3 <- function(...){
  	D=c(150);sigma=c(1:10)*0.1;loop.time=100;c(1,10,15)
  	m1=readRDS(system.file("extdata/REg.rds",package="PathCpClust"))
		t1=readRDS(system.file("extdata/REm.rds",package="PathCpClust")) 
		m2=readRDS(system.file("extdata/REg2.rds",package="PathCpClust"))
		t2=readRDS(system.file("extdata/REm2.rds",package="PathCpClust")) 
		m3=readRDS(system.file("extdata/REg3.rds",package="PathCpClust"))
		t3=readRDS(system.file("extdata/REm3.rds",package="PathCpClust")) 
		##########end of simulation3

		#Figure 3
		plotchar=c(0,3)
		linetype=c(1,3)
		COL=c("blue","red")

		par(mfrow=c(1,3))
		 
		par(mar=c(4,5,5,1)+0.1,fig=c(0,0.333,0,1) )

		plot(sigma,m1[1,],ylab=expression(paste("S"^"*")),xlab=expression(sigma), ylim=c(-1,1),type="b", lwd=1.5,lty=linetype[1],
		 pch=plotchar[1], col=COL[1],main=expression(paste('n'[1]==1, ", ",'n'[2]==10, ", ",'n'[3]==15) ) ) 
		lines(sigma,t1[1,],type="b", lwd=1.5,lty=linetype[2],
		 pch=plotchar[2], col=COL[2]) 

		par(mar=c(4,5,5,1)+0.1,fig=c(0.333,0.666,0,1), new=TRUE)

		plot(sigma,m2[1,],
		main=expression(paste('n'[1]==5,", ",'n'[2]==20,", ", 'n'[3]==30)),ylab=expression(paste("S"^"*")),xlab=expression(sigma), ylim=c(-1,1),type="b", lwd=1.5,
		lty=linetype[1],
		 pch=plotchar[1], col=COL[1]) 

		lines(sigma,t2[1,],type="b", lwd=1.5,lty=linetype[2],
		 pch=plotchar[2], col=COL[2]) 

		par(mar=c(4,5,5,1)+0.1,fig=c(0.666,0.999,0,1), new=TRUE)
		plot(sigma,m3[1,],
		main=expression(paste('n'[1]==20,", ",'n'[2]==100, ", ",'n'[3]==200)),ylab=expression(paste("S"^"*")),xlab=expression(sigma), ylim=c(-1,1),type="b", lwd=1.5,
		lty=linetype[1],
		 pch=plotchar[1], col=COL[1]) 
		 
		lines(sigma,t3[1,],type="b", lwd=1.5,lty=linetype[2],
		 pch=plotchar[2], col=COL[2]) 

		add_legend("top", legend=c("Graph-based clustering    ","Model-based clustering    "), 
		cex=1, bty='n', horiz=TRUE, lty=linetype,lwd=1.5, pch=plotchar, col=COL  ) 
  }
  # par(mar=c(2,2,2,1),oma = c(0, 0, 3, 0))
  # layout(matrix(c(1,2),2,1,byrow=F))
  fig2()
  fig3()
  # return(list(fig2=fig2(),fig3=fig3()))
}

