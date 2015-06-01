
#' MSE function
#' @param x
#' @param y
#' @export
MSE<-function(x,y){
  N<-dim(x)
  if (is.null(N)){ N<-c(1,length(x))}
  return( sum((x-y)^2)/(N[1]*N[2])   )
  
}

#' ARE function
#' @param x
#' @param y
#' @export
ARE<-function(x,y){
  N<-dim(x)
  if (is.null(N)){ N<-c(1,length(x))}
  return( sum(abs(x-y))/(N[1]*N[2])   )
}

#' create folds for validation
#' 
#' @param data a data.frame
#' @param k positive integer number of folds
#' @param random, logical, random order?
#' @return partitions in folds
#' @export
create_folds<-function(data,k=10,random=F){
  D<-list()
  N<-dim(data)[1]
  avoid<-c(N+2)
  k<-min(k,N)
  D$folds<-list()
  for (i in 1:k){
    if (k==N){ tst<-i}
    else {
      if (random) {tst<-sample((1:N)[-avoid],size=floor(N/k),replace=F)}
      else{tst<-(floor(N/k)*(i-1)+1):(floor(N/k)*i)  }
    }
    D$folds[[i]]<-list()
    D$folds[[i]]$training<-(1:N)[-tst]
    D$folds[[i]]$test<-tst
    D$folds[[i]]$ix<-i
    avoid<-c(avoid,tst)
  }
  D$info$k<-k
  D$info$random<-random
  D$info$N<-N
  return(D)
}



