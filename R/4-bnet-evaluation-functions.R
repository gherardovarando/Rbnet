
#' validate regression bnet models
#' 
#' @param object optional value, a bnet object
#' @param modelfun function to compute a bnet regression model, optional 
#' @param data data.frame a dataset of observation
#' @param targets names of target variables
#' @param predictors names of predictor variables
#' @param folds folds for validation 
#' @param numfolds number of folds 
#' @param randomfolds logical
#' @param train indx of train set
#' @param test  indx of test set
#' @param evalfun function, must works as MSE or ARE
#' @param mle logica
#' @param mc.cores positive integer
#' @return mean and standard deviation of the \code{numfolds} evaluations,
#'        each one of the evaluations is computed using \code{evalfun} 
#'        (usally \code{MSE}) over the prediction and the true value for the 
#'        given fold (ot train/test set).  
#' @export
validation<-function(x=naive_bayes.bnet,data=NULL,targets=NULL,
                        predictors=NULL,folds=NULL,numfolds=10,randomfolds=F,
                        train=NULL,test=NULL,evalfun=MSE,mc.cores=1){
  prl<-F
  if (is.null(data)){
    warning("provide dataset")
    return(NULL)
  }
  if (is.null(folds)){
    if ((is.null(train))&(!(is.null(test)))){
      train<-(1:dim(data)[1])[-test]    
    }
    if ((is.null(test))&(!(is.null(train)))){
      test<-(1:dim(data)[1])[-train]    
    }
    if (!(is.null(train))&(!(is.null(test)))){
      folds<-list()
      folds$folds$test<-test
      folds$folds$training<-train
    }
    if (is.null(folds$folds)){
      folds<-create_folds(data = data,k = numfolds,random = randomfolds)
    }
  }
  
  if (is.function(x)){
  modelfun<-x
      eval<-parallel::mclapply(mc.cores = mc.cores,X = folds$folds,FUN = function(ff){
        model<-modelfun(targets=targets,predictors=predictors,
                       data=data[ff$training, ])
        prd<-predict(object = model,newdata = data[ff$test, ])
        return(evalfun(prd,data[ff$test, targets]))
      })
      eval<-as.numeric(eval)
    return(list(mean=mean(eval),sd=sd(eval)))
    
}
  
  if (is.bnet(x)){
    object<-x
  eval<-parallel::mclapply(mc.cores = mc.cores,X = folds$folds,FUN = function(ff){
      fitted<-fit.bnet(object = object,data = data[ff$training, ],mle=mle,
                       search=search)
      prd<-predict(object = fitted,newdata = data[ff$test, ])
      return(evalfun(prd,data[ff$test, object$targets]))
    })
  eval<-as.numeric(eval)
  return(list(mean=mean(eval),sd=sd(eval)))
  }
  
  
  return(NULL)
  
  
}




#' zero one accuracy
#' 
#' @param x a data.frame or a matrix
#' @param y a data.frame or a matrix
#' @return positive numbers, the fraction of matching entries between \code{x}
#' and \code{y}
#' @export
zero_one_accuracy<-function(x,y){
  I<-array(dim = dim(x),data = 1)
  I[x!=y]=0
  return(sum(I)/prod(dim(I)))
}



