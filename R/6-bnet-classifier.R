

#' Multinet classifier using bnet
#' 
#' Learn a multi net classifier unsing bnet, in particular learn a model
#' where for every class value a bnet object is estimated. 
#' As particular case naive Bayes (\code{bnets=NULL}).
#' @param targets
#' @param predictors
#' @param prior string, NULL, or vector of positive numbers
#' @param data data set of observations
#' @param bnets optionally a list of bnet object, one for each class value
#' @return a classifier_bnet object
#' @export
classifier_bnet<-function(formula=NULL,targets,predictors,prior=NULL,data=NULL
                          ,bnets=NULL,...){
  
  if (inherits(formula,"formula")){
    targets<-as.character(formula[[2]])
    targets<-targets[! (targets %in% c("+","~","-","*") )]
    predictors<-all.vars(formula)[! (all.vars(formula) %in% targets)]
  }
  
  if (is.null(bnets)){
    levs<-levels(as.factor(data[,targets]))
    bnets <- lapply(levs,FUN = function(x){
      new_bnet(variables = predictors)})
    names(bnets)<-levs
  }
  Maxs<-lapply(predictors,FUN = function(a){ max(data[,a])})
  names(Maxs)<-predictors
  Mins<-lapply(predictors,FUN = function(a){ min(data[,a])})
  names(Mins)<-predictors
  nb<-list()
  if (!is.null(prior)){
    if (prior=="uniform"){
      prior<-rep(1,length(levels(data[,targets])))
    }
    nb$prior<-prior/sum(prior) 
  }
  
  if (!is.null(data)){
    data[, targets]<-as.factor(data[, targets])
    nb$levels<-levels(as.factor(x = data[,targets]))
    if (is.null(prior)){
      nb$prior<-sapply(X = nb$levels,FUN = function(l){
        return((1+length(data[data[,targets]==l,targets])))
      })
    }
    nb$prior<-nb$prior/sum(nb$prior)
    
    nb$bnets<-lapply(nb$levels,FUN = function(l){
      return(fit.bnet(object = bnets[[l]],
                      data =   data[ data[,targets]==l, predictors]
                      ,Maxs=Maxs,Mins=Mins))
    })
    names(nb$bnets)<-nb$levels
  }
  nb$targets<-targets
  nb$predictors<-predictors
  class(nb)<-"bnetClassifier"
  return(nb)
}



#' prediction of \code{classifier_bnet} object
#' @export
predict.bnetClassifier<-function(object,newdata,log.odds=F,...){
  targets<-object$targets
  newdata<-as.data.frame(newdata[,object$predictors])
  prd<-array(dim = c(dim(newdata)[1],1))
  for (i in 1:(dim(newdata)[1])){
    
    prob<-sapply(X = object$levels,FUN = function(l){
      return(probability.bnet(x = newdata[i,],object = object$bnets[[l]]
                              ,log = T))
    })
    prob[is.null(prob)] <- -Inf
    prob<-prob+log(object$prior)
    sorted<-sort(decreasing = T,x = prob,index.return=T)
    if (log.odds){
      prd[i]<-prob[1]-prob[2]
    }
    else{
    prd[i]<-object$levels[sorted$ix[1]]
    }
  }
  return(prd)
  
}




#' Plot various densities for bnet classifier
#'
#' @param x a bnetClassifier object
#' @param predictor the name of one of the predictor
#' @param ... other parameters for compatibility
#' @return invisible()
#' @export 
densplot.bnetClassifier<-function(x,predictor=NULL,...){
  if (is.null(predictor)){
    predictor<-x$bnets[[1]]$variables[[1]]$name
  }
  l<-list()
    for (i in x$levels){
      l[[i]]<-x$bnets[[i]]$variables[[predictor]]$prob
    }
  m<-lower(x$bnets[[1]]$variables[[predictor]]$prob)
  M<-upper(x$bnets[[1]]$variables[[predictor]]$prob)
    Rbmop::comparison_plot(l,dtrue = NULL,
                           names.bmop = x$levels,
                           colors = rainbow(length(l)) ,ylim=c(0,3/(M-m)))
  return(invisible())
}


#' Plot a bnet classifier
#'
#' @param x a bnetClassifier object
#' @param ... other parameters for compatibility
#' @return invisible()
#' @export 
plot.bnetClassifier<-function(x,...){
    plot(naive_bayes.bnet(bnet = x$bnets[[1]],targets =x$targets,
                          predictors=variables(x$bnets[[1]]) ))
return(invisible())
}


#' Classifier One vs All 
#' 
#' @param target
#' @param predictors
#' @param data data set of observations
#' @return a multiclass_bnet object
#' @export
onevsall_bnet<-function(targets,predictors,data=NULL){
  if (length(levels(data[[targets]]))<=2){warning("Use this functions with 
                                                 multi-class problems")
                                         return(classifier_bnet(targets,
                                                                predictors, 
                                                                data=data))
  }
  for (l in levels(data[[targets]])){
    data[[l]]<-rep(0,dim(data)[1])
    data[[l]][data$Layer==l]<-1
    data[[l]]<-as.factor(data[[l]])
  }
  mod<-lapply(levels(data[[targets]]),FUN = function(l){
    return(classifier_bnet(targets = l,predictors = predictors,prior = NULL,
                           data=data,bnets=NULL))
  })
  names(mod)<-levels(data[[targets]])
  class(mod)<-"onevsall_bnet"
  return(mod)
  
}


#' predict function for onevsall_bnet object
#' @export
predict.onevsall_bnet<-function(object,newdata,...){
  odds<-lapply(object,FUN = predict,newdata=newdata,log.odds=T)
  odds<-as.data.frame(odds)
  odds$unknow<-rep(0,dim(odds)[1])
  res<-apply(odds,MARGIN = 2,FUN = function(x){
    x<-sort(x,index.return=T)
    return(names(odds)[x$ix])
  })
  return(res)
}

#' variables for bnetClassifier
#' 
#' @param object a bnetClassiifer object
#' @return vector of variables name
#' @export
variables.bnetClassifier<-function(object){
  return(c(variables(object$bnets[[1]])))
}


