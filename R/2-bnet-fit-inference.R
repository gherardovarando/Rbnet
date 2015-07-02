.bnetenv<-new.env()
lockBinding(sym = ".bnetenv",env = environment())
assign(".bnetPars",value = list(searchBmop=FALSE,
                                intMethod="adaptative"), envir=.bnetenv  )
lockEnvironment(env = .bnetenv,bindings = ".bnetPars")

#' Set fitting parameters
#' 
#' Set appropriate parameters to be used by fitting
#' functions.
#' 
#'  @param ... see details.
#'  @details This function set parameters used for estimation.
#'           \code{N}:  The number of knots in every dimensions.
#'           \code{order}: The order of the B-spline in every dimensions.
#'           \code{alpha}: The exponent to compute the number of knots.
#'          \code{knotsMethod}: \code{"uniform"} or \code{"quantiles"} knots. 
#'           
#'           @export
bnetPar<-function(...){
  
  l <- list(...)
  
  old <- get(".bnetPars",envir = .bnetenv)
  if (length(l) == 0){
    return(old)
  }
  unlockBinding(sym = ".bnetPars",env = .bnetenv)
  for (a in names(l)){
    if (a %in% names(old)){
      if (a=="intmethod"){
        if (!(l[[a]] %in% c("adaptative","sampling"))){
          l[[a]]<-"adaptative"
          warning("intMethod parameter can be only adaptative or sampling, 
                  set to adaptative by default")
        }
      }
      old[[a]] <- l[[a]]
    }
    else{
      warning(paste(a,"is not a Rbnet parameter"))
    }
  }
  assign(".bnetPars",old,envir = .bnetenv)
  lockBinding(sym = ".bnetPars",env = .bnetenv)
}




# Evaluate an object over points (INTERNAL FUNCTION)
# 
# This is a generic functions
# @param x points over perform evaluations
# @param object R object to be evaluated \code{evaluate.foo} must be available 
# for object of class \code{foo}
# @return array of evaluations
evaluate<-function(x,object,...){
  UseMethod("evaluate",object)
}


# (INTERNAL FUNCTION)
evaluate.evidence<-function(x,object,MIN=10^(-10)){
  if (is.null(dim(x))){
    dim(x)<-c(1,length(x))
  }
  if (dim(x)[1]>1){ return(apply(x,MARGIN = 1,FUN =
                                   evaluate.evidence,object=object,MIN=MIN))}
  if (x[1]==object){ return(1)}
  else{return(max(MIN,0))}
}


# Upper and Lower values for an object (INTERNAL FUNCTION)
# 
# This is a generic function, used by package bnet to compute lower and upper
# boundaries for the \code{$prob} values of a node.
# @param object an R object, e.g. a bmop object
# @return numeric value
lower<-function(object){ UseMethod("lower",object)}

# Upper and Lower values for an object (INTERNAL FUNCTION)
# 
# This is a generic function, used by package bnet to compute lower and upper
# boundaries for the \code{$prob} values of a node.
# @param object an R object, e.g. a bmop object
# @return numeric value
upper<-function(object){ UseMethod("upper",object)}


# internal lower for evidence (INTERNAL FUNCTION)
lower.evidence<-function(object){
  return(object)
}

# internal upper for evidence (INTERNAL FUNCTION)
upper.evidence<-function(object){
  return(object)
}


# fast fit, internal function (INTERNAL FUNCTION)
bmop_fit.bnet<-function(object,data,nodes,unif.variables=NULL,Mins=NULL,
                        Maxs=NULL,...){
  nodes<-nodes[nodes %in% variables.bnet(object)]
  for (a in nodes){
    var<-c(a,object$variables[[a]]$parents)
    D<-data[,var]
    if (!is.element(el = a,set = unif.variables)){
   
          object$variables[[a]]$prob<-Rbmop::bmop_fit(data = D,conditional = T,
                                                           Min = Mins[var],
                                                           Max=Maxs[var],
                                                           ...)
    

    }
    else{
      object$variables[[a]]$prob<-
        Rbmop::new_bmop(knots = Rbmop::generate_knots(data=D,N = 1),
                        order = 1,ctrpoints = 1)
      object$variables[[a]]$prob<-
        Rbmop::normalize.bmop(object$variables[[a]]$prob)  
      object$info$fit$unif.variables<-unif.variables
    }
  }
  object$info$fitted<-TRUE
  return(object)
}



#' Fit a bnet object
#' 
#' @param object a bnet obejct
#' @param data data.frame with names that include the variables of \code{object}
#' @param search logical, use penalized loglik search or not
#' @param ... additional parameters to be passed to the bmop fitting functions
#' @return a bnet object with fitted densities and conditional densities
#' @export
#' @examples
#' \dontrun{require(bnlearn)
#' data(gaussian.test)
#' bn<-hc(gaussian.test,maxp=2)
#' plot(bn)
#' bnet<-as.bnet(bn)
#' plot(bnet)
#' is.fitted.bnet(bnet)
#' bnet<-fit.bnet(bnet,gaussian.test[1:100,]) # just 100 observation
#'  are used to make computation lighter
#' is.fitted.bnet(bnet)
#' plot(bnet$variables$A$prob)
#' plot(bnet$variables$B$prob)}
#' X<-rnorm(100)
#' Y<-rnorm(100)
#' U<-rnorm(100,mean=X+Y)
#' V<-rnorm(100,mean=X)
#' Z<-rnorm(100,mean=U)
#' data<-data.frame(X,Y,U,V,Z)
#' mat<-matrix(nrow=5, c(0,0,1,1,0,
#'                      0,0,1,0,0,
#'                      0,0,0,0,1,
#'                      0,0,0,0,0))
#' bnet<-new_bnet(names(data),mat)
#' bnet<-fit.bnet(bnet,data)
fit.bnet<-function(object,data,nodes=NULL,Mins=NULL,Maxs=NULL,...){
  if (is.null(nodes)){
    nodes<-variables.bnet(object)
  }
    
  if (is.null(Mins)){
    Mins<-lapply(data,min)
    names(Mins)<-variables(object)
  }
  if (is.null(Maxs)){
    Maxs<-lapply(data,max)
    names(Maxs)<-variables(object)
  }
    return(bmop_fit.bnet(object = 
                                  object,data = data,nodes=nodes,
                         Mins=Mins,Maxs=Maxs,...))

}


#' put evidence on a bnet object
#' 
#' @param object a bnet object
#' @param evidence numeric with names or list, the value of evidence
#' @param store.old logical
#' @param propagate logical, see details
#' @return a bnet object with evidence
#' @details If  \code{propagate} is set to \code{TRUE}, the evidence will be 
#' naively propagate, that is propagate just to the direct child.
#' @export
put_evidence.bnet<-function(object,evidence,store.old=T,propagate=T){
  evidence<-as.list(evidence)
  for (a in names(evidence)){
    if (is.null(object$variables[[a]]$name)){
      object$variables[[a]]$name<-a
      
    }
    if (store.old){
        object$variables[[a]]$prob.old<-object$variables[[a]]$prob
      if (propagate){
      for (n in child.bnet(object = object,nodes = a)){
        object$variables[[n]]$prob.old<-object$variables[[n]]$prob
      }  
      }
    }
    object$variables[[a]]$prob<-evidence[[a]]
    class(object$variables[[a]]$prob)<-"evidence"
  }
  
  if (propagate){
    for (a in names(evidence)){
    for (n in child.bnet(object = object,nodes = a)){
      par<-object$variables[[n]]$parents
      if (!class(object$variables[[n]]$prob)=="evidence"){
        
      object$variables[[n]]$prob<-Rbmop::put_evidence.bmop(
        evd.pos = (2:(length(par)+1))[par==a] ,
        evidence = evidence[[a]]  ,
        object = object$variables[[n]]$prob  )
      }
    }  
  }
  }
  object$store.old<-store.old
  return(object)
}


#' remove all evidence from a bnet object
#' 
#' @param object a bnet object
#' @return a bnet object, the same as 
#' \code{object} if old probabilities were not stored when evidence was set, 
#' otherwise the restored bnet object
#' @export
clear_evidence.bnet<-function(object){
  if (is.null(object$store.old)){
    warning("No old probability, sure evidence was set?...")
    return(object)
  }
  if (!object$store.old){
    warning("Old probability was lost putting evidence...
            Next time try saving it")
    return(object)
  }
  else{
    for (a in variables(object)){
      object$variables[[a]]$prob<-object$variables[[a]]$prob.old
      object$variables[[a]]$prob.old<-NULL
    }
  }
  object$store.old<-NULL
}


#' compute probability for a bnet 
#' 
#' @param x vector of evidence
#' @param object 
#' @param log logical, if \code{TRUE} the logarithm of the probabilities 
#' is returned
#' @return numeric or \code{NULL} if partial evidence is provided
#' @export
probability.bnet<-function(x,object,log=F){
  if (!object$info$fitted){ 
    warning("not fitted bnet object, use fit.bnet to fit the 
            Bayesian network, NULL is returned")
    return(NULL)}
  if (is.null(dim(x))){ dim(x)<-c(1,length(x))}
  if (is.null(names(x))){ names(x)<-variables.bnet(object)}
  x<-as.data.frame(x)
  vars<-variables.bnet(object)
  if (all(vars %in% names(x))){
    return(apply(x,MARGIN = 1,FUN = function(d){
      if (!log) {
      return(prod(sapply(object$variables,FUN = function(vari){
        return(evaluate(object = vari$prob,x = d[c(vari$name,vari$parents)], 
                        MIN = 0))
      })))
      
      }
      else {
        return(sum(log((sapply(object$variables,FUN = function(vari){
          return(evaluate(object = vari$prob,x = d[c(vari$name,vari$parents)]))
        })))))
      }
    }))
  }
else {
  to.integrate<-var[!(var %in%  names(x))]
  warning("partial evidence is not yet implemented")
  return(NULL)
}
}

#' predict posterior values for a bnet object
#' 
#' @param object a bnet object
#' @param newdata new data, evidence
#' @param targets targets variables
#' @param predictors predictors variables
#' @param method string "posteriorMean" or "posteriorMode"
#' @param optmethod \link{optim}
#' @param intmethod string "adaptative" or "sampling"
#' @return the values predicted
#' @export
predict.bnet<-function(object,newdata,targets=NULL,predictors=NULL,
                       method="posteriorMean",optmethod="Nelder-Mead",
                       intMethod=bnetPar()$intMethod,...){
  if (!object$info$fitted){ 
    warning("not fitted bnet object, use fit.bnet to fit the Bayesian
            network, NULL is returned")
    return(NULL)}
  if (is.null(targets)){
    if (is.null(object$targets)){
      targets<-variables.bnet(object)[1]
    }
    else{
    targets<-object$targets
    }
  }
  if (is.null(predictors)){
    if (is.null(object$predictors)){
      predictors<-variables.bnet(object)
      predictors<-predictors[!(predictors %in% targets) ]
    }
    else{
    predictors<-object$predictors
    }
  }
  nam<-names(newdata)
  for (a in targets){
    if (!(a %in% nam)){
      newdata[,a]<-rep(0,dim(newdata)[1])
    }
  }
  mb<-markovblanket.bnet(object = object,nodes=targets)  
  if (!(all(mb %in% nam))){
    warning("please provide a dataset with the needed predictor variables,
            that is the markov blanket of the targets")
    return(NULL)
  }
    if (dim(newdata)[1]==1){  
      part<-sapply(X = object$variables[targets],FUN=function(n){
        return(lower(n$prob))
      } )
      fin<-sapply(X = object$variables[targets],FUN=function(n){
        return(upper(n$prob))
      } )
      if (method=="posteriorMode"){
        if (is.null(dim(newdata))){
          dim(newdata)=c(1,length(newdata))
        }
        if (length(targets)==1){ optmethod="Brent" }
        ctr=list()
        ctr$fnscale=-1
        return(optim(lower =part ,upper =fin ,
                     method=optmethod,control=ctr,par=part,fn = function(x){ 
          newdata[,targets]<-x
          return(
            sum(sapply(object$variables,FUN = 
                         function(a){
                           log(evaluate(
                             x = as.numeric(newdata[,c(a$name,a$parents)]) ,
                             object = a$prob))}))
          )
        })$par)     
      }
      if (method=="posteriorMean"){  
        child<-child.bnet(object,targets)
        #cc<-prod(sapply(object$variables[!(variables.bnet(object)
        #%in% child) ],FUN =function(a){ evaluate(MIN = 10^(-10),x =
          #             as.numeric(newdata[,c(a$name,a$parents)]) ,
        #object = a$prob)} ))
        densityi<-function(x){
          newdata[,targets]<-x
          return(prod(sapply(object$variables[c(targets,child)],
                             FUN = function(a){
                               evaluate(MIN = 10^(-10),
                                        x = as.numeric(
                                          newdata[,c(a$name,a$parents)]) ,
                                        object = a$prob)})))
        }
        if (intMethod=="adaptative"){
          R<-rep(0,times = length(targets))
          for (i in 1:length(targets)){
            konst<-
              cubature::adaptIntegrate(f =densityi,lowerLimit = part,
                                       upperLimit = fin )$integral
            R[i]<-cubature::adaptIntegrate(f =function(x){x[i]*densityi(x)},
                                           lowerLimit = part,upperLimit = fin )$
              integral/konst
          }
        }
        if (intMethod=="sampling"){
          D<-length(targets)
          xstart<-(part+fin)/2
          sample<-Rbmop::sampler_MH(N = 500,d = D,densit = densityi,
                                    h = 3,M=500,xstart=xstart,max=fin,min=part)
          if (D>1){colnames(sample)<-targets}
          R<-colMeans(sample)
        }
        return(R)
      }
    }
    else { 
      l<-1:dim(newdata)[1]
      dim(l)<-dim(newdata)[1]
      return(t(apply(l,MARGIN = 1,FUN = function(i){ return(
        predict.bnet(method=method,object=object,targets=targets,
                     predictors=predictors,newdata=newdata[i,],
                     intmethod=intmethod))})))
    
    }
  
}


