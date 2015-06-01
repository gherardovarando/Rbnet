
#' New bnet object
#' 
#' Create a bnet object with prescribed variables and optionally arcs selected 
#' via adjacency matrix.
#' 
#' @param variables vector of variables names
#' @param adj.mat  adjacency matrix, not symmetric, acyclicity will be checked 
#' for every possible arcs.
#' @return a bnet object, that is a list with components $variables, 
#' \code{bnet$variables} is also a list with component named
#'         as variables, each element represents a node and is composed of 
#'         three elements: \code{name}, \code{parents} and
#'          \code{prob} (initially set to NULL until fitting of the bnet)
#' @export
#' @examples
#' bnet1<-new_bnet(c("a","b"))
#' summary(bnet1)   
#' plot(bnet1)
new_bnet<-function(variables,adj.mat=NULL){
  bn<-list()
  bn$info$model<-"Bayesian Network"
  bn$variables<-list()
  if (!is.null(variables)){
  for (a in variables){
    bn$variables[[a]]<-list()
    bn$variables[[a]]$name<-a
    bn$variables[[a]]$parents<-c()
  }
  }
  bn$info$fitted<-F
  class(bn)<-"bnet"
  if (!is.null(adj.mat)){
    for (i in 1:(dim(adj.mat)[1])){
      for (j in 1:(dim(adj.mat)[2])){
        if ((adj.mat[i,j]!=0)&(i!=j)){
          bn<-add_arc.bnet(bn,from=variables[i],to=variables[j],check = T)
        }
      }
    }
  }
  return(bn)
}

#' Variables of an object
#' 
#' This is a generic functions that returns the variables/nodes of an R object.
#' @param object an R object
#' @return vector of string
#' @export
#' @examples
#' bnet1<-new_bnet(c("a","b","c","d"))
#' variables(bnet1)
variables<-function(object) UseMethod("variables",object)


#' variables bnet
#' 
#' @param object an object \code{bnet}
#' @return vector of string
#' @export
variables.bnet<-function(object){
  return(attributes(object$variables)$names)
}


#' Convert a bnlearn object to bnet object
#' 
#' @param bn a bn object from package bnlearn
#' @return a bnet object
#' @seealso bnet_to_bn
#' @export
#' @examples
#' \dontrun{data(asia,package= "bnlearn")
#' bn<-bnlearn::hc(asia)
#' bnet<-bn_to_bnet(bn)
#' plot(bnet)
#' summary(bnet)}
bn_to_bnet<-function(bn){
  mopbn<-list()
  mopbn$variables<-list()
  for (a in attributes(bn$nodes)$names){
    mopbn$variables[[a]]<-list()
    mopbn$variables[[a]]$name<-a
    mopbn$variables[[a]]$parents<-bn$arcs[,1][bn$arc[,2]==a]
  }
  class(mopbn)<-"bnet"
  mopbn$info$bnlearn<-T
  return(mopbn)
}

#' Convert a bnet object to bnlearn object
#' 
#' @param bnet a bnet object 
#' @return a bn object from package bnlearn
#' @seealso bn_to_bnet
#' @export
bnet_to_bn<-function(bnet){
  maxp<-max(sapply(bnet$variables,FUN = function(x){length(x$parents)}))
  bn<-list()
  class(bn)<-"bn"
  bn$nodes<-lapply(bnet$variables,FUN = function(x){
    l<-list()
    l$parents<-x$parents
  })
  attributes(bn$nodes)$names<-variables.bnet(bnet)
  bn$arcs<-arcs.bnet(bnet)
  return(bn)
}

#' adjacency matrix for a bnet object
#' 
#' @param object a bnet object
#' @return a matrix, the adjacency matrix relative to \code{object}
#' @export
adj_matrix.bnet<-function(object){
  adjM<-matrix(nrow=length(variables.bnet(object)),ncol=length(variables.bnet(object)),0)
  colnames(adjM)<-variables.bnet(object)
  rownames(adjM)<-colnames(adjM)
  for (a in variables.bnet(object)){
    if (!is.null(object$variables[[a]]$parents)){
      adjM[object$variables[[a]]$parents,a]<-1}
  }
  return(adjM)
}

#' Plot a bnet object
#' 
#' @param x a bnet object
#' @param ... compatibility with plot.default
#' @return a plot of a bnet structure, with targets (red) and predictors (green) 
#' @seealso points.bmop
#' @export
#' @examples
#' bnet<-new_bnet(c("A","B"),matrix(nrow=2,c(0,1,0,0)))
#' plot(bnet)
plot.bnet<-function(x,y=NULL,...){
  attrs<-list()
  attrs$node$fontsize<-"24"
  attrs$edge$arrowsize<-"0.5"
  attrs$edge$style<-"bold"
  g<-graph::graphAM(adjMat = adj_matrix.bnet(x),edgemode = "directed")
  if (is.null(x$targets)){ return(Rgraphviz::plot(g,attrs=attrs ) )}
  nAttrs<-list()
  if (is.null(x$predictors)){
    predictors<-variables.bnet(object = x)
    predictors<-predictors[!(predictors %in% x$targets)]
  }
  else{
    predictors<-x$predictors
  }
  nAttrs$fillcolor<-c(rep(x = "red",times = length(x$targets)),
                      rep("green",times = length(predictors)))
  names(nAttrs$fillcolor)<-c(x$targets,predictors)
  Rgraphviz::plot(g,nodeAttrs=nAttrs,attrs=attrs,... )
}



#' Arcs of an object
#' 
#' @param object 
#' @return a matrix with two columns, each row represent an arc from 
#' the first element to the second.
#' @export
arcs<-function(object){ UseMethod("arcs",object) }

#' Arcs of an object
#' 
#' @param object 
#' @return a matrix with two columns, each row represent an arc from 
#' the first element to the second.
#' @export
#' @examples
#' bnet<-new_bnet(c("A","B"),matrix(nrow=2,c(0,1,0,0)))
#' plot(bnet)
#' arcs(bnet)
arcs.bnet<-function(object){
  maxp<-max(sapply(object$variables,FUN = function(x){length(x$parents)}))
  arcs<-matrix(ncol=2,nrow=length(variables.bnet(object))*maxp)
  i<-1
  for (v in variables.bnet(object)){
    for (b in object$variables[[v]]$parents){
      arcs[i,1]<-b
      arcs[i,2]<-v
      i<-i+1
    }
  }
  return(arcs[!is.na(arcs[,1]),])
}

#' leaf, internal function
leaf<-function(arcs,nodes=NULL){
  if (is.null(nodes)){
    nodes==unique(arcs[,1])
  }
  ix<-!(nodes %in% arcs[,2])
  return((nodes[ix])[1])
}

#' isacyclic, internal function
isacyclic.arcs<-function(arcs,nodes){
  if (is.null(arcs)){return(TRUE)}  
  if (is.null(dim(arcs))){return(TRUE)}  
  if (dim(arcs)[1]==0){ return(TRUE)}
  if (is.null(nodes)){ return(TRUE)}
  l<-leaf(arcs,nodes)
  if (is.null(l)){ return(FALSE)}
  if (is.na(l)){return(FALSE)}
  return(isacyclic.arcs(arcs[arcs[,1]!=l,],nodes[nodes!=l]))
}

#' check if an object fullfil the acyclic property
#' 
#' @param object
#' @return logical
#' @export
is.acyclic<-function(object){UseMethod("is.acyclic",object)}

#' isacyclic bnet
is.acyclic.bnet<-function(object){
  return(isacyclic.arcs(arcs = arcs.bnet(object),
                        nodes = variables.bnet(object)))
}

#' markovblanket of a bnet object
#' 
#' @param object a bnet object
#' @param nodes a vector of nodes of \code{object}
#' @return a vector with the markov blanket of the nodes in \code{nodes}
#' @export
#' @examples
#' \dontrun{data(asia,package="bnlearn")
#' bn<-bnlearn::hc(asia)
#' bnet<-bn_to_bnet(bn)
#' plot(bnet)
#' summary(bnet)
#' markovblanket.bnet(bnet,names(asia)[1])}
markovblanket.bnet<-function(object,nodes){
  mbtotal<-c()
  for (node in nodes){
    child<-c()
    for (a in variables.bnet(object)){
      if (is.element(node,object$variables[[a]]$parents)){
        child<-unique(c(child,a))
      }  
    }
    mb<-child
    for (a in child){
      mb<-unique(c(mb,object$variables[[a]]$parents))
    }
    mbtotal<-unique(c(mb,mbtotal,object$variables[[node]]$parents))
  }
  return(unique(c(mbtotal,nodes))) 
}

#' Summary of a bnet object
#' 
#' @param object a bnet object
#' @export
#' @examples
#' \dontrun{data(gaussian.test,package="bnlearn")
#' bn<-bnlearn::hc(gaussian.test)
#' bnet<-bn_to_bnet(bn)
#' summary(bnet)}
summary.bnet<-function(object,...){
  summ<-list()
  str<-object$variables[[1]]$name
  if (!(length(object$variables[[1]]$parents)==0)){
    
    str<-paste(str,"|",paste(object$variables[[1]]$parents,collapse=":"),sep="")
  }
  for (a in object$variables[-1]){
    str<-paste(str,"][",a$name,sep="")
    if (!(length(a$parents)==0)){
      str<-paste(str,"|",paste(a$parents,collapse=":"),sep="")
    }
  }
  summ$structure<-paste("[",str,"]",sep="")
  summ$nodes<-length(variables.bnet(object))
  summ$arcs<-dim(arcs.bnet(object))[1]
  summ$averagemb<-mean(sapply(object$variables,FUN = function(x){
    return(length(markovblanket.bnet(object = object,nodes = x$name))-1)
  }))
  summ$averageparents<-mean(sapply(object$variables,FUN = function(x){
    return(length(x$parents))
  }))
  summ$info<-object$info
  summ$targets<-object$targets
  summ$predictors<-object$predictors
  class(summ)<-"summary.bnet"
  return(summ)
}

#' summary internal function
print.summary.bnet<-function(x,...){
  cat(x$info$model)
  cat("\n")
  if (!is.null(x$targets)){
    cat("Targets:    ",paste(x$targets,collapse = ", "),"\n")
  }
  if (!is.null(x$predictors)){
    cat("Predictors: ",paste(x$predictors,collapse=", "),"\n")
  }
  cat("\n")
  cat("Model structure: \n")
  cat("  ",x$structure)
  cat("\n")
  cat(paste("  nodes:                       ",x$nodes,sep="\t \t"))
  cat("\n")
  cat(paste("  arcs:                        ",x$arcs,sep="\t \t"))
  cat("\n")
  cat(paste("  average markov blanket size: ",x$averagemb,sep="\t \t"))
  cat("\n")
  cat(paste("  average parents per node:    ",x$averageparents,sep="\t \t"))
  cat("\n")
  cat("\n")
  if (!is.null(x$info$wrapper)){
  if (x$info$wrapper$state){
  cat("Wrapper informations: \n")
  cat("  evaluation function: ",x$info$wrapper$evalfun,"\n")
  cat("  evaluation:          ",x$info$wrapper$eval,"\n")
  }
  }
}

#' Print a bnet object
#' 
#' @param x a bnet object
#' @export
print.bnet<-function(x,...){
  print(summary.bnet(x))
  #   cat("BMoP Bayesian network \n \n")
  #   for (a in x$variables){
  #     cat(a$name)
  #     if (!(length(a$parents)==0)){
  #       cat("\t parents: ")
  #       cat(paste(a$parents,collapse=","))
  #       cat("\n")
  #     }
  #     cat(a$mop$knots)
  #     cat(a$mop$ctrpoints)
  #     cat("\n")
  #   }
}

#' child of selected nodes in a bnet object
#' 
#' @param object a bnet object
#' @param nodes a vector of nodes of bnet
#' @return a vector with the names of the child of \code{nodes}
#' @export
child.bnet<-function(object,nodes){
  child<-c()
  for (node in nodes){
   for (n in object$variables){
     if (node %in% n$parents){
       child<-c(child,n$name)
     }
   }
  }
  return(unique(child))
}

#' check if a bnet is fitted
#' 
#' @param object a bnet object
#' @return logical
#' @export
is.fitted.bnet<-function(object){
  if (is.null(object$info$fitted)){ return(FALSE)}
  return(object$info$fitted)
}


#' Check if an object is abnet object
#' 
#' @param x an R object
#' @return logical
#' @export
is.bnet <- function(x){
  return(inherits(x = x,what = "bnet"))
}

#' Convert object to bnet 
#' 
#' ..when possible (matrix, bnlearn, graph)
#' @param x a R object
#' @param ... parameter for compatibility reason
#' @return a bnet object
#' @export
#' @examples
#' mat <- rbind(c(0, 0, 1, 1),
#' c(0, 0, 1, 1),
#' c(0, 0, 0, 1),
#' c(0, 0, 0, 0))
#' g1 <- graph::graphAM(adjMat=mat,edgemode = "directed")
#' formula <- Y ~ X1 + X2 + X3 + X4 
#' bnet1 <- as.bnet(g1)
#' bnet2 <- as.bnet(mat)
#' bnet3 <- as.bnet(3*mat)
#' bnet4 <- as.bnet(formula)
as.bnet<-function(x,...){
  
  if (inherits(x,what= "bnet")){
    return(x)
  }
  
  if (inherits(x,what = "bn")){
    return(bn_to_bnet(bn = x))
  }
  
  if (is.matrix(x)){
    if (is.null(colnames(x))){
      var<-as.character(1:(dim(x)[2]))
    }
    else{
      var<-colnames(x)
    }
    return(new_bnet(variables =var,adj.mat = x ))
  }
  
  if (is.data.frame(x)){
    return(new_bnet(variables = names(x),adj.mat = as.matrix(x)))
  }
  
  if (is.character(x)){
   return(new_bnet(variables = x))  
  }
  
  
  
  if (class(x)=="graphNEL"){
    x<-as(x,"graphAM")}
  
  if (class(x)=="graphAM"){
    return(new_bnet(variables= graph::nodes(x),adj.mat =as(x,"matrix") ) )
  }
  
  if (inherits(x,"formula")){
    return(naive_bayes.bnet(formula=x))
  }
  
  if (is.na(x)){
    warning("NA introduced by coercion")
    return(NA)
  }
  
  return(as.bnet(as.matrix(x)))
}


AIC.bnet <- function(object, data, ... , k=2 ){
  aa <- 0
  for (a in variables(object)){
    aa <- aa + AIC(object$variables[[a]]$prob , k = k , 
                   data = data[, c(a, object$variables[[a]]$parents )] ) 
  }
  return(aa)
}

BIC.bnet <- function(object, data, ... ){
 return(AIC.bnet(object,data,k=log(dim(data)[1])))
}