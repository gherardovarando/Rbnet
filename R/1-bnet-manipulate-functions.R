
#' add nodes to a bnet object
#' 
#' @param object a bnet obejct
#' @param nodes a vector of names for the nodes to be added
#' @param data optionally a data.frame for fitting the new nodes
#' @param ... additional parameters to be passed to bmop fitting functions
#'  @return a bnet object
#' @export
#' @examples
#' mat <- rbind(c(0, 0, 1, 1),
#' c(0, 0, 1, 1),
#' c(0, 0, 0, 1),
#' c(0, 0, 0, 0))
#' bnet<-as.bnet(mat)
#' plot(bnet)
#' bnet<-add_nodes.bnet(bnet,"5")
#' plot(bnet)
add_nodes.bnet<-function(object,nodes,data=NULL,...){
  for (node in nodes){
  if (node %in% variables.bnet(object)){
    warning("node is already in the network, original bnet is returned")
  }
  else{
    object$variables[[node]]$names<-node
  }
  if (!is.null(data)){
    object<-fit.bnet(object,data,nodes=nodes,,...)
  }
  }
  return(object)
}

#' add arc to a bnet object
#' 
#' @param object a bnet obejct
#' @param from node to start the arc from
#' @param to node to end the arc
#' @param data optionally a data.frame for fitting the \code{to} node that 
#' now has an additional parent
#' @param check logical, to check or not for acyclic property
#' @param ... additional parameters to be passed to bmop fitting functions
#'  @return a bnet object
#' @export
#' @examples
#' mat <- rbind(c(0, 0, 1, 1),
#' c(0, 0, 1, 1),
#' c(0, 0, 0, 1),
#' c(0, 0, 0, 0))
#' bnet<-as.bnet(mat)
#' plot(bnet)
#' bnet<-add_arc.bnet(bnet,"1","2")
#' plot(bnet)
add_arc.bnet<-function(object,from,to,data=NULL,check=T,...){
  old<-object
  variables<-variables.bnet(object)
  if (from==to){
    warning("Start and end node are the same, original object is returned")
    return(object)}
  if (!(from %in% variables)){ object<-add_nodes.bnet(object,
                                                      from,data=data[,from])}
  if (!(to %in% variables)){ object<-add_nodes.bnet(object,to)}
  object$variables[[to]]$parents<-unique(c(object$variables[[to]]$parents,from))
  object$info$fitted<-F
  if (!is.acyclic.bnet(object = object)){
   warning("a cycle is formed by adding the arc")
   if (check){
     warning("original network is returned")
     return(old)
   }
  }
  if (!is.null(data)){
    object<-fit.bnet(object,data,nodes=to,...)
  }
  return(object)
}

#'delete arc from a bnet object
#' 
#' @param object a bnet obejct
#' @param from node to start the arc from
#' @param to node to end the arc
#' @param data optionally a data.frame for fitting 
#' @param ... additional parameters to be passed to bmop fitting functions
#'  @return a bnet object
#' @export
#' @examples
#' mat <- rbind(c(0, 0, 1, 1),
#' c(0, 0, 1, 1),
#' c(0, 0, 0, 1),
#' c(0, 0, 0, 0))
#' bnet<-as.bnet(mat)
#' plot(bnet)
#' bnet<-delete_arc.bnet(bnet,"1","3")
#' plot(bnet)
delete_arc.bnet<-function(object,from,to,data=NULL,...){
  variables<-variables.bnet(object)
  arcs <- arcs.bnet(object = object)
  tos <- arcs[arcs[, 1]==from, 2]
  if (!(to %in% tos)){
       warning("the arc is not present in the network")
       return(object)
  }
  pp<-object$variables[[to]]$parents
  object$variables[[to]]$parents<-pp[!(pp==from)]
  object$info$fitted<-F
  if (!is.null(data)){
    object<-fit.bnet(object,data,nodes=to,...)
  }
  return(object)
}

#' delete node (one!) from a bnet object
#' 
#' @param object a bnet obejct
#' @param node node to be deleted
#' @param data optionally a data.frame for fitting 
#' @param ... additional parameters to be passed to bmop fitting functions
#' @return a bnet object
#' @export
#' @examples
#' mat <- rbind(c(0, 0, 1, 1),
#' c(0, 0, 1, 1),
#' c(0, 0, 0, 1),
#' c(0, 0, 0, 0))
#' bnet<-as.bnet(mat)
#' plot(bnet)
#' bnet<-delete_node.bnet(bnet,"1")
#' plot(bnet)
delete_node.bnet<-function(object,node,data=NULL,...){
  var<-variables.bnet(object)
  if (!(node %in% var)){
    warning("node is note present in the network, original object is returned")
    return(object)
  }
  child<-child.bnet(object = object,node)
  object$variables[[node]]<-NULL
  for (a in child){
    object$variables[[a]]$parents<-object$variables[[a]]$
      parents[object$variables[[a]]$parents!=node]
  }
  object$info$fitted<-F
  if (!is.null(data)){
   object<-fit.bnet(object = object,nodes=child,data,...)
  }
  return(object)
  
}

#' Bnet variables ordering 
#' 
#' @param object a bnet object
#' @return a bnet object with variables ordered
#' @export
orden.bnet<-function(object){
  var<-c()
  old<-variables(object)
  while (length(old)>0){
    if (all(object$variables[[old[1]]]$parents %in% var )){
      var<-c(var,old[1])
      old<-old[-1]
    }
    else{
      old<-c(old[-1],old[1])
    }
  }
  object$variables<-object$variables[var]
  return(object)
}


#' Marginal Computation
#' 
#' @param object an bnet object
#' @return an bnet object with marginal computed by integration
#' @export
compute_marginal.bnet<- function(object){
  object<-orden.bnet(object)
  for (a in variables(object)){
    
    mopC<- object$variables[[a]]$prob
    mop<-new_bmop(knots = mopC$knots[MARGIN],order = mopC$order[MARGIN])
    
    C<-integration_constants(bmop = mopC)
    c<-integration_constants(bmop=mop)
    mop$ctrpoints<-apply(C*mopC$ctrpoints,MARGIN = MARGIN,FUN = sum)
    mop$ctrpoints<-mop$ctrpoints/c
    
  }
  
}