

#' Hill-Climbing structure search
#' 
#' Function that perform a stocastic hill-climbing search in the space 
#' of possible Bayesian networks, danities and conditional densities are 
#' estimated with a bmop penalized-logLik search, package \code{Rbmop}
#' 
#'  @param data data frame of observation
#'  @param bnet optional, a bnet object to start the search from
#'  @param whitelist optional, a matrix or data frame with two columns 
#'         each row indicates an arc
#'  @param blacklist optional, not implemented
#'  @param score a score function to be minimized in the search, usually 
#'  \code{BIC} or\code{AIC}
#'  @param maxp positive integer, the maximum number of parents for each node
#'  @param maxrep positive integer, the maximum number of iterations or jump
#'         in the search space
#'  @param ... additional parameters to be passed to bmop fitting functions
#'  @return A bnet object, the bnet among the ones visited that minimized 
#'         the score. 
#'         @export
hillClimb.bnet<-function(data, bnet=NULL, whitelist=NULL, blacklist=NULL, score=BIC,
                  maxp=dim(data)[2], maxrep=1000, ...  ){
  
  cl<-match.call()
  if (is.null(bnet)){
    bnet <- new_bnet(variables = names(data))
  }
  
  if (!is.null(whitelist)){
   for (i in dim(whitelist)[1]){
     bnet<-add_arc.bnet(object = bnet,from = whitelist[i,1],to = whitelist[i,2]
                        , check = TRUE)
   }
  }
  
  bnet <- fit.bnet(object = bnet,data = data,...)
  finish<-FALSE
  rep <- 0
  #actualscore <- score(object=bnet, data=data) #slower version
  actualscore<-sum(sapply(bnet$variables,FUN = function(nod){
    return(nod$prob$logLik)
  }))
  arcCandidates <- expand.grid(variables(bnet),variables(bnet),
                               stringsAsFactors = FALSE)
  names(arcCandidates)<-c("from","to")
  while (!finish){
   
    rep <- rep + 1
    print(rep)
    
    #fromto <- sample(variables(bnet) , size = 2, replace = F ) #naive version
    idx <- sample(x = 1:(dim(arcCandidates)[1]),size = 1)
    fromto <- arcCandidates[idx,]
    arcCandidates<-arcCandidates[-idx,]
    
    candidateAdd <- add_arc.bnet(object = bnet, from = fromto$from, to=fromto$to, 
                              check = T, data = data,
                              ... )
    
    
    candscoreAdd<-sum(sapply(candidateAdd$variables,FUN = function(nod){
      return(nod$prob$logLik)
    }))
    #candscore1 <- score(object=candidate1 , data=data) #slower version
    
    if (candscoreAdd < actualscore){
      diff <- actualscore - candscoreAdd
      actualscore <- candscoreAdd
      bnet <- candidateAdd
      
    }
    
    if (rep>=length(variables(bnet))^2){ finish <- TRUE}
    if (rep>=maxrep){ finish <- TRUE}
    
   }
  
  bnet$info$model <- "Bayesian Network, Stocastic Hill-Climbing "
  bnet$info$call<- cl
  
  return(bnet)
  
}



