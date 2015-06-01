

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
#'  @param mle logical, if use mle estimation of bmop densities       
#'  @param search logical, if use greedy penalized logLik search
#'  @param ... additional parameters to be passed to bmop fitting functions
#'  @return A bnet object, the bnet among the ones visited that minimized 
#'         the score. 
hc.bnet<-function(data, bnet=NULL, whitelist=NULL, blacklist=NULL, score=BIC,
                  maxp=dim(data)[2], maxrep=1000, mle=F, search=T, ...  ){
  
  cl<-match.call()
  if (is.null(bnet)){
    bnet <- new_bnet(variables = names(data))
  }
  
  if (!is.null(whitelist)){
   for (i in dim(whitelist)[1]){
     bnet<-add_arc.bnet(object = bnet,from = whitelist[i,1],to = whitelist[i,2]
                        , check = TRUE, mle=FALSE, search=TRUE)
   }
  }
  
  bnet <- fit.bnet(object = bnet,data = data,mle = mle,search=search,...)
  finish<-FALSE
  rep <- 0
  actualscore <- score(object=bnet, data=data)
  #actualarcs <- arcs.bnet(bnet)
  while (!finish){
   
    rep <- rep + 1
    print(rep)
    
    fromto <- sample(variables(bnet) , size = 2, replace = F )
    
    candidate <- add_arc.bnet(object = bnet, from = fromto[1], to=fromto[2], 
                              check = T, data = data, mle=mle, search=search,
                              ... )
    
    
    candscore <- score(object=candidate , data=data)
    
    if (candscore < actualscore){
      diff <- actualscore - candscore
      
      actualscore <- candscore
      bnet <- candidate
    }
    
    
    if (rep>maxrep){ finish <- TRUE}
    
   }
  
  bnet$info$model <- "Bayesian Network, Stocastic Hill-Climbing "
  bnet$info$call<- cl
  
  return(bnet)
  
}



