
#' Naive Bayes model or BAN model
#' 
#' Naive Bayes structure generator for bnet object
#' @param formula a formula object that specify the model
#' @param targets vector of names of the target variables
#' @param predictors vector of names of the predictor variables
#' @param data optionally a data.frame for bnet fitting
#' @param bnet a bnet object, optionally a bnet over the predictor variables
#' @return a Naive Bayes structure, or a BAN model with \code{bnet} as predictor
#' subgraph
#' @export
naive_bayes.bnet<-function(formula=NULL,targets=NULL,predictors=NULL,data=NULL,
                           bnet=NULL,...){
  cl<-match.call()
  if (inherits(formula,"formula")){
    targets<-as.character(formula[[2]])
    targets<-targets[! (targets %in% c("+","~","-","*") )]
    predictors<-all.vars(formula)[! (all.vars(formula) %in% targets)]
  }
  
  if ( (is.null(targets)) | (is.null(predictors)) ){ return(NULL) }
  
  if (is.null(bnet)){
  NB<-new_bnet(variables = c(targets,predictors))}
  else{
    NB<-bnet
    NB<-add_nodes.bnet(object = NB,nodes = targets)
  }
  for (p in predictors){
    NB$variables[[p]]$parents<-c(NB$variables[[p]]$parents,targets)
  }
  NB$targets<-targets
  NB$predictors<-predictors
  NB$info$model<-"Naive Bayes"
  NB$info$call <- cl
  if (!is.null(data)){
    return(fit.bnet(NB,data,...))
  }
  else{
    return(NB)
  }
}

#' Discriminative model
#' 
#' Discriminative model structure generator (inverso of naive bayes)
#' @param targets vector of names for targets variables
#' @param predictors vector of names for predictor variables
#' @param data
#' @export
discriminative.bnet<-function(targets,predictors,data=NULL,...){
  cl<-match.call()
  BN<-new_bnet(variables = c(targets,predictors))
  for (a in targets){
    BN$variables[[a]]$parents<-predictors
  }
  BN$targets<-targets
  BN$predictors<-predictors
  BN$info$model<-"Discriminative regression"
  BN$info$call<-cl
  if (!is.null(data)){
    return(fit.bnet(BN,data))
  }
  else{
    return(BN)
  }
}

#' wrapper greedy search for structure
#' 
#' wrapper greedy search based on evaluation of training-error for general funciotn that define a structure
#' @param targets vector of names for target variables
#' @param predictors vector of names for target variables
#' @param data data.frame of observations
#' @param modelfun function, function that define bnet 
#'          must accept data, targets and predictors as input and output a bnet
#' @param bnet bnet object, starting strucure or NULL     
#' @param evalfun string, one of "MSE" or "MAE", or other string 
#'        naming a scoring function
#' @param exit logica, if \code{TRUE} teh search stops at first not improvement
#'         of the score    
wrapper<-function(targets,predictors,data,modelfun=naive_bayes.bnet,bnet=NULL,
                  evalfun='MSE',exit=F){
  if (is.null(bnet)){
   bnet<-modelfun(targets,predictors=NULL,data=NULL)
   pred<-predictors
  }
  else{
  pred<-predictors[!(predictors %in% bnet$predictors)] 
  bnet$targets<-targets
  }
  bnet<-fit.bnet(object = bnet,data = data)
  value<-predict.bnet(object = bnet,newdata = data)
  arg <- list(x = value, y = data[,targets])
  eval <- do.call(evalfun, arg)
  for (p in predictors){
    bnetnew<-modelfun(targets,predictors<-unique(c(bnet$predictors,p)),
                      data=data)
    valuenew<-predict.bnet(object = bnetnew,newdata = data)
    args <- list(x=valuenew, y=data[,targets])
    evalnew <- do.call(evalfun, args)    
    if (evalnew<eval){
      bnet<-bnetnew
      eval<-evalnew
    }
    else{
      if (exit){
        break
      }
    }
    
  }
  bnet$info$wrapper$state<-T
  bnet$info$wrapper$eval<-eval
  bnet$info$wrapper$evalfun<-evalfun
  bnet$info$wrapper$exit<-exit
  return(bnet)
}

