% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/6-bnet-classifier.R
\name{classifier_bnet}
\alias{classifier_bnet}
\title{Multinet classifier using bnet}
\usage{
classifier_bnet(formula = NULL, targets, predictors, prior = NULL,
  data = NULL, bnets = NULL, ...)
}
\arguments{
\item{targets}{}

\item{predictors}{}

\item{prior}{string, NULL, or vector of positive numbers}

\item{data}{data set of observations}

\item{bnets}{optionally a list of bnet object, one for each class value}
}
\value{
a classifier_bnet object
}
\description{
Learn a multi net classifier unsing bnet, in particular learn a model
where for every class value a bnet object is estimated.
As particular case naive Bayes (\code{bnets=NULL}).
}

