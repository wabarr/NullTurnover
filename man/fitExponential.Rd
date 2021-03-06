% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitExponential.R
\name{fitExponential}
\alias{fitExponential}
\title{Internal Function to Fit exponential model to FADs
Function is to recreate images from Vrba (1995) classic turnover pulse fig}
\usage{
fitExponential(myTree = NULL, FADs = NULL, showExpected = FALSE,
  showTree = FALSE)
}
\arguments{
\item{myTree}{an object of class "phylo". You must include EITHER a value for myTree OR a value for FADs}

\item{FADs}{Vector of taxon first appearances.  You must include EITHER a value for myTree OR a value for FADs}

\item{showExpected}{whether or not to show the plot of the expected versus actual.  Default is FALSE}

\item{showTree}{Whether or not to plot the tree. Default is FALSE}
}
\value{
a Nonlinear Least Squares object of class "nls" representing the exponential fit
}
\description{
Internal Function to Fit exponential model to FADs
Function is to recreate images from Vrba (1995) classic turnover pulse fig
}

