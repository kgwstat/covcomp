% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covcomp.R
\name{ComplnFVE}
\alias{ComplnFVE}
\title{Completion}
\usage{
ComplnFVE(Cov, Omega, FVE = 0.95)
}
\arguments{
\item{Cov}{The partial covariance matrix with same dimensions as Omega.}

\item{Omega}{A matrix of 0s and 1s; indicating the domain of partial covariance.}

\item{FVE}{Fraction of variance explained; 0.95 by default.}
}
\value{
A matrix of the same dimensions: the covariance completion of Cov
}
\description{
Calculates the canonical completion of partial covariance assuming the domain is given by Omega.
}
