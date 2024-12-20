% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcFrC.R
\name{calcFrC}
\alias{calcFrC}
\title{Calculate Fractional Contribution (FrC)}
\usage{
calcFrC(FinalInt)
}
\arguments{
\item{finalres}{Dataframe resulting from \emph{autoQ} function}

\item{val.to.use}{Type of quantification to transform. Either "Area" or "Maxo"}
}
\value{
A data.frame containing a FrC value for every sample and compound.
}
\description{
Fractional contribution quantifies the contribution of a labeled nutrient to the metabolite of interest.
}
\references{
doi: 10.1016/j.copbio.2015.02.003
}
