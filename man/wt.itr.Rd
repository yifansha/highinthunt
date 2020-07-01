\name{wt.itr}
\alias{wt.itr}
\title{Iteratively Grows Weighted Random Forests}
\author{
Yifan Sha and Min Lu
}
\usage{
wt.itr(formula, data,
       pmd.initial,
       obj.initial = NULL,
       wt = function(btpmd,digpmd){log(1/btpmd/digpmd)},
       iteration = 3)
}
\description{
This function iteratively grows Random Forests using pairwise minimal depth (PMD) weights.
}
\arguments{
\item{formula}{A symbolic description of the model to be fit.}
\item{data}{Data frame containing the y-outcome and x-variables.}
\item{pmd.initial}{Inital PMD matrix. Note that variables should be arranged as the same order
as data.}
\item{obj.initial}{Inital object of class \emph{(rfsrc, grow)}. New object will be created using the data and formula if it is set to null.}
\item{wt}{A function calculating variable weights using PMD matrix where \emph{btpmd}[i] is the
average of \emph{i}th row in the PMD matrix where only off diagonal elements that have smaller values are used and \emph{digpmd}[i] is the \emph{i}th diagonal element in the PMD matrix.}
\item{iteration}{Number of iterations.}
}

\value{
A list of each iteration's output:
\item{pmd}{PMD matrix after \emph{j}th iteration (first \emph{p} columns in the dataframe. When there is only
one iteration, this equals to \emph{pmd.initial}.}
\item{digpmd}{Diagonal elements of the PMD matrix used in the \emph{j}th iteration.}
\item{btpmd}{The average of \emph{i}th row in the PMD matrix used in the \emph{j}th iteration,
where only off diagonal elements that have smaller values are used.}
\item{wt}{Variable weights used in the \emph{j}th iteration.}
}


\references{
Ishwaran H. (2007). Variable importance in binary regression trees and forests, \emph{Electronic J. Statist.}, 1:519-537.

Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and Lauer M.S. (2010). High-dimensional variable selection for survival data. \emph{J. Amer. Statist. Assoc.}, 105:205-217.

Ishwaran H., Kogalur U.B., Chen X. and Minn A.J. (2011). Random survival forests for high-dimensional data. \emph{Statist. Anal. Data Mining}, 4:115-132.
}

\examples{
data(express)
o.0 <- rfsrc(y~., data = express[,1:11])
pmd.0 <- find.interaction(o.0, method="maxsubtree", sorted = FALSE, verbose = FALSE)

o <- wt.itr(y~., data = express[,1:11],
            pmd.initial = pmd.0,
            wt = function(btpmd,digpmd){log(1/btpmd/digpmd)},
            iteration = 2)
o[[2]]
}
