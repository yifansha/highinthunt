\name{subRF}
\alias{subRF}
\title{Grows Random Forests with Subset of Features}
\author{
Yifan Sha and Min Lu
}
\usage{
subRF(formula, data, w0,
subvars = function(data){ceiling((ncol(data)-1)/5)},
n.RF = 5,
wtRF = T)
}
\description{
This function grows Random Forests with subset of features. Features are selected according to a probability vector.
}
\arguments{
\item{formula}{A symbolic description of the model to be fit.}
\item{data}{Data frame containing the y-outcome and x-variables.}
\item{w0}{A probability vector, according to which features are selected.}
\item{subvars}{Number of variables selected for fitting each random forest.}
\item{n.RF}{Number of random forests to grow.}
\item{wtRF}{logical. Should weighted random forests grown?}
}

\value{
A list of each random forest's output:
\item{subdata}{Dataframe used for growing the \emph{i}th random forest.}
\item{pmd}{PMD matrix from the \emph{i}th random forest.}
\item{subRF.o}{An object of class \emph{(rfsrc, grow)} of the \emph{i}th random forest.}
}


\references{
Ishwaran H. (2007). Variable importance in binary regression trees and forests, \emph{Electronic J. Statist.}, 1:519-537.

Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and Lauer M.S. (2010). High-dimensional variable selection for survival data. \emph{J. Amer. Statist. Assoc.}, 105:205-217.

Ishwaran H., Kogalur U.B., Chen X. and Minn A.J. (2011). Random survival forests for high-dimensional data. \emph{Statist. Anal. Data Mining}, 4:115-132.
}

\examples{
data(express)
o.0 <- rfsrc(y~., data = express[,1:11])
w.0 <- vimp(o.0,o.0$xvar.names,importance = "permute")
o <- subRF(y~., data = express[,1:11],
           w0 = w.0$importance,
           subvars = function(data){ceiling((ncol(data)-1)/5)},
           n.RF = 5,
           wtRF = TRUE)
o[[2]]
}
