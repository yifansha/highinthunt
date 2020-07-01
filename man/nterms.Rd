\name{nterms}
\alias{nterms}
\title{Highest order interaction terms up to search}
\author{
Yifan Sha and Min Lu
}
\usage{
nterms(obj)
}
\description{
This function decides the highest order interaction terms that a random forest object could potentially provide.
}
\arguments{
\item{obj}{An object of class \emph{(rfsrc, grow)}. Note that when \emph{rfsrc} function is used, set \emph{statistics = TRUE}.}
}

\value{
\item{tree.depth}{Average tree depth. This number will be returned only when n/p>10 that the distribution of minimal depth statistics is less skewed. }
\item{min.depth}{Average minimal depth from all variables. This number will be returned only when n/p>10 that the distribution of minimal depth statistics is less skewed.}
\item{terms}{Highest way of interaction terms to detect. When n/p>10, this equals to the rounded value of tree.depth minus min.depth, which is the average depth of maximal Xi-subtrees. Otherwise, this equals to floor(log2((n)/nodesize)).}
}


\references{
Ishwaran H. (2007). Variable importance in binary regression trees and forests, \emph{Electronic J. Statist.}, 1:519-537.

Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and Lauer M.S. (2010). High-dimensional variable selection for survival data. \emph{J. Amer. Statist. Assoc.}, 105:205-217.

Ishwaran H., Kogalur U.B., Chen X. and Minn A.J. (2011). Random survival forests for high-dimensional data. \emph{Statist. Anal. Data Mining}, 4:115-132.
}

\examples{
data(express)
obj <- rfsrc(y~., data = express[,1:6], statistics = TRUE)
nterms(obj)
}
