\name{hvp}
\alias{hvp}
\title{Interaction term ranking with normalized high order interaction variable importance and deviance criteria.}
\author{
Yifan Sha and Min Lu
}
\usage{
hvp(joint, obj, importance = "permute",
    block = 1, choice = "hivimp")}
\description{
This function ranks the interaction terms in a \emph{rfsrc} object according to normalized high order interaction variable importance and compute other interaction variable importance.
}
\arguments{
\item{joint}{A matrix that each row stores the variable combination in each interaction term.}
\item{obj}{An object of class \emph{(rfsrc, grow)}.}
\item{importance}{Method for computing variable importance (VIMP). It is the same as ``importance" in \emph{rfsrc} function in \emph{randomForestSRC} package. }
\item{block}{Specifies number of trees in a block when calculating VIMP. It is the same as ``block.size" in \emph{vimp} function in \emph{randomForestSRC} package.}
\item{choice}{Method(s) used for ranking interaction terms. Choose ``hivimp" for high order interaction variable importance or/and choose ``acuvimp" for exact high-order interaction deviance.}
}

\value{
A dataframe with ranking criteria for each interaction term in the row and several methods in the column. ``joinvimp" for joint vimp; ``HIvimp" for high order interaction variable importance; ``nmHIvimp" for normalized high order interaction variable importance.
}


\references{
Ishwaran H. (2007). Variable importance in binary regression trees and forests, \emph{Electronic J. Statist.}, 1:519-537.

Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and Lauer M.S. (2010). High-dimensional variable selection for survival data. \emph{J. Amer. Statist. Assoc.}, 105:205-217.

Ishwaran H., Kogalur U.B., Chen X. and Minn A.J. (2011). Random survival forests for high-dimensional data. \emph{Statist. Anal. Data Mining}, 4:115-132.
}

\examples{
data(express)
o.0 <- rfsrc(y~., data = express[,1:11])
joint <- t(combn(1:length(o.0$xvar.names),3))
o <- hvp(joint, obj = o.0, importance = "permute",
         block = 1, choice = c("hivimp", "acuvimp"))
head(o)
}
