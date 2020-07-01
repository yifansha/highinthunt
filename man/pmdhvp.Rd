\name{pmdhvp}
\alias{pmdhvp}
\title{Interaction Term Ranking with PMD criteria}
\author{
Yifan Sha and Min Lu
}
\usage{
pmdhvp(obj, inter, verbose = T)
}
\description{
This function ranks the interaction terms in a \emph{rfsrc} object according to pairwise minimal depth (PMD) matrix.
}
\arguments{
\item{obj}{An object of class \emph{(rfsrc, grow)}.}
\item{inter}{A list of vectors that store combinations in interaction terms.}
\item{verbose}{Set to TRUE for verbose output.}
}

\value{
\item{pmdvp}{PMD variable importance for each interaction terms in \emph{inter}.}
\item{pmd}{PMD matrix calculated from input \emph{obj}.}
}


\references{
Ishwaran H. (2007). Variable importance in binary regression trees and forests, \emph{Electronic J. Statist.}, 1:519-537.

Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and Lauer M.S. (2010). High-dimensional variable selection for survival data. \emph{J. Amer. Statist. Assoc.}, 105:205-217.

Ishwaran H., Kogalur U.B., Chen X. and Minn A.J. (2011). Random survival forests for high-dimensional data. \emph{Statist. Anal. Data Mining}, 4:115-132.
}

\examples{
data(express)
o.0 <- rfsrc(y~., data = express[,1:11])
## computing combinations up to 5 way interaction terms
cmbn <- lapply(2:5,function(i){t(combn(1:length(o.0$xvar.names),i))})
inter <- unlist(lapply(1:length(cmbn), function(i){
     lapply(1:nrow(cmbn[[i]]),function(j){
     matrix(cmbn[[i]][j,],1,length(cmbn[[i]][j,])) })
                      }), recursive = FALSE )

o <- pmdhvp(obj = o.0, inter, verbose = TRUE)
o$pmdvp
## rank(o$pmdvp)[which(rownames(o$pmdvp) == "x1578_x1430_x692_x1223")]
}
