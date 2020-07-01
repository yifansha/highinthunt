\name{var.sel}
\alias{var.sel}
\title{Variable Selection for High order Interaction Detection}
\author{
Yifan Sha and Min Lu
}
\usage{
var.sel(formula, data,
subvars = function(data){ceiling((ncol(data)-1)/5)},
w.initial = "vimp",
wt = function(btpmd,digpmd){log(1/btpmd/digpmd)},
itrSub = 5, wtSub = T,
itrWt = 1,
subQtl = function(btpmd,digpmd){
     which(btpmd < quantile(btpmd, probs = 0.1))},
verbose = TRUE,
obj = NULL)
}
\description{
For high order interaction detection in high dimentional data, firstly, random forests using subset of features are grown (subRFs, see subRF function) with output of Pairwise Minimal Depth matrix (subPMD); secondly, features with smaller values (according to the rule in \emph{subQtl}) in their rows of subPMDs are selected. Iteratively growing weighted Random Forests using Pairwise Minimal Depth weights is optional.
}
\arguments{
\item{formula}{A symbolic description of the model to be fit.}
\item{data}{Data frame containing the y-outcome and x-variables.}
\item{subvars}{Number of variables selected for fitting each random forest.}
\item{w.initial}{A probability vector, according to which features are selected (see w0 in subRF function) for each subRF. When equals to ``vimp", then variable importance is used. When equals to ``md", then minimal depth from maximal subtree is used.}
\item{wt}{A function calculating variable weights using PMD matrix where \emph{btpmd}[i] is the
average of \emph{i}th row in the PMD matrix where only off diagonal elements that have smaller values are used and \emph{digpmd}[i] is the \emph{i}th diagonal element in the PMD matrix.}
\item{itrSub}{Number of subRFs to grow.}
\item{wtSub}{logical. Should weighted random forests grown?}
\item{itrWt}{Number of iterations of each weighted subRF.}
\item{subQtl}{Which features will be selected from each subRF. ``digpmd” is the diagonal elements of the PMD matrix and
``btpmd” is the average of \emph{i}th row in the PMD matrix used
where only off diagonal elements that have smaller values are used.}
\item{verbose}{Set to TRUE for verbose output.}
\item{obj}{Inital object of class \emph{(rfsrc, grow)}. New object will be created using the data and formula if it is set to null.}
}

\value{
\item{var.sl}{Names of variables selected.}
\item{var.sl.list}{A list of each subRF's output: the \emph{j}th object is the output of the \emph{j}th subRF with two elements: \emph{sl.w} contains the weights (calculated using input \emph{wt}) of selected variable and \emph{sl.var} contains the names of selected variable. Final result (\emph{var.sl}) is the union of each set \emph{sl.var}.}
\item{itrwt}{A list of each weighted subRF's weights: the \emph{j}th object is the weight path of the \emph{j}th subRF with \emph{itrWt} elements. See the output of \emph{wt.itr} function. }
}


\references{
Ishwaran H. (2007). Variable importance in binary regression trees and forests, \emph{Electronic J. Statist.}, 1:519-537.

Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and Lauer M.S. (2010). High-dimensional variable selection for survival data. \emph{J. Amer. Statist. Assoc.}, 105:205-217.

Ishwaran H., Kogalur U.B., Chen X. and Minn A.J. (2011). Random survival forests for high-dimensional data. \emph{Statist. Anal. Data Mining}, 4:115-132.
}

\examples{
data(express)
o <- var.sel(y~., data = express[,1:200],
             subvars = function(data){ceiling((ncol(data)-1)/2)},
             w.initial = "vimp",
             wt = function(btpmd,digpmd){log(1/btpmd/digpmd)},
             itrSub = 3, wtSub = TRUE,
             itrWt = 2,
             subQtl = function(btpmd,digpmd){
                    which(btpmd < quantile(btpmd, probs = 0.1))},
             verbose = TRUE)
o$var.sl
o$itrwt[[3]][[2]] # variable weights in the second iteration of the 5th subRF
}
