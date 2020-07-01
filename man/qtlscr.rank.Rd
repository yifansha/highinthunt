\name{qtlscr.rank}
\alias{qtlscr.rank}
\title{Interactions by 95th percentile of normalized high order interaction variable importance density estimation.}
\author{
Yifan Sha and Min Lu
}
\usage{
qtlscr.rank(o,data.nm,
            qtlcut=0.95,plot=TRUE)}
\description{
This function output the interaction terms of \emph{HIvimp} object by 95th percentile of normalized high order interaction variable importance density estimation. 
}
\arguments{
\item{o}{An object of class \emph{HIvimp}.}
\item{data.nm}{Data name provided by user for output and plot presentation.}
\item{qtlcut}{Quantile score cut to decide final interactions to be presented.}
\item{plot}{Option to generate quantile score plots for output interactions shown in groups, where all other group elements are strict subset of each group identifier in red color. If this option is turned off, only output interactions with quantile scores will be generated.}
}

\value{
A list of interactions with their quantile scores. Only interactions with quantile scores less than (1-qtlcut) will be shown. Plots of output interactions in groups with group elements being strict subset of each group identifier in red color.
}

\examples{
data(express)
formula <- y~.
o <- HIvimp(y~., data = data, method = "hivimp", terms = NULL, PMDorder=TRUE, 
             importance = "permute", block = 1,
             subRF = T, subvars = NULL,
             w.initial = "vimp",
             wt = function(btpmd,digpmd){log(1/digpmd/btpmd)},
             itrSub = NULL, wtSub = TRUE,
             itrWt = 1, subQtl = NULL, max.var = 50,
             nHIvimp = function(all.cmbn){min(length(all.cmbn),150)},
             qtlscr = TRUE,
             verbose = TRUE)
qtlscr.rank(o,data.nm="express")
}