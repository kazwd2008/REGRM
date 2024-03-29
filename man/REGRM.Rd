% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/REGRM.r
\name{REGRM}
\alias{REGRM}
\title{REGRM: Robust estimator for a generalized ratio model
by iteratively re-weighted least squares (IRLS) algorithm for M-estimation}
\usage{
REGRM(x1, y1, gm = "b", wf = "T", scale = "AAD", c2 = 8,
  rp.max = 100, cg.rt = 0.01)
}
\arguments{
\item{x1}{single explanatory variable}

\item{y1}{objective variable to be imputed}

\item{gm}{indication of gamma value as follows:  \cr
gm="a": gamma=1 \cr
gm="b": gamma=1/2    (conventional ratio model) \cr
gm="c"; gamma=0    (regression model without intercept) \cr}

\item{wf}{weight function (wf=T : Tukey, wf=H : Huber) \cr}

\item{scale}{scale parameter (scale=AAD; scale=MAD)}

\item{c2}{tuning parameter (c2=4, 6, 8) for weight function. Smaller figure is more robust.}

\item{rp.max}{maximum number of iteration}

\item{cg.rt}{convergence condition to stop iteration (default: cg1=0.001)}
}
\value{
a list with the following elements
\describe{
\item{\code{par}}{robustly estimated ratio of y1 to x1 (beta)}
\item{\code{res}}{homoscedastic quasi-residuals}
\item{\code{wt}}{robust weights}
\item{\code{rp}}{total number of iteration}
\item{\code{s1}}{changes of the scale (AAD or MAD)}
\item{\code{efg}}{error flag. 1: acalculia (all weights become zero)  0: successful termination}
}
}
\description{
REGRAM function integrates 12 functions below included in RrT.r and RrH.r for easy
use for users.  Please note that the values for the tuning parameter c2
allowed in this function is standardized and limited.  If you prefer to use
other values, please use the functions contained in RrT.r and/or RrH.r directly.
The turning parameters of those functions are not standardized.
}
\details{
Child functions included RrT.r using Tukey's biweight function and RrH.r of Huber
weight function.  AAD stands for Average Absolute Deviation and MAD, Median Absolute Deviation.
\itemize{
\item RrTa.aad: Tukey's biweight function & AAD, gamma=1
\item RrTb.aad: Tukey's biweight function & AAD, gamma=1/2
\item RrTc.aad: Tukey's biweight function & AAD, gamma=0
\item RrTa.mad: Tukey's biweight function & MAD, gamma=1
\item RrTb.mad: Tukey's biweight function & MAD, gamma=1/2
\item RrTc.mad: Tukey's biweight function & MAD, gamma=0
\item RrHa.aad: Huber's weight function & AAD, gamma=1
\item RrHb.aad: Huber's weight function & AAD, gamma=1/2
\item RrHc.aad: Huber's weight function & AAD, gamma=0
\item RrHa.mad: Huber's weight function & MAD, gamma=1
\item RrHb.mad: Huber's weight function & MAD, gamma=1/2
\item RrHc.mad: Huber's weight function & MAD, gamma=0
}
}
\examples{
\dontrun{
require(REGRM)

x <- seq(1, 10, by=0.1)
#e <- rnorm(length(x))
e <- rt(length(x), df=3)   # error term following t distribution

b <- 2		# true value of slope

y1 <- b*x + x*e			# example 1: gamma=1
y2 <- b*x + sqrt(x)*e   # example 2: gamma=1/2

o1 <- REGRM(x, y1, gm="a")
o2 <- REGRM(x, y2, gm="b")

o1$par;  o2$par     # estimated slope

require(RColorBrewer)
cols = brewer.pal(11, "PiYG")
cl1 <- round((o1$wt)*10+1)
cl2 <- round((o2$wt)*10+1)

par(mfrow=c(1,2))
plot(x, y1, col=cols[cl1], pch=20)
plot(x, y2, col=cols[cl2], pch=20)
}

}
