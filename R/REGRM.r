###############################################################################################
# REGRM:  Robust estimator for a generalized ratio model
#          	by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#   Weight function:  		Tukey's biweight function & Huber weight function
#   Scale of quasi-residuals:
#               Average absolute deviation (AAD) or median absolute deviation (MAD)
#---------------------------------------------------------------------------------------------#
#       Implemented by K. Wada (NSTAC, Japan)
#---------------------------------------------------------------------------------------------#
#       Ver. 0.0 [2019.11.12]
#---------------------------------------------------------------------------------------------#
#  Parameters
#   x1       single explanatory variable
#   y1       objective variable to be imputed
#   gm	     gamma value
#              gm="a": gamma=1
#              gm="b": gamma=1/2	(conventional ratio model)
#              gm="c"; gamma=0    (regression model without intercept)
#   wf       weight function (wf=T : Tukey, wf=H : Huber)
#   scale    scale parameter (scale=AAD; scale=MAD)
#   c2       tuning parameter (c1=4, 6, 8) see the table of "Tuning parameter c1" below
#   rp.max   maximum number of iteration (default: rp.max=100)
#   cg.rt    convergence condition to stop iteration (default: cg1=0.001)
#---------------------------------------------------------------------------------------------#
#  Tuning parameter c1
#            More robust <----> Less robust | default setting
#  ------------+--------+--------+----------+------------------
#         c1   |  4     |  6     |   8      |
#  ------------+--------+--------+----------+------------------
#  Tukey  SD   |  5.01  |  7.52  |  10.03   |  10.03
#        AAD   |  4     |  6     |   8      |   8
#        MAD   |  3.38  |  5.07  |   6.76   | * Use the value for SD for mad function in R
#  ------------+--------+--------+----------+------------------
#  Huber  SD   |  1.44  |  2.16  |   2.88   |   2.88
#        AAD   |  1.15  |  1.72  |   2.30   |   2.30
#        MAD   |  0.97  |  1.46  |   1.94   | * Use the value for SD for mad function in R
#---------------------------------------------------------------------------------------------#
#  Return values
#   par      robustly estimated ratio of y1 to x1 (beta)
#   res	     homoscedastic quasi-residuals
#   wt       robust weights
#   rp       total number of iteration
#   s1       changes of the scale (AAD or MAD)
#   efg	     error flag, 1: acalculia (all weights become zero)  0: successful termination
################################################################################################
#
#' @title REGRM: Robust estimator for a generalized ratio model
#'               by iteratively re-weighted least squares (IRLS) algorithm for M-estimation
#'
#' @description REGRAM function integrates 12 functions below included in RrT.r and RrH.r for easy
#'              use for users.  Please note that the values for the tuning parameter c2
#'              allowed in this function is standardized and limited.  If you prefer to use
#'              other values, please use the functions contained in RrT.r and/or RrH.r directly.
#'              The turning parameters of those functions are not standardized.
#'
#' @details Child functions included RrT.r using Tukey's biweight function and RrH.r of Huber
#' weight function.  AAD stands for Average Absolute Deviation and MAD, Median Absolute Deviation.
#'
#' * RrTa.aad: Tukey's biweight function & AAD, gamma=1
#' * RrTb.aad: Tukey's biweight function & AAD, gamma=1/2
#' * RrTc.aad: Tukey's biweight function & AAD, gamma=0
#' * RrTa.mad: Tukey's biweight function & MAD, gamma=1
#' * RrTb.mad: Tukey's biweight function & MAD, gamma=1/2
#' * RrTc.mad: Tukey's biweight function & MAD, gamma=0
#' * RrHa.aad: Huber's weight function & AAD, gamma=1
#' * RrHb.aad: Huber's weight function & AAD, gamma=1/2
#' * RrHc.aad: Huber's weight function & AAD, gamma=0
#' * RrHa.mad: Huber's weight function & MAD, gamma=1
#' * RrHb.mad: Huber's weight function & MAD, gamma=1/2
#' * RrHc.mad: Huber's weight function & MAD, gamma=0
#' @md
#'
#' @param x1 single explanatory variable
#' @param y1 objective variable to be imputed
#' @param gm indication of gamma value as follows:  \cr
#'           gm="a": gamma=1 \cr
#'           gm="b": gamma=1/2	(conventional ratio model) \cr
#'           gm="c"; gamma=0    (regression model without intercept) \cr
#' @param wf weight function (wf=T : Tukey, wf=H : Huber) \cr
#' @param scale scale parameter (scale=AAD; scale=MAD)
#' @param c2 tuning parameter (c2=4, 6, 8) for weight function. Smaller figure is more robust.
#' @param rp.max maximum number of iteration
#' @param cg.rt convergence condition to stop iteration (default: cg1=0.001)
#' @return a list with the following elements
#' \describe{
#'   \item{\code{par}}{robustly estimated ratio of y1 to x1 (beta)}
#'   \item{\code{res}}{homoscedastic quasi-residuals}
#'   \item{\code{wt}}{robust weights}
#'   \item{\code{rp}}{total number of iteration}
#'   \item{\code{s1}}{changes of the scale (AAD or MAD)}
#'   \item{\code{efg}}{error flag. 1: acalculia (all weights become zero)  0: successful termination}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' require(REGRM)
#'
#' x <- seq(1, 10, by=0.1)
#' #e <- rnorm(length(x))
#' e <- rt(length(x), df=3)   # error term following t distribution
#'
#' b <- 2		# true value of slope
#'
#' y1 <- b*x + x*e			# example 1: gamma=1
#' y2 <- b*x + sqrt(x)*e   # example 2: gamma=1/2
#'
#' o1 <- REGRM(x, y1, gm="a")
#' o2 <- REGRM(x, y2, gm="b")
#'
#' o1$par;  o2$par     # estimated slope
#'
#' require(RColorBrewer)
#' cols = brewer.pal(11, "PiYG")
#' cl1 <- round((o1$wt)*10+1)
#' cl2 <- round((o2$wt)*10+1)
#'
#' par(mfrow=c(1,2))
#' plot(x, y1, col=cols[cl1], pch=20)
#' plot(x, y2, col=cols[cl2], pch=20)
#' }
#'
REGRM <- function(x1, y1, gm="b", wf="T", scale="AAD", c2=8, rp.max=100, cg.rt=0.01){

#----------------------------------------------- check arguments
  stopifnot(gm=="a"| gm=="b"| gm=="c")
  stopifnot(wf=="T" | wf=="H")
  stopifnot(scale=="AAD"| scale=="MAD")
  stopifnot(c2==4 | c2==6 | c2==8)

#----------------------------------------------- select an appropriate function
  c3 <- c2/2-1  # (4, 6, 8) => (1, 2, 3)

  if (wf=="T") {
     # source("RrT.r")
     if (scale=="AAD"){

        c1 <- switch(c3, 4, 6, 8)
        if (gm=="a") ot <- RrTa.aad(x1, y1, c1, rp.max, cg.rt)
        if (gm=="b") ot <- RrTb.aad(x1, y1, c1, rp.max, cg.rt)
        if (gm=="c") ot <- RrTc.aad(x1, y1, c1, rp.max, cg.rt)

     } else {

        c1 <- switch(c3, 5.01, 7.52, 10.03)  # tuning constant for SD (MAD in R)
        if (gm=="a") ot <- RrTa.mad(x1, y1, c1, rp.max, cg.rt)
        if (gm=="b") ot <- RrTb.mad(x1, y1, c1, rp.max, cg.rt)
        if (gm=="c") ot <- RrTc.mad(x1, y1, c1, rp.max, cg.rt)

     }
     return(list(par=ot$par, res=ot$res, wt=ot$wt, rp=ot$rp, s1=ot$s1, efg=ot$efg))
  }

  if (wf=="H") {
     # source("RrH.r")
     if (scale=="AAD"){

        c1 <- switch(c3, 1.15, 1.72, 2.30)
        if (gm=="a") ot <- RrHa.aad(x1, y1, c1, rp.max, cg.rt)
        if (gm=="b") ot <- RrHb.aad(x1, y1, c1, rp.max, cg.rt)
        if (gm=="c") ot <- RrHc.aad(x1, y1, c1, rp.max, cg.rt)

     } else {

        c1 <- switch(c3, 1.44, 2.16, 2.88)  # tuning constant for SD (MAD in R)
        if (gm=="a") ot <- RrHa.mad(x1, y1, c1, rp.max, cg.rt)
        if (gm=="b") ot <- RrHb.mad(x1, y1, c1, rp.max, cg.rt)
        if (gm=="c") ot <- RrHc.mad(x1, y1, c1, rp.max, cg.rt)

     }
     return(list(par=ot$par, res=ot$res, wt=ot$wt, rp=ot$rp, s1=ot$s1, efg=ot$efg))
  }

  return(list(par=NA, res=NA, wt=NA, rp=0, s1=NA, efg=1))

}
