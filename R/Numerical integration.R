#' @title Numerical integration with integrate using R
#' @description Numerical integration with integrate using R
#' @param x samples of one distribution
#' @param y samples of one distribution
#' @return the value of the Numerical integration with integrate \code{n}
#' @examples
#' \dontrun{
#' p <- 0
#' q <- 1
#' x<-seq(-10,10,0.1)
#' f<-function(x,p,q){      #density function
#'  1/(q*pi*(1+((x-p)/q)^2))
#'  }
#' n<-length(x)
#' res<-numeric(n)
#' plot(x,cdf(x,p,q))
#' @export
cdf <- function(x,p,q) {
  n<-length(x)
  res<-numeric(n)
  for(i in 1:n){
    res[i]<-integrate(f, lower=-Inf, upper=x[i],
                      rel.tol=.Machine$double.eps^0.25,
                      p=p,q=q)$value
  }
  return(res)
}
