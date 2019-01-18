#' @title Cramer-von Mises statistic using R
#' @description Cramer-von Mises statistic using R
#' @param x samples of one distribution
#' @param y samples of one distribution
#' @return the value of the Cramer-von Mises statistic \code{n}
#' @examples
#' \dontrun{
#' attach(chickwts)
#' x <- sort(as.vector(weight[feed == "soybean"]))
#' y <- sort(as.vector(weight[feed == "linseed"]))
#' detach(chickwts)
#' w2.0 <- cramer(x,y)
#' print(w2.0)
#' }
#' @export
cramer <-function(x,y){ #compute the Cramer-von Mises statistic
  n <- length(x)
  m <- length(y)
  Fn <- ecdf(x)
  Gm <- ecdf(y)
  W2 <- ((m*n)/(m+n)^2)*
    (sum((Fn(x)-Gm(x))^2)+sum((Fn(y)-Gm(y))^2))
  return(W2)
}
