#' Saddlepoint approximation to the Tweedie Compound Poisson-gamma likelihood
#' 
#' @param y response
#' @param p index parameter
#' @param mu mean vector
#' @param phi dispersion vector
#' @keywords saddle_likelihood
#' @import tweedie
#' @examples 
#' \dontrun{
#' require(tweedie)
#' N <- 10
#' mu <- rnorm(N, 10, 2)
#' phi <- rgamma(N, 2, 2)
#' y <- rtweedie(N, 1.5, mu, phi)
#' saddle_likelihood(y=y,p=1.5,mu=mu,phi=phi)
#' }
saddle_likelihood <- function(y = NULL,
                              p = NULL,
                              mu = NULL,
                              phi = NULL,
                              epsilon = 1/6){
  lik_vec <- 1 / phi * ((y^(2 - p) - y * mu^(1 - p)) / (1 - p)-(y^(2 - p) - mu^(2 - p)) / (2 - p)) +
    log(2 * pi * phi * (y + epsilon)^p) / 2
  lik_vec[y == 0] <- mu[y == 0]^(2 - p) / (phi[y == 0] * (2 - p))
  
  return(sum(lik_vec))
}
