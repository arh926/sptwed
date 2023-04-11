#' Internal function
#' 
#' Not for external use
#' 
#' @param y response
#' @param phi dispersion vector
#' @param p index parameter
#' @param xi index parameter
#' @param index.y.0 index for 0 values in response
#' @keywords del2afun
del2afun <- function(y=NULL,
                     phi=NULL,
                     p=NULL,
                     xi=NULL,
                     index.y.0=NULL,
                     epsilon=1e-18){
  res <- rep(0,length(y)) # 1
  res[!index.y.0] <- apply(cbind(y[!index.y.0],phi[!index.y.0]),1, function(x){
    k=x[1]^(2-p)/((2-p)*x[2])
    w_max <- k^2*(x[1]/(p-1))^(k*xi)/(2-p)^k/factorial(k)/factorial(k*xi-1)/x[2]^(k*(1+xi))
    k_min <- 1
    while(k_min^2*(x[1]/(p-1))^(k_min*xi)/(2-p)^k_min/factorial(k_min)/factorial(k_min*xi-1)/x[2]^(k_min*(1+xi))<epsilon*w_max) k_min <- k_min+1
    k_max <- k+1
    while(k_max^2*(x[1]/(p-1))^(k_max*xi)/(2-p)^k_max/factorial(k_max)/factorial(k_max*xi-1)/x[2]^(k_max*(1+xi))>epsilon*w_max) k_max <- k_max+1
    return(sum((1+xi)^2*(k_min:k_max)^2*(x[1]/(p-1))^((k_min:k_max)*xi)/(2-p)^(k_min:k_max)/factorial((k_min:k_max))/factorial((k_min:k_max)*xi-1)/x[2]^((k_min:k_max)*(1+xi)))/x[1])
  })
  res
}


