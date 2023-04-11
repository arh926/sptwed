#' Tweedie Compound Poisson-Gamma (CP-g) Double Generalized Linear Spatial Model
#'
#' Fits a Double Generalized Linear Model: \eqn{\log(\mu({\bf s}))=x(\{\bf s})^T\beta+f({\bf s})}^Tw({\bf s}) and \eqn{\log(\phi)=z^T\gamma}. Parameters not listed below are optional. This relies on a co-ordinate descent approach.
#'
#' @param y observed response
#' @param X covariates for the mean model
#' @param Z covariates for the dispersion model
#' @param index incidence matrix
#' @param index.y.0 index for 0 entries
#' @param beta.init initial values for mean coefficients
#' @param gamma.init initial values for dispersion coefficients
#' @param alpha.init initial values for spatial coefficients
#' @param pen_mat penalty matrix---the graph Laplacian
#' @keywords spatial_tweedie
#' @import stats tweedie Matrix
spatial_tweedie <- function(y = NULL,
                            X = NULL,
                            Z = NULL,
                            index = NULL,
                            index.y.0 = NULL,
                            beta.init = NULL,
                            gamma.init = NULL,
                            alpha.init = NULL,
                            pen_mat=NULL,
                            p = 1.5,
                            tol = 1e-3,
                            miter = 1e4,
                            inf = FALSE,
                            saddle.approx = TRUE,
                            verbose = TRUE,
                            p.update = TRUE){
  obj_track <- matrix(0,nrow=miter,ncol=2)
  iter <- 1

  # Theta initial
  beta <- beta.init
  alpha <- alpha.init
  gamma <- gamma.init
  ngamma <- length(gamma)
  sp_id <- match(index, names(alpha))

  Xtb <- crossprod(t(X),beta)+matrix(alpha[sp_id],ncol=1)
  Ztg <- crossprod(t(Z),gamma)

  if(saddle.approx){
    if(sum(pen_mat!=0)==0){
      obj_init <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha,pen_mat)),alpha))
      if(verbose) cat("Objective Fn:: ", obj_init,"\t")

      while(iter <= miter){

        #######################
        # alpha and beta step #
        #######################
        t1 <- exp(((2-p)*Xtb)-Ztg)
        t2 <- y*exp(((1-p)*Xtb)-Ztg)

        # construct combined design matrix
        Xtm <- X*as.vector((2-p)*t1-(1-p)*t2)
        XmX <- crossprod(Xtm,X)
        XmL <- t(apply(t(Xtm),1,function(x) aggregate(x,list(index),sum)[,2]))
        m.a <- aggregate(((2-p)*t1-(1-p)*t2),list(index),sum)[,2]

        scale.mat1 <- as.matrix(rbind(cbind(XmX,XmL),
                                      cbind(t(XmL),diag(as.vector(m.a))+0.01)),
                                nrow=(length(beta)+length(alpha)),
                                ncol=(length(beta)+length(alpha)))

        scale.mat2 <- as.matrix(rbind(cbind(XmX,XmL),
                                      cbind(t(XmL),diag(as.vector(m.a)))),
                                nrow=(length(beta)+length(alpha)),
                                ncol=(length(beta)+length(alpha)))

        # update beta and alpha simultaneously
        update.next <-  crossprod(chol2inv(chol(scale.mat1)),
                                  crossprod(scale.mat2,matrix(c(beta,alpha),ncol=1))-matrix(c(crossprod(X,matrix((t1-t2),ncol=1)),aggregate((t1-t2), list(index),sum)[,2]), ncol=1))
        beta.next <- update.next[1:length(beta),]
        alpha.next <-   update.next[-(1:length(beta)),]; alpha.next <- alpha.next-mean(alpha.next)
        names(alpha.next) <- names(alpha)

        ##############
        # gamma-step #
        ##############
        Xtb <- crossprod(t(X),beta.next)+matrix(alpha.next[sp_id],ncol=1)
        obj_track[iter,1] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
        d <- (y^(2-p)-y*exp((1-p)*Xtb))/(1-p)-(y^(2-p)-exp((2-p)*Xtb))/(2-p); d[y==0] <- exp((2-p)*Xtb[y==0])/(2-p)

        # update gamma
        grad1.gamma <- crossprod(Z,-exp(-Ztg)*d+0.5)
        grad2.gamma <- crossprod(t(crossprod(Z,diag(as.vector(exp(-Ztg)*d)))),Z)
        step <- crossprod(solve(grad2.gamma),grad1.gamma)
        check <- crossprod(t(crossprod(grad1.gamma,grad2.gamma)),grad1.gamma)
        if(check>=0){
          c <- 1; times <- 1
          gamma.next <- gamma-c*step
          Ztg <- crossprod(t(Z),gamma.next)

          # objective function value at Theta (t)
          obj_track[iter,2] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))


          if(iter==1){
            while(obj_track[iter,2]>obj_init){
              c <- c/2
              gamma.next <- gamma-c*step
              Ztg <- crossprod(t(Z),gamma.next)
              obj_track[iter,2] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
              times <- times+1
              if(times>1e2) break
            }
          }else{
            while(obj_track[iter,2]>obj_track[(iter-1),2]){
              c <- c/2
              gamma.next <- gamma-c*step
              Ztg <- crossprod(t(Z),gamma.next)
              obj_track[iter,2] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
              times <- times+1
              if(times>1e2) break
            }
          }
        }else{
          c <- 1; times <- 1
          gamma.next <- gamma+c*step
          Ztg <- crossprod(t(Z),gamma.next)

          # objective function value at Theta (t)
          obj_track[iter,2] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))


          if(iter==1){
            while(obj_track[iter,2]>obj_init){
              c <- c/2
              gamma.next <- gamma+c*step
              Ztg <- crossprod(t(Z),gamma.next)
              obj_track[iter,2] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
              times <- times+1
              if(times>1e2) break
            }
          }else{
            while(obj_track[iter,2]>obj_track[(iter-1),2]){
              c <- c/2
              gamma.next <- gamma+c*step
              Ztg <- crossprod(t(Z),gamma.next)
              obj_track[iter,2] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
              times <- times+1
              if(times>1e2) break
            }
          }
        }
        if((obj_track[iter,2]>obj_track[iter,1])){
          gamma.next <- gamma
          Ztg <- crossprod(t(Z),gamma.next)
          obj_track[iter,2] <- obj_track[iter,1]
        }
        # update p
        if(p.update){
          p <- seq(1.2,1.9,by=1e-3)[which.min(sapply(seq(1.2,1.9,by=1e-3), function(x) saddle_likelihood(y=y,p=x,mu=exp(Xtb),phi=exp(Ztg))))]
          cat("(",p,")","\t")
          obj_track[iter,2] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
        }

        # convergence criteria
        if(iter>1) dif <- (obj_track[(iter-1),2]-obj_track[iter,2])/obj_track[(iter-1),2] else dif <- (obj_init-obj_track[iter,2])/obj_init


        if(dif < tol){
          if(verbose)  cat("Done::",round(dif,6),"\n")
          break
        }else{
          if(verbose) cat(obj_track[iter,2],"(", obj_track[iter,1],")","\t")
          beta <- beta.next
          gamma <- gamma.next
          alpha <- alpha.next
          iter <- iter+1
        }
      }
    }else{
      obj_init <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha,pen_mat)),alpha))
      if(verbose) cat("Objective Fn.:: ", obj_init,"\t")

      while(iter <= miter){
        #######################
        # alpha and beta step #
        #######################
        t1 <- exp(((2-p)*Xtb)-Ztg)
        t2 <- y*exp(((1-p)*Xtb)-Ztg)

        # construct combined design matrix
        Xtm <- X*as.vector(((2-p)*t1)-((1-p)*t2))
        XmX <- crossprod(Xtm,X)
        XmL <- t(apply(t(Xtm),1,function(x) aggregate(x,list(index),sum)[,2]))
        m.a <- aggregate(((2-p)*t1-(1-p)*t2),list(index),sum)[,2]

        scale.mat1 <- as.matrix(rbind(cbind(XmX,XmL),
                                      cbind(t(XmL),diag(as.vector(m.a))+pen_mat)))

        scale.mat2 <- as.matrix(rbind(cbind(XmX,XmL),
                                      cbind(t(XmL),diag(as.vector(m.a)))),
                                nrow=(length(beta)+length(alpha)),
                                ncol=(length(beta)+length(alpha)))

        # update beta and alpha simultaneously
        update.next <-  crossprod(solve(scale.mat1),crossprod(scale.mat2,matrix(c(beta,alpha),ncol=1))-matrix(c(crossprod(X,matrix((t1-t2),ncol=1)),aggregate((t1-t2), list(index),sum)[,2]), ncol=1))
        beta.next <- update.next[1:length(beta),]
        alpha.next <-   update.next[-(1:length(beta)),]; alpha.next <- alpha.next-mean(alpha.next)
        names(alpha.next) <- names(alpha)

        ##############
        # gamma-step #
        ##############
        Xtb <- crossprod(t(X),beta.next)+matrix(alpha.next[sp_id],ncol=1)
        obj_track[iter,1] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
        d <- (y^(2-p)-y*exp((1-p)*Xtb))/(1-p)-(y^(2-p)-exp((2-p)*Xtb))/(2-p); d[y==0] <- exp((2-p)*Xtb[y==0])/(2-p)

        # update gamma
        grad1.gamma <- crossprod(Z,-exp(-Ztg)*d+0.5)
        grad2.gamma <- crossprod(t(crossprod(Z,diag(as.vector(exp(-Ztg)*d)))),Z)
        step <- crossprod(solve(grad2.gamma),grad1.gamma)
        check <- crossprod(t(crossprod(grad1.gamma,grad2.gamma)),grad1.gamma)
        if(check>=0){
          c <- 1; times <- 1
          gamma.next <- gamma-c*step
          Ztg <- crossprod(t(Z),gamma.next)

          # objective function value at Theta (t)
          obj_track[iter,2] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))


          if(iter==1){
            while(obj_track[iter,2]>obj_init){
              c <- c/2
              gamma.next <- gamma-c*step
              Ztg <- crossprod(t(Z),gamma.next)
              obj_track[iter,2] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
              times <- times+1
              if(times>1e2) break
            }
          }else{
            while(obj_track[iter,2]>obj_track[(iter-1),2]){
              c <- c/2
              gamma.next <- gamma-c*step
              Ztg <- crossprod(t(Z),gamma.next)
              obj_track[iter,2] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
              times <- times+1
              if(times>1e2) break
            }
          }
        }else{
          c <- 1; times <- 1
          gamma.next <- gamma+c*step
          Ztg <- crossprod(t(Z),gamma.next)

          # objective function value at Theta (t)
          obj_track[iter,2] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))


          if(iter==1){
            while(obj_track[iter,2]>obj_init){
              c <- c/2
              gamma.next <- gamma+c*step
              Ztg <- crossprod(t(Z),gamma.next)
              obj_track[iter,2] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
              times <- times+1
              if(times>1e2) break
            }
          }else{
            while(obj_track[iter,2]>obj_track[(iter-1),2]){
              c <- c/2
              gamma.next <- gamma+c*step
              Ztg <- crossprod(t(Z),gamma.next)
              obj_track[iter,2] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
              times <- times+1
              if(times>1e2) break
            }
          }
        }
        if((obj_track[iter,2]>obj_track[iter,1])){
          gamma.next <- gamma
          Ztg <- crossprod(t(Z),gamma.next)
          obj_track[iter,2] <- obj_track[iter,1]
        }

        # update p
        if(p.update){
          p <- seq(1.2,1.9,by=1e-3)[which.min(sapply(seq(1.2,1.9,by=1e-3), function(x) saddle_likelihood(y=y,p=x,mu=exp(Xtb),phi=exp(Ztg))))]
          cat("(",p,")","\t")
          obj_track[iter,2] <- saddle_likelihood(y=y,p=p,mu=exp(Xtb),phi=exp(Ztg)) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
        }

        # convergence criteria
        if(iter>1) dif <- (obj_track[(iter-1),2]-obj_track[iter,2])/obj_track[(iter-1),2] else dif <- (obj_init-obj_track[iter,2])/obj_init


        if(dif < tol){
          if(verbose)  cat("Done::",round(dif,6),"\n")
          break
        }else{
          if(verbose) cat(obj_track[iter,2],"(", obj_track[iter,1],")","\t")
          beta <- beta.next
          gamma <- gamma.next
          alpha <- alpha.next
          iter <- iter+1
        }
      }
    }
  }else{
    obj_init <- -sum(log(dtweedie.series(y=y, power = p, mu=exp(Xtb), phi=exp(Ztg)))) + 0.5*as.vector(crossprod(t(crossprod(alpha,pen_mat)),alpha))
    if(verbose) cat("Objective Fn.:: ", obj_init,"\t")

    while(iter <= miter){
      #######################
      # alpha and beta step #
      #######################
      t1 <- exp(((2-p)*Xtb)-Ztg)
      t2 <- y*exp(((1-p)*Xtb)-Ztg)

      # construct combined design matrix
      Xtm <- X*as.vector(((2-p)*t1)-((1-p)*t2))
      XmX <- crossprod(Xtm,X)
      XmL <- t(apply(t(Xtm),1,function(x) aggregate(x,list(index),sum)[,2]))
      m.a <- aggregate(((2-p)*t1-(1-p)*t2),list(index),sum)[,2]

      scale.mat1 <- as.matrix(rbind(cbind(XmX,XmL),
                                    cbind(t(XmL),diag(as.vector(m.a))+pen_mat)),
                              nrow=(length(beta)+length(alpha)),
                              ncol=(length(beta)+length(alpha)))

      scale.mat2 <- as.matrix(rbind(cbind(XmX,XmL),
                                    cbind(t(XmL),diag(as.vector(m.a)))),
                              nrow=(length(beta)+length(alpha)),
                              ncol=(length(beta)+length(alpha)))

      # update beta and alpha simultaneously
      update.next <-  crossprod(solve(scale.mat1),crossprod(scale.mat2,matrix(c(beta,alpha),ncol=1))-matrix(c(crossprod(X,matrix((t1-t2),ncol=1)),aggregate((t1-t2), list(index),sum)[,2]), ncol=1))
      beta.next <- update.next[1:length(beta),]
      alpha.next <-   update.next[-(1:length(beta)),]; alpha.next <- alpha.next-mean(alpha.next)
      names(alpha.next) <- names(alpha)

      ##############
      # gamma-step #
      ##############
      Xtb <- crossprod(t(X),beta.next)+matrix(alpha.next[sp_id],ncol=1)
      obj_track[iter,1] <- -sum(log(dtweedie.series(y=y, power = p, mu=exp(Xtb), phi=exp(Ztg)))) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
      d <- (y*exp((1-p)*Xtb)/(1-p)-exp((2-p)*Xtb)/(2-p)); xi <- (2-p)/(p-1)
      # update gamma
      grad1.gamma <- crossprod(Z,exp(-Ztg)*d+del1afun(y=y,phi=exp(Ztg),p=p, xi=xi,index.y.0 =index.y.0)/afun(y=y,phi=exp(Ztg),p=p, xi=xi,index.y.0 =index.y.0))
      grad2.gamma <- -crossprod(t(crossprod(Z,diag(as.vector(exp(-Ztg)*d+del1afun(y=y,phi=exp(Ztg),p=p, xi=xi,index.y.0 =index.y.0)^2/afun(y=y,phi=exp(Ztg),p=p, xi=xi,index.y.0 =index.y.0)^2+del2afun(y=y,phi=exp(Ztg),p=p, xi=xi,index.y.0 =index.y.0)/afun(y=y,phi=exp(Ztg),p=p, xi=xi,index.y.0 =index.y.0))))),Z)
      step <- crossprod(solve(grad2.gamma),grad1.gamma)
      check <- crossprod(t(crossprod(grad1.gamma,grad2.gamma)),grad1.gamma)
      if(check>=0){
        c <- 1; times <- 1
        gamma.next <- gamma-c*step
        Ztg <- crossprod(t(Z),gamma.next)

        # objective function value at Theta (t)
        obj_track[iter,2] <- -sum(log(dtweedie.series(y=y, power = p, mu=exp(Xtb), phi=exp(Ztg)))) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))

        if(iter==1){
          while(obj_track[iter,2]>obj_init){
            c <- c/2
            gamma.next <- gamma-c*step
            Ztg <- crossprod(t(Z),gamma.next)
            obj_track[iter,2] <- -sum(log(dtweedie.series(y=y, power = p, mu=exp(Xtb), phi=exp(Ztg)))) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
            times <- times+1
            if(times>1e2) break
          }
        }else{
          while(obj_track[iter,2]>obj_track[(iter-1),2]){
            c <- c/2
            gamma.next <- gamma-c*step
            Ztg <- crossprod(t(Z),gamma.next)
            obj_track[iter,2] <- -sum(log(dtweedie.series(y=y, power = p, mu=exp(Xtb), phi=exp(Ztg)))) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
            times <- times+1
            if(times>1e2) break
          }
        }
      }else{
        c <- 1; times <- 1
        gamma.next <- gamma+c*step
        Ztg <- crossprod(t(Z),gamma.next)

        # objective function value at Theta (t)
        obj_track[iter,2] <- -sum(log(dtweedie.series(y=y, power = p, mu=exp(Xtb), phi=exp(Ztg)))) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))

        if(iter==1){
          while(obj_track[iter,2]>obj_init){
            c <- c/2
            gamma.next <- gamma+c*step
            Ztg <- crossprod(t(Z),gamma.next)
            obj_track[iter,2] <- -sum(log(dtweedie.series(y=y, power = p, mu=exp(Xtb), phi=exp(Ztg)))) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
            times <- times+1
            if(times>1e2) break
          }
        }else{
          while(obj_track[iter,2]>obj_track[(iter-1),2]){
            c <- c/2
            gamma.next <- gamma+c*step
            Ztg <- crossprod(t(Z),gamma.next)
            obj_track[iter,2] <- -sum(log(dtweedie.series(y=y, power = p, mu=exp(Xtb), phi=exp(Ztg)))) + 0.5*as.vector(crossprod(t(crossprod(alpha.next,pen_mat)),alpha.next))
            times <- times+1
            if(times>1e2) break
          }
        }
      }
      if((obj_track[iter,2]>obj_track[iter,1])){
        gamma.next <- gamma
        Ztg <- crossprod(t(Z),gamma.next)
        obj_track[iter,2] <- obj_track[iter,1]
      }


      # convergence criteria
      if(iter>1) dif <- (obj_track[(iter-1),2]-obj_track[iter,2])/obj_track[(iter-1),2] else dif <- (obj_init-obj_track[iter,2])/obj_init


      if(dif < tol){
        if(verbose) cat("Done::",round(dif,6),"\n")
        break
      }else{
        if(verbose) cat(obj_track[iter,2],"(", obj_track[iter,1],")","\t")
        beta <- beta.next
        gamma <- gamma.next
        alpha <- alpha.next
        iter <- iter+1
      }
    }
  }

  if(inf){
    se.mean <- sqrt(diag(chol2inv(chol(crossprod(scale.mat1,scale.mat1))))) [1:length(beta)]
    se.disp <- sqrt(diag(chol2inv(chol(crossprod(Z,Z)))))

    # Prepare model object
    mean.summ <- round(cbind(beta.next,se.mean,beta.next/se.mean,pnorm(abs(beta.next/se.mean),lower.tail=F)),4)
    colnames(mean.summ) <- c("Estimates", "Std. Error","Wald z-val","p-val")
    disp.summ <- round(cbind(gamma.next,se.disp,gamma.next/se.disp,pnorm(abs(gamma.next/se.disp),lower.tail=F)),4)
    colnames(disp.summ) <- c("Estimates", "Std. Error","Wald z-val","p-val")

    return(list(optim_pars = list(alpha = alpha.next,
                                  beta = beta.next,
                                  gamma = as.vector(gamma.next),
                                  p=p,
                                  objFn = obj_track[1:iter,]),
                model_summary = list(mean.model = mean.summ,
                                     disp.model = disp.summ)))
  }else{
    return(list(alpha = alpha.next,
                beta = beta.next,
                gamma = as.vector(gamma.next),
                p=p,
                objFn = obj_track[1:iter,]))
  }

}
