# ============================================================
# Basic Functions for OCRL and Kernel
#
# This file contains the core functions used in the proposed Online Change Rate Learning (OCRL) framework, including:
# 1. Online local polynomial estimators;
# 2. Offline kernel estimators;
# 3. Bandwidth-related pilot estimation functions;
# 4. Data generation utilities;
# 5. Kernel and auxiliary utility functions.
#
# The functions support both:
#   d = 1 : the mean derivative estimation;
#   d = 2 : the derivative estimation of the covariance .
# ============================================================
# FUNCTION: online_LL
  # Purpose: Perform the online local linear regression for streaming functional data.
  # Inputs:
  #   x: Numeric vector (if d=1) or matrix of dimension n x d (if d=2). 
  #   y: Numeric vector of responses with length n. 
  #   eval: Numeric vector (length EV) or matrix (EV x 2) of evaluation points.
  #   h: Numeric. Initial bandwidth parameter.
  #   L: Integer. Number of candidate bandwidths.
  #   res_list: List. Contain the online sufficient statistics and centroids.
  #   N: Integer. The total number of measurements for all the subjects.
  #   n: Integer. The number of subjects.
  #   d: Integer (1 or 2). Dimensionality of the covariate x.
  # Output: Updated res_list containing the updated sufficient statistics 
  #         and the updated candidate bandwidth centroids for online estimation based on local linear regression.
# ============================================================
{
  online_LL <- function(x, y, eval, h, L, res_list, N, n, d){
    
    eta <- sapply(1:L, function(l){ ((L-l+1) / L) ^ (1/(d+4)) * h})
    
    if(K>1){
      idx <- sapply(1:L,function(l){which.min(abs(eta[l] - res_list$centroids))})
    }else{
      idx <- 1:L
    }
    
    res_list$centroids <- (res_list$centroids[idx] * (N-n) + eta * n) / N
    
    if(d==1){
      
      EV <- length(eval)
      
      for(l in 1:L){
        
        Pnew <- array(0, dim = c(2,2,EV)); qnew <- matrix(0,2,EV)
        
        for(i in 1:EV){
          side <- cbind(1, x - eval[i])
          K_vec <- Epan((x - eval[i])/eta[l])/eta[l]
          Pnew[,,i] <- matrix(c(
            sum(K_vec*side[,1]^2), sum(K_vec*side[,1]*side[,2]),
            sum(K_vec*side[,2]*side[,1]),  sum(K_vec*side[,2]^2)
          ),2,2) / n
          qnew[,i] <- matrix(c(
            sum(K_vec*side[,1]*y), sum(K_vec*side[,2]*y)
          ),2,1) / n
        }
        
        res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
        res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
        
      }
    }else{
      
      EV <- nrow(eval)
      
      for(l in 1:L){
        
        Pnew <- array(0, dim = c(3,3,EV)); qnew <- matrix(0,3,EV)
        
        for(i in 1:EV){
          side <- cbind(1, x[,1] - eval[i,1], x[,2] - eval[i,2])
          K_vec <- (Epan((x[,1] - eval[i,1])/eta[l])/eta[l]
                    *Epan((x[,2] - eval[i,2])/eta[l])/eta[l])
          Pnew[,,i] <- matrix(c(
            sum(K_vec*side[,1]^2), sum(K_vec*side[,1]*side[,2]), sum(K_vec*side[,1]*side[,3]),
            sum(K_vec*side[,2]*side[,1]), sum(K_vec*side[,2]^2), sum(K_vec*side[,2]*side[,3]),
            sum(K_vec*side[,3]*side[,1]), sum(K_vec*side[,3]*side[,2]), sum(K_vec*side[,3]^2)
          ),3,3) / n
          qnew[,i] <- matrix(c(
            sum(K_vec*side[,1]*y), sum(K_vec*side[,2]*y),
            sum(K_vec*side[,3]*y)
          ),3,1) / n
        }
        
        res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
        res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
      }
    }
    
    return(res_list)
  }
# ============================================================
# FUNCTION: online_LQuad2
#
# Purpose: Perform the online local quadratic regression for streaming functional data. 
# Inputs: x, y, eval, h, L, res_list, N, n, d: Same as online_LL.
# Output: Updated res_list containing the updated sufficient statistics 
#         and the updated candidate bandwidth centroids for online estimation based on local quadratic regression.
# ============================================================
  online_LQuad2 <- function(x, y, eval, h, L, res_list, N, n, d){
    
    eta <- sapply(1:L, function(l){ ((L-l+1) / L) ^ (1/(6+d)) * h}) 
    
    if(K>1){
      idx <- sapply(1:L,function(l){which.min(abs(eta[l] - res_list$centroids))})
    }else{
      idx <- 1:L
    }
    
    res_list$centroids <- (res_list$centroids[idx] * (N-n) + eta * n) / N
    
    if(d==1){ 
      
      EV <- length(eval)
      
      for(l in 1:L){
        
        Pnew <- array(0, dim = c(3,3,EV)); qnew <- matrix(0,3,EV)
        
        for(i in 1:EV){
          
          side <- cbind(1, x - eval[i], (x - eval[i])^2)
          K_vec <- Epan((x - eval[i])/eta[l])/eta[l]
          for(nr in 1:3){
            for(nc in 1:3){
              Pnew[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/n
            }
            qnew[nr,i] <-  sum(K_vec*side[,nr]*y)/n
          }
        }
        
        res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
        res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
      }
    }else{
      
      EV <- nrow(eval)
      
      for(l in 1:L){
        Pnew <- array(0, dim = c(6,6,EV)); qnew <- matrix(0,6,EV)
        for(i in 1:EV){
          side <- cbind(1, x[,1] - eval[i,1], x[,2] - eval[i,2],
                        (x[,1] - eval[i,1])^2,
                        (x[,1] - eval[i,1])*(x[,2] - eval[i,2]), 
                        (x[,2] - eval[i,2])^2)# dim: N*6
          K_vec <- (Epan((x[,1] - eval[i,1])/eta[l])/eta[l]
                    *Epan((x[,2] - eval[i,2])/eta[l])/eta[l])
          for(nr in 1:6){
            for(nc in 1:6){
              Pnew[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/n
            }
            qnew[nr,i] <-  sum(K_vec*side[,nr]*y)/n
          }
        }
        
        res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
        res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
      }
    }
    
    return(res_list)
  }
# ============================================================
# FUNCTION: online_LCub3
#
# Purpose: Perform the online local cubic regression for streaming functional data. 
# Inputs: x, y, eval, h, L, res_list, N, n, d: Same as online_LL.
# Output: Updated res_list containing the updated sufficient statistics 
#         and the updated candidate bandwidth centroids for online estimation based on local cubic regression.
# ============================================================
  online_LCub3 <- function(x, y, eval, h, L, res_list, N, n, d){
    
    eta <- sapply(1:L, function(l){ ((L-l+1) / L) ^ (1/(d+8)) * h})
    
    if(K>1){
      idx <- sapply(1:L,function(l){which.min(abs(eta[l] - res_list$centroids))})
    }else{
      idx <- 1:L
    }
    
    res_list$centroids <- (res_list$centroids[idx] * (N-n) + eta * n) / N
    
    if(d==1){
      
      EV <- length(eval)
      
      for(l in 1:L){
        
        Pnew <- array(0, dim = c(4,4,EV)); qnew <- matrix(0,4,EV)
        
        for(i in 1:EV){
          side <- cbind(1, x - eval[i], 
                        (x - eval[i])^2,
                        (x - eval[i])^3)
          K_vec <- Epan((x - eval[i])/eta[l])/eta[l]
          for(nr in 1:4){
            for(nc in 1:4){
              Pnew[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/n
            }
            qnew[nr,i] <-  sum(K_vec*side[,nr]*y)/n
          }
        }
        
        res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
        res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
        
      }
    }else{
      
      EV <- nrow(eval)
      
      for(l in 1:L){
        
        Pnew <- array(0, dim = c(10,10,EV)); qnew <- matrix(0,10,EV)
        
        for(i in 1:EV){
          side <- cbind(1, x[,1] - eval[i,1], x[,2] - eval[i,2],
                        (x[,1] - eval[i,1])^2,
                        (x[,1] - eval[i,1])*(x[,2] - eval[i,2]), 
                        (x[,2] - eval[i,2])^2,
                        (x[,1] - eval[i,1])^3, 
                        (x[,1] - eval[i,1])^2*(x[,2] - eval[i,2]),
                        (x[,1] - eval[i,1])*(x[,2] - eval[i,2])^2,
                        (x[,2] - eval[i,2])^3)
          K_vec <- (Epan((x[,1] - eval[i,1])/eta[l])/eta[l]
                    *Epan((x[,2] - eval[i,2])/eta[l])/eta[l])
          for(nr in 1:10){
            for(nc in 1:10){
              Pnew[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/n
            }
            qnew[nr,i] <-  sum(K_vec*side[,nr]*y)/n
          }
        }
        
        res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
        res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
      }
    }
    
    return(res_list)
  }
# ============================================================
# FUNCTION: online_LCub3
#
# Purpose: Perform the online local quartic regression for streaming functional data. 
# Inputs: x, y, eval, h, L, res_list, N, n, d: Same as online_LL.
# Output: Updated res_list containing the updated sufficient statistics 
#         and the updated candidate bandwidth centroids for online estimation based on local quartic regression.
# ============================================================
  online_LQuar4 <- function(x, y, eval, h, L, res_list, N, n, d){
    
    eta <- sapply(1:L, function(l){ ((L-l+1) / L) ^ (1/(d+8)) * h})
    
    if(K>1){
      idx <- sapply(1:L,function(l){which.min(abs(eta[l] - res_list$centroids))})
    }else{
      idx <- 1:L
    }
    
    res_list$centroids <- (res_list$centroids[idx] * (N-n) + eta * n) / N
    
    if(d==1){
      
      EV <- length(eval)
      
      for(l in 1:L){
        
        Pnew <- array(0, dim = c(5,5,EV)); qnew <- matrix(0,5,EV)
        
        for(i in 1:EV){
          side <- cbind(1, x - eval[i], 
                        (x - eval[i])^2,
                        (x - eval[i])^3, 
                        (x - eval[i])^4)
          K_vec <- Epan((x - eval[i])/eta[l])/eta[l]
          for(nr in 1:5){
            for(nc in 1:5){
              Pnew[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/n
            }
            qnew[nr,i] <-  sum(K_vec*side[,nr]*y)/n
          }
        }
        
        res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
        res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
        
      }
    }else{
      
      EV <- nrow(eval)
      
      for(l in 1:L){
        
        Pnew <- array(0, dim = c(15,15,EV)); qnew <- matrix(0,15,EV)
        
        for(i in 1:EV){
          side <- cbind(1, x[,1] - eval[i,1], x[,2] - eval[i,2],
                        (x[,1] - eval[i,1])^2,
                        (x[,1] - eval[i,1])*(x[,2] - eval[i,2]), 
                        (x[,2] - eval[i,2])^2,
                        (x[,1] - eval[i,1])^3, 
                        (x[,1] - eval[i,1])*(x[,2] - eval[i,2])^2,
                        (x[,1] - eval[i,1])^2*(x[,2] - eval[i,2]),
                        (x[,2] - eval[i,2])^3,
                        (x[,1] - eval[i,1])^4, 
                        (x[,1] - eval[i,1])^3*(x[,2] - eval[i,2]),
                        (x[,1] - eval[i,1])^2*(x[,2] - eval[i,2])^2,
                        (x[,1] - eval[i,1])*(x[,2] - eval[i,2])^3,
                        (x[,2] - eval[i,2])^4)
          K_vec <- (Epan((x[,1] - eval[i,1])/eta[l])/eta[l]
                    *Epan((x[,2] - eval[i,2])/eta[l])/eta[l])
          for(nr in 1:15){
            for(nc in 1:15){
              Pnew[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/n
            }
            qnew[nr,i] <-  sum(K_vec*side[,nr]*y)/n
          }
        }
        
        res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
        res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
      }
    }
    
    return(res_list)
  }
# ============================================================
# FUNCTION: batch_LL
# Purpose: Perform the offline local linear regression for streaming functional data.
# Inputs: x, y, eval, h, L, res_list, N, n, d: Same as online_LL.
# Output: A list containing the density estimation and the mean/covariance function estimation at the evaluation points.
# ============================================================
  batch_LL <- function(x, y, eval, h, N, d){
    
    if(d==1){
      
      EV <- length(eval)
      
      P <- array(0, dim = c(2,2,EV)); q <- matrix(0,2,EV)
      
      for(i in 1:EV){
        
        side <- cbind(1, x - eval[i])
        K_vec <- Epan((x - eval[i])/h)/h
        P[,,i] <- matrix(c(
          sum(K_vec*side[,1]^2), sum(K_vec*side[,1]*side[,2]),
          sum(K_vec*side[,2]*side[,1]),  sum(K_vec*side[,2]^2)
        ),2,2) / N
        q[,i] <-  matrix(c(
          sum(K_vec*side[,1]*y), sum(K_vec*side[,2]*y)
        ),2,1) / N
        
      }
      
      den <- sapply(1:EV, function(i){P[1,1,i]})
      est <- sapply(1:EV, function(i){
        (solve(P[,,i]+diag(1e-12,2)) %*% matrix(q[,i],2,1))[1]
      })
      estder <- sapply(1:EV, function(i){
        (solve(P[,,i]+diag(1e-12,2)) %*% matrix(q[,i],2,1))[2]
      })
      
    }else{
      
      EV <- nrow(eval)
      
      P <- array(0, dim = c(3,3,EV)); q <- matrix(0,3,EV)
      
      for(i in 1:EV){
        side <- cbind(1, x[,1] - eval[i,1], x[,2] - eval[i,2])
        K_vec <- (Epan((x[,1] - eval[i,1])/h)/h
                  *Epan((x[,2] - eval[i,2])/h)/h)
        P[,,i] <- matrix(c(
          sum(K_vec*side[,1]^2), sum(K_vec*side[,1]*side[,2]), sum(K_vec*side[,1]*side[,3]),
          sum(K_vec*side[,2]*side[,1]), sum(K_vec*side[,2]^2), sum(K_vec*side[,2]*side[,3]),
          sum(K_vec*side[,3]*side[,1]), sum(K_vec*side[,3]*side[,2]), sum(K_vec*side[,3]^2)
        ),3,3) / N
        q[,i] <- matrix(c(
          sum(K_vec*side[,1]*y), sum(K_vec*side[,2]*y),
          sum(K_vec*side[,3]*y)
        ),3,1) / N
      }
      
      den <- sapply(1:EV, function(i){P[1,1,i]})
      est <- sapply(1:EV, function(i){
        (solve(P[,,i]+diag(1e-12,3)) %*% matrix(q[,i],3,1))[1]
      }) 
      estder <- sapply(1:EV, function(i){
        (solve(P[,,i]+diag(1e-12,3)) %*% matrix(q[,i],3,1))[2]
      }) 
    }
    
    res<- list(den,est,estder)
    names(res) <- c('den','est','estder')
    return(res)
  }
# ============================================================
# FUNCTION: batch_LQuad2
# Purpose: Perform the offline local quadratic regression for streaming functional data.
# Inputs: x, y, eval, h, N, d: Same as batch_LL.
# Output: A list containing the density estimation and the first-order derivative estimation at the evaluation points.
# ============================================================
  batch_LQuad2 <- function(x, y, eval, h, N, d){
    
    if(d==1){
      
      EV <- length(eval)
      
      P <- array(0, dim = c(3,3,EV)); q <- matrix(0,3,EV)
      
      for(i in 1:EV){
        
        side <- cbind(1, x - eval[i],(x - eval[i])^2)
        K_vec <- Epan((x - eval[i])/h)/h
        for(nr in 1:3){
          for(nc in 1:3){
            P[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/N
          }
          q[nr,i] <-  sum(K_vec*side[,nr]*y)/N
        }
      }
      
      den <- sapply(1:EV, function(i){P[1,1,i]})
      est <- sapply(1:EV, function(i){
        (solve(P[,,i]+diag(1e-12,3)) %*% matrix(q[,i],3,1))[2]
      })
      #sec_der <- sapply(1:EV, function(i){
      # 2*(solve(P[,,i]+diag(1e-12,3)) %*% matrix(q[,i],3,1))[3]
      #})
      
    }else{
      
      EV <- nrow(eval)
      
      P <- array(0, dim = c(6,6,EV)); q <- matrix(0,6,EV)
      
      for(i in 1:EV){
        side <- cbind(1, x[,1] - eval[i,1], x[,2] - eval[i,2], (x[,1] - eval[i,1])^2,
                      (x[,1] - eval[i,1])*(x[,2] - eval[i,2]), 
                      (x[,2] - eval[i,2])^2)
        K_vec <- (Epan((x[,1] - eval[i,1])/h)/h
                  *Epan((x[,2] - eval[i,2])/h)/h)
        for(nr in 1:6){
          for(nc in 1:6){
            P[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/N
          }
          q[nr,i] <-  sum(K_vec*side[,nr]*y)/N
        }
      }
      
      den <- sapply(1:EV, function(i){P[1,1,i]})
      est <- sapply(1:EV, function(i){
        (solve(P[,,i]+diag(1e-12,6)) %*% matrix(q[,i],6,1))[2]
      }) 
      #sec_der <- sapply(1:EV, function(i){
      #  2*(solve(P[,,i]+diag(1e-12,6)) %*% matrix(q[,i],6,1))[c(4,6)]
      #}) 
      #}
    }
    res<- list(den,est)#,sec_der)
    names(res) <- c('den','est')#,'sec_der')
    return(res)
  }
# ============================================================
# FUNCTION: batch_LCub3
# Purpose: Perform the offline local cubic regression for streaming functional data.
# Inputs: x, y, eval, h, N, d: Same as batch_LL.
# Output: A list containing the density estimation, the second-order derivative estimation, 
#         and the pilot estimation used in the plug-in bandwidth selection procedure at the evaluation points. 
# ============================================================
  batch_LCub3 <- function(x, y, eval, h, N, d){
    
    if(d==1){ 
      
      EV <- length(eval)
      
      P <- array(0, dim = c(4,4,EV)); q <- matrix(0,4,EV)
      
      for(i in 1:EV){
        side <- cbind(1, x - eval[i], (x - eval[i])^2,(x - eval[i])^3)
        K_vec <- Epan((x - eval[i])/h)/h
        for(nr in 1:4){
          for(nc in 1:4){
            P[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/N
          }
          q[nr,i] <-  sum(K_vec*side[,nr]*y)/N
        }
      }
      den <- sapply(1:EV, function(i){P[1,1,i]})
      sec_der <- sapply(1:EV, function(i){
        2*(solve(P[,,i]+diag(1e-12,4)) %*% matrix(q[,i],4,1))[3]
      })
      theta <- mean(P[1,1,6:95] * sec_der[6:95]^2)
      #thi_der <- sapply(1:EV, function(i){
      #  6*(solve(P[,,i]+diag(1e-12,4)) %*% matrix(q[,i],4,1))[4]
      #})
    }else{
      
      EV <- nrow(eval)
      
      P <- array(0, dim = c(10,10,EV)); q <- matrix(0,10,EV)
      for(i in 1:EV){
        side <- cbind(1, x[,1] - eval[i,1], x[,2] - eval[i,2],
                      (x[,1] - eval[i,1])^2,
                      (x[,1] - eval[i,1])*(x[,2] - eval[i,2]), 
                      (x[,2] - eval[i,2])^2,
                      (x[,1] - eval[i,1])^3, 
                      (x[,1] - eval[i,1])^2*(x[,2] - eval[i,2]),
                      (x[,1] - eval[i,1])*(x[,2] - eval[i,2])^2,
                      (x[,2] - eval[i,2])^3)
        K_vec <- (Epan((x[,1] - eval[i,1])/h)/h
                  *Epan((x[,2] - eval[i,2])/h)/h)
        for(nr in 1:10){
          for(nc in 1:10){
            P[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/N
          }
          q[nr,i] <-  sum(K_vec*side[,nr]*y)/N
        }
      }
      den <- sapply(1:EV, function(i){P[1,1,i]})
      sec_der <- sapply(1:EV, function(i){
        2*(sum((solve(P[,,i]+diag(1e-12,10)) %*% matrix(q[,i],10,1))[c(4,6)]))
      })
      sec_der <- matrix(sec_der, EV2, EV2)
      den <- matrix(P[1,1,], EV2, EV2)
      theta <-  mean(den[3:48,3:48] * sec_der[3:48,3:48]^2)
      #thi_der <- sapply(1:EV, function(i){
      #  6*(sum((solve(P[,,i]+diag(1e-12,10)) %*% matrix(q[,i],10,1))[c(7,10)]))
      #})
    }
    
    #return(theta)
    res<- list(den,sec_der,theta)
    names(res) <- c('den','sec_der','theta')
    return(res)
  }
# ============================================================
# FUNCTION: batch_LQuar4
# Purpose: Perform the offline local quartic regression for streaming functional data.
# Inputs: x, y, eval, h, N, d: Same as batch_LL.
# Output: A list containing the density estimation, the third-order derivative estimation,
#         and the pilot estimation used in the plug-in bandwidth selection procedure at the evaluation points. 
# ============================================================
  batch_LQuar4 <- function(x, y, eval, h, N, d){
    
    if(d==1){ 
      
      EV <- length(eval)
      
      P <- array(0, dim = c(5,5,EV)); q <- matrix(0,5,EV)
      
      for(i in 1:EV){
        side <- cbind(1, x - eval[i], (x - eval[i])^2,(x - eval[i])^3,(x - eval[i])^4)
        K_vec <- Epan((x - eval[i])/h)/h
        for(nr in 1:5){
          for(nc in 1:5){
            P[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/N
          }
          q[nr,i] <-  sum(K_vec*side[,nr]*y)/N
        }
      }
      
      thi_der <- sapply(1:EV, function(i){
        6*(solve(P[,,i]+diag(1e-12,5)) %*% matrix(q[,i],5,1))[4]
      })
      
      index1 <- 6 * (EV==EV1) + 3 * (EV==EV2) 
      index2 <- 95 * (EV==EV1) + 48 * (EV==EV2) 
      
      theta <- mean(P[1,1,index1:index2] * thi_der[index1:index2]^2)
      den <- sapply(1:EV, function(i){P[1,1,i]})
      
    }else{
      
      EV <- nrow(eval)
      
      P <- array(0, dim = c(15,15,EV)); q <- matrix(0,15,EV)
      for(i in 1:EV){
        side <- cbind(1, x[,1] - eval[i,1], x[,2] - eval[i,2],
                      (x[,1] - eval[i,1])^2,
                      (x[,1] - eval[i,1])*(x[,2] - eval[i,2]), 
                      (x[,2] - eval[i,2])^2,
                      (x[,1] - eval[i,1])^3, 
                      (x[,1] - eval[i,1])*(x[,2] - eval[i,2])^2,
                      (x[,1] - eval[i,1])^2*(x[,2] - eval[i,2]),
                      (x[,2] - eval[i,2])^3,
                      (x[,1] - eval[i,1])^4, 
                      (x[,1] - eval[i,1])^3*(x[,2] - eval[i,2]),
                      (x[,1] - eval[i,1])^2*(x[,2] - eval[i,2])^2,
                      (x[,1] - eval[i,1])*(x[,2] - eval[i,2])^3,
                      (x[,2] - eval[i,2])^4)
        K_vec <- (Epan((x[,1] - eval[i,1])/h)/h
                  *Epan((x[,2] - eval[i,2])/h)/h)
        for(nr in 1:15){
          for(nc in 1:15){
            P[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/N
          }
          q[nr,i] <-  sum(K_vec*side[,nr]*y)/N
        }
      }
      
      thi_der <- sapply(1:EV, function(i){
        6*(sum((solve(P[,,i]+diag(1e-12,15)) %*% matrix(q[,i],15,1))[c(7,8)]))
      })
      thi_der <- matrix(thi_der, EV2, EV2)
      den <- matrix(P[1,1,], EV2, EV2)
      theta <-  mean(den[3:48,3:48] * thi_der[3:48,3:48]^2)
    }
    
    #return(theta)
    res<- list(den,thi_der,theta)
    names(res) <- c('den','thi_der','theta')
    return(res)
  }
# ============================================================
# FUNCTION: Epan
# Purpose: Compute the Epanechnikov kernel function.
# z: Numeric vector or scalar.
# Output: Numeric vector of kernel weights evaluated at z.
# ============================================================
  Epan <- function(z){
    return( 3/4 * (1-z^2) * (abs(z)<1) )
  } 
# ============================================================
# FUNCTION: gene_data
# Purpose: Generate the simulated functional data under the model considered in the manuscript.
# Inputs:
# mk: Integer. Number of subjects.
# njk: Integer vector of length mk containing the number of measurements for each subject.
# Output: A list containing the time points and the observed responses for each subject.
# ============================================================
  gene_data <- function(mk, njk){
    
    t <- lapply(1:mk, function(i) runif(njk[i],a,b))
    e <- lapply(1:mk, function(i) rnorm(njk[i],0,0.5))
    mu <- lapply(t, fun_mu)
    phi <- lapply(t, fun_phi)
    kesi <- c()
    for(i in 1:Mpc){
      kesi <- cbind(kesi,c(rnorm(mk,0,sqrt(lam[i]))))
    }
    y <- lapply(1:mk, function(i) mu[[i]] + 
                  rowSums(matrix(rep(kesi[i,],each=njk[i]),ncol=Mpc)
                          * phi[[i]]) + e[[i]])
    #z <- mu_z + as.vector(t(beta0)%*%t(kesi)) + rnorm(mk,0,0.1) 
    mylist <- list(t,y)
    names(mylist) <- c('t','y')
    return(mylist)
  }
# ============================================================
# FUNCTION: gene_gam_data
# Purpose: Construct pseudo observations for the covariance estimation based on the residual products.
# Inputs:
# NK: Integer. The total number of observations for all subjects.
# njk: Integer vector of length mk containing the number of measurements for each subject.
# mk: The number of subjects.
# x: Numeric vector of the time points with length NK.
# y: Numeric vector of the observed responses with length NK.
# Output: A list containing the numeric matrix of dimension M x 2 paired time points for the covariance function estimation 
#         and the numeric vector of length M residual cross-products corresponding to u.
# ============================================================
  gene_gam_data<-function(NK, njk, mk, x, y, mu_est){
    
    coor1 <- unlist(sapply(1:NK, function(i) 
      rep((1:NK)[i],njk[min(which(i<=cumsum(njk)))]-1)))
    II <- c(0,cumsum(njk))
    coor2 <- unlist(sapply(1:mk,function(i){
      a <- (1:NK)[(II[i]+1):II[i+1]]
      b <- as.vector(sapply(1:length(a), function(j) a[-j]))
    }))
    rm(II)
    u <- t(sapply(1:length(coor1), function(i) c(x[coor1[i]],x[coor2[i]])))
    v <- sapply(1:length(coor1), function(i) 
      (y[coor1[i]]-mu_est[coor1[i]])*(y[coor2[i]]-mu_est[coor2[i]]))
    rm(coor1,coor2)
    
    res <- list(u,v)
    names(res) <- c('u','v')
    return(res)
  }  
# ============================================================
# FUNCTION: sign_complex
# Purpose: Compute the sign of real or complex inputs.
# For the complex numbers, the sign is defined using the modulus.
# Inputs:
# x: Numeric or complex vector/scalar.
# Output: Numeric sign (-1, 0, or 1).
# ============================================================
  sign_complex <- function(x) {
    if (is.complex(x)) {
      return(sign(Mod(x)))
    } else {
      return(sign(x))
    }
  }
}
