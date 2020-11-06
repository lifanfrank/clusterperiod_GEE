################################################################################################
# Simulate correlated binary responses
# in a cohort stepped wedge CRT with a nested exchangeable correlation structure
# using the conditional linear family method (Qaqish, 2003)

# Ref: Qaqish, B. F. (2003). A family of multivariate binary distributions 
# for simulating correlated binary variables. Biometrika 90, 455-463.

# INPUT
# n: Number of clusters
# m: Cluster size a matrix of n * t
# t: Number of periods
# delta: Effect size in log odds ratio
# beta: Vector of period effects
# alpha: Vector of correlations: (alpha_0, alpha_1)
#        alpha_0 = within period correlation
#        alpha_1 = inter-period correlation
##############################################################################################

binGEN_nex <- function(n, m, t, delta, beta, alpha){
  
  ########################################################
  # Create block exchangeable correlation matrix.
  ########################################################
  
  bxch<-function(alpha, n_i){
    alpha_0 = alpha[1]
    alpha_1 = alpha[2]

    n_i = c(n_i)
    n_length = length(n_i)
    cor_matrix = matrix(alpha_1, sum(n_i), sum(n_i))

    loc1 = 0
    loc2 = 0 
    for(t in 1:n_length){
        n_t = n_i[t]
        loc1 = loc2 + 1
        loc2 = loc1 + n_t - 1
        for(i in loc1:loc2){
            for(j in loc1:loc2){
                if(i != j){
                    cor_matrix[i, j] = alpha_0
                }else{
                    cor_matrix[i, j] = 1
                }
            }
        }
    }
    return(cor_matrix)
  }
  
  ########################################################
  # a[1:n, 1:n] is the input covariance matrix of Y[1:n].
  # Returns  b[1:n,1:n] such that b[1:t-1, t] are the 
  # slopes for regression of y[t] on y[1:t-1], for t=2:n.
  # Diagonals and lower half of b[,] are copied from a[,].
  # a[,] is assumed +ve definite symmetric, not checked.
  ########################################################
  
  allreg<-function(a){

    n <- nrow(a)
    b <- a
    for(t in 2:n){
      t1 <- t - 1
      gt <- a[1:t1, 1:t1]
      st <- a[1:t1, t]
      bt <- solve(gt,st)
      b[1:t1,t] <- bt
    }

    return(b)
  }
  
  ########################################################
  # returns variance matrix of binary variables with mean
  # vector u[] and corr matrix r[,].
  ########################################################
  
  cor2var <- function(r, u){
    s <- diag(sqrt(u * (1 - u)))
    return(s %*% r %*% s)
  }
  
  ########################################################
  # r[1:n, 1:n] is the corr mtx
  # u[1:n] is the mean of a binary vector
  # checks that pairwise corrs are in-range for the given u[]
  # only upper half of r[,] is checked.
  # return 0 if ok
  # return 1 if out of range
  ########################################################
  
  chkbinc <- function(r, u){
    n <- length(u)
    s <- sqrt(u * (1 - u))
    for(i in 1:(n - 1)){
      for(j in (i + 1):n){
       uij <- u[i] * u[j] + r[i,j] * s[i] * s[j]
       ok <- ((uij <= min(u[i], u[j])) & (uij >= max(0, u[i]+u[j]-1)))
       if(!ok) {return(1)}
      }
    }
    return(0)
  }
  
  ########################################################
  # Multivariate Binary Simulation by Linear Regression.
  # Simulate a single vector.
  # Returns a simulated binary random vector y[1:n] with mean 
  # u[1:n] and regression coefs matrix b[1:n,1:n] (obtained 
  # by calling allreg() above).
  # y[] and u[] are column vectors.
  # Returns -1 if the cond. linear family not reproducible
  ########################################################
  
  mbslr1 <- function(b, u){
    n <- nrow(b)
    y <- rep(-1, n)
    y[1] <- rbinom(1, 1, u[1])
    for(i in 2:n){
      i1 <- i - 1
      r <- y[1:i1] - u[1:i1]              # residuals
      ci <- u[i] + sum(r * b[1:i1, i])       # cond.mean
      if(ci < 0 | ci > 1){
        stop(paste("mbslr1: ERROR:", ci))
        return(-1)
      }
      y[i] <- rbinom(1, 1, ci)
    }
    return(y)
  }
  

    
  # Create treatment sequences
  trtSeq <- matrix(0, t - 1, t)
  trtSeq[upper.tri(trtSeq)] <- 1
  g <- n / (t - 1)                         
   # number of clusters per step
  
  # Simulate correlated binary outcomes
  y <- NULL
  for(i in 1:n){
      n_i = m[i, ]
      r = bxch(alpha, n_i)
      j = (i-1) %/% g + 1
      u_c <- c(plogis(beta + delta*trtSeq[j,]))
      u <- rep(u_c, c(n_i))
      v <- cor2var(r,u) 
      oor <- chkbinc(r,u)  
      if(oor){
        stop("ERROR: Corr out of range for given mean")
        }
      b <- allreg(v) 
      y <- c(y, mbslr1(b,u)) 
  }


  # Return simulated data vector
  return(y)
}




