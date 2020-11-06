
library(purrr)
library(MASS)
library(rootSolve)

# data generation functions
source('binGEN_nex.R')
source('binGEN_ed.R')

# cluster-period GEE functions
source('cpgee_nex.R')
source('cpgee_ed.R')

# SW-CRT design
n <- 12           # number of clusters
t <- 5            # number of periods
n_lower <- 20     # cluster-period size lower bound
n_upper <- 30     # cluster-period size upper bound
maxiter <- 500    # maximum number of iterations for Fisher scoring updates
epsilon <- 0.001  # tolerance for convergence

# true parameter values
delta <- log(0.5)
beta <- cumsum(c(log(0.35/(1-0.35)),-0.1/2,-0.1/(2^2),-0.1/(2^3), -0.1/(2^4)))

set.seed(1234)
m <- matrix(rdunif(n*t, n_lower, n_upper), n, t) # matrix of cluster-period sizes


# individual-level marginal mean design matrix including period and treatment indicators
X <- NULL
trtSeq <- matrix(0, t - 1, t)
trtSeq[upper.tri(trtSeq)] <- 1
g <- n/(t - 1) # number of clusters per step

for(i in 1:n){
    j = (i-1) %/% g + 1
    base_mat <- cbind(diag(t), trtSeq[j,])
    for(k in 1:t){X <- rbind(X, t(replicate(m[i,k] , base_mat[k, ])))}
}


# period label 
period <- c()
for(i in 1:n){
    for(j in 1:t){
        period <- c(period, rep(j, m[i, j]))
    }
}


#####################################################################
### (1) Nested exchangeable correlation structure ###################
#####################################################################
# alpha values
alpha <-  c(0.1, 0.05)


## 1.1 Generate cross sectional SW-CRT correlated binary outcome data
## with the nested exchangeable correlation structure
set.seed(1234)
y <- binGEN_nex(n = n, m = m, t = t, delta = delta, beta = beta, alpha = alpha)


# Create id
id <- rep(1:n, apply(m, 1, sum))         # create cluster id
#rep(1:n, each=t*m)
clsize <- apply(m, 1, sum)


# cluster-period outcome 
clp_mu <- tapply(y, list(id, period), FUN = mean)
y_cp <- c(t(clp_mu))


# cluster-period design matrix, id, cluster sample size vector and cluster-period size vector
trt <- tapply(X[, t + 1], list(id, period), FUN=mean)
trt <- c(t(trt))

time <- tapply(period,list(id,period),FUN=mean)
time <- c(t(time))


X_cp <- matrix(0, n * t, t)

s = 1
for(i in 1:n){for(j in 1:t){X_cp[s, time[s]] <- 1; s = s + 1}}
X_cp <- cbind(X_cp, trt); id_cp <- rep(1:n, each= t); n_cp <- rep(t, n); m_cp <-  c(t(m))


### MAEE
est_maee_nex <- cpgee_nex(y = y_cp, X = X_cp, id = id_cp, n = n_cp, m = m_cp, maxiter = maxiter, epsilon = epsilon, printrange = TRUE, alpadj = TRUE)
#####################################################################################
# GEE and MAEE for Cluster-Period Summaries 
# Number of Clusters: 12 
# Maximum Cluster-Period Size: 30 
# Minimum Cluster-Period Size: 20 
# Number of Iterations: 4 
# Correlation Structure: Nested exchangeable 
# Results for marginal mean parameters 
# Beta   Estimate MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 -0.7063758 0.2190786  0.1786765  0.1866340  0.1949465  0.1875385
# [2,]    1 -0.9776584 0.2419961  0.2925541  0.3092297  0.3269758  0.3077233
# [3,]    2 -0.9749175 0.2661480  0.2816550  0.3006636  0.3212165  0.2994892
# [4,]    3 -1.1629618 0.3201781  0.2701559  0.2869750  0.3050911  0.2874543
# [5,]    4 -1.1932122 0.3715894  0.3397117  0.3650041  0.3924780  0.3673190
# [6,]    5 -0.3387314 0.2688767  0.2519009  0.2750049  0.3002816  0.2720360
# 
# Results for correlation parameters 
# Alpha   Estimate BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]     0 0.09140198 0.02773135 0.02895367 0.03023050 0.02892813
# [2,]     1 0.06534831 0.02668092 0.02783518 0.02903978 0.02782306
#####################################################################################

### UEE
est_uee_nex <- cpgee_nex(y = y_cp, X = X_cp, id = id_cp, n = n_cp, m = m_cp, maxiter = maxiter, epsilon = epsilon, printrange = TRUE, alpadj = FALSE)
#####################################################################################
# GEE and MAEE for Cluster-Period Summaries 
# Number of Clusters: 12 
# Maximum Cluster-Period Size: 30 
# Minimum Cluster-Period Size: 20 
# Number of Iterations: 3 
# Correlation Structure: Nested exchangeable 
# Results for marginal mean parameters 
# Beta   Estimate MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 -0.7064384 0.2087681  0.1787601  0.1867233  0.1950422  0.1876938
# [2,]    1 -0.9786578 0.2306638  0.2925892  0.3092723  0.3270288  0.3078128
# [3,]    2 -0.9750309 0.2532540  0.2827734  0.3019083  0.3226023  0.3009627
# [4,]    3 -1.1638213 0.3048610  0.2705990  0.2874016  0.3054902  0.2881891
# [5,]    4 -1.1936248 0.3530481  0.3409917  0.3664460  0.3941003  0.3690849
# [6,]    5 -0.3395339 0.2550355  0.2523288  0.2755175  0.3008924  0.2726413
# 
# Results for correlation parameters 
# Alpha   Estimate BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]     0 0.07922396 0.02544368 0.02656716 0.02774086 0.02654176
# [2,]     1 0.06015348 0.02424914 0.02529966 0.02639609 0.02528767
#####################################################################################




#####################################################################
### (2) Exponential decay correlation structure #####################
#####################################################################
# alpha values
alpha <- c(0.03, 0.8)


## 2.1 Generate cross sectional SW-CRT correated binary outcome data
## with the exponential decay correlation structure
set.seed(1234)
y <- binGEN_ed(n = n, m = m, t = t, delta = delta, beta = beta, alpha = alpha)

# Create id
id <- rep(1:n, apply(m, 1, sum))         # create cluster id
#rep(1:n, each=t*m)
clsize <- apply(m, 1, sum)


# cluster-period outcome 
clp_mu <- tapply(y, list(id, period), FUN = mean)
y_cp <- c(t(clp_mu))


### MAEE
est_maee_ed <- cpgee_ed(y = y_cp, X = X_cp, id = id_cp, n = n_cp, m = m_cp, 
                        maxiter = maxiter, epsilon = epsilon, printrange = TRUE, 
                        alpadj = TRUE, rho.init = NULL)
#####################################################################################
# GEE and MAEE for Cluster-Period Summaries 
# Number of Clusters: 12 
# Maximum Cluster-Period Size: 30 
# Minimum Cluster-Period Size: 20 
# Number of Iterations: 6 
# Correlation Structure: Exponential decay 
# Results for marginal mean parameters 
# Beta   Estimate MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 -0.6491236 0.1620461  0.1027036  0.1072544  0.1120075  0.1073747
# [2,]    1 -0.8100334 0.1787675  0.2341744  0.2479077  0.2625393  0.2465087
# [3,]    2 -0.7726590 0.2009338  0.1981924  0.2120216  0.2270366  0.2081359
# [4,]    3 -0.8649966 0.2494992  0.1974919  0.2103508  0.2241981  0.2115468
# [5,]    4 -0.7834666 0.2991502  0.3038734  0.3270794  0.3521934  0.3259170
# [6,]    5 -0.6557009 0.2345900  0.2329788  0.2535574  0.2760114  0.2503363
# 
# Results for correlation parameters 
# Alpha   Estimate BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]     0 0.03273627 0.01102476 0.01152433 0.01204858 0.01150553
# [2,]     1 0.71912460 0.13216743 0.13841914 0.14497231 0.13868140
#####################################################################################


### UEE
est_uee_ed <- cpgee_ed(y = y_cp, X = X_cp, id = id_cp, n = n_cp, m = m_cp, 
                       maxiter = maxiter, epsilon = epsilon, printrange = TRUE, 
                       alpadj = FALSE, rho.init = NULL)
#####################################################################################
# GEE and MAEE for Cluster-Period Summaries 
# Number of Clusters: 12 
# Maximum Cluster-Period Size: 30 
# Minimum Cluster-Period Size: 20 
# Number of Iterations: 5 
# Correlation Structure: Exponential decay 
# Results for marginal mean parameters 
# Beta   Estimate MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 -0.6493925 0.1546994  0.1027768  0.1073365  0.1120993  0.1074349
# [2,]    1 -0.8122632 0.1709768  0.2339503  0.2476299  0.2622017  0.2461086
# [3,]    2 -0.7774128 0.1923282  0.1970108  0.2108822  0.2259501  0.2067343
# [4,]    3 -0.8722992 0.2394327  0.1954702  0.2083336  0.2222106  0.2093342
# [5,]    4 -0.7944567 0.2866587  0.3034874  0.3267139  0.3518572  0.3256770
# [6,]    5 -0.6458147 0.2250621  0.2323677  0.2529255  0.2753647  0.2498266
# 
# Results for correlation parameters 
# Alpha   Estimate BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]     0 0.02619168 0.01004136 0.01050376 0.01098955 0.01048291
# [2,]     1 0.75537815 0.13876265 0.14524387 0.15203232 0.14540603
#####################################################################################
