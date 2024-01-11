######################################
# Geometric Typicality Test Function #
######################################

## Based on code by B. Le Roux 
## https://helios2.mi.parisdescartes.fr/~lerb/livres/CIGDA/DataR_fr.html

gt.test <- function(base, wt, notable_D = 0.4, 
                    alpha = 0.05, max_number = 1e+06, n_dir = 200, 
                    seed = 2022, centre = TRUE, mean = FALSE) {
  n     <- dim(base)[1]
  K     <- dim(base)[2]
  X.IK  <- as.matrix(base, nrow= n, ncol= K)
  out <- list()
  if(centre == TRUE){
    X_P <- rep(0,K)
  } else {
    X_P <- mean
  }
  if(!missing(wt)){
    Mcov.KK <- cov.wt(X.IK, wt = wt, method = "ML")$cov
    mean.KK <- cov.wt(X.IK, wt = wt, method = "ML")$center
  } else {
    Mcov.KK <- cov.wt(X.IK, method = "ML")$cov
    mean.KK <- cov.wt(X.IK, method = "ML")$center
  }
  eig <- eigen(Mcov.KK, symmetric = TRUE)
  L <- sum(eig$values > 1.5e-8) ; lambda.L <- eig$values[1:L] 
  BasisChange.KL <- eig$vectors[ ,1:L] %*% diag(1/sqrt(lambda.L), nrow=L)
  Z_GM.IL  <- sweep(X.IK, 2, mean.KK, "-") %*% BasisChange.KL
  Z_GP.L   <-  t(BasisChange.KL) %*% t(t(mean.KK - X_P))
  norm2_PG <- sum(Z_GP.L^2)
  if(norm2_PG < notable_D^2){
    print("Descriptively, the deviation is small: there is no point in doing the test.")
    out$norm2_PG <- norm2_PG
    out$mdist <- round(sqrt(norm2_PG),3)
    return(out)
  } else if (2^(n-1) <= max_number){
    cardJ <- 2^(n-1)                    
    Epsilon.JI <- matrix(1L, nrow= cardJ, ncol= n)
    for(i in 1:(n-1))
      Epsilon.JI[ , i+1] <- rep(c(1L,-1L), each=2^(n-i-1), times=2^(i-1))
    Epsilon.J <- t(t(rowSums(Epsilon.JI)))/n
    Z_GGj.JL  <- Epsilon.JI %*% Z_GM.IL/n
    rm(Epsilon.JI)
  } else {
    cardJ     <-  max_number; set.seed(seed)
    Epsilon.J <- matrix(0L, nrow= cardJ, ncol= 1)
    Z_GGj.JL  <- matrix(0,  nrow= cardJ, ncol= L)
    for (j in 1:cardJ) {
      epsilon_j.I   <- t(sample(c(-1L, 1L), size= n, replace= TRUE))
      Epsilon.J[j]  <- sum(epsilon_j.I)/n
      Z_GGj.JL[j, ] <- epsilon_j.I %*% Z_GM.IL/n
    }
  }
  if (norm2_PG >= notable_D^2) {
    Z_PPj.JL <- Z_GGj.JL + Epsilon.J %*% t(Z_GP.L)
    norm2_PPj.J <- rowSums(Z_PPj.JL^2)
    d2P.J <- norm2_PPj.J - (Z_PPj.JL %*% Z_GP.L)^2 / (1 + norm2_PG)
    rm(Z_PPj.JL)
    d2P_obs     <- norm2_PG/(1 + norm2_PG)
    n_sup       <- sum(d2P.J >= d2P_obs * (1 - 1e-12))
    p_value     <- n_sup/cardJ
    cat(" M-distance D = ", round(sqrt(norm2_PG),3), "\n",
        "Descriptively, the deviation is notable.", sep= "")
    if(K>1){
      cat(" p-value = ", 2*n_sup, "/", 2*cardJ, " = ",
          format(round(p_value, 3), nsmall= 3),"\n", sep= "")
      p_value <- round(p_value, 3)
    } else {
      cat(" p-value = ", n_sup, "/", 2*cardJ, " = ",
          format(round(p_value/2, 3), nsmall= 3)," (one-sided) \n", sep= "") 
      p_value <- round(p_value/2, 3)
    }
  }
  out$d2P_obs <- d2P_obs
  out$p_value <- p_value
  out$m_dist <- round(sqrt(norm2_PG),3)
  #out$mean <- mean.KK
  out$Mcov <- Mcov.KK
  
  set.seed(seed); rank_inf <- trunc(alpha * cardJ) + 1
  # cat(rank_inf)
  C.J <- t(t(rowSums(Z_GGj.JL^2)))                 
  Anul <- which(abs(C.J + Epsilon.J^2 - 1)  < 1e-12)
  if (K==1) n_dir <- 1
  
  kappa.D <- rep(0, 2 * n_dir)
  for (d in 1: n_dir){
    xy.L <- runif(L, -1, 1) ; beta.L <- t(t(xy.L / sqrt(sum(xy.L^2))))
    U.J <- Z_GGj.JL %*% beta.L
    A.J <- C.J - U.J^2 + Epsilon.J^2- 1
    B.J <- -Epsilon.J * U.J
    DeltaP.J <- abs(B.J^2 - A.J*C.J) 
    X1.J <- (-B.J + sqrt(DeltaP.J))/A.J
    X2.J <- (-B.J - sqrt(DeltaP.J))/A.J
    X1.J[Anul] <- -Inf ; X2.J[Anul] <- Inf
    
    X1_sorted  <- sort(X1.J); X2_sorted  <- sort(X2.J)
    kappa.D[2*d - 1] <- abs(X1_sorted[rank_inf])
    kappa.D[2*d] <- X2_sorted[cardJ + 1 - rank_inf]
  }
  if (max(kappa.D) == Inf) {
    warning("Finite compatibility region is not accessible")
  } else {
    if (L > 1){
      cat("\n Adjusted ", 100*(1 - alpha), "% compatibility region:\n", 
          "  principal ellipsoid of the cloud with scale parameter", 
          " kappa = ", round(mean(kappa.D), 3), " \n", 
          "  (mean of ", 2*max(L,n_dir)," kappa values in the range ", 
          floor(min(kappa.D)*1000)/1000, " to ", 
          ceiling(max(kappa.D)*1000)/1000, "). \n", sep= "")
      kappa.D <- round(mean(kappa.D), 3)
    } else { 
      x1 <- mean.KK -  beta.L*kappa.D[1]*sqrt(Mcov.KK[1,1])
      x2 <- mean.KK +  beta.L*kappa.D[2]*sqrt(Mcov.KK[1,1])
      cat("\n ", 100*(1 - alpha), "% compatibility interval [",
          round( min(x1,x2),3)," ; ", round( max(x1,x2), 3),"]", sep= "")
      c_int <- c(round(min(x1,x2),3),round( max(x1,x2), 3))
    }
    out$kappa.D <- kappa.D
    if(L <= 1) out$c_int <- c_int
    return(out)
  }
}