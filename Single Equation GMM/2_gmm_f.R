efficient2step_gmm_f<- function(y, x, z, flag_print) 
{  
  nobs <- nrow(x)
  
  e_hat2 <- resid( line_3 )[[1]]
  S_hat2 <- matrix(0, ncol=length(x[1,]), nrow=length(x[1,]))
  for (i in 1:nobs)
    S_hat2 <- S_hat2 + e_hat2[i]^2 * x[i,] %*% t(x[i,])/nobs
  
  dim_xz2 <- dim(x[1,] %*% t(z[1,]))
  
  S_xz2 <- matrix(0, nrow=dim_xz2[1], ncol=dim_xz2[2] )
  for (i in 1:nobs)
    S_xz2 <- S_xz2 +  x[i,] %*% t(z[i,])/nobs
  
  s_xy2 <- matrix(0, nrow=length(x[1,]), ncol=1 )
  for (i in 1:nobs)
    s_xy2 <- s_xy2 +  x[i,] * y[i]/nobs
  
  W_hat2 <- solve(S_hat2)
  
  d_gmm2 <- solve( t(S_xz2) %*% W_hat2 %*% S_xz2 ) %*% t(S_xz2) %*% W_hat2 %*% s_xy2
  d_gmm2
  
  g_n2 <- s_xy2 - S_xz2 %*% d_gmm2
  
  J_stat2 <- nobs * t(g_n2) %*% W_hat2 %*% g_n2
  J_stat2
  
  W_hat1 <- solve(S_hat2[1:15,1:15])
  d_gmm1 <- solve( t(S_xz2[1:15,]) %*% W_hat1 %*% S_xz2[1:15,] ) %*% t(S_xz2[1:15,]) %*% W_hat1 %*% s_xy2[1:15,]
  print(d_gmm1)
  g_n1 <- s_xy2[1:15] - S_xz2[1:15,] %*% d_gmm1
  
  
  J_stat1 <- nobs * t(g_n1) %*% W_hat1 %*% g_n1
  J_stat1
  
  C_stat <- J_stat2 - J_stat1
  C_stat
  
  df <- nobs-ncol(x)
  mean_y <- mean(y)
  s_yy <- t(y)%*%y/nobs
  var_y <- s_yy-mean_y^2
  print(d_gmm2)
  # Print results if desired
  if (flag_print==1){
    cat("********************* 2 Step GMM function ***********************\n")
    cat(" Number of Observations:",nobs);cat("\n")
    cat(" Degrees of freedom :", df);cat("\n")
    cat(" Mean of the Dependent Variable:",mean_y);cat("\n")
    cat(" Variance of the Dependent Variable:",var_y);cat("\n")
    cat(" Hansen	J-statistic1 :",J_stat1);cat("\n")   
    cat(" Hansen	J-statistic2 :",J_stat2);cat("\n")
    cat(" C_stat :",C_stat);cat("\n")
  }
}