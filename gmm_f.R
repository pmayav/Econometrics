gmm_f<- function(y, x , z, W_hat, flag_print)
{
  
  nobs <- nrow(x)
  s_yy <- t(y)%*%y/nobs
  S_xx <- t(x)%*%x/nobs
  
  dim_xz <- dim(x[1,] %*% t(z[1,]))
  
  S_hat <- matrix(0, ncol=length(x[1,]), nrow=length(x[1,]))
  for (i in 1:nobs)
    S_hat <- S_hat + e_hat[i]^2 * x[i,] %*% t(x[i,])/nobs
  
  S_xz <- matrix(0, nrow=dim_xz[1], ncol=dim_xz[2] )
  for (i in 1:nobs)
    S_xz <- S_xz +  x[i,] %*% t(z[i,])/nobs
  
  s_xy <- matrix(0, nrow=length(x[1,]), ncol=1 )
  for (i in 1:nobs)
    s_xy <- s_xy +  x[i,] * y[i]/nobs
  
  S_hat <- matrix(0, ncol=length(x[1,]), nrow=length(x[1,]))
  for (i in 1:nobs)
    S_hat <- S_hat + e_hat[i]^2 * x[i,] %*% t(x[i,])/nobs

  W_hat <- solve(S_hat)
  d_gmm <- solve( t(S_xz) %*% W_hat %*% S_xz ) %*% t(S_xz) %*% W_hat %*% s_xy
  
  Avar_d <- solve(t(S_xz)%*%W_hat%*%S_xz)
  sqrt( diag(Avar_d)/nobs )
  sqrt(sum((y - z %*% d_gmm)^2)/nobs)
  
  g_n <- s_xy - S_xz %*% d_gmm
  
  J_stat <- nobs * t(g_n) %*% W_hat %*% g_n
  J_stat
  chi_sqr_stat <- 1-pchisq(J_stat, 2)
  print(d_gmm)
  df <- nobs-ncol(x)  # degrees of freedom
  mean_y <- mean(y)
  var_y <- s_yy-mean_y^2
  
  # Print results if desired
  if (flag_print==1){
    cat("********************* GMM function ***********************\n")
    cat(" Number of Observations:",nobs);cat("\n")
    cat(" Degrees of freedom :", df);cat("\n")
    cat(" Mean of the Dependent Variable:",mean_y);cat("\n")
    cat(" Variance of the Dependent Variable:",var_y);cat("\n")
    cat(" Hansen	J-statistic :",J_stat);cat("\n")
    cat(" Chi-squared statistics :",chi_sqr_stat);cat("\n")}
}

