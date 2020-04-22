twostageLS_f	<- function(y, X, Z, flag_print)
{
  nobs <- nrow(X)
  # number of observations
  df <- nobs-ncol(X)
  mean_y <- mean(y)
  
  P <- X%*%(solve(t(X)%*%X))%*%t(X) #projection matrix
  b <- solve((t(Z)%*%P%*%Z))%*%t(Z)%*%P%*%y

  e_hat<- y-Z%*%b
  s_sqr<- t(e_hat)%*%e_hat/ nobs
  var <- (solve(t(Z)%*%P%*% Z))
  J <- t(y-Z%*%b) %*% P %*% (y-Z%*%b)/ s_sqr
  Sargans <- t(e_hat) %*% P %*% (e_hat) / s_sqr
  p_value <- (1-pchisq(Sargans, 2))
  
  print(b)
  
  s_yy <- t(y)%*%y/nobs
  var_y <- s_yy-mean_y^2
  
  # Print results if desired
  if (flag_print==1){
    cat("********************* 2SLS ***********************\n")
    cat(" Number of Observations:",nobs);cat("\n")
    cat(" Degrees of freedom:",df);cat("\n")
    cat(" Mean of the Dependent Variable:",mean_y);cat("\n")
    cat(" Variance of the Dependent Variable::",var_y);cat("\n")
    cat(" s_sqr:",s_sqr);cat("\n")
    cat(" J:",J);cat("\n")
    cat(" Sargans:",Sargans);cat("\n")
    cat(" P_value:",p_value);cat("\n")
  }
  
}
