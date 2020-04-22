# Function gmm_f.R
gmm_f <- function(y, X, Z,	W, flag_print)
{
  
  nobs <- nrow(Z)  # number of observations
  df <- nobs-ncol(Z)  # degrees of freedom
  S_xx <- t(X)%*%X/nobs  # sample cross moments
  s_xy <- t(X)%*%y/nobs  # x'y/n
  s_yy <- t(y)%*%y/nobs  # y'y/n
  s_xz <- t(X)%*%Z/nobs  #x'z/n
  s_zy <- t(Z)%*%y/nobs  #x'z/n
  mean_y <- mean(y)     # mean of y
  
  b1 <- t(s_xz)%*%W%*%s_xz   # first part of GMM estimator
  b1_inv <- solve(b1)            
  
  b2 <- t(s_xz)%*%W%*%s_xy  # second part of GMM estimator

  
  beta <- b1_inv%*%b2  # GMM estimator
  #cat("Beta:",beta)
  #cat("Sxy:",s_xy)
  

  
  # Calculate standard errors
  
  e <- y-Z%*%beta
  e <- c(e)   # convert e to a vector
  temp <- e*X
  S_hat <- t(temp)%*%temp/nobs
  S_xx_inv <- solve(S_xx)

  ssr <- t(e) %*% e   # sum of squared residuals
  var_y <- s_yy-mean_y^2  # variance of y
  R_sq <- 1-(ssr/nobs)/var_y   # R-squared
  see <- sqrt(ssr/df)     # standard error of the equation
    
  v1 <- solve(t(s_xz)%*%W%*%s_xz)    # first and third part of Avar(b)
  v2 <- t(s_xz)%*%W%*%S_hat%*%W%*%s_xz # second part of Avar(b)
  
  v <- v1%*%v2%*%v1    # calculate Avar(b)
  
  se <- sqrt(diag(v/nobs))  # standard errors
  t <- beta/se   # t-value
  
  # Hansen J statistics
  g_a <- colSums(e*X)/nobs
  J <- nobs*t(g_a)%*%solve(S_hat)%*%g_a
  #chi <- nobs*g_a%*%t(beta)%*%solve(S_hat)%*%g_a%*%beta
  #uu <- chisq.test(chi)
  sgs <- nobs^2*t(s_xy-s_xz%*%beta)%*%solve(S_xx)%*%(s_xy-s_xz%*%beta)/ssr
  sgs1 <- t(e)%*%X%*%W%*%t(X)%*%e/nobs
  
  K_L <- dim(X)[2] - dim(Z)[2]
  p <- pchisq(sgs,K_L,lower.tail = FALSE)
  
  # Print results if desired
  if (flag_print==1){
    cat("************* GMM Regression ***************\n")
    cat(" Number of Observations:",nobs);cat("\n")
    cat(" Degree of Freedom (K-L):",K_L);cat("\n")
    cat(" Degree of Freedom (Z):",df);cat("\n")
    cat(" Centered R-squared:",R_sq);cat("\n")
    cat(" Standard Error of the Equation:",see);cat("\n")
    cat(" Sum of Squared Residuals:",ssr);cat("\n")
    cat(" Hansen J statistics:",J);cat("\n")
    cat(" Significance level P:",p);cat("\n")
    cat(" Sargan's Statistic:",sgs);cat("\n")
    cat("--------------------------------------------------------------\n")
    cat("variable no.   estimate       s.e.        t\n")
    cat("--------------------------------------------------------------\n")
    for (i in 1:length(beta)){
#      cat(b[i],"   ",se[i],"   ",t[i],"\n")
      cat(sprintf('%2.0f  %20.6f %12.6f %12.6f\n', i, beta[i], se[i], t[i]))
    }
    cat("--------------------------------------------------------------\n")
  }
    
    
    
  return(list
         (b=beta,se=se,res=e,t=t,v=v,nobs=nobs,mean_y=mean_y,R_sq=R_sq,see=see,ssr=ssr,J=J,df=K_L,sgs=sgs,sgs1=sgs1))
}