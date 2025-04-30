LinMod_sim <- function(X,Beta = rep(0,p),  errors_gen = rnorm, B =10000){
  n <- nrow(X)
  p <- ncol(X)
  
  #Aha, so we multiply our XB value - which gives us a vector of len n (sample size), then because 
  #we're adding on a matrix of dimension $n*B$ this will add XB to each col (as in R - addition in this way cycles around!)
  #So, we end up with $B$ different instances!
  y_mat <- as.vector(X %*% Beta) + matrix(errors_gen(n*B) , n ,B)
  #Next, invert the matrix X^TX
  inv_gram <- solve(t(X) %*% X )
  P <- X %*% inv_gram %*% t(X)
  #Again, Beta hat is going to be a matrix with B cols as we're taking B samples
  Beta_hat_mat <- inv_gram %*% t(X) %*% y_mat
  residue_mat <- y_mat - P%*% y_mat
  
  #Now - for each instance, we'll get a sigma hat which estimatest the variance,
  #We convert this to the form useful to the t stat by mutliplying by (n/n-p)
    
  sigma_tilde_vec <- colSums(residue_mat^2) * 1/ (n-p)
  #And finally - we can calculate our t-statistic estimators
  t_stat_mat <- Beta_hat_mat / sqrt(diag(inv_gram))
  t_stat_mat <- t_stat_mat / rep(sqrt(sigma_tilde_vec), each =p)
  #Note that because we're getting a statistic for each index of beta from 1 to p
  #we've replicated the sigma of that instance $p$ times!
  #Now can return our output!
  return(t_stat_mat)
}


qqplot_F <- function(f_stat, df1, df2){
  #dfs are the degrees of freedom, not some dataframes!
  f_stat <- sort(f_stat)
  B <- length(f_stat)
  #print(B)
  theoretical_quantiles <- qf( (1:B) /(B+1) , df1, df2)
  #print(theoretical_quantiles)
  plot(theoretical_quantiles, f_stat)
  #Then we add a line (a,b being y = bx +a)
  abline(0,1,col = 'red')
  return(mean(f_stat<= qf(0.95,df1,df2)))
  #We're returning the percentage of f_stat that lies below the 95% point of the f distrib.
}



LinMod_sim_modified <- function(X,Beta = rep(0,p),  errors_gen = rnorm, B =10000){
  #In this function we modify the function so that the error variances are unequal,
  #the way he mentions is to consider multiplying the error matrix to scale each row
  #by some random amount (taken from a normal) We have written this below. 

  n <- nrow(X)
  p <- ncol(X)
  
  y_mat <- as.vector(X %*% Beta) + matrix(errors_gen(n*B) , n ,B) * rnorm(n, )
  #Next, invert the matrix X^TX
  inv_gram <- solve(t(X) %*% X )
  P <- X %*% inv_gram %*% t(X)
  #Again, Beta hat is going to be a matrix with B cols as we're taking B samples
  Beta_hat_mat <- inv_gram %*% t(X) %*% y_mat
  residue_mat <- y_mat - P%*% y_mat
  
  #Now - for each instance, we'll get a sigma hat which estimatest the variance,
  #We convert this to the form useful to the t stat by mutliplying by (n/n-p)
  
  sigma_tilde_vec <- colSums(residue_mat^2) * 1/ (n-p)
  #And finally - we can calculate our t-statistic estimators
  t_stat_mat <- Beta_hat_mat / sqrt(diag(inv_gram))
  t_stat_mat <- t_stat_mat / rep(sqrt(sigma_tilde_vec), each =p)
  #Note that because we're getting a statistic for each index of beta from 1 to p
  #we've replicated the sigma of that instance $p$ times!
  #Now can return our output!
  return(t_stat_mat)
}


LinMod_sim_mod_Q2 <- function(X,Beta = rep(0,p),  errors_gen = rnorm, B =10000,
                              G_0 = rep(1,p)){
  n <- nrow(X)
  p <- ncol(X)
  
  #So G_0 is supposed to select the variables that we want to choose, we shall represent
  #This as a vector of length $p$ where the variables to be ignored are given a zero
  #Default is that no variables are ignored. 
  
  p0 = sum(G_0)
  X0 = X[G_0]
  Beta0 = Beta[G_0]
  
  
  y_mat <- as.vector(X %*% Beta) + matrix(errors_gen(n*B) , n ,B)
  #Next, invert the matrix X^TX
  inv_gram <- solve(t(X) %*% X )
  P <- X %*% inv_gram %*% t(X)
  #Again, Beta hat is going to be a matrix with B cols as we're taking B samples
  Beta_hat_mat <- inv_gram %*% t(X) %*% y_mat
  residue_mat <- y_mat - P%*% y_mat
  
  #Now - for each instance, we'll get a sigma hat which estimatest the variance,
  #We convert this to the form useful to the t stat by mutliplying by (n/n-p)
  
  sigma_tilde_vec <- colSums(residue_mat^2) * 1/ (n-p)
  #And finally - we can calculate our t-statistic estimators
  t_stat_mat <- Beta_hat_mat / sqrt(diag(inv_gram))
  t_stat_mat <- t_stat_mat / rep(sqrt(sigma_tilde_vec), each =p)
  #Note that because we're getting a statistic for each index of beta from 1 to p
  #we've replicated the sigma of that instance $p$ times!
  #Now can return our output!
  return(t_stat_mat)
}






