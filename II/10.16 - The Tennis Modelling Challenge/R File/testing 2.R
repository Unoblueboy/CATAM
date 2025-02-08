test_coefs = array(0, dim = c(100,59))

for (j in 1:100) {
  Training_Y = round(runif(Training_length,0,1))
  
  X <- matrix(data=0, nrow=Training_length, ncol = length(tennis_names)-1)
  colnames(X) <- tennis_names[-1]
  
  for (i in 1:Training_length) {
    winner = Training_winners[i]
    loser = Training_losers[i]
    y_factor = Training_Y[i]*2-1
    if (winner != tennis_names[1]) {
      X[i, winner] = 1 * y_factor
    }
    if (loser != tennis_names[1]) {
      X[i, loser] = -1 * y_factor
    }
  }
  
  coefs <- as.vector(coef(TennisGLM1, s=0))
  coefs <- coefs[2: length(coefs)]
  names(coefs) <- colnames(X)
  test_coefs[j, ] = coefs
}
