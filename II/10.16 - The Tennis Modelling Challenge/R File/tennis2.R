library(glmnet)
library(quadprog)
setwd('E:/Google Drive/School/Cambridge/Maths/II/04 - CATAM/10.16 - The Tennis Modelling Challenge/R File')
Tennis <- read.csv("II-10-16-2019-mensResults.csv")
Tennis$Date <- as.Date(Tennis$Date,format='%d/%m/%y')
lam_seq = c(exp(seq(log(1), log(1e-10), length.out = 100)), 0)

invlogit <- function(x) exp(x) / (1 + exp(x))
K <- function(theta) log(1 + exp(theta))

Training_Data <- Tennis[(Tennis$Date <= as.Date("2014/12/31")) & (!is.na(Tennis$Date)), ]
Training_length = dim(Training_Data)[1]
Test_Data <- Tennis[(Tennis$Date >= as.Date("2015/01/01")) & (!is.na(Tennis$Date)), ]
Test_length = dim(Test_Data)[1]

tennis_names = levels(Tennis$Winner)
Training_winners <- as.vector(Training_Data$Winner)
Training_losers <- as.vector(Training_Data$Loser)
Test_winners <- as.vector(Test_Data$Winner)
Test_losers <- as.vector(Test_Data$Loser)

# QUESTION 2

bin_blank_data = matrix(0, ncol=4, nrow=length(tennis_names)*(length(tennis_names)-1)/2)

for (i in 1:(length(tennis_names)-1)) {
  beg_ind = (2*length(tennis_names) - i)*(i-1)/2
  for (j in (i+1):length(tennis_names)) {
    bin_blank_data[beg_ind + j - i, ] = c(tennis_names[i], tennis_names[j], 0, 0)
  }
}

bin_training_data <- bin_blank_data

for (i in 1:Training_length) {
  winner = Training_winners[i]
  loser = Training_losers[i]
  win_ind = which(tennis_names == winner)
  los_ind = which(tennis_names == loser)
  
  if (win_ind<los_ind) {
    # ind = which(apply(bin_training_data[, c(1,2)], 1, function(x) {all(x==c(winner, loser))}))
    ind = (2*length(tennis_names) - win_ind)*(win_ind-1)/2 + los_ind - win_ind
    bin_training_data[ind, 3] = as.numeric(bin_training_data[ind, 3]) + 1
  } else {
    #ind = which(apply(bin_training_data[, c(1,2)], 1, function(x) {all(x==c(loser, winner))}))
    ind = (2*length(tennis_names) - los_ind)*(los_ind-1)/2 + win_ind - los_ind
    bin_training_data[ind, 4] = as.numeric(bin_training_data[ind, 4]) + 1
  }
}
colnames(bin_training_data) <- c("player1", "player2", "win1", "win2")
# Remove rows that have no observations
bin_training_data = bin_training_data[(bin_training_data[, 3]!= 0) | (bin_training_data[, 4]!= 0), ]

X = matrix(0, ncol = length(tennis_names), nrow=dim(bin_training_data)[1])
colnames(X) = tennis_names

for (i in 1:dim(bin_training_data)[1]) {
  # if (l[i, 1] != tennis_names[1]) {
    X[i, bin_training_data[i, 1]] = 1
  # }
  X[i, bin_training_data[i, 2]] = -1
}
# Make sure X is sparse and remove columns with 0 parameter
X <- X[, -1]
X <- as(X, "sparseMatrix")
Training_Y = cbind(as.numeric(bin_training_data[, 4]), as.numeric(bin_training_data[, 3]))

TennisGLM1 <- glmnet(X, Training_Y, lambda=0, intercept=FALSE, family="binomial", standardize = FALSE, thresh = 1e-14)
summary(TennisGLM1)


coefs <- as.numeric(coef(TennisGLM1, s=0))
coefs <- coefs[-c(1)]
names(coefs) <- tennis_names[-1]

LL1 = 0
for (i in 1:dim(X)[1]) {
  row_X = X[i, ]
  Xb = row_X %*% coefs
  y = Training_Y[i, 2]
  n = Training_Y[i, 1] + Training_Y[i, 2]
  mu = invlogit(Xb)
  term = log(choose(n, y)) + y*log(mu) + (n-y)*log(1-mu)
  LL1 = LL1 + term
}
LL1 = LL1*-1/dim(X)[1]
print(paste("TennisGLM1 Log Loss for Training Set:", LL1))

bin_Test_data <- bin_blank_data

for (i in 1:Test_length) {
  winner = Test_winners[i]
  loser = Test_losers[i]
  win_ind = which(tennis_names == winner)
  los_ind = which(tennis_names == loser)
  
  if (win_ind<los_ind) {
    # ind = which(apply(bin_training_data[, c(1,2)], 1, function(x) {all(x==c(winner, loser))}))
    ind = (2*length(tennis_names) - win_ind)*(win_ind-1)/2 + los_ind - win_ind
    bin_Test_data[ind, 3] = as.numeric(bin_Test_data[ind, 3]) + 1
  } else {
    #ind = which(apply(bin_training_data[, c(1,2)], 1, function(x) {all(x==c(loser, winner))}))
    ind = (2*length(tennis_names) - los_ind)*(los_ind-1)/2 + win_ind - los_ind
    bin_Test_data[ind, 4] = as.numeric(bin_Test_data[ind, 4]) + 1
  }
}
colnames(bin_Test_data) <- c("player1", "player2", "win1", "win2")
# Remove rows that have no observations
bin_Test_data = bin_Test_data[(bin_Test_data[, 3]!= 0) | (bin_Test_data[, 4]!= 0), ]

Test_X = matrix(0, ncol = length(tennis_names), nrow=dim(bin_Test_data)[1])
colnames(Test_X) = tennis_names

for (i in 1:dim(bin_Test_data)[1]) {
  # if (l[i, 1] != tennis_names[1]) {
  Test_X[i, bin_Test_data[i, 1]] = 1
  # }
  Test_X[i, bin_Test_data[i, 2]] = -1
}
# Make sure Test_X is sparse
Test_X <- Test_X[, -1]
Test_X <- as(Test_X, "sparseMatrix")
Test_Y = cbind(as.numeric(bin_Test_data[, 4]), as.numeric(bin_Test_data[, 3]))

LL2 = 0
for (i in 1:dim(Test_X)[1]) {
  row_X = Test_X[i, ]
  Xb = row_X %*% coefs
  y = Test_Y[i, 2]
  n = Test_Y[i, 1] + Test_Y[i, 2]
  mu = invlogit(Xb)
  term = log(choose(n, y)) + y*log(mu) + (n-y)*log(1-mu)
  LL2 = LL2 + term
}
LL2 = LL2*-1/dim(Test_X)[1]
print(paste("TennisGLM1 Log Loss for Test Set:", LL2))


# QUESTION 3

d <- rep(0, 59)
names(d) <- tennis_names[-1]
d["Federer R."] = 1
d["Murray A."] = -1
pred <- t(d) %*% coefs
critval <- qnorm(0.5 + 0.68/2)

W_diag = rep(0, dim(Training_Y)[1])
for (i in 1:dim(Training_Y)[1]) {
  theta = t(X[i, ]) %*% coefs
  W_diag[i] = (Training_Y[i,1]+Training_Y[i,2])*invlogit(theta)*(1-invlogit(theta))
}
W = diag(W_diag)
Covar_mat = solve(t(X) %*% W %*% X)
se = sqrt(as.numeric(t(d) %*% Covar_mat %*% d))
CI_prob_Bounds <- invlogit(c(pred - critval*se, pred + critval*se))
print(paste0("Confidence interval for probability that Federer beats Murray is [",
             CI_prob_Bounds[1],", ", CI_prob_Bounds[2], "]"))


# QUESTION 4

Training_surfaces <- as.vector(Training_Data$Surface)
Test_surfaces <- as.vector(Test_Data$Surface)

tennis_surfaces = levels(Tennis$Surface)
column_names = rep("",(length(tennis_names))*(length(tennis_surfaces)))
for (i in 1:length(tennis_names)) {
  column_names[5*(i-1)+1] <- tennis_names[i]
  #for (j in 2:length(tennis_surfaces))
  #column_names[4*(i-2)+(j)] <- paste(tennis_names[i], tennis_surfaces[j-1])
  for (j in 1:length(tennis_surfaces)) {
    column_names[5*(i-1)+(j+1)] <- paste(tennis_names[i], tennis_surfaces[j])
  }
}

bin_blank_data_2 = matrix(0, ncol=4, nrow=length(column_names)*(length(column_names)-1)/2)
colnames(bin_blank_data_2) <- c("player1", "player2", "win1", "win2")
blank_X2 = matrix(0, ncol = length(column_names), nrow=length(column_names)*(length(column_names)-1)/2)
colnames(blank_X2) = column_names

for (i in 1:(length(column_names)-1)) {
  beg_ind = (2*length(column_names) - i)*(i-1)/2
  for (j in (i+1):length(column_names)) {
    tot_ind = beg_ind + j - i
    bin_blank_data_2[tot_ind, ] = c(column_names[i], column_names[j], 0, 0)
    player_1 = sub(paste("", tennis_surfaces, collapse="|"), "", column_names[i])
    player_2 = sub(paste("", tennis_surfaces, collapse="|"), "", column_names[j])
    blank_X2[tot_ind, column_names[i]] = 1
    blank_X2[tot_ind, player_1] = 1
    if (player_1 != player_2){
      blank_X2[tot_ind, column_names[j]] = -1
      blank_X2[tot_ind, player_2] = -1
    }
    # The ceil basically makes it so the base parameter is always 1
  }
}

bin_training_data_2 <- bin_blank_data_2



for (i in 1:Training_length) {
  winner = Training_winners[i]
  loser = Training_losers[i]
  surface = Training_surfaces[i]
  win_surf_ind = which(column_names == paste(winner, surface))
  los_surf_ind = which(column_names == paste(loser, surface))
  win_ind = which(column_names == winner)
  los_ind = which(column_names == loser)
  
  if (win_ind<los_ind) {
    # ind = which(apply(bin_training_data[, c(1,2)], 1, function(x) {all(x==c(winner, loser))}))
    surf_ind = (2*length(column_names) - win_surf_ind)*(win_surf_ind-1)/2 + los_surf_ind - win_surf_ind
    ind = (2*length(column_names) - win_ind)*(win_ind-1)/2 + los_ind - win_ind
    bin_training_data_2[surf_ind, 3] = as.numeric(bin_training_data_2[ind, 3]) + 1
    bin_training_data_2[ind, 3] = as.numeric(bin_training_data_2[ind, 3]) + 1
  } else {
    #ind = which(apply(bin_training_data[, c(1,2)], 1, function(x) {all(x==c(loser, winner))}))
    surf_ind = (2*length(column_names) - los_surf_ind)*(los_surf_ind-1)/2 + win_surf_ind - los_surf_ind
    ind = (2*length(column_names) - los_ind)*(los_ind-1)/2 + win_ind - los_ind
    bin_training_data_2[ind, 4] = as.numeric(bin_training_data_2[ind, 4]) + 1
    bin_training_data_2[surf_ind, 4] = as.numeric(bin_training_data_2[ind, 4]) + 1
  }
}

permitted_rows = (bin_training_data_2[, 3]!= 0) | (bin_training_data_2[, 4]!= 0)
bin_training_data_2 = bin_training_data_2[permitted_rows, ]
X2 = blank_X2[permitted_rows, ]

# # rename the columns of X2 so parameters have nice names
# colnames(X2)[which(colnames(X2) %in% paste(tennis_names, "Hard"))] = tennis_names

# Make sure X2 is sparse and remove unnecessary columns
X2 <- X2[, -(1:5)]
X2 <- as(X2, "sparseMatrix")
Training_Y2 = cbind(as.numeric(bin_training_data_2[, 4]), as.numeric(bin_training_data_2[, 3]))

TennisGLM2 <- glmnet(X2, Training_Y2, lambda=lam_seq, intercept=FALSE, family="binomial", standardize = FALSE, thresh = 1e-10)
summary(TennisGLM2)


coefs2 <- coef(TennisGLM2, s=0)
new_names <- rownames(coefs2)
coefs2 <- as.numeric(coefs2[-c(1)])
names(coefs2) <- new_names[-1]

LL3 = 0
for (i in 1:dim(X2)[1]) {
  row_X2 = X2[i, ]
  X2b = row_X2 %*% coefs2
  y = Training_Y2[i, 2]
  n = Training_Y2[i, 1] + Training_Y2[i, 2]
  mu = invlogit(X2b)
  term = log(choose(n, y)) + y*log(mu) + (n-y)*log(1-mu)
  LL3 = LL3 + term
}
LL3 = LL3*-1/dim(X2)[1]
print(paste("TennisGLM2 Log Loss for Training Set:", LL3))


bin_test_data_2 <- bin_blank_data_2

for (i in 1:Test_length) {
  winner = Test_winners[i]
  loser = Test_losers[i]
  surface = Test_surfaces[i]
  win_surf_ind = which(column_names == paste(winner, surface))
  los_surf_ind = which(column_names == paste(loser, surface))
  win_ind = which(column_names == winner)
  los_ind = which(column_names == loser)
  
  if (win_ind<los_ind) {
    # ind = which(apply(bin_training_data[, c(1,2)], 1, function(x) {all(x==c(winner, loser))}))
    surf_ind = (2*length(column_names) - win_surf_ind)*(win_surf_ind-1)/2 + los_surf_ind - win_surf_ind
    ind = (2*length(column_names) - win_ind)*(win_ind-1)/2 + los_ind - win_ind
    bin_test_data_2[surf_ind, 3] = as.numeric(bin_test_data_2[ind, 3]) + 1
    bin_test_data_2[ind, 3] = as.numeric(bin_test_data_2[ind, 3]) + 1
  } else {
    #ind = which(apply(bin_training_data[, c(1,2)], 1, function(x) {all(x==c(loser, winner))}))
    surf_ind = (2*length(column_names) - los_surf_ind)*(los_surf_ind-1)/2 + win_surf_ind - los_surf_ind
    ind = (2*length(column_names) - los_ind)*(los_ind-1)/2 + win_ind - los_ind
    bin_test_data_2[ind, 4] = as.numeric(bin_test_data_2[ind, 4]) + 1
    bin_test_data_2[surf_ind, 4] = as.numeric(bin_test_data_2[ind, 4]) + 1
  }
}
permitted_rows = (bin_test_data_2[, 3]!= 0) | (bin_test_data_2[, 4]!= 0)
bin_test_data_2 = bin_test_data_2[permitted_rows, ]
Test_X2 = blank_X2[permitted_rows, ]
Test_Y2 = cbind(as.numeric(bin_test_data_2[, 4]), as.numeric(bin_test_data_2[, 3]))

# Make sure Test_X2 is sparse and remove unnecessary columns
Test_X2 <- Test_X2[, -(1:5)]
Test_X2 <- as(Test_X2, "sparseMatrix")

LL4 = 0
for (i in 1:dim(Test_X2)[1]) {
  row_Test_X2 = Test_X2[i, ]
  Test_X2b = row_Test_X2 %*% coefs2
  y = Training_Y2[i, 2]
  n = Training_Y2[i, 1] + Training_Y2[i, 2]
  mu = invlogit(Test_X2b)
  term = log(choose(n, y)) + y*log(mu) + (n-y)*log(1-mu)
  LL4 = LL4 + term
}
LL4 = LL4*-1/dim(Test_X2)[1]
print(paste("TennisGLM2 Log Loss for Test Set:", LL4))


# QUESTION 5
# Time to do anova test

coef0 = rep(0, length(coefs2))
coef0[which(names(coefs2) %in% tennis_names)] = coefs
names(coef0) <- names(coefs2)

LLa = 0
for (i in 1:dim(X2)[1]) {
  row_X2 = X2[i, ]
  X2b = row_X2 %*% coef0
  y = Training_Y2[i, 2]
  n = Training_Y2[i, 1] + Training_Y2[i, 2]
  mu = invlogit(X2b)
  term = log(choose(n, y)) + y*log(mu) + (n-y)*log(1-mu)
  LLa = LLa + term
}


test_stat = (-LLa*2-LL3*dim(X2)[1]*2)
p_val = 1-pchisq(test_stat, length(coefs2)-length(coefs))
print(paste("P-value for the likelihood ratio test is:",p_val))


# QUESTION 7

w <- rep(c(0,1,1,1,1), length(tennis_names)-1)
names(w) <- column_names[-(1:5)]
w <- w[names(w) %in% colnames(X2)]

# Variables don't need to be standardised as we have already designed the matrix so they are in the
# correct units

TennisLassoGLM1 <- cv.glmnet(X2, Training_Y2, family="binomial", alpha=1, penalty.factor=w, intercept=FALSE,
                             standardize = FALSE, thresh = 1e-14, maxit = 1e5)

min_lambda = TennisLassoGLM1$lambda.min
temp_coefs = coef(TennisLassoGLM1)
coefs3 <- as.vector(temp_coefs)
names(coefs3) <- rownames(temp_coefs)
coefs3 <- coefs3[-1]
print(paste("Minimum Lambda is:", min_lambda))


# QUESTION 8

non_zero_surface_terms = 0
for (surface_param in rownames(coefs3)[!rownames(coefs3) %in% tennis_names]) {
  print(coefs3[surface_param, ])
  print(surface_param)
  if (coefs3[surface_param, ] != 0) {
    non_zero_surface_terms = non_zero_surface_terms + 1
  }
}
print(paste("The number of non zero surface terms are:", non_zero_surface_terms))

LL5 = 0
for (i in 1:dim(X2)[1]) {
  row_X2 = X2[i, ]
  X2b = row_X2 %*% coefs3
  y = Training_Y2[i, 2]
  n = Training_Y2[i, 1] + Training_Y2[i, 2]
  mu = invlogit(X2b)
  term = log(choose(n, y)) + y*log(mu) + (n-y)*log(1-mu)
  LL5 = LL5 + term
}
LL5 = LL5*-1/dim(X2)[1]
print(paste("TennisLassoGLM1 Log Loss for Training Set:", LL5))

LL6 = 0
for (i in 1:dim(Test_X2)[1]) {
  row_Test_X2 = Test_X2[i, ]
  Test_X2b = row_Test_X2 %*% coefs3
  y = Training_Y2[i, 2]
  n = Training_Y2[i, 1] + Training_Y2[i, 2]
  mu = invlogit(Test_X2b)
  term = log(choose(n, y)) + y*log(mu) + (n-y)*log(1-mu)
  LL6 = LL6 + term
}
LL6 = LL6*-1/dim(Test_X2)[1]
print(paste("TennisGLM2 Log Loss for Test Set:", LL6))


# QUESTION 9

# QUESTION 11

portfolio <- function(v, coeffs) {
  Bet_Data_2015 = Tennis[(Tennis$Date >= as.Date("2015/01/01")) &
                           (Tennis$Date <= as.Date("2015/12/31")) & 
                           (!is.na(Tennis$Date)) & 
                           (!is.na(Tennis$B365W)) & 
                           (!is.na(Tennis$B365L)), ]
  num_bets = dim(Bet_Data_2015)[1]
  q = rep(0, 2*num_bets)
  Q = matrix(0, 2*num_bets, 2*num_bets)
  for (i in 1:num_bets) {
    match = Bet_Data_2015[i, ]
    winner = match$Winner
    loser = match$Loser
    bet_win = as.numeric(match$B365W)
    bet_loss = as.numeric(levels(match$B365L))[match$B365L]
    surface = match$Surface
    match_vec = rep(0, dim(X2)[2])
    names(match_vec) <- colnames(X2)
    if (winner %in% names(match_vec)) {
      match_vec[winner] = 1
      if (paste(winner, surface) %in% names(match_vec)) {
        match_vec[paste(winner, surface)] = 1
      }
    }
    if (loser %in% names(match_vec)) {
      match_vec[loser] = -1
      if (paste(loser, surface) %in% names(match_vec)) {
        match_vec[paste(loser, surface)] = -1
      }
    }
    mu_i = as.numeric(invlogit(match_vec %*% coeffs))
    
    q[2*i-1] = (bet_win - 1)*mu_i
    q[2*i] = (bet_loss - 1)*(1 - mu_i)
    
    A_i = matrix(c((bet_win - 1)^2, -(bet_win - 1)*(bet_loss - 1),
                   -(bet_win - 1)*(bet_loss - 1), (bet_loss - 1)^2),2,2, byrow = TRUE)
    Q[c(2*i-1, 2*i), c(2*i-1, 2*i)] = A_i*mu_i*(1-mu_i)
  }
  A_mat = rbind(matrix(rep(1, 2*num_bets), nrow=1), diag(2*num_bets))
  b_vec = c(1000, rep(0, 2*num_bets))
  # Note Q will always be singular so will find approximate solution using
  # Q + 10^-8
  sol = solve.QP(2*v*(Q)+diag(x=1e-8, 2*num_bets), q, t(A_mat), b_vec, meq=1)
  # Note Solve.QP isn't great as it may still have -ve amounts in the solution
  # To fix this we will set all values s.t they are >= 0, then rescale all values accordingly
  # Tried max, isn'n working well
  solu = sol$solution - min(sol$solution)
  solu = solu*1000/sum(solu)
  # print(paste("The Value of Solution for v =",v,"is:",v*t(solu)%*%Q%*%solu-t(q)%*%solu))
  # print(v*t(solu)%*%Q%*%solu)
  # print(t(q)%*%solu)
  return(solu)
}

# As could have probably been guessed the value of v which maximises the potential profit
# minimises the impact of the risk term
calc_portfolio_benefits <- function(pf) {
  Bet_Data_2015 = Tennis[(Tennis$Date >= as.Date("2015/01/01")) &
                           (Tennis$Date <= as.Date("2015/12/31")) & 
                           (!is.na(Tennis$Date)) & 
                           (!is.na(Tennis$B365W)) & 
                           (!is.na(Tennis$B365L)), ]
  num_bets = dim(Bet_Data_2015)[1]
  total = -1000
  for (i in 1:num_bets) {
    match = Bet_Data_2015[i, ]
    bet_win = as.numeric(match$B365W)
    total = total + bet_win*pf[2*i-1]
  }
  return(total)
}

v = 10^seq(-3, 3, length.out = 7) # 31)
profits1 = rep(0,length(v))
profits2 = rep(0,length(v))
profits3 = rep(0,length(v))
for (i in 1:length(v)) {
  profits1[i] <- calc_portfolio_benefits(portfolio(v[i], coefs2))
  profits2[i] <- calc_portfolio_benefits(portfolio(v[i], coefs3))
  profits3[i] <- calc_portfolio_benefits(portfolio(v[i], coef0))
}

plot(v, profits1,
     xlab="v", ylab="Total Profit",
     ylim=c(-1000, 0))
lines(v, profits1, col="black")
par(new=TRUE)
# points(v, profits2)
lines(v, profits2, col="red")
lines(v, profits3, col="blue")
legend("right", legend = c("TennisGLM2", "TennisLassoGLM1", "Basic"), 
       col=c("black", "red", "blue"), lty=c(1, 1), inset=0.05)


# QUESTION 13

kelly_fraction <- function(match, coeffs, rho, grid, debug) {
  winner = match$Winner
  loser = match$Loser
  bet_win = as.numeric(match$B365W)
  bet_loss = as.numeric(levels(match$B365L))[match$B365L]
  surface = match$Surface
  match_vec = rep(0, dim(X2)[2])
  names(match_vec) <- colnames(X2)
  if (winner %in% names(match_vec)) {
    match_vec[winner] = 1
    if (paste(winner, surface) %in% names(match_vec)) {
      match_vec[paste(winner, surface)] = 1
    }
  }
  if (loser %in% names(match_vec)) {
    match_vec[loser] = -1
    if (paste(loser, surface) %in% names(match_vec)) {
      match_vec[paste(loser, surface)] = -1
    }
  }
  mu_i = as.numeric(invlogit(match_vec %*% coeffs))
  opt_win_loss = c(0,0)
  opt_win_loss_val = -Inf
  for (Win in seq(0, 1, length.out=grid + 1)) {
    for (Loss in seq(0, 1, length.out=grid + 1)) {
      if (Win+Loss > 1) {
        break
      } else {
        val = log((1-Loss) + Win*(bet_win-1))*(mu_i) + log((1-Win) + Loss*(bet_loss-1))*(1-mu_i)
        if (!is.na(val) && val != -Inf && val > opt_win_loss_val) {
          opt_win_loss = c(Win, Loss)
          opt_win_loss_val = val
        }
      }
    }
  }
  return(rho*opt_win_loss)
}

eval_strat <- function(coeffs, rho, grid, bankroll){
  Bet_Data_2015 = Tennis[(Tennis$Date >= as.Date("2015/01/01")) &
                           (Tennis$Date <= as.Date("2015/12/31")) & 
                           (!is.na(Tennis$Date)) & 
                           (!is.na(Tennis$B365W)) & 
                           (!is.na(Tennis$B365L)), ]
  Bet_Data_2015 = Bet_Data_2015[order(Bet_Data_2015$Date), ]
  num_bets = dim(Bet_Data_2015)[1]
  b = rep(0, num_bets+1)
  b[1]=bankroll
  for (i in 1:num_bets) {
    match = Bet_Data_2015[i, ]
    fractions = kelly_fraction(match, coeffs, rho, grid, FALSE)
    bet_win = as.numeric(match$B365W)
    b[i+1] = b[i]*(1 - fractions[2] + fractions[1]*(bet_win-1))
  }
  return(c(b, as.Date("2015/1/1"), Bet_Data_2015$Date))
}

strat_results = eval_strat(coefs2, 0.1, 1000, 1000)
par(new = FALSE)
plot(as.Date(strat_results[-(1:length(strat_results)/2)], origin = "1970-01-01"), 
     strat_results[1:(length(strat_results)/2)], type="l",
     xlab = "Date", ylab = "Bankroll")
lines(as.Date(strat_results[-(1:length(strat_results)/2)], origin = "1970-01-01"), 
      strat_results[1:(length(strat_results)/2)])
