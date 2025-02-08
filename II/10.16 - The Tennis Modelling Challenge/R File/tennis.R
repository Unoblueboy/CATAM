library(glmnet)
library(quadprog)
library(latex2exp)
setwd('E:/Google Drive/School/Cambridge/Maths/II/04 - CATAM/10.16 - The Tennis Modelling Challenge/R File')
Tennis <- read.csv("http://www.damtp.cam.ac.uk/user/catam/data/II-10-16-2019-mensResults.csv")
Tennis$Date <- as.Date(Tennis$Date,format='%d/%m/%y')
lam_seq = c(exp(seq(log(1), log(1e-10), length.out = 1000)), 0)

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
Training_Y <- rep(c(1,0), length.out = Training_length)
Test_Y <- rep(c(1,0), length.out = Test_length)

# QUESTION 2

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

# Make X a sparse Matrix
X <- as(X, "sparseMatrix")

TennisGLM1 <- glmnet(X, Training_Y, lambda=0, intercept=FALSE, family="binomial", standardize = FALSE, thresh = 1e-14)
summary(TennisGLM1)

coefs <- as.vector(coef(TennisGLM1, s=0))
coefs <- coefs[2: length(coefs)]
names(coefs) <- colnames(X)
write.table(coefs, "question 2 coefs.txt", col.names = F, sep = " : ", quote = F)

LL1 = 0
for (i in 1:Training_length) {
  row_X = X[i, ]
  Xb = row_X %*% coefs
  term = Training_Y[i]*log(invlogit(Xb)) - (1-Training_Y[i])*K(Xb)
  LL1 = LL1 + term
}
LL1 = LL1*-1/Training_length
print(paste("TennisGLM1 Log Loss for Training Set:", LL1))


Test_X <- matrix(data=0, nrow=Test_length, ncol = length(tennis_names)-1)
colnames(Test_X) <- tennis_names[-1]

for (i in 1:Test_length) {
  winner = Test_winners[i]
  loser = Test_losers[i]
  y_factor = Test_Y[i]*2-1
  if (winner != tennis_names[1]) {
    Test_X[i, winner] = 1 * y_factor
  }
  if (loser != tennis_names[1]) {
    Test_X[i, loser] = -1 * y_factor
  }
}

LL2 = 0
for (i in 1:Test_length) {
  row_Test_X = Test_X[i, ]
  Test_Xb = row_Test_X %*% coefs
  term = Test_Y[i]*log(invlogit(Test_Xb)) - (1-Test_Y[i])*K(Test_Xb)
  LL2 = LL2 + term
}
LL2 = LL2*-1/Test_length
print(paste("TennisGLM1 Log Loss for Test Set:", LL2))

write(sprintf("Training Data Logistic Loss %.15f
Test Data Logistic Loss %.15f", LL1, LL2), file = "question 2 LogLoss.txt")

# QUESTION 3

d <- rep(0, 59)
names(d) <- tennis_names[-1]
d["Federer R."] = 1
d["Murray A."] = -1
pred <- t(d) %*% coefs
critval <- qnorm(0.5 + 0.68/2)

# W = matrix(0, nrow = Training_length, ncol = Training_length)
W_diag = rep(0, Training_length)
for (i in 1:Training_length) {
  theta = t(X[i, ]) %*% coefs
  W_diag[i] = invlogit(theta)*(1-invlogit(theta))
}
W = diag(W_diag)
Covar_mat = solve(t(X) %*% W %*% X)
se = sqrt(as.numeric(t(d) %*% Covar_mat %*% d))
CI_prob_Bounds <- invlogit(c(pred - critval*se, pred + critval*se))
CI_string = paste0("Confidence interval the for probability that Roger Federer beats Andy Murray is 
[",
                   CI_prob_Bounds[1],", ", CI_prob_Bounds[2], "]")
print(CI_string)

write(CI_string, file = "question 3 CI.txt")

# QUESTION 4

tennis_surfaces = levels(Tennis$Surface)
column_names = rep("",(length(tennis_names)-1)*(length(tennis_surfaces)))
X2 <- matrix(data=0, nrow=Training_length, ncol = length(column_names))
for (i in 2:length(tennis_names)) {
  column_names[4*(i-2)+1] <- tennis_names[i]
  for (j in 2:length(tennis_surfaces))
    column_names[4*(i-2)+(j)] <- paste(tennis_names[i], tennis_surfaces[j-1])
}
colnames(X2) <- column_names

Training_surfaces <- as.vector(Training_Data$Surface)
Test_surfaces <- as.vector(Test_Data$Surface)

for (i in 1:Training_length) {
  winner = Training_winners[i]
  loser = Training_losers[i]
  surface = Training_surfaces[i]
  y_factor = 2*Training_Y[i]-1
  if (winner != tennis_names[1]) {
    X2[i, winner] <- 1*y_factor
    if (surface != tennis_surfaces[4]) {
      X2[i, paste(winner, surface)] <- 1*y_factor
    }
  }
  if (loser != tennis_names[1]) {
    X2[i, loser] <- -1*y_factor
    if (surface != tennis_surfaces[4]) {
      X2[i, paste(loser, surface)] <- -1*y_factor
    }
  }
}
# Delete empty columns
X2 <- X2[, colSums(X2==0) != nrow(X2)]
# Make X2 a sparse Matrix
X2 <- as(X2, "sparseMatrix")

TennisGLM2 <- glmnet(X2, Training_Y, lambda=lam_seq, intercept=FALSE, 
                     family="binomial", standardize = FALSE, thresh = 1e-12, maxit = 1e5)
summary(TennisGLM2)


coefs2 <- as.vector(coef(TennisGLM2, s=0))
coefs2 <- coefs2[2: length(coefs2)]
names(coefs2) <- colnames(X2)

write.table(coefs2, "question 4 coefs.txt", col.names = F, sep = ": ", quote = F)

LL3 = 0
for (i in 1:Training_length) {
  row_X2 = X2[i, ]
  X2b = row_X2 %*% coefs2
  term = Training_Y[i]*log(invlogit(X2b)) - (1-Training_Y[i])*K(X2b)
  LL3 = LL3 + term
}
LL3 = LL3*-1/Training_length
print(paste("TennisGLM2 Log Loss for Training Set:", LL3))

Test_X2 <- matrix(data=0, nrow=Test_length, ncol = (length(tennis_names)-1)*(length(tennis_surfaces)))
colnames(Test_X2) <- column_names

for (i in 1:Test_length) {
  winner = Test_winners[i]
  loser = Test_losers[i]
  surface = Test_surfaces[i]
  y_factor = 2*Test_Y[i]-1
  if (winner != tennis_names[1]) {
    Test_X2[i, winner] <- 1*y_factor
    if (surface != tennis_surfaces[4]) {
      Test_X2[i, paste(winner, surface)] <- 1*y_factor
    }
  }
  if (loser != tennis_names[1]) {
    Test_X2[i, loser] <- -1*y_factor
    if (surface != tennis_surfaces[4]) {
      Test_X2[i, paste(loser, surface)] <- -1*y_factor
    }
  }
}
# Delete columns not in X2
Test_X2 <- Test_X2[, colnames(Test_X2) %in% colnames(X2)]

LL4 = 0
for (i in 1:Test_length) {
  row_Test_X2 = Test_X2[i, ]
  Test_X2b = row_Test_X2 %*% coefs2
  term = Test_Y[i]*log(invlogit(Test_X2b)) - (1-Test_Y[i])*K(Test_X2b)
  LL4 = LL4 + term
}
LL4 = LL4*-1/Test_length
print(paste("TennisGLM2 Log Loss for Test Set:", LL4))

write(sprintf("Training Data Logistic Loss %.15f
Test Data Logistic Loss %.15f", LL3, LL4), file = "question 4 LogLoss.txt")

# Question 5

test_stat = (LL1*Training_length*2-LL3*Training_length*2)
# Variance known to be 1 so use chi squared dist
p_val = 1-pchisq(test_stat, dim(X2)[2]-dim(X)[2])# length(coefs2)-length(coefs))
print(paste("P-value for the likelihood ratio test is:",p_val))
write(sprintf("Test Statistic: %.15f
P-value for the likelihood ratio test: %.15f", test_stat, p_val), file = "question 5 results.txt")

# Question 7

w <- rep(c(0,1,1,1), length(tennis_names)-1)
names(w) <- column_names
w <- w[names(w) %in% colnames(X2)]

TennisLassoGLM1 <- cv.glmnet(X2, Training_Y, family="binomial", alpha=1, penalty.factor=w, intercept=FALSE,
                             standardize = FALSE, thresh = 1e-14, maxit = 1e5, nfolds = 10)
min_lambda = TennisLassoGLM1$lambda.min
temp_coefs = coef(TennisLassoGLM1, s=min_lambda)
coefs3 <- as.vector(temp_coefs)
names(coefs3) <- rownames(temp_coefs)
coefs3 <- coefs3[-1]
print(paste("Minimum Lambda is:", min_lambda))

write.table(coefs3, "question 7 coefs.txt", col.names = F, sep = ": ", quote = F)
write(paste("Lambda which minimises the mean cross-validated error is:", min_lambda), file = "question 7 minLam.txt")

# Question 8

non_zero_surface_terms = 0
for (surface_param in names(coefs3)[!names(coefs3) %in% tennis_names]) {
  if (coefs3[surface_param] != 0) {
    non_zero_surface_terms = non_zero_surface_terms + 1
  }
}
print(paste("The number of non zero surface terms are:", non_zero_surface_terms))
write(paste("The number of non zero surface terms are:", non_zero_surface_terms), file="question 8 nonzero.txt")


LL5 = 0
for (i in 1:Training_length) {
  row_X2 = X2[i, ]
  X2b = row_X2 %*% coefs3
  term = Training_Y[i]*log(invlogit(X2b)) - (1-Training_Y[i])*K(X2b)
  LL5 = LL5 + term
}
LL5 = LL5*-1/Training_length
print(paste("TennisLassoGLM1 Log Loss for Training Set:", LL5))

LL6 = 0
for (i in 1:Test_length) {
  row_Test_X2 = Test_X2[i, ]
  Test_X2b = row_Test_X2 %*% coefs3
  term = Test_Y[i]*log(invlogit(Test_X2b)) - (1-Test_Y[i])*K(Test_X2b)
  LL6 = LL6 + term
}
LL6 = LL6*-1/Test_length
print(paste("TennisLassoGLM1 Log Loss for Test Set:", LL6))

write(sprintf("Training Data Logistic Loss %.15f
Test Data Logistic Loss %.15f", LL5, LL6), file = "question 8 LogLoss.txt")


# QUESTION 9

# Create the Design Matrix

column_names = rep("", 2*(length(tennis_names)-1))
X3 <- matrix(data=0, nrow=Training_length, ncol = length(column_names))
Test_X3 <- matrix(data=0, nrow=Test_length, ncol = length(column_names))
for (i in 2:length(tennis_names)) {
  column_names[2*i-3] <- tennis_names[i]
  column_names[2*i-2] <- paste(tennis_names[i], "Yrly.")
}
colnames(X3) <- column_names
colnames(Test_X3) <- column_names

Training_Date <- Training_Data$Date
Test_Date <- Test_Data$Date

for (i in 1:Training_length) {
  winner = Training_winners[i]
  loser = Training_losers[i]
  year = as.numeric(format(Training_Date[i],'%Y'))
  t = (year-2000)/(2014-2000)
  y_factor = 2*Training_Y[i]-1
  if (winner != tennis_names[1]) {
    X3[i, winner] <- 1*y_factor
    X3[i, paste(winner, "Yrly.")] <- t*y_factor
  }
  if (loser != tennis_names[1]) {
    X3[i, loser] <- -1*y_factor
    X3[i, paste(loser, "Yrly.")] <- -t*y_factor
  }
}

# Delete empty columns
X3 <- X3[, colSums(X3==0) != nrow(X3)]
# Make X3 into a sparse Matrix
X3 <- as(X3, "sparseMatrix")

TennisLassoGLM2 <- glmnet(X3, Training_Y, family="binomial", lambda = 0, 
                          intercept=FALSE, standardize = FALSE, thresh = 1e-9, maxit = 1e5)
temp_coefs = coef(TennisLassoGLM2, s = 0)
coefs4 <- as.vector(temp_coefs)
names(coefs4) <- rownames(temp_coefs)
coefs4 <- coefs4[-1]

write.table(coefs4, "question 9 coefs 1.txt", col.names = F, sep = ": ", quote = F)

for (i in 1:Test_length) {
  winner = Test_winners[i]
  loser = Test_losers[i]
  year = as.numeric(format(Test_Date[i],'%Y'))
  t = (year-2000)/(2014-2000)
  y_factor = 2*Test_Y[i]-1
  if (winner != tennis_names[1]) {
    Test_X3[i, winner] <- 1*y_factor
    Test_X3[i, paste(winner, "Yrly.")] <- t*y_factor
  }
  if (loser != tennis_names[1]) {
    Test_X3[i, loser] <- -1*y_factor
    Test_X3[i, paste(loser, "Yrly.")] <- -t*y_factor
  }
}
# Delete columns not in X3
Test_X3 <- Test_X3[, colnames(Test_X3) %in% colnames(X3)]

LL7 = 0
for (i in 1:Training_length) {
  row_X3 = X3[i, ]
  X3b = row_X3 %*% coefs4
  term = Training_Y[i]*log(invlogit(X3b)) - (1-Training_Y[i])*K(X3b)
  LL7 = LL7 + term
}
LL7 = LL7*-1/Training_length
print(paste("TennisLassoGLM2 Log Loss for Training Set:", LL7))

LL8 = 0
for (i in 1:Test_length) {
  row_X3 = X3[i, ]
  X3b = row_X3 %*% coefs4
  term = Test_Y[i]*log(invlogit(X3b)) - (1-Test_Y[i])*K(X3b)
  LL8 = LL8 + term
}
LL8 = LL8*-1/Test_length
print(paste("TennisLassoGLM2 Log Loss for Test Set:", LL8))

write(sprintf("Training Data Logistic Loss %.15f
Test Data Logistic Loss %.15f", LL7, LL8), file = "question 9 LogLoss noWeights.txt")


#Weights
w2 <- rep(c(0,1), length(tennis_names)-1)
names(w2) <- column_names
w2 <- w2[names(w2) %in% colnames(X3)]

TennisLassoGLM3 <- cv.glmnet(X3, Training_Y, family="binomial", alpha=1, penalty.factor = w2, 
                             intercept=FALSE, standardize = FALSE, thresh = 1e-9, maxit = 1e5)

temp_coefs = coef(TennisLassoGLM3, s = "lambda.min")
coefs5 <- as.vector(temp_coefs)
names(coefs5) <- rownames(temp_coefs)
coefs5 <- coefs5[-1]

write.table(coefs5, "question 9 coefs 2.txt", col.names = F, sep = ": ", quote = F)

LL9 = 0
for (i in 1:Training_length) {
  row_X3 = X3[i, ]
  X3b = row_X3 %*% coefs5
  term = Training_Y[i]*log(invlogit(X3b)) - (1-Training_Y[i])*K(X3b)
  LL9 = LL9 + term
}
LL9 = LL9*-1/Training_length
print(paste("TennisLassoGLM3 Log Loss for Training Set:", LL9))

LL10 = 0
for (i in 1:Test_length) {
  row_X3 = X3[i, ]
  X3b = row_X3 %*% coefs5
  term = Test_Y[i]*log(invlogit(X3b)) - (1-Test_Y[i])*K(X3b)
  LL10 = LL10 + term
}
LL10 = LL10*-1/Test_length
print(paste("TennisLassoGLM2 Log Loss for Test Set:", LL10))

write(sprintf("Training Data Logistic Loss %.15f
Test Data Logistic Loss %.15f", LL9, LL10), file = "question 9 LogLoss Weights.txt")





Fed_Nad_Mat = matrix(0, nrow = 17, ncol=dim(X3)[2])
colnames(Fed_Nad_Mat) <- colnames(X3)
for (year in 2000:2016) {
  Fed_Nad_Mat[year - 1999, "Federer R."] = 1
  Fed_Nad_Mat[year - 1999, "Federer R. Yrly."] = (year-2000)/(2014-2000)
  Fed_Nad_Mat[year - 1999, "Nadal R."] = -1
  Fed_Nad_Mat[year - 1999, "Nadal R. Yrly."] = -(year-2000)/(2014-2000)
}
pdf("federer_nadal_probs.pdf")
plot(2000:2016, invlogit(Fed_Nad_Mat %*% coefs4), col = "blue", type="p", pch=19,
     xlab="Year", ylab="Probability Federer beats Nadal", ylim=c(0, 1))
lines(2000:2016, invlogit(Fed_Nad_Mat %*% coefs4), col = "blue", lty=1)
par(new=T)
plot(2000:2016, invlogit(Fed_Nad_Mat %*% coefs5), col = "red", type="p", pch=18,
     xlab="Year", ylab="Probability Federer beats Nadal", ylim=c(0, 1))
lines(2000:2016, invlogit(Fed_Nad_Mat %*% coefs5), col = "red", lty=2)
legend("topright", legend=c("No Lasso penalty", "Lasso Penalty"), col=c("blue", "red"),
       lty=1:2, pch = 19:18)
dev.off()

# We will compare all the tests by looking at both their log loss as well as their AIC for both the training and Test Data
# We will divide the AIC by the number of observations to get an approximation of the 2 * Kullback-Liebler Divergence
# So KB ~ AIC / 2n = Logloss*-1/n + p/n, where p is the number of parameters
KB1 = LL1 + length(coefs)/Training_length
KB2 = LL2 + length(coefs)/Test_length
KB3 = LL3 + length(coefs2)/Training_length
KB4 = LL4 + length(coefs2)/Test_length
KB5 = LL5 + length(coefs3)/Training_length
KB6 = LL6 + length(coefs3)/Test_length
KB7 = LL7 + length(coefs4)/Training_length
KB8 = LL8 + length(coefs4)/Test_length
KB9 = LL9 + length(coefs5)/Training_length
KB10 = LL10 + length(coefs5)/Test_length
print("The following are the log-loss and approximations of the Kullback Liebler Divergence of the models with the true distribution")
print("TennisGLM1")
print(paste("KB: Training Set:", KB1, "Test Set:", KB2))
print(paste("LL: Training Set:", LL1, "Test Set:", LL2))
print("TennisGLM2")
print(paste("KB: Training Set:", KB3, "Test Set:", KB4))
print(paste("LL: Training Set:", LL3, "Test Set:", LL4))
print("TennisLassoGLM1")
print(paste("KB: Training Set:", KB5, "Test Set:", KB6))
print(paste("LL: Training Set:", LL5, "Test Set:", LL6))
print("TennisLassoGLM2")
print(paste("KB: Training Set:", KB7, "Test Set:", KB8))
print(paste("LL: Training Set:", LL7, "Test Set:", LL8))
print("TennisLassoGLM3")
print(paste("KB: Training Set:", KB9, "Test Set:", KB10))
print(paste("LL: Training Set:", LL9, "Test Set:", LL10))

write.table(c("Model 1",
              paste("KB: Training Set:", KB1),
              paste("LL: Training Set:", LL1, "Test Set:", LL2),
              "Model 2",
              paste("KB: Training Set:", KB3),
              paste("LL: Training Set:", LL3, "Test Set:", LL4),
              "Model 3",
              paste("KB: Training Set:", KB5),
              paste("LL: Training Set:", LL5, "Test Set:", LL6),
              "Model 4",
              paste("KB: Training Set:", KB7),
              paste("LL: Training Set:", LL7, "Test Set:", LL8),
              "Model 5",
              paste("KB: Training Set:", KB9),
              paste("LL: Training Set:", LL9, "Test Set:", LL10)), file = "question 9 results.txt", sep = "\n", 
            row.names = F, col.names = F, quote = F)


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
    A_i = matrix(c((bet_win - 1)^2              , -(bet_win - 1)*(bet_loss - 1),
                   -(bet_win - 1)*(bet_loss - 1), (bet_loss - 1)^2             ),2,2, byrow = TRUE)
    Q[c(2*i-1, 2*i), c(2*i-1, 2*i)] = A_i*mu_i*(1-mu_i)
  }
  A_mat = rbind(matrix(rep(1, 2*num_bets), nrow=1), diag(2*num_bets))
  b_vec = c(1000, rep(0, 2*num_bets))
  # Note Q will always be singular so will find approximate solution using Q + 10^-8
  sol = solve.QP(2*v*(Q)+diag(x=1e-8, 2*num_bets), q, t(A_mat), b_vec, meq=1)
  # Note that Solve.QP isn't great as it may still have -ve entries in the solution
  # To fix this we will translate all values s.t they are >= 0, then rescale all values accordingly
  solu = sol$solution - min(sol$solution)
  solu = solu*1000/sum(solu)
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

v = 10^seq(-3, 3, length.out =7)
profits1 = rep(0,length(v))
profits2 = rep(0,length(v))
for (i in 1:length(v)) {
  profits1[i] <- calc_portfolio_benefits(portfolio(v[i], coefs2))
  profits2[i] <- calc_portfolio_benefits(portfolio(v[i], coefs3))
}

pdf("markowitz_portfolio.pdf")
plot(log10(v), profits1, type = "p", col="blue", pch=19,
     xlab=TeX("$\\log_{10}\\nu$"), ylab="Total Profit",
     ylim=c(-1000, 0))
lines(log10(v), profits1, col="blue", lty=1)
par(new=TRUE)
# points(v, profits2)
plot(log10(v), profits2, type = "p", col="red", pch=18,
     xlab=TeX("$\\log_{10}\\nu$"), ylab="Total Profit",
     ylim=c(-1000, 0))
lines(log10(v), profits2, col="red", lty=2)
legend("topright", legend = c("Model 2", "Model 3"), col=c("blue", "red"), lty=c(1, 2), pch=c(19, 18), inset=0.1)
dev.off()


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
pdf("question 13 bankroll.pdf")
plot(as.Date(strat_results[-(1:length(strat_results)/2)], origin = "1970-01-01"), 
     strat_results[1:(length(strat_results)/2)], type="p",
     xlab = "Date", ylab = "Bankroll", pch=20, col="black")
lines(as.Date(strat_results[-(1:length(strat_results)/2)], origin = "1970-01-01"), 
      strat_results[1:(length(strat_results)/2)], lty=1, col="black")
dev.off()
