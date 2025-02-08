install.packages('latex2exp')
library(latex2exp)

setwd('E:/Google Drive/School/Cambridge/Maths/II/04 - CATAM/10.09 - Markov Chain Monte Carlo/R File')
Football <- read.csv("http://www.damtp.cam.ac.uk/user/catam/data/II-10-9-2019football.csv")
rownames(Football) <- Football[, 1]
Football <- Football[, -1]

# QUESTION 3
print("Question 3")


params <- rep(0, 5)
names(params) <- c("theta_0", "alpha_0", "beta_0", "mu_0", "tau_0")
params["sigma_0"] = 10
params["alpha_0"] = 1e-5
params["beta_0"] = 1e-3
params["mu_0"] = 60
params["tau_0"] = 20

gibbs_sampler <- function(data, start, iter, d_params) {
  # note that in the data the columns are different times
  # while the rows are different teams
  num_teams = dim(data)[1]
  num_years = dim(data)[2]
  
  sum_data_over_time = rowSums(data)
  
  result = matrix(c(start, rep(0,iter*(2*num_teams+1))),  
                  ncol=2*num_teams+1, byrow = TRUE)
  for (i in 1:iter) {
    cur_sol = result[i, ]
    # Calculate all mu's
    for (k1 in 1:num_teams) {
      sigma_k_squared = cur_sol[num_teams+k1]
      dist_var = 1/(num_years*sigma_k_squared^(-1)+d_params["sigma_0"]^(-2))
      dist_mu = (sum_data_over_time[k1]*sigma_k_squared^(-1) + 
                   d_params["sigma_0"]^(-2)*cur_sol[2*num_teams+1])*dist_var
      cur_sol[k1] = rnorm(1, mean = dist_mu, sd = sqrt(dist_var))
    }
    # Calculate all sigma
    faux_var = rowSums((data-cur_sol[1:num_teams])^2)
    for (k2 in 1:num_teams) {
      alpha_i = d_params["alpha_0"] + num_years/2
      beta_i = d_params["beta_0"] + (1/2)*faux_var[k2]
      cur_sol[num_teams + k2] = 1/rgamma(1, shape = alpha_i, rate = beta_i)
    }
    # calculate theta
    dist_var = 1/(num_teams*d_params["sigma_0"]^(-2)+d_params["tau_0"]^(-2))
    dist_mu = (d_params["sigma_0"]^(-2)*sum(cur_sol[1:num_teams])-
                 d_params["mu_0"]*d_params["tau_0"]^(-2))*dist_var
    cur_sol[2*num_teams + 1] = rnorm(1, mean = dist_mu, sd = sqrt(dist_var))
    result[i+1, ] = cur_sol
  }
  # Just for my own readability
  result_colnames = c(
    unname(sapply(rownames(data), function(x){paste(x, "mean")})),
    unname(sapply(rownames(data), function(x){paste(x, "variances")})),
    "theta"
  )
  colnames(result) <- result_colnames
  return(result)
}


# We made the starting point estimations of each of the parameters
# i.e starting means is the mean score of team 
starting_means = rowMeans(Football)
# starting vars is the vars of each team according to the data
starting_vars = rowMeans(Football^2)-starting_means^2
# starting theta is the mean of all scores
starting_theta = mean(starting_means) # valid as all rows have same number of entries
gibbs_start = c(starting_means, starting_vars, starting_theta)
gibbs_start = unname(gibbs_start)

iterations = 10000

result = gibbs_sampler(Football, gibbs_start, iterations, params)
K = dim(Football)[1]

# QUESTION 4
print("Question 4")



post_mean = colMeans(result)
mean_strings = paste("The posterior mean of",names(post_mean), ":", post_mean)
write.table(mean_strings, file = "means.txt", sep="\n", row.names = FALSE, col.names = FALSE, quote=FALSE)
post_theta_samples = result[, 2*K+1]
pdf("post_theta_dist.pdf")
post_theta_hist_info = hist(post_theta_samples, freq = F, xlab=TeX("$\\theta$"), main = TeX("Histogram of the posterior distribution of $\\theta$"), ylim=c(0,0.14))
dev.off()

# QUESTION 5
print("Question 5")

apply(result, 1, function(x) {as.numeric(x[1]>x[2*K+1])})
prob_mu_gt_theta = colMeans(t(apply(result, 1, function(x) {x[1:K]>x[2*K+1]})))
names(prob_mu_gt_theta) <- rownames(Football)
prob_strings = paste("Probability mean score of", names(prob_mu_gt_theta), 
                     "is greater than average is:", prob_mu_gt_theta)
write.table(prob_strings, file = "probs.txt", sep="\n", row.names = FALSE, col.names = FALSE, quote=FALSE)

# QUESTION 6
print("Question 6")

N = 3000
trials = 100
trial_results = array(0, dim = c(N+1, 2*K+1, trials))
for (i in 1:trials) {
  trial_results[,,i] = gibbs_sampler(Football, gibbs_start, N, params)
}

# Points represents the number of points I want to sample e.g. If I want to sample at
# 5 different values of N from N = 2 to 1001, then points = 5
points = 40
var_mat_probs = matrix(0, nrow = points, ncol = K)
var_mat_mus = matrix(0, nrow = points, ncol = K)
colnames(var_mat_probs) <- rownames(Football)
# n_seq is a list representing the stopping point of several trials that is to simulate 5 
# trials with N = 5, 10, 15, 20, can do 1 trial of length 20, then just look at the data points
# up until N = 5, 10, 15, 20
n_seq = floor(seq(2,N+1, length.out = points))

for (i in 1:points) {
  n = n_seq[i]
  # Here we are approximating P(mu_k > theta) with the given data (i.e iterations 1 to n)
  probs_mu_gt_theta = apply(trial_results[1:n,,], 3, function(x) {
  colMeans(t(apply(x, 1, function(z) {
    z[1:K]>z[2*K+1]
    })))
  })
  # Here we are approximating mu_k with the given data (i.e iterations 1 to n)
  mus = apply(trial_results[1:n,,], 3, function(x) {
    colMeans(x[, 1:K])
  })
  
  var_probs_mu_gt_theta = rowMeans(probs_mu_gt_theta^2)-rowMeans(probs_mu_gt_theta)^2
  var_mat_probs[i,] = var_probs_mu_gt_theta
  var_mus = rowMeans(mus^2)-rowMeans(mus)^2
  var_mat_mus[i,] = var_mus
}

# We could look at indvidual team k, but we could also look at all teams then consider the 
# maximum variance, thus looking at some kind of uniform convergence
pdf("var_probs.pdf")
max_var_probs = apply(var_mat_probs, 1, max)
plot(n_seq, max_var_probs, type="l",
     xlab="n", ylab="Maximum Variance over all teams (probabilities)")
dev.off()

pdf("var_means.pdf")
max_var_mus = apply(var_mat_mus, 1, max)
plot(n_seq, max_var_mus, type="l",
     xlab="n", ylab="Maximum Variance over all teams (means)")
dev.off()


# Now to investigate the relation between the variance and n
pdf("log_var_probs.pdf")
log_max_var_probs = log(max_var_probs)
plot(n_seq, log_max_var_probs, type="l",
     xlab="n", ylab="Log Maximum Variance over all teams (probabilities)")
dev.off()

pdf("log_var_means.pdf")
log_max_var_mus = log(max_var_mus)
plot(n_seq, log_max_var_mus, type="l",
     xlab="n", ylab="Log Maximum Variance over all teams (means)")
dev.off()


# QUESTION 7
print("Question 7")


min_M = 100
max_M = 3000
num_M = 5
M_seq = floor(seq(min_M, max_M, length.out = num_M))
N = 1000
trials = 100
trial_results = array(0, dim = c(N+max_M+1, 2*K+1, trials))
for (i in 1:trials) {
  trial_results[,,i] = gibbs_sampler(Football, gibbs_start, N+max_M, params)
}

points = 40
var_mat_probs = array(0, dim = c(points, num_M, K))
var_mat_mus = array(0, dim = c(points, num_M, K))
n_seq = floor(seq(2,N+1, length.out = points))
for (i in 1:points) {
  for (j in 1:num_M) {
    n = n_seq[i]
    m = M_seq[j]
    probs_mu_gt_theta = apply(trial_results[(m+1):(m+n),,], 3, function(x) {
      colMeans(t(apply(x, 1, function(z) {
        z[1:K]>z[2*K+1]
      })))
    })
    mus = apply(trial_results[(m+1):(m+n),,], 3, function(x) {
      colMeans(x[, 1:K])
    })
    var_probs_mu_gt_theta = rowMeans(probs_mu_gt_theta^2)-rowMeans(probs_mu_gt_theta)^2
    var_mat_probs[i, j, ] = var_probs_mu_gt_theta
    var_mus = rowMeans(mus^2)-rowMeans(mus)^2
    var_mat_mus[i, j, ] = var_mus
  }
}

# We could look at indvidual team k, but we could also look at all teams then consider the 
# maximum variance, thus looking at some kind of uniform convergence
cols = rainbow(num_M)
pdf("var_probs_2.pdf")
matplot(t(matrix(rep(n_seq,num_M), nrow=num_M, byrow = T)), 
        apply(var_mat_probs, 2, function(x) {apply(x, 1, max)}), type="l",
        col = cols, xlab="n", ylab="Maximum Variance over all teams (probabilities)")
legend("topright", legend = M_seq, col=cols, lwd=2, ncol=2, cex=0.8, title = "Number, M,  of pre-iterations")
dev.off()

pdf("var_means_2.pdf")
matplot(t(matrix(rep(n_seq,num_M), nrow=num_M, byrow = T)), 
        apply(var_mat_mus, 2, function(x) {apply(x, 1, max)}), type="l",
     xlab="n", ylab="Maximum Variance over all teams (means)")
legend("topright", legend = M_seq, col=cols, lwd=2, ncol=2, cex=0.8, title = "Number, M,  of pre-iterations")
dev.off()


min_M = 1
max_M = 1000
num_M = 100
M_seq = floor(seq(min_M, max_M, length.out = num_M))
N = 100
trials = 100
trial_results = array(0, dim = c(N+max_M+1, 2*K+1, trials))
for (i in 1:trials) {
  trial_results[,,i] = gibbs_sampler(Football, gibbs_start, N+max_M, params)
}

points = 4
var_mat_probs = array(0, dim = c(points, num_M, K))
var_mat_mus = array(0, dim = c(points, num_M, K))
# colnames(var_mat_probs) <- rownames(Football)
n_seq = round(exp(seq(log(2),log(N), length.out=points)))
for (i in 1:points) {
  for (j in 1:num_M) {
    n = n_seq[i]
    m = M_seq[j]
    probs_mu_gt_theta = apply(trial_results[(m+1):(m+n),,], 3, function(x) {
      colMeans(t(apply(x, 1, function(z) {
        z[1:K]>z[2*K+1]
      })))
    })
    mus = apply(trial_results[(m+1):(m+n),,], 3, function(x) {
      colMeans(x[, 1:K])
    })
    var_probs_mu_gt_theta = rowMeans(probs_mu_gt_theta^2)-rowMeans(probs_mu_gt_theta)^2
    var_mat_probs[i, j, ] = var_probs_mu_gt_theta
    var_mus = rowMeans(mus^2)-rowMeans(mus)^2
    var_mat_mus[i, j, ] = var_mus
  }
}

cols = rainbow(points)
pdf("var_probs_3.pdf")
matplot(M_seq, 
        apply(var_mat_probs, 1, function(x) {apply(x, 1, max)}), type=rep("p", points),
        col = cols, xlab="m", ylab="Maximum Variance over all teams (probabilities)", pch=rep(1, points))
legend("right", bg="white", legend = n_seq, col=cols, pch=rep(1, points), ncol=2, cex=0.8, title = "Number, n, of iterations")
dev.off()

pdf("var_means_3.pdf")
matplot(M_seq, 
        apply(var_mat_mus, 1, function(x) {apply(x, 1, max)}), type=rep("p", points),
        col = cols, xlab="m", ylab="Maximum Variance over all teams (means)", pch=rep(1, points))
legend("topright", bg="white", legend = n_seq, col=cols, pch=rep(1, points), ncol=2, cex=0.8, title = "Number, n, of iterations")
dev.off()



#QUESTION 8
print("Question 8")
num_start_points = 200
iterations=200
max_pre_iterations = 500
num_theta_analysis = 5

starting_points = matrix(
  c(runif(K*num_start_points, 0, 120),
    runif(K*num_start_points, 0.001,40),
    runif(num_start_points, 0, 120)), nrow = 2*K+1, byrow = TRUE)

all_data = array(0, dim=c(num_start_points, max_pre_iterations+ iterations+1, 2*K+1))

for (i in 1:num_start_points) {
  data = gibbs_sampler(Football, starting_points[, i], max_pre_iterations + iterations, params)
  all_data[i, , ]=data
}

write.table(all_data[1:5,1:5 ,2*K+1 ], file = "post_theta_iters_1.txt", row.names=FALSE, col.names=FALSE)
write.table(all_data[1:5,(iterations - 3):(1+iterations) ,2*K+1 ], file = "post_theta_iters_2.txt", row.names=FALSE, col.names=FALSE)

pdf("theta_dist_conv.pdf")
plot(post_theta_hist_info, freq=FALSE, ylim = c(0, 0.14), 
     xlab=TeX("$\\theta$"), main = TeX("Histogram of the posterior distribution of $\\theta$"))
cols = rainbow(num_theta_analysis)
for (i in 1:num_theta_analysis) {
  break_points = seq(floor(min(all_data[i,1:(1+iterations) , 2*K+1])) - (floor(min(all_data[i,1:(1+iterations) , 2*K+1]))%%2),
                     ceiling(max(all_data[i,1:(1+iterations) , 2*K+1])) + (ceiling(max(all_data[i,1:(1+iterations) , 2*K+1]))%%2),
                     by=2)
  hist_data = hist(all_data[i,1:(1+iterations), 2*K+1], plot = FALSE, breaks = break_points)
  lines(hist_data$mids, hist_data$density, col=cols[i])
}
dev.off()

means = t(apply(all_data[,1:(1+iterations),], 1, function(x) {colMeans(x)}))
probs = t(apply(all_data[,1:(1+iterations),], 1, function(x) {
  rowMeans(
    apply(x, 1, function(z) {
      z[1:K] > z[2*K + 1]
    })
  )
}))

all_estimates = cbind(means, probs)
all_variance = colMeans(all_estimates^2)-colMeans(all_estimates)^2
max_sds = c(max(all_variance[1:K]),
             max(all_variance[(K+1):(2*K)]),
             max(all_variance[2*K+1]),
             max(all_variance[(2*K+2):(3*K+1)]))^0.5

write.table(c(paste("Maximum standard deviation for posterior mean of mu_k is", max_sds[1]),
              paste("Maximum standard deviation for posterior mean of sigma_k^2 is", max_sds[2]),
              paste("Maximum standard deviation for posterior mean of theta is", max_sds[3]),
              paste("Maximum standard deviation for posterior probability mu_k > theta is", max_sds[4]))
          ,file = "max_estimate_sds.txt", sep="\n", row.names = FALSE, col.names = FALSE, quote=FALSE)

write.table(t(rbind(colMeans(all_estimates)[(K+1):(2*K)], all_variance[(K+1):(2*K)]^0.5)),
            file = "sigma_sds.txt", sep=" ", row.names = FALSE, col.names = c("Expected value", "Standard Deviation"), quote=FALSE)

pre_iter_max_vars = matrix(rep(c(0,0,0,0), max_pre_iterations+1), nrow=max_pre_iterations+1)

for (i in 0:max_pre_iterations) {
  means = t(apply(all_data[,(i+1):(i+1+iterations),], 1, function(x) {colMeans(x)}))
  probs = t(apply(all_data[,(i+1):(i+1+iterations),], 1, function(x) {
    rowMeans(
      apply(x, 1, function(z) {
        z[1:K] > z[2*K + 1]
      })
    )
  }))
  
  all_estimates = cbind(means, probs)
  all_variance = colMeans(all_estimates^2)-colMeans(all_estimates)^2
  max_vars = c(max(all_variance[1:K]),
              max(all_variance[(K+1):(2*K)]),
              max(all_variance[2*K+1]),
              max(all_variance[(2*K+2):(3*K+1)]))
  pre_iter_max_vars[i+1, ] = max_vars
}

pdf("max_var_estimates_1.pdf")
plot(0:max_pre_iterations, pre_iter_max_vars[, 1], ylim=c(0,0.7), pch=1, xlab="m", ylab="Maximum Variance over all teams (means)")
dev.off()
pdf("max_var_estimates_2.pdf")
plot(0:max_pre_iterations, pre_iter_max_vars[, 2], ylim=c(0,5000), pch=1, xlab="m", ylab="Maximum Variance over all teams (variance)")
dev.off()
pdf("max_var_estimates_3.pdf")
plot(0:max_pre_iterations, pre_iter_max_vars[, 3], ylim=c(0,0.1), pch=1, xlab="m", ylab="Variance over all teams (theta)")
dev.off()
pdf("max_var_estimates_4.pdf")
plot(0:max_pre_iterations, pre_iter_max_vars[, 4], ylim=c(0,0.002), pch=1, xlab="m", ylab="Maximum Variance over all teams (probs)")
dev.off()
