l = matrix(0, ncol=4, nrow=length(tennis_names)*(length(tennis_names)-1)/2)

index = 1
for (i in 1:(length(tennis_names)-1)) {
  for (j in (i+1):length(tennis_names)) {
    l[index, ] = c(tennis_names[i], tennis_names[j], 0, 0)
    index = index + 1
  }
}

for (i in 1:Training_length) {
  winner = Training_winners[i]
  loser = Training_losers[i]
  win_ind = which(tennis_names == winner)
  los_ind = which(tennis_names == loser)
  
  # if (any(apply(l[, c(1,2)], 1, function(x) {all(x==c(winner, loser))}))) {
  if (win_ind<los_ind) {
    ind = which(apply(l[, c(1,2)], 1, function(x) {all(x==c(winner, loser))}))
    l[ind, 3] = as.numeric(l[ind, 3]) + 1
  } else {
    ind = which(apply(l[, c(1,2)], 1, function(x) {all(x==c(loser, winner))}))
    l[ind, 4] = as.numeric(l[ind, 4]) + 1
  }
}
colnames(l) <- c("player1", "player2", "win1", "win2")
l = l[(l[, 3]!= 0) | (l[, 4]!= 0), ]

X_oof = matrix(0, ncol = length(tennis_names), nrow=dim(l)[1])
colnames(X_oof) = tennis_names # [-1]

for (i in 1:dim(l)[1]) {
  # if (l[i, 1] != tennis_names[1]) {
    X_oof[i, l[i, 1]] = 1
  # }
  X_oof[i, l[i, 2]] = -1
}

Y_oof = cbind(as.numeric(l[, 3]), as.numeric(l[, 4]))
