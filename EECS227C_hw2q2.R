library(rootSolve)
library(ggplot2)

# Part i
set.seed(1)
A <- matrix(data = rnorm(400, 0, 1), nrow = 20, ncol = 20)
f <- function(x) {
  return(sum((A %*% x) ^2))
}

x1 <- rep(1, 20)

gd <- function(f, x1, iters, eta) {
  xn <- x1
  obj <- rep(0, iters)
  obj[1] <- f(x1)
  for (t in 2:iters) {
    g <- gradient(f, xn)
    xn <- xn - eta * c(g)
    obj[t] <- min(f(xn), obj[t-1])
  }
  return(obj)
}

agd <- function(f, x1, iters, eta) {
  xn1 <- x1
  xn2 <- x1
  obj <- rep(0, iters)
  obj[1] <- f(x1)
  for (t in 2:iters) {
    y <- xn2 + t/(t+3)*(xn2 - xn1)
    g <- gradient(f, y)
    xn1 <- xn2
    xn2 <- y - eta * c(g)
    obj[t] <- min(f(xn2), obj[t-1])
  }
  return(obj)
}

iters <- 5000
res_gd <- gd(f, x1, iters, 0.005)
res_agd <- agd(f, x1, iters, 0.005)
f_min <- min(res_agd)
ggplot(data = data.frame(`Log_suboptimality_gap` = -c(log(res_gd),
                                                     log(res_agd)),
                         iter = log(c(1:5000, 1:5000)),
                         algo = c(rep('gd', 5000), rep('agd', 5000)))) + 
  geom_line(aes_string(y = "-log_suboptimality_gap", x = "log-time", col = "algo"))


# Part ii
A <- matrix(nrow = 20, ncol = 20)
for(i in 1:20) {
  for(j in 1:20) {
    if (i == j) {
      A[i, j] <- 2
    } else if(abs(i-j) == 1) {
      A[i, j] <- -1
    } else {
      A[i, j] <- 0
    }
  }
}

f <- function(x) {
  return(sum((x %*% (A %*% x)) ^2))
}

iters <- 5000
res_gd <- gd(f, x1, iters, 0.05)
res_agd <- agd(f, x1, iters, 0.05)
f_min <- min(res_agd)
ggplot(data = data.frame(`Log_suboptimality_gap` = -c(log(res_gd),
                                                     log(res_agd)),
                         iter = log(c(1:5000, 1:5000)),
                         algo = c(rep('gd', 5000), rep('agd', 5000)))) + 
  geom_line(aes_string(y = "-log_suboptimality_gap", x = "log-time", col = "algo"))
