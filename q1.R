library(rootSolve)
library(ggplot2)

# Get matrix A
getA <- function(k) {
  A <- matrix(nrow = 10, ncol = 10)
  for (i in 1:10) {
    if (i < 10){
      for (j in (i+1):10) A[i, j] <- exp(i/j)*cos(i*j)*sin(k)
    }
    if (i > 1) {
      for (j in 1:(i-1)) A[i, j] <- A[j, i] 
    }
    A[i, i] <- i*abs(sin(k))/10 + sum(abs(A[i, ]), na.rm = T)
  }
  return(A)
}

# Get vector b
getb <- function(k) {
  b <- sapply(1:10, function(i) exp(i/k)*sin(i*k))
  return(b)
}

A <- lapply(1:5, function(k) getA(k))
b <- lapply(1:5, function(k) getb(k))

# Objective function
f <- function(x, A, b) {
  fk <- sapply(1:5, function(k) x %*% (A[[k]] %*% x) - b[[k]] %*% x)
  return(max(fk))
}

# Init
x1 <- rep(1, 10)
f(x1, A, b)

# subgrad descent
subgrad <- function(f, x1, iters, C, A, b) {
  xn <- x1
  obj <- rep(0, iters)
  obj[1] <- f(x1, A, b)
  for (t in 2:iters) {
    eta  <- C/sqrt(t)
    g <- gradient(f, xn, A = A, b = b)
    xn <- xn - c(eta * g / sqrt(sum(g^2)))
    obj[t] <- min(f(xn, A, b), obj[t-1])
  }
  return(obj)
}

iters <- 1000
res_sg <- subgrad(f, x1, iters, 1, A = A, b = b)

f_min <- -0.83

# subgrad Polyak stepsize
subgrad_polyak <- function(f, x1, iters, f_min, A, b) {
  xn <- x1
  obj <- rep(0, iters)
  obj[1] <- f(x1, A, b)
  for (t in 2:iters) {
    g <- gradient(f, xn, A = A, b = b)
    norm <- sqrt(sum(g^2))
    eta <- (f(xn, A, b) - f_min) / norm
    xn <- xn - c(eta * g / norm)
    obj[t] <- min(f(xn, A, b), obj[t-1])
  }
  return(obj)
}
res_polyak <- subgrad_polyak(f, x1, iters, f_min, A = A, b = b)

ggplot(data = data.frame(`Log_suboptimality_gap` = c(log(res_sg - f_min),
                                                     log(res_polyak - f_min)),
                         iter = c(1:1000, 1:1000),
                         algo = c(rep('subgrad', 1000), rep('polyak', 1000)))) + 
  geom_line(aes_string(y = "Log_suboptimality_gap", x = "iter", col = "algo"))
