
params <- list(sigma = 10, b = 8/3, r = 28)
y <- c(x = 1, y = 1, z = 1)


Lorenz <- function(t, y, parms,...) {
  x_dot <- params$sigma*(y[2] - y[1])
  y_dot <- -y[3]*y[1] + params$r*y[1] - y[2]
  z_dot <- y[1]*y[2] - params$b*y[3]

  list(c(x_dot, y_dot, z_dot))
}


sol <- deSolve::rk4(y, seq(0, 8300, by = 1/24), Lorenz, params)

g <- as.matrix(sol)[, -1]
g <- scale(g, scale = FALSE)

plot((g %*% r)[1:500], type = "l")


a <- acf((g %*% r)[, 1], lag.max = 15*24)

plot(a$lag/24, a$acf, type = "l")


eof <- svd(scale(g))


a <- acf(eof$u[, 1], lag.max = 15*24)

plot(a$lag/24, a$acf, type = "l")
