OPP <- function(data, Tau, method = c("t1", "t2")) {
  fun <- switch (method[1],
                 t1 = opp_t1,
                 t2 = opp_t2,
                 stop("Unrecognized method.")
  )
}


opp_t1 <- function(data, Tau) {
  C <- function(tau) {
    cov(g[(1 + t):nrow(g), ], g[1:(nrow(g) - t), ])
  }

  C0 <- C(0)

  M <- trap_integral(C, 0:Tau)


  X <- solve(C0) %*% (M %*% t(M))
  eig <- eigen(X)

  data_prime <- as.data.frame(data %*% eig$vectors)
  colnames(data_prime) <- paste0("V", seq_len(ncol(data_prime)))

  r <- C0 %*% eig$vectors

  return(list(T1 = eig$values,
              data = data_prime,
              OPP = r))
}


opp_t2 <- function(g, Tau) {
  force(g)
  force(Tau)
  C <- function(tau) {
    cov(g[(1 + t):nrow(g), ], g[1:(nrow(g) - t), ])
  }

  # Define T2(X) and its gradient wrt X
  C0 <- C(0)
  C0_eigen <- eigen(C0)
  C0_12 <- with(C0_eigen, vectors %*% diag(sqrt(values), nrow = length(values)) %*% t(vectors) )
  C0_12_inv <- solve(C0_12)

  T2 <- function(X) {
    denom <- sum(X^2)
    left_num <-  t(X) %*% C0_12_inv
    right_num <- C0_12_inv %*% X

    numer.fun <- function(tau) {
      ((left_num %*% C(tau) %*% right_num)^2)[1, 1]
    }
    integral <- trap_integral(numer.fun, 0:Tau)
    t <- 2 * integral / denom^2
    return(t)
  }


  DT2 <- function(X) {
    denom <- sum(X^2)

    numer.fun <- function(tau) {
      C_t <-  C(tau)
      (t(X)  %*% (C0_12_inv %*% C_t %*% C0_12_inv) %*% X)[1, 1] *
        ((C0_12_inv %*% (C_t + t(C_t)) %*% C0_12_inv) %*% X)
    }

    integral <- trap_integral(numer.fun, 0:Tau)

    4*(integral/denom^2  - (T2(X) * X)/denom)
  }


  # Optimisation.
  prev_sol <- list(value = NA, par =  runif(3)) # Random start
  for (i in 1:10) {
    new_sol <- optim(prev_sol$par, T2, gr = DT2, control = list(fnscale = -1),
                     method = "CG")
    meta_convergence <- isTRUE(all.equal(new_sol$par, prev_sol$par))

    if (meta_convergence) {
      break
    }
    prev_sol <-  new_sol
    prev_sol$par <- prev_sol$par / sqrt(sum(prev_sol$par^2))
  }

  e <- C0_12_inv %*% new_sol$par

  r <- C0 %*% e
}



trap_integral <- function(fun, values) {
  integral <- lapply(values, fun)
  integral[[1]] <- integral[[1]]/2
  integral[[length(integral)]] <- integral[[length(integral)]]/2
  Reduce("+", integral)
}

microbenchmark::microbenchmark(as.numeric(t(X) %*% X), sum(X^2))


T2_vector <- function(X) {
  print(X)
  denom <- sum(X^2)
  left_num <-  t(X) %*% C0_12_inv
  right_num <- C0_12_inv %*% X

  numer.fun <- function(tau) {
    left <- left_num %*% C(tau)
    dim(left) <- c(3, length(tau))
    left <- t(left)

    ((left %*% right_num)^2)[, 1]
  }
  numer <- trap_integral(numer.fun, 1:Tau)
  t <- 2 * numer / denom
  print(t)
  return(t)
}
