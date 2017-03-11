library("plot3Drgl")
doublesweep <- function(A, B, C, F, a, b) {
  # Solves the equation 
  # A[i]*v[i-1] - C[i]*v[i] + B[i]*v[i+1] = F[i]
  # for v[i], i = 1,...,N-1, with boundary conditions
  # v[0]=a and v[N]=b
  
  # Check the lengths of the vectors
  N <- length(C) + 1
  if ((length(B) != N-1) || (length(A) != N-1) || (length(F) != N-1)) {
    stop("The lengths of the vector arguments need to be equal")
  }
  
  alpha <- rep(0, N)
  beta <- rep(0, N)
  beta[1] <- a
  
  #sweep up to get alphas and betas
  for (i in 1:(N-1)) {
    alpha[i+1] <- B[i] / (C[i]-alpha[i]*A[i])
    beta[i+1] <- (beta[i]*A[i] - F[i]) / (C[i] - alpha[i]*A[i])
  }
  
  v <- rep(0, N-1 )
  v[N-1] <- alpha[N]*b + beta[N]
  
  #sweep down to find v's
  for (i in (N-1):2) {
    v[i-1] <- alpha[i]*v[i] + beta[i]    
  }
  
  return(v)
}

#---------------------Backward Difference 1------------------------#
# First version of backward difference from lab 8.
backwardDifference1 <- function(u0, K, f, L=1, N=30, T=1, M=30) {
  # set up space grid
  h <- L/N
  x <- h*(1:(N-1))
  xFull <- h*(0:N)
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  gamma <- tau/h^2
  
  w <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  w[ , 1] <- u0(x)  # Initial value
  # Loop over time steps
  for (j in 1:M) {
    Kv <- K(xFull, t[j], c(0, w[, j], 0)) 
    chi <- (Kv[1:N] + Kv[2:(N+1)])/2
    a <- 1 + gamma*(chi[2:N]+chi[1:(N-1)])
    b <- -gamma*chi[1:(N-1)]
    w[, j+1] <- doublesweep(b[1:(N-1)], c(b[2:(N-1)],0), -a, 
                            w[, j]+tau*f(x, t[j], w[, j]), 0, 0)
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=x, t=t, w=w))
}

#------------------------Exercise 1-----------------------------#
sol <- backwardDifference1(u0=function(x) pmax(0, 1-10*abs(x-0.5)), 
                           K=function(x, t, u) x*t*u^2,
                           f=function(x, t, u) 8*x*sin(8*pi*x*t))

persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
plotrgl(lighting = TRUE)

#-----------------------Backward difference 2-------------------------#
backwardDifference2 <- function(u0=function(x) pmax(0,1-10*abs(x-0.5)), 
                                K=function(x, t, u) u^2/2,
                                f=function(x, t, u) -u^2,
                                L=1, N=30, T=1, M=30,
                                max_iteration=20, tolerance=0.000001) {
  # set up space grid
  h <- L/N
  x <- h*(1:(N-1))
  xFull <- h*(0:N)
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  gamma <- tau/h^2
  
  w <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  iterations <- rep(0, M)  # to hold number of iterations needed at each step
  precision <- rep(0, M)  # to hold precision achieved
  
  w[ , 1] <- u0(x)  # Initial value
  # Loop over time steps
  for (j in 1:M) {
    wn <- w[, j]  # initial guess for w_{j+1}
    for (s in 1:max_iteration) {
      ws <- wn
      Kv <- K(xFull, t[j], c(0, ws, 0))
      chi <- (Kv[1:N] + Kv[2:(N+1)])/2
      a <- 1 + gamma*(chi[2:N]+chi[1:(N-1)])
      b <- -gamma*chi[1:(N-1)]
      # update guess for w_{j+1}
      wn <- doublesweep(b[1:(N-1)], c(b[2:(N-1)],0), -a, 
                        w[, j]+tau*f(x, t[j], ws), 0, 0)
      # break if tolerance limit is reached
      if (max(abs(wn-ws)) < tolerance) {
        break
      }
    }
    iterations[j] <- s
    precision[j] <- max(abs(wn-ws))
    w[, j+1] <- wn
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=x, t=t, w=w, iterations=iterations, precision=precision))
}

#---------------------------Exercise 2--------------------------------#
sol3 <- backwardDifference2(u0=function(x) pmax(0, 1-10*abs(x-0.5)), 
                            K=function(x, t, u) u^2/2,
                            f=function(x, t, u) 8*x*sin(8*pi*x*t),
                            max_iteration = 20,tolerance=0.000001)

persp3D(sol3$x, sol3$t, sol3$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
plotrgl()

mean(sol3$precision) # Calculates the mean precision attained which is 5.92803e-07
sum(sol3$iterations) # Calculates the total number of iterations performed which is 352
mean(sol3$iterations) # Calculates the avaerage number of iterations performed per time step which is 11

sol4 <- backwardDifference2(u0=function(x) pmax(0, 1-10*abs(x-0.5)), 
                            K=function(x, t, u) u^2/2,
                            f=function(x, t, u) 8*x*sin(8*pi*x*t),
                            max_iteration = 4,tolerance=0.000001)
mean(sol4$precision) # A mean precision of 0.00524723 is achieved when max iteration is restriced to 4
