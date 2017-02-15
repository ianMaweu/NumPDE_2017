library(plot3Drgl)

backwardDifference <- function(u0=function(x) 2*sin(2*pi*x), 
                               K=1, L=1, N=30, T=0.1, M=30) {
  # set up space grid
  h <- L/N
  x <- h*(1:(N-1))
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  # set up vectors with initial condition
  w <- u0(x)
  
  gamma <- K*tau/(h^2)
  
  Temperature <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  Temperature[ , 1] <- w  # Initial value
  # Loop over time steps
  for (j in 1:M) {
    w <- doublesweep(rep(gamma, N-1), rep(gamma, N-1), 
                     rep(1 + 2* gamma, N-1), -w, 0, 0)
    Temperature[ , j+1] <- w
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=x, t=t, w=Temperature))
}

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

maxError <- function(N, M, method) {
  # numerical solution
  numSol <- method(M=M, N=N)
  # exact solution
  x <- numSol$x
  t <- numSol$t
  xy <- mesh(x, t)
  u <- with(xy, 2*sin(2*pi*x)*exp(-4*pi^2*y))
  
  return(max(abs(u - numSol$w)))
}

N <- 15*2^(0:7)
errbd <- sapply(N, function(N) maxError(N, N, backwardDifference))
#sapply applies the function maxError(N,N) for each N[i] in N then outputs the return values as a vector
errCN <- sapply(N, function(N) maxError(N, N, CrankNicolson))
plot(errbd ~ N, type="b", log="xy", ylim=c(0.0000001, 0.1), ylab="Error")
lines(errCN ~ N, type="b", col="blue")


#--------------------------------------------------------------#
#Exercise 6.1
CrankNicolson <- function(u0=function(x) 2*sin(2*pi*x), 
                               K=1, L=1, N=30, T=0.1, M=30) {
  # set up space grid
  h <- L/N
  x <- h*(1:(N-1))
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  # set up vectors with initial condition
  w <- u0(x)
  
  gamma <- K*tau/(h^2)
  
  Temperature <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  Temperature[ , 1] <- w  # Initial value
  
  #Construct B
  B <- diag(1-gamma, N-1,N-1)
  
  for (i in 1:(N-2))
  {
    B[i, i+1] <- gamma/2
    B[i+1, i] <- gamma/2
  }
  
  # Loop over time steps
  for (j in 1:M) {
    w <- doublesweep(rep(gamma/2, N-1), rep(gamma/2, N-1), 
                     rep(1 + gamma, N-1), B%*%-w, 0, 0)
    Temperature[ , j+1] <- w
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=x, t=t, w=Temperature))
}

sol <- backwardDifference()
persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks
plotrgl()

#Exercise 6.2
N <- 15*2^(0:7)
errbd <- sapply(N, function(N) maxError(N,60, backwardDifference))
#sapply applies the function maxError(N,N) for each N[i] in N then outputs the return values as a vector
errCN <- sapply(N, function(N) maxError(N, 60, CrankNicolson))
plot(errbd ~ N, type="b", log="xy", ylim=c(0.0000001, 0.1), ylab="Error")
lines(errCN ~ N, type="b", col="blue")

#Exercise 6.3
N <- 15*2^(0:7)
errbd <- sapply(N, function(N) maxError(60,N, backwardDifference))
#sapply applies the function maxError(N,N) for each N[i] in N then outputs the return values as a vector
errCN <- sapply(N, function(N) maxError(60, N, CrankNicolson))
plot(errbd ~ N, type="b", log="xy", ylim=c(0.0000001, 0.1), ylab="Error")
lines(errCN ~ N, type="b", col="blue")

