library("plot3Drgl")

backwardDifference1 <- function(u0=function(x) 2*sin(2*pi*x), 
                               K=1, L=1, N=30, T=0.1, M=30) {
  # set up space grid
  h <- L/N
  x <- h*(1:(N-1))
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  # set up vectors with initial condition
  w <- u0(x)
  
  # Set up evolution matrix
  gamma <- K*tau/(h^2)
  A <- diag(1+2*gamma, N-1)
  for (k in 1:(N-2)) {
    A[k,k+1] <- -gamma
    A[k+1,k] <- -gamma
  }
  
  Ainv <- solve(A)
  
  # Compute soln
  Temperature <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  Temperature[ , 1] <- w  # Initial value
  # Loop over time steps
  for (j in 1:M) {
    w <- Ainv %*% w 
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




#Exercise 5.1
#N = L/h
numSol <- backwardDifference1(N=500)
persp3D(numSol$x, numSol$t, numSol$w,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks
plotrgl()


#Exercise 5.2
A <- c(0,2,3,1)
B <- c(2,1,4,0)
C <- c(-1,-1,-1,-7)
F <- c(1,2,6,-6)

x <- doublesweep(A,B,-C,F,0,0)
D <- t(matrix(c(-1,2,0,0, 2,-1,1,0, 0,3,-1,4, 0,0,1,-7), nrow = 4, ncol = 4))
D%*%x

#-------------------------------------------------------------------------------#
#Exercise 5.3
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
  
  # Set up evolution matrix
  Ci <- rep(1+2*gamma,N-1)
  Bi <- rep(gamma,N-1)
  Ai <- rep(gamma,N-1)
  a <- 0
  b <- 0
  
  Temperature <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  Temperature[ , 1] <- w  # Initial value
  
  for (j in 1:M)
  {
    w <- doublesweep(Ai,Bi,Ci,-w,a,b)
    Temperature[ ,j+1] <- w
  }
  
  # Return a list consisting of time intervals, length intervals and solution to heat eqn at the intervals.
  return(list(x=x, t=t, w=Temperature))
  
}

sol <- backwardDifference()
persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks
plotrgl()

#-------------------------------------------------------------------#

maxError <- function(N, M) {
  # numerical solution
  numSol <- backwardDifference(M=M, N=N)
  # exact solution
  x <- numSol$x
  t <- numSol$t
  xy <- mesh(x, t)
  u <- with(xy, 2*sin(2*pi*x)*exp(-4*pi^2*y))
  
  return(max(abs(u - numSol$w)))
}

#Exercise 5.4
varyStepSize <- function(n=7, m =10)
{
  N <- 20*(1:n)+10
  M <- 20*(1:m)+10
  error <- matrix(0,nrow = m, ncol = n)
  
  for (j in 1:m)
  {
    for (i in 1:n)
    {
      error[j,i] <- maxError(N[i],M[j])
    }
  }
  
  return(list(N=N, M=M, E=error))
}


err <- varyStepSize()
persp3D(err$M, err$N, err$E, xlab = "M", ylab = "N", zlab = "Error",ticktype="detailed")
plotrgl()
# The error appears to decrease faster with an decrease in length step size (quadratic) compared to 
# a decrease in time step size (linearly)