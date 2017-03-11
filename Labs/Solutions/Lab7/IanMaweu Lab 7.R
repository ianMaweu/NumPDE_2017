library(plot3Drgl)

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

#-----------------Backward Difference--------------------#

backwardDifference <- function(f=function(x,t) 0, 
                               u0=function(x) 2*sin(2*pi*x),
                               a= function(t) 0, b=function(t) 0,
                               K=1, L=1, N=30, T=0.1, M=30) {
  # set up space grid
  h <- L/N
  x <- h*(0:(N))
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  # set up vectors with initial condition
  w <- u0(x)
  
  gamma <- K*tau/(h^2)
  
  Temperature <- matrix(0, N+1, M+1)  # Matrix to hold the solution
  Temperature[ , 1] <- w  # Initial value
  # Loop over time steps
  for (j in 1:M) {
    F <- f(x,t[j])
    A <- a(t[j])
    B <- b(t[j])
    w <- doublesweep(rep(gamma, N+1), rep(gamma, N+1), 
                     rep(1 + 2* gamma, N+1), -w-tau*F, A, B)
    Temperature[ , j+1] <- w
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=x, t=t, w=Temperature))
}

#-----------------Solution plot--------------------#
sol <- backwardDifference(u0=function(x) x, 
                         f=function(x,t) 10*sin(10*pi*t)*exp(-10*(x-0.5)^2),
                         a=function(t) sin(20*pi*t), b=function(t) cos(20*pi*t),
                         K=1, L=1, T=0.4, N=30, M=720)
persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
plotrgl()

#Changes made to backward difference formula:
# 1) Increased x by 2 entries, done by setting the sequance from 0 to N instead of 1 to N-1
# 2) Increased the sizes of vectors gamma and (1+2*gamma) by 2 by changing thier value from N-1 to N+1
# 3) Added vectors F,A and B which calculate the non-homogeneous term f(x,t) and the functions
#    a(t) and b(t) respectively. They reside within the for loop since they are now time
#    dependent functions
# 4) Increased the number of rows of Temperature by 2 by changing N-1 to N+1.

# Note that the backward difference method does not require that you manually edit
#  the boundary conditions since the doublesweep() method already takes them into account.