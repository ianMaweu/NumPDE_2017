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

#------------------ADI Algorithm Exercise 1-----------------------#
ADI <- function(u0, K=1, f, L1=1, N1=30, L2=1, N2=30, T=1, M=30) {
  # set up space grids
  h1 <- L1/N1
  x <- h1*(0:N1)
  h2 <- L2/N2
  y <- h2*(0:N2)
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  gamma1 <- K*tau/h1^2
  gamma2 <- K*tau/h2^2
  # Vectors to be later used in double sweep method
  A1 <- rep(gamma1/2, N1-1)
  C1 <- rep(1+gamma1, N1-1)
  A2 <- rep(gamma2/2, N2-1)
  C2 <- rep(1+gamma2, N2-1)
  
  w <- array(0, dim=c(N1+1, N2+1, M+1))  # Array to hold the solution
  w[, , 1] <- outer(x, y, u0)  # Initial value
  
  # Loop over time steps
  for (n in 1:M) {
    
    # incorparate non-homogenoues boundary conditions on the (0,y,t) and (L1,y,t) boundaries
    w[1,,n] <- t[n]*sin(2*pi*y) 
    w[N1+1,,n] <- t[n]*y*(1-y)
    
    # Matrix with contributions from inhomogeneous term $\tau/2 f^{n+1/2}$
    # f^{n+(1/2)} = (1/2)*[f^{n} + f^{n+1}]
    Fh = tau*(outer(x, y, f, t=t[n]) + outer(x, y, f, t=t[n+1]))/4
    # first half step
    wh <- matrix(0, nrow=N1+1, ncol=N2+1)  # matrix to hold w^{n+1/2}
    
    for (j in 2:N2) {
      F1 <- gamma2/2*(w[2:N1, j+1, n] + w[2:N1, j-1, n]) + 
        (1-gamma2)*w[2:N1, j, n] + Fh[2:N1, j]
      wh[2:N1, j] <- doublesweep(A1, A1, C1, -F1, w[1,j,n], w[N1+1,j,n]) #Added '-' sign to F1
      #Also changed the argumnets a and b of the doublesweep() function from a=0 & b=0 to a=w[1,j,n] & b=w[N1+1,j,n]) 
    }
    
    # second half step
    for (k in 2:N1) {
      F2 <- gamma1/2*(wh[k+1, 2:N2] + wh[k-1, 2:N2]) + 
        (1-gamma1)*wh[k, 2:N2] + Fh[k, 2:N2]
      w[k, 2:N2, n+1] <- doublesweep(A2, A2, C2, -F2, 0, 0) #Added '-' sign to F2
    }
  }
  
  # Return a list consisting of grid and solution
  return(list(x=x, y=y, t=t, w=w))
}

#---------------Solution with ADI------------------#
sol <- ADI(u0=function(x, y) sin(pi*x)*sin(pi*y), K=1, 
           f=function(x, y, t) 20*pi^2*sin(2*pi*x)*sin(2*pi*y),
           L1=1, L2=1, N1=100, N2=100, T=0.2, M=100)

persp3D(sol$x, sol$y, sol$w[, , 1],
        xlab="x", ylab="y", zlab="w",
        ticktype="detailed", nticks=4, phi=10, theta=90)

persp3D(sol$x, sol$y, sol$w[, , 21],
        xlab="x", ylab="y", zlab="w",
        ticktype="detailed", nticks=4, phi=10, theta=90)

for (n in 1:21) {
  persp3D(sol$x, sol$y, sol$w[, , n],
          xlab="x", ylab="y", zlab="w", zlim=c(-0.7, 1), clim=c(-0.7, 1),
          ticktype="detailed", nticks=4, phi=12, theta=90)
}
