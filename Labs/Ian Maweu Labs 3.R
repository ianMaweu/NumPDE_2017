N <- 150            #Number of steps in interval
L <- 1             #Length of interval
h <- L/N           #Step size in space
x <- h*(1:(N-1))   #grid points for x
w <- 2*sin(2*pi*x) #Initial condition (u(x,0) = w^(0))

#w <- 2*cos(5*pi*x)
#w <- (x-0.5)^3+0.5 #Change T to 1 and M to 2000
#w <- 2*cos(9*pi*x)*exp(-x)


T <- 0.1       #Last Time value
M <- 1000       #Number of time steps
tau <- T/M     #size of time steps
t <- tau*(0:M) 
K <- 1         #Thermal diffusivity
gamma <- K*tau/(h^2)


# Constructs the matrix A which evolves the statw W^(j) to w^(j+1) (the next ime step)
A <- diag(1-2*gamma, N-1)
for(k in 1:N-2)
{
  A[k,k+1] <- gamma
  A[k+1,k] <- gamma
}

Temperature <- matrix(0,N-1,M+1) #Matrix to hold the soln
Temperature[ , 1] <- w

# loop over time steps to generate soln
for(j in 0:(M))
{
  w <- A%*%w
  Temperature[ , j+1] <- w
}

library("plot3Drgl")

#plot soln
persp3D(x,t,Temperature, xlab = "x", ylab = "t", zlab = "Temperature", ticktype ="detailed", nticks =4)
plotrgl()
#plotrgl(smooth = TRUE, lighting = TRUE)


# Check stability
if (gamma < 0.0 || gamma >0.5)
{
  print("The forward difference method is unstable ")
}else
{
  print("The finite difference method is stable ")
}