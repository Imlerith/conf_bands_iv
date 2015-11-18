library("foreach")
library("MASS")
library("quantreg")
library("KernSmooth")
library("doParallel")

### Functions

kernelq = function(u){
  dnorm(u, mean = 0, sd = 1)
}

# M-type smoother
lnrob = function(x, y, h, x0 = seq(0,1, length.out = 100)){
  
  xx = sort(x0)
  xx = (xx - min(xx))/(max(xx) - min(xx)) # EB
  fv = xx
  dv = xx
  x  = (x - min(x))/(max(x)-min(x))
  
  for (i in 1:length(xx)) {
    z     = x - xx[i]
    wx    = dnorm(z/h)
    r     =  rlm(y ~ z, weights = wx, method = "M")
    u     =  r$wresid
    fv[i] =  r$coefficients[[1]]
    dv[i] =  r$coefficients[[2]]
  }
  
  "psi1" = r$psi
  
  return( list(xx = xx, fv = fv, dv = dv, "psi1" = psi1) )
}

# Quantile regression with specific tau
lprq2 = function(x, y, h, tau, x0)
{
  # xx = x0
  xx = sort(x0) # EB
  xx = (xx - min(xx))/(max(xx) - min(xx)) # EB
  fv = xx
  dv = xx
  x  = (x - min(x)) / (max(x) - min(x)) # EB
  for(i in 1:length(xx)){
    z     = x - xx[i]
    wx    = dnorm(z/h)
    r     = rq(y ~ z, tau = tau, weights = wx, method = "br")
    fv[i] = r$coef[1.]
    dv[i] = r$coef[2.]
  }
  list(xx = xx, fv = fv, dv = dv)
}

# Quantile regression with random tau
# EB: Should'nt it be the same tau over whole support?
lprq3 = function(x, y, h, x0){
  # xx  = x0
  xx  = sort(x0) # EB
  xx  = (xx - min(xx)) / (max(xx) - min(xx))
  fv  = xx
  dv  = xx
  x   = (x - min(x)) / (max(x) - min(x)) # EB
  tau = runif(1) # EB
  for(i in 1:length(xx)) {
    z     = x - xx[i]
    wx    = dnorm(z/h)
    r     = rq(y ~ z, weights = wx, tau = runif(1), ci = FALSE) # EB
    # r     = rq(y ~ z, weights = wx, tau = tau, ci = FALSE) # EB
    fv[i] = r$coef[1.]
    dv[i] = r$coef[2.]
  }
  list(xx = xx, fv = fv, dv = dv)
}


## Setup
cl    = 7     # Number of cores to use
n     = 500
B     = 1000
alpha = 0.05
gridn = 100

h     = 0.06
hg    = 0.15
g     = 0.20
# Generate variables

x            = runif(n, 0, 1)
prob         = 0.9
aa           = rbern(n, prob)
xxa          = rnorm(length(aa[aa==1]), mean = 0, sd = 1)
xxb          = rnorm(length(aa[aa==0]), mean = 0, sd = 10)
aa[aa==1]    = xxa
aa[aa==0]    = xxb
y            = 8 * sin(x * pi) + aa

rm(prob, aa, xxa, xxb)

# Sort
y   = y[order(x)]
x   = x[order(x)]

# True curve
ytrue = 8 * sin(x * pi)

plot(x, y)
lines(x, ytrue, col = "red")

# Scale x to [0, 1]
xmin = min(x)
xmax = max(x)
x = (x - xmin) / (xmax - xmin)


# Initial fit
yhat.h = lnrob(x, y, h = h, x0 = x)
yhat.g = lnrob(x, y, h = g, x0 = x)

yhat.grid.h = lnrob(x, y, h = h, x0 = seq(0, 1, length.out = gridn))
yhat.grid.g = lnrob(x, y, h = g, x0 = seq(0, 1, length.out = gridn))


# Empirical pdf of x at gridpoints
fxd = bkde(x, gridsize = gridn, range.x = c(yhat.grid.h$xx[1], yhat.grid.h$xx[gridn]))

# Conditional pdf f(e|x) and E(psi^2(e)) at gridpoints
fl  = vector(length = gridn, mode = "numeric")
fll = vector(length = gridn, mode = "numeric")

for (k in 1: gridn){
  # Conditional pdf f(e|x)at gridpoints
  nom   = sum((kernelq((x - yhat.grid.h$xx[k]) / (hg)) *
                 yhat.grid.h$psi1((y - yhat.grid.h$fv[k]), deriv = 1)))
  denom = sum(kernelq((x - yhat.grid.h$xx[k])/(hg)))
  fl[k] = nom / denom
  
  # Conditional E(psi^2(e))
  nom    = sum((kernelq((x - yhat.grid.h$xx[k])/(h)) *
                  yhat.grid.h$psi1((y - yhat.grid.h$fv[k]), deriv = 0)^2))
  denom  = sum(kernelq((x -  yhat.grid.h$xx[k])/(hg)))
  fll[k] = nom / denom
}

bandt = (fxd$y)^(1/2) * abs(fl / sqrt(fll))


# Bootstrapping

pack = c("MASS", "KernSmooth", "Rlab", "quantreg")
cl   = makeCluster(cl)
registerDoParallel(cl)

d    = vector(length = B, mode = "numeric")

d = foreach(i = 1:B, .packages = pack)%dopar%{
  ehat    = lprq3( yhat.h$xx, (y - yhat.h$fv), h = hg, x0 = yhat.grid.h$xx )
  ystar   = yhat.grid.g$fv + ehat$fv
  fitstar = lnrob(yhat.grid.h$xx, ystar, h = h, x0 = yhat.grid.h$xx )
  d.m     = max(abs(bandt*abs(fitstar$fv - yhat.grid.g$fv)))
  d.m
}

stopCluster(cl)

d = unlist(d)


dstar = quantile(d[d!=0], probs = 1 - alpha)
dstar = dstar * {bandt}^(-1)

# Asymptotic Bands
cc     = 1 / 4
lambda = 1 / 2 / sqrt(pi)
delta  = - log(h) / log(n)
dd     = sqrt(2 * delta * log(n)) + (2 * delta * log(n))^(-1/2) * log(cc / 2 / pi)
cn     = log(2) - log(abs(log(1 - alpha)))
band   = (n * h)^(- 1/2) * bandt^{-1} * (dd + cn * (2 * delta * log(n))^(-1/2)) * sqrt(lambda)

# Scale back
x      = ( x              * (xmax - xmin) ) + xmin
x.grid = ( yhat.grid.h$xx * (xmax - xmin) ) + xmin


plot(x, y, xlab = "", ylab = "")
lines(x, ytrue, col = "green", lwd = 2)
lines(x.grid, yhat.grid.h$fv, lwd = "3", col = "blue")


lines(x.grid, (yhat.grid.h$fv - dstar), col = "red",   lty = 3, lwd = 3)
lines(x.grid, (yhat.grid.h$fv + dstar), col = "red",   lty = 3, lwd = 3)
lines(x.grid, (yhat.grid.h$fv - band),  col = "yellow", lty = 3, lwd = 3)
lines(x.grid, (yhat.grid.h$fv + band),  col = "yellow", lty = 3, lwd = 3)