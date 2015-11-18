rm(list = ls(all=TRUE))

#Load libraries

library("foreach")
library("MASS")
library("quantreg")
library("KernSmooth")
library("doParallel")
source("cbands_funcs.r")


## Setup
cl    = 7     # Number of cores to use
B     = 1000
alpha = 0.05
gridn = 100


# Load data
# 1:1032
monivdata = read.table('monivdata.csv',sep=',')
x         = as.matrix(monivdata[,1])
y         = as.matrix(monivdata[,2])

#Choose bandwidths
n  = nrow(x)
hg    = 0.15
g     = n^(-1/9)

# x            = runif(n, 0, 1)
# prob         = 0.9
# aa           = rbern(n, prob)
# xxa          = rnorm(length(aa[aa==1]), mean = 0, sd = 1)
# xxb          = rnorm(length(aa[aa==0]), mean = 0, sd = 10)
# aa[aa==1]    = xxa
# aa[aa==0]    = xxb
# y            = 8 * sin(x * pi) + aa

# rm(prob, aa, xxa, xxb)

# Sort
y   = y[order(x)]
x   = x[order(x)]

# True curve
# ytrue = 8 * sin(x * pi)

plot(x, y)
# lines(x, ytrue, col = "red")

# Scale x to [0, 1]
xmin = min(x)
xmax = max(x)
x = (x - xmin) / (xmax - xmin)
h = median(abs(x-median(x)))/0.6745*(4/3/n)^0.2


# Initial fit
#one can choose ANY of these two: gridded or non-gridded
yhat.h = lnrob(x, y, h = h, maxiter = 100, x0 = x)
yhat.g = lnrob(x, y, h = g, maxiter = 100, x0 = x)

yhat.grid.h = lnrob(x, y, h = h, maxiter = 100, x0 = seq(0, 1, length.out = gridn))
yhat.grid.g = lnrob(x, y, h = g, maxiter = 100, x0 = seq(0, 1, length.out = gridn))


# Empirical pdf of x at gridpoints
fxd = bkde(x, gridsize = gridn, range.x = c(yhat.grid.h$xx[1], yhat.grid.h$xx[gridn]))

# Bootstrapping

pack = c("MASS", "KernSmooth", "Rlab", "quantreg")
cl   = makeCluster(cl)
registerDoParallel(cl)

d    = vector(length = B, mode = "numeric")

d = foreach(i = 1:B, .packages = pack)%dopar%{
  estar = lprq3( yhat.h$xx, (y - yhat.h$fv), h = hg, x0 = yhat.grid.h$xx )
  eh    = median(abs(estar$fv-median(estar$fv)))/0.6745*(4/3/n)^0.2
  
  #THE WHOLE BANDT THING IN THE BOOTSTRAP PROCEDURE NOW!
  # Conditional pdf f(e|x) and E(psi^2(e)) at gridpoints
  fl  = vector(length = gridn, mode = "numeric")
  fll = vector(length = gridn, mode = "numeric")
  
  for (k in 1: gridn){
    # Conditional pdf f(e|x)at gridpoints
    #HERE THE BANDWIDTH CHANGED TO BANDWIDTH FOR EHAT!
    #also ehat inserted instead of y - yhat.grid.h$fv; makes sense, because 
    #in the paper the density is estimated at ehat
    nom   = sum((kernelq((x - yhat.grid.h$xx[k]) / (eh)) *
                   yhat.grid.h$psi1((estar$fv[k]), deriv = 1)))
    denom = sum(kernelq((x - yhat.grid.h$xx[k])/(eh)))
    fl[k] = nom / denom
    
    # Conditional E(psi^2(e))
    nom    = sum((kernelq((x - yhat.grid.h$xx[k])/(eh)) *
                    yhat.grid.h$psi1((estar$fv[k]), deriv = 0)^2))
    denom  = sum(kernelq((x -  yhat.grid.h$xx[k])/(eh)))
    fll[k] = nom / denom
  }
  bandt = (fxd$y)^(1/2) * abs(fl / sqrt(fll))
  
  ystar   = yhat.grid.g$fv + estar$fv
  fitstar = lnrob(yhat.grid.h$xx, ystar, h = h, maxiter = 50, x0 = yhat.grid.h$xx )
  d.m     = max(abs(bandt*abs(fitstar$fv - yhat.grid.g$fv)))
}

stopCluster(cl)

d = unlist(d)

# Conditional pdf f(e|x) and E(psi^2(e)) at gridpoints
fl  = vector(length = gridn, mode = "numeric")
fll = vector(length = gridn, mode = "numeric")

for (k in 1: gridn){
  # Conditional pdf f(e|x)at gridpoints
  #HERE THE BANDWIDTH CHANGED TO H!
  nom   = sum((kernelq((x - yhat.grid.h$xx[k]) / (h)) *
                 yhat.grid.h$psi1((y - yhat.grid.h$fv[k]), deriv = 1)))
  denom = sum(kernelq((x - yhat.grid.h$xx[k])/(h)))
  fl[k] = nom / denom
  
  # Conditional E(psi^2(e))
  nom    = sum((kernelq((x - yhat.grid.h$xx[k])/(h)) *
                  yhat.grid.h$psi1((y - yhat.grid.h$fv[k]), deriv = 0)^2))
  denom  = sum(kernelq((x -  yhat.grid.h$xx[k])/(hg)))
  fll[k] = nom / denom
}

bandt = (fxd$y)^(1/2) * abs(fl / sqrt(fll))


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
# lines(x, ytrue, col = "green", lwd = 2)
lines(x.grid, yhat.grid.h$fv, lwd = "3", col = "blue")


lines(x.grid, (yhat.grid.h$fv - dstar), col = "red",   lty = 3, lwd = 3)
lines(x.grid, (yhat.grid.h$fv + dstar), col = "red",   lty = 3, lwd = 3)
lines(x.grid, (yhat.grid.h$fv - band),  col = "yellow", lty = 3, lwd = 3)
lines(x.grid, (yhat.grid.h$fv + band),  col = "yellow", lty = 3, lwd = 3)
