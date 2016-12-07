# SOS surface example to include with paper
# (c) Stephen R. Carpenter 2016-11-30

rm(list = ls())
graphics.off()

library(fields)

# Functions ********************************************************************************

dx.roots = function(I,r,K,p,h,H) { # Function to find equilibria of x
  parvec = c( I*h*h, 
              (h*h*r - h*h*H),
              (I - p - h*h*r/K),
              (r-H),
              -(r/K) )
  xeq = polyroot(parvec)    
  return(xeq)
}

N.real.roots = function(I,r,K,p,h,H) { # Function to count real roots of polynomial
  parvec = c( I*h*h, 
              (h*h*r - h*h*H),
              (I - p - h*h*r/K),
              (r-H),
              -(r/K) )
  xeq = polyroot(parvec) 
  Iparts = abs(Im(xeq))
  isreal = ifelse(Iparts<1.e-8,1,0)
  nreal = sum(isreal)
  return(nreal)
}

N.realpos.roots = function(I,r,K,p,h,H) { # Function to count real positive roots of polynomial
  parvec = c( I*h*h, 
              (h*h*r - h*h*H),
              (I - p - h*h*r/K),
              (r-H),
              -(r/K) )
  xeq = polyroot(parvec) 
  Iparts = abs(Im(xeq))
  isreal = ifelse(Iparts<1.e-8,1,0)
  nreal = sum(isreal)
  isrealpos = ifelse(isreal>=0.5 & Re(xeq)>0,1,0)
  nrealpos = sum(isrealpos)
  return(nrealpos)
}

# Function for finding x as a function of tau using largest real root of polynomial
x.maxreal = function(tau.e) {
  H.temp = q*N*tau.e
  x.poly = dx.roots(I,r,K,p,h,H.temp)
  Iparts = abs(Im(x.poly))
  x.real = ifelse(Iparts<1.e-8,x.poly,0)
  x.Rmax = max(Re(x.real))
  return(x.Rmax)
}

# Function for angler utility as a function of tau
U.tau = function(tau){
  x.tau = x.maxreal(tau)
  # positive utility
  U = ( (fish.mult*q*tau*x.tau)^a + S*Tall - S*tau)
  return(U)
}

# Function for angler utility given expected tau
U.tau.tau.e = function(tau,tau.e){
  x.tau.e = x.maxreal(tau.e)
  # positive utility
  U = ( (fish.mult*q*tau*x.tau.e)^a + S*Tall - S*tau)
  return(U)
}

# Function for finding the rational expectations tau - best way to do it
ftau.RE2 = function(taubar) {
  x.tau.e = x.maxreal(taubar)
  # function for positive utility
  U = function(tau) {U = ( (fish.mult*q*tau*x.tau.e)^a + S*Tall - S*tau)}
  taustar.tau.e = optimize(U,interval=c(0,2),maximum=T)
  dev = (taubar - taustar.tau.e$maximum)
  return(dev)
}

dxdt.fun = function(x,I,p,K,H) {  # Function to return rate given a vector of x
  rate = I + r*x - (r/K)*x*x - (p*x*x/(h*h + x*x)) - H*x
  return(rate)
}

findcrit = function(par,p,K)   {  # function that is zero at critical x and H, given p, K
  x = exp(par[1])
  H = exp(par[2])
  denom = h*h + x*x
  dxdt = I + r*x - (r/K)*x*x - (p*x*x/denom) - H*x
  dfdx = r - (2*r*x/K) - ( 2*p*h*h*x/(denom*denom) ) - H
  dev = dxdt*dxdt + dfdx*dfdx 
  return(dev)
}

# End of functions *************************************************************************

# Main Program <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Biological parameters
I = 0.1 # stocking
r=0.5
p=0.5
K=8  
h=0.2

# Add economics; start with solving the angler's problem
# At this stage, assume a value of catchability q
q = 0.2 # catchability
N = 1  # number of anglers
a = 0.5 # exponent of Cobb-Douglas function
S = 0.5 # value of non-fishing activity
Tall = 1 # maximum time available
fish.mult = 1.5 # multiplier for fun of fishing

# Set up loops over K and q
NK = 15
Kvec = seq(3,8,length.out=NK)
Nq = 15
qvec = seq(0.05,0.3,length.out=Nq)

# Make a vector for habitat loss relative to max K, sort, recalculate Kvec
Kloss = sort( 1-(Kvec/max(Kvec)) )
Kvec = max(Kvec)*(1-Kloss)

# matrices to save output
Hcrit.mat = matrix(0,nr=Nq,nc=NK)
Tau.mat = matrix(0,nr=Nq,nc=NK)
Hopt.mat = matrix(0,nr=Nq,nc=NK)
xopt.mat = matrix(0,nr=Nq,nc=NK)
SOS.mat = matrix(0,nr=Nq,nc=NK)
prod.mat = matrix(0,nr=Nq,nc=NK)

# Start loop
for(j in 1:NK) {  # loop over columns
  K = Kvec[j]
  for(i in 1:Nq)   { # loop over rows
    q=qvec[i]
    print('',quote=F)
    print(c('row and col',i,j),quote=F)
    print(c('q and S',q,S),quote=F)
    
    # Find taubar
    tau2 = uniroot(f=ftau.RE2,interval=c(0,4),extendInt='upX')
    taubar = tau2$root
    Tau.mat[i,j] = taubar
    print(c('Rational expectations estimate of opt. tau = ',taubar),quote=F)
    
    # Calculate H that corresponds to taubar -- optimal H
    H.opt = q*N*taubar
    print(c('H.opt = ',H.opt),quote=F)
    Hopt.mat[i,j] = H.opt
    
    # Find equilibria given taubar (i.e. at H.opt)
    x.taubar = x.maxreal(taubar)
    print(c('x eq at H.opt ',round(x.taubar,3)),quote=F)
    
    # Hold x at H.opt
    x.Hopt = x.taubar
    
    # Save equilibrium fish pop and production at equilibrium
    xopt.mat[i,j] = x.taubar
    production = r*x.taubar - (r/K)*x.taubar^2
    prod.mat[i,j] = production
    print(c('Equilibrium biomass & production',round(x.taubar,3),round(production,3)),quote=F)
    
    # Find critical x and H
    xeq = sort( dx.roots(I,r,K,p,h,H.opt) )
    print('polynomial roots given H.opt',quote=F)
    print(round(xeq,3),quote=F)
    
    # Get number of real roots and compute Hcrit if there are 4 real roots
    #NRR = N.real.roots(I,r,K,p,h,H.opt)
    # Get number of real positive roots and compute Hcrit if there are 3 real positive roots
    NRR = N.realpos.roots(I,r,K,p,h,H.opt)
    print(c('Number of real positive roots = ',NRR),quote=F)
    if(NRR==3) {
      xguess = x.Hopt
      Hguess = 3*H.opt
      guess = log( c(xguess,Hguess) )
      # Nelder-Mead
      CritEst = optim(guess,findcrit,gr=NULL,p,K,method='Nelder-Mead') # Nelder Mead option
      # Gradient with bounds; sometimes faster, sometimes fails
      #LL.xH = log(c(max(Re(xeq[3]),1.e-6),1.e-6))  # gradient option with bounds
      #UL.xH = log(c(x.Hopt,10*H.opt))
      #CritEst = optim(guess,findcrit,gr=NULL,p,K,method='L-BFGS-B',lower=LL.xH,upper=UL.xH)
      # Critical point
      CritxH = exp(CritEst$par)
      print(c('Critical x and H ',round(CritxH,3)),quote=F)
      
      # Save results for plotting
      Hcrit.mat[i,j] = CritxH[2]
      SOS.mat[i,j] = Hcrit.mat[i,j] - Hopt.mat[i,j]
    } # end critical point calculation for case with 4 real roots
  }  # end inner loop over rows
}  # end outer loop over columns

# Save results if needed
#save(qvec,Kvec,Kloss,Hopt.mat,SOS.mat,xopt.mat,prod.mat,Tau.mat,file='q+habitat-gradients_V8.Rdata')

# Two-panel plot of SOS and biomass ------------------------------------------------------------------
windows(width=12,height=5)
par(mfrow=c(1,2),mar=c(4,5,3,6) + 0.1, cex.axis=1.6,cex.lab=1.6)

# SOS
# Color image plot for SOS -------------------------------------------------------

# Increase the data density by interpolation
OBJ = list(x=Kloss,y=qvec,z=t(SOS.mat))
Sdense = seq(min(Kloss),max(Kloss),length.out=50)
qdense = seq(min(qvec),max(qvec),length.out=50)
grid.list = list(x=Sdense,y=qdense)
SOSdense = interp.surface.grid(OBJ,grid.list)

# Generate plots

# Plot SOS
image.plot(SOSdense,
           col=rainbow(128,s = 1, v = 1, start = 0, end = max(1, 56)/64, alpha = 0.5),
           #col=topo.colors(128),
           xlab='Habitat Loss',
           ylab='Catchability',
           main='A. SOS')
# overlay contours
zlevels = c(0.02)
contour(x = Kloss, y = qvec, z = t(SOS.mat), #log='y',
        levels=zlevels, labcex=1, col='black',lwd=2, 
        xlab=' ', ylab=' ', add=T)

# Plot biomass with SOS overlaid
# Increase the data density by interpolation
# Increase the data density by interpolation
OBJ = list(x=Kloss,y=qvec,z=t(xopt.mat))
Sdense = seq(min(Kloss),max(Kloss),length.out=50)
qdense = seq(min(qvec),max(qvec),length.out=50)
grid.list = list(x=Sdense,y=qdense)
xdense = interp.surface.grid(OBJ,grid.list)

# Plot equilibrium fish population

image.plot(xdense,
           col=rainbow(128,s = 1, v = 1, start = 0, end = max(1, 56)/64, alpha = 0.5),
           #col=topo.colors(128),
           xlab='Habitat Loss',
           ylab='Catchability',
           main='B. Fish Biomass')
# overlay contours
zlevels = c(1)
contour(x = Kloss, y = qvec, z = t(xopt.mat), #log='y',
        levels=zlevels, labcex=1, col='black',lwd=2, 
        xlab=' ', ylab=' ', add=T)

# overlay SOS contour
zlevels = c(0.02)
contour(x = Kloss, y = qvec, z = t(SOS.mat), #log='y',
        levels=zlevels, labels='SOS boundary', labcex=1.3, method='flattest',
        col='black',lwd=2,lty=3, 
        xlab=' ', ylab=' ', add=T)