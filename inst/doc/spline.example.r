
# example of how to do least squares approximation with splines
# maintained assumption: #data points = #basis functions (=> no OLS needed)
# approximate f(x,y,z): tensor product of 3 univariate splines. this uses the RcppSimpleTensor library to deal with the high (>2) dimension

rm(list=ls(all=TRUE))
# load libraries
library(splines)
library(RcppSimpleTensor)

# Usage examples
# -------------

a <- array(1:8,dim=c(2,2,2))
b <- array(1,dim=c(2,1))
B <- array(1,dim=c(2,2))

### matrix mults

# multiply elements of a with b along common index i (col of a), keeping indices j and k 
a.col <- TI( a[i,j,k] * b[i], j + k)

# multiply elements of a with b along common index j (row of a), keeping indices i and k 
a.row <- TI( a[i,j,k] * b[j], i + k)

# multiply elements of a with B along common indices i and j
aB.k <- TI( a[i,j,k] * B[i,j], k)

# multiply elements of a with B along common indices j and k
aB.i <- TI( a[i,j,k] * B[j,k], i)


### other operations

# sum over powers along index k
a.pow <- TI( a[i,j,k]^b[k], i + j )

# sum over cosines between elements of a and B
a.cos <- TI( cos(a[i,j,k],B[i,j]), k )


# Spline Example
# ---------

### Data Setup

# 3 dimensions x,y,z

# integration setup: choose number of integration nodes and method
library(statmod)
num.int.z <- 50
int.z <- gauss.quad(n=num.int.z,kind="hermite")
 
num.int.y <- 40
int.y <- gauss.quad(n=num.int.y,kind="hermite")

# number of evaluation points in each dimension
num.x <- 10
num.y <- 8
num.z <- 4

# select degree of splines
degree <- 3

# get spline knot vector defined on nodes
xknots <- c(rep(0,times=degree),seq(0,10,le=num.x),rep(10,times=degree)) 
zknots <- c(rep(min(int.z$nodes),times=degree),seq(int.z$nodes[1],int.z$nodes[num.int.z],le=num.z),rep(max(int.z$nodes),times=degree))
yknots <- c(rep(min(int.y$nodes),times=degree),seq(int.y$nodes[1],int.y$nodes[num.int.y],le=num.y),rep(max(int.y$nodes),times=degree))

# get grid points where to evaluate function

xdata <- as.array(seq(0,10,length=length(xknots)-degree-1))
ydata <- as.array(seq(int.y$nodes[1],int.y$nodes[num.int.y],length=length(yknots)-degree-1))
zdata <- as.array(seq(int.z$nodes[1],int.z$nodes[num.int.z],length=length(zknots)-degree-1))

# design matrices
X <- splineDesign(knots=xknots,x=xdata)
Y <- splineDesign(knots=yknots,x=ydata)
Z <- splineDesign(knots=zknots,x=zdata)


# define the function that will give us data (data generating process)
DGP <- function( x,y,z ) {
    return((x + y - 5)^2 + (z-5)^2)
}

# Fast Multidimensional Array Operations
# --------------------------------------

# evaluate function in traditional way: either loop over or 
# do an apply on a grid of values
data.grid <- expand.grid(x=xdata,y=ydata,z=zdata)
trad.time <- system.time(traditional <- array( with(data.grid, mapply(DGP,x,y,z)), dim=c(length(xdata),length(ydata),length(zdata))))

# evaluate function with RcppSimpleTensor
# define a tensor function to calculate function values:
RcppVals <- tensorFunction( R[i,j,k] ~ (X[i] + Y[j] - 5)^2 + (Z[k] - 5)^2 )
# read: return array indexed by [i,j,k], defined as (x[i] + y[j] - 5)^2 + (z[k]-5)^2
Rcpp.time <- system.time( RcppArray <- RcppVals(xdata,ydata,zdata) )
sum(abs(traditional - RcppArray))
print(rbind(trad.time,Rcpp.time))



# use RcppSimpleTensor to evaluate 3 Dimensional Spline
# -----------------------------------------------------


##### plot original data at 2 values of z

zind <- 1:length(zdata)
persp(x=xdata,y=ydata,z=RcppArray[ , ,head(zind,1)],main=paste("original at z =",head(zind,1)),ticktype="detailed",xlab="x",ylab="y",zlab="value",theta=300,phi=30)
persp(x=xdata,y=ydata,z=RcppArray[ , ,tail(zind,1)],main=paste("original at z =",tail(zind,1)),ticktype="detailed",xlab="x",ylab="y",zlab="value",theta=300,phi=30)

# calculate spline coeffs by solving system of equations exactly:
b <- solve(kronecker(Z,kronecker(Y,X)), as.vector(RcppArray))

# put coefficients into a 3-dimensional array
bb <- array(b,dim=c(length(xdata),length(ydata),length(zdata)))

##### predict function at new, random values

# new random data
new_xdata <- sort( runif(n=50,min=min(xdata),max=max(xdata)) )
new_ydata <- sort( runif(n=55,min=min(int.y$nodes),max=max(int.y$nodes) ) )
new_zdata <- sort( runif(n=31,min=min(int.z$nodes),max=max(int.z$nodes) ) )

# basis for new values
new_X <- splineDesign(knots=xknots,x=new_xdata)
new_Y <- splineDesign(knots=yknots,x=new_ydata)
new_Z <- splineDesign(knots=zknots,x=new_zdata)

##### use Rcpp function TI() to get approximation: sum over certain dimensions

# define RcppSimpleTensor function
spline.eval <- tensorFunction( R[nx,ny,nz] ~ coeffs[mx,my,mz] * Xbase[nx,mx] * Ybase[ny,my] * Zbase[nz,mz] )
pred.vals <- spline.eval( bb, new_X, new_Y, new_Z )
      
# or simply with inline:
TIpred.vals <- TI( bb[m1,m2,m3] * new_X[n1,m1] * new_Y[n2,m2] * new_Z[n3,m3], n1+n2+n3)
sum(pred.vals - TIpred.vals)

zind <- 1:length(new_zdata)
# test by plotting at different values of Z
persp(x=new_xdata,y=new_ydata,z=pred.vals[ , ,head(zind,1)],ticktype="detailed",xlab="new_x",ylab="new_y",zlab="approx",theta=300,phi=30)
persp(x=new_xdata,y=new_ydata,z=pred.vals[ , ,tail(zind,1)],ticktype="detailed",xlab="new_x",ylab="new_y",zlab="approx",theta=300,phi=30)


# integrate w.r.t. z and y
# ------------------------

##### obtain data on integration nodes

# calculate basis on integration nodes
Int.base.z <- splineDesign(knots=zknots,x=int.z$nodes)
Int.base.y <- splineDesign(knots=yknots,x=int.y$nodes)
Intdata <- spline.eval( bb, X, Int.base.y, Int.base.z )

yweights <- int.y$weights
zweights <- int.z$weights


#### integrate

#i.e. weighted sum over corresponding dimensions
RcppIntFun <- tensorFunction( R[nx] ~ Data[nx,ny,nz] * yweight[ny] * zweight[nz] )
Int.y.z <- RcppIntFun( Intdata, yweights, zweights )

# or inline

TI.Int.y.z <- TI( Intdata[nx,ny,nz] * yweights[ny] * zweights[nz], nx )
sum(Int.y.z - TI.Int.y.z)

dim(Int.y.z)


# integrate w.r.t. z only 
# ------------------------

Int.z <- TI( Intdata[nx,ny,nz] * zweights[nz], nx + ny)
dim(Int.z)
persp(x=xdata,y=int.y$nodes,z=Int.z,main="Integral w.r.t z",ticktype="detailed",xlab="xdata",ylab="int.y$nodes",zlab="Value",theta=300,phi=30)


# integrate w.r.t. y only 
# ------------------------

Int.y <- TI( Intdata[nx,ny,nz] * yweights[ny], nx + nz )
persp(x=xdata,y=int.z$nodes,z=Int.y,main="Integral w.r.t y",ticktype="detailed",xlab="xdata",ylab="int.z$nodes",zlab="Value",theta=300,phi=30)

