
## R-code for Evolutionary Design of Experiments (EDoE)
## Copyright N. Packard, ProtoLife, May 2008


# for a gaussian data object
# with cliffs on some fraction of the Gaussian peaks.
# object is a list, with
# ctr = matrix of gaussian centers, nrow=#centers, ncol=dimension
# stddev = vector of stddevs of the gaussians
# weight = vector of weights (height of each gaussian)
# cliff = vector, if component nonzero, that gaussians has a cliff in direction of that componenent's axis

## cliffGauss = class for data generator using Gaussians with cliffs.
## essentially a container for all the information to
## specify Gaussians and cliffs.
## number of Gaussians = nrow(ctr)
##
cliffGauss = function(ctr = matrix(c(0.5,0.5),nrow=1),
  stddev = 0.25,
  weight=1,
  cliff = 1)                            # vec of which coords have cliffs
{
  if(!is.matrix(ctr))
    ctr = t(ctr)                        #change vec to matrix
  if(length(stddev)!=length(weight)
     || length(stddev)!= nrow(ctr)
     || length(weight) != nrow(ctr)
     || length(cliff) != nrow(ctr))
    stop("incompatible arguments...")
  if(nrow(ctr)!= length(stddev))
    stop("ctr incompatible with stddev")
  ## collect class member objects in a list
  res = list()
  res$d = ncol(ctr)
  res$ctr = ctr
  res$stddev = stddev
  res$weight = weight
  res$cliff = cliff
  class(res) = "cliffGauss"
  res
}

## method sampleData for the cliffGauss class
## samples cliffdat generator
##
sampleData = function(ob,...) UseMethod("sampleData",ob)	#GGG: will apply to object 'ob' method sampleData for class 'class(ob)',
															#GGG which implies there should be a method called sampleData.[classOfOb]
sampleData.cliffGauss = function(       # use R dispatch to create method for cliffGauss class
    cg,                                 # cliffGauss object
    x)                                  # matrix of points to sample the cliffGauss generator
{
  if(!is.matrix(x))
    x=t(x)                              # turn vec into 1-row matrix
  if(ncol(x) != ncol(cg$ctr))
    stop("dimensions not compatible...")
  gg = numeric(0)
  for(j in 1:nrow(x)){                  # for each row of x, compute the superposition of cliffed Gaussians
    g = 0
    for(i in 1:nrow(cg$ctr)){
      ctr = cg$ctr[i,]
      r2 = sum((ctr-x[j,])*(ctr-x[j,]))
                                        #		r2 = sum((rev(ctr)-x)*(rev(ctr)-x)) # why rev ctr???
                                        #		g = g + exp(-r2/(cg$stddev[i]*cg$stddev[i])) / sqrt(2*pi)
      cliff = cg$cliff
      if(cliff[i] == 0)               # no cliff; add gaussian normally
        g = g + cg$weight[i]*exp(-r2/(cg$stddev[i]*cg$stddev[i])) # max of weight[i] at peak...
      else{
        coord = cliff[i]
        if(coord< -nrow(cg$ctr) || coord > nrow(cg$ctr) || coord==0)
          stop("bad coord spec for cliff")
        if(coord <0){
          coord = - coord
          if((ctr-x[j,])[coord]<0)  # add cliff in neg direction
            g = g + cg$weight[i]*exp(-r2/(cg$stddev[i]*cg$stddev[i])) # max of weight[i] at peak...
        }
        else{
          if((ctr-x[j,])[coord]>0)  # add cliff in pos direction
            g = g + cg$weight[i]*exp(-r2/(cg$stddev[i]*cg$stddev[i])) # max of weight[i] at peak...
        }
      }
    }
    g = g/nrow(cg$ctr)           # normalization for max...
    gg = rbind(gg,g)
  }
  gg
}


## rancliff creates a cliffGauss object:
## chooses N Gaussians with random centers,
## random widths, with a fraction ncliff of them with cliffs in
## random directions
## returns a cliffGauss object.
##
rancliff = function(N=4,                # number of Gaussians
  d=2,                                  # dimension
  mnmn = 0.5,                           # mean of dist of means
  mnsdev = 0.2,                         # sdev of dist of means
  sdevmn = 1/N,                          # mn of dist of sdevs
  sdevsdev = 1/N,                        # sdev of dist of sdevs
  mnlim = c(0,1),                       # limit of location of mn
  ncliff = 0.5                          # fraction of gaussians that are cliffs
  )
{
    ctr = matrix(nrow=N,ncol=d)
# means first
   tstmn = function(m,lm){
    if(m<lm[1]) return(1) else if(m>lm[2]) return(1) else return(0)}
  for(i in 1:N){
    ctr[i,] =rnorm(d,mean=mnmn,sd=mnsdev)
    for(j in 1:d)
      while(tstmn(ctr[i,j],mnlim))
        ctr[i,j] =rnorm(1,mean=mnmn,sd=mnsdev)
  }
# now sdevs
  sdev = rnorm(N,mean=mnsdev,sd=sdevsdev)
  sdev = abs(sdev)                      # to make sure all sdevs>0
# now weights
  wt = 0.5+runif(N)/10
  wt = rep(0.5,N)
# now cliffs
  ncliff = ncliff*N
  cliff = 1+floor(runif(ncliff)*d)                            # choose random coords
  cliff = cliff * (2*as.numeric(rnorm(length(cliff))<0) -  1) # randomize sign
  cliff = c(cliff,rep(0,N-length(cliff)))
  dat = cliffGauss(ctr,sdev,wt,cliff)
  dat
}
