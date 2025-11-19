
## R-code for Evolutionary Design of Experiments (EDoE)
## Copyright N. Packard, ProtoLife, May 2008


## object for a data generator to sample a distribution of Gaussians
# 
# returns datGauss object, which is a list with:
# ctr = matrix of gaussian centers, nrow=#centers, ncol=dimension
# stddev = vector of stddevs of the gaussians
# weight = vector of weights (height of each gaussian)

########################################################
datGauss = function(ctr = matrix(c(0.5,0.5),nrow=1),
		stddev = 0.25,
		weight=1)
{
	if(nrow(ctr)!= length(stddev))
		stop("ctr incompatible with stddev")
        if(!is.matrix(ctr))
          ctr = t(ctr)                  # turn vector to matrix
	res = list()
        res$d = ncol(ctr)
	res$ctr = ctr
	res$stddev = stddev
	res$weight = weight
	class(res) = "datGauss"             # turn list into a class
	res
}

## method sampleData for the datGauss class
## samples cliffdat generator
## 
## 
sampleData = function(ob,...) UseMethod("sampleData",ob) #
sampleData.datGauss = function(gd,x)
{
  if(!is.matrix(x))
    x=t(x)                              # turn vec into 1-row matrix
  if(ncol(x) != ncol(gd$ctr))
    stop("dimensions not compatible...")
  gg = numeric(0)
  for(j in 1:nrow(x)){
    g=0
    for(i in 1:nrow(gd$ctr)){
      ctr = gd$ctr[i,]
      r2 = sum((ctr-x[j,])*(ctr-x[j,]))   #
                                        # r2 = sum((rev(ctr)-x)*(rev(ctr)-x)) # why rev ctr???
                                        # g = g + exp(-r2/(gd$stddev[i]*gd$stddev[i])) / sqrt(2*pi)
      g = g + gd$weight[i]*exp(-r2/(gd$stddev[i]*gd$stddev[i])) # max of weight[i] at peak...
    }
    g = g/nrow(gd$ctr)           # normalization for max...
    gg = rbind(gg,g)
  }
gg
}

## ranGauss creates a datGauss object:
## chooses N Gaussians with random centers,
## random widths, with a fraction ncliff of them with cliffs in
## random directions
## returns a cliffdat object.
## 
rangauss = function(N=4,                # number of Gaussians
  d=2,                                  # dimension
  mnmn = 0.2,                           # mean of dist of means
  mnsdev = 0.2,                         # sdev of dist of means
  sdevmn = 1/N,                          # mn of dist of sdevs
  sdevsdev = 1/N,                        # sdev of dist of sdevs
  wtsdev = 0.1,							# sdev of weight (height) of peaks
  mnlim = c(0,1)                        # limit of location of mn
  )
{
# means first
  ctr = matrix(nrow=N,ncol=d)
  tstmn = function(m,lm){
    if(m<lm[1]) return(1) else if(m>lm[2]) return(1) else return(0)}
  for(i in 1:N){
    ctr[i,] =rnorm(d,mean=mnmn,sd=mnsdev)
    for(j in 1:d)
      while(tstmn(ctr[i,j],mnlim))
        ctr[i,j] =rnorm(1,mean=mnmn,sd=mnsdev)
  }
# now sdevs
  sdev = rnorm(N,mean=sdevmn,sd=sdevsdev)
  sdev = abs(sdev)                      # to make sure all sdevs>0
# now weights
  wt = (1-wtsdev)+wtsdev*runif(N)
  dat = datGauss(ctr,sdev,wt)
  dat
}



