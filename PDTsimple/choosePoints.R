
## R-code for Evolutionary Design of Experiments (EDoE)
## Copyright N. Packard, ProtoLife, May 2008

## point chooser objects
# all must have method choosem(ob,...)
# [choose already taken...

choosePoints = function (object, ...) 
UseMethod("choosePoints",object)

###########################################################################
# choose random points, no constraints
chooseRand = function(N=100)  {
res = list(N=N)                            #vestigal
class(res)="chooseRand"
res
}

## choosePoints method for chooseRand class
## pass model param for consistency of format with other choosers 
## 
choosePoints.chooseRand = function(chooser,model,d=2)
{
  N=chooser$N
  matrix(runif(N*d),ncol=d)
}

############################################################################
## chooseRandFit class: # GGG container of parameters for response-proportional sampling
##  choose random points weighted by fitness according to model
##  make amp > 0 (bug?)
chooseRandFit = function(N=100,         # number to be chosen
  d=2,                                  # dimension of chosen
  fac=20,                               # population of pts = fac*N from which top N will be chosen
  amp=1                                 # amplification of more fit
                                        # amp = 1 => linear dependence on fitness
                                        # amp > 1 => super-lin dependence on fitness
  )  {          #always passed an arg..
res = list(N=N,d=d,fac=fac,amp=amp)
class(res)="chooseRandFit"
res
}


## choosePoints method for chooseRandFit class
## ### GG: actual function that runs response-proportional sampling, based on parameters in 'chooser' list 
choosePoints.chooseRandFit = function(chooser,model)
{
  cat("chooseRandFit choosing... ")
  N=chooser$N
  d=chooser$d
  fac=chooser$fac
  amp=chooser$amp

  raw = matrix(runif(fac*N*d),ncol=d)
  fit = predict(model,raw)
  mx = max(fit); mn = min(fit)
  fit = (fit-mn)/(mx-mn)
# following formula for amplification doesn't seem to work...
#  fit = fit^(1/amp)                     # note 1/amp because 0<fit<1
  idx = sample(1:length(fit),N,prob=fit)
  res = raw[idx,]
  browser()
  cat("done.\n")
  res
}

##########################################################################
# choose top
#	choosetop = function(model,train,N,res = 0.01){
  ## do the coarse grain...
#  trans = function(xx,r=res){round(((xx-mean(xx))/(max(xx)-min(xx)))/r)}
#  xx = trans(train[,1])
#  yy = trans(train[,2])
#  for(i in 1:nrow(train))
    ## unfinished...
#}    

###########################################################################
# choose smart random ...
# later...


