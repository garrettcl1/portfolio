    library('AMORE')

## R-code for Evolutionary Design of Experiments (EDoE)
## Copyright N. Packard, ProtoLife, May 2008

# package amorenet more or less like nnet...
# with predict method

amorenet = function(x, y, weights, size=10, layers=NA, Wts, mask,
     linout = FALSE, entropy = FALSE, softmax = FALSE,
     censored = FALSE, skip = FALSE, rang = 0.7, decay = 0,
     maxit = 100, Hess = FALSE, trace = TRUE, MaxNWts = 1000,
     abstol = 1.0e-4, reltol = 1.0e-8, ...){
  if(length(layers)==1) neur = c(ncol(x),size,ncol(y)) # check for NA
  else neur=layers
  net <- newff(	n.neurons= neur,
               learning.rate.global=1e-2,
               momentum.global=0.5,
                        error.criterium="LMS",
               Stao=NA,
               hidden.layer="tansig", 
#             output.layer="purelin",    # default for amore
               output.layer="tansig",
             method="ADAPTgdwm")
  result <- train(net, x, y, error.criterium="LMS",
                  report=TRUE, show.step=maxit/10, n.shows=10 )
  class(result) = "amorenet"
  result
}



## for non-release amore with RBF-nets...
xxxamorenet = function(x, y,net.type="MLPnet", weights, size=10, layers=NA, Wts, mask,
     linout = FALSE, entropy = FALSE, softmax = FALSE,
     censored = FALSE, skip = FALSE, rang = 0.7, decay = 0,
     maxit = 100, Hess = FALSE, trace = TRUE, MaxNWts = 1000,
     abstol = 1.0e-4, reltol = 1.0e-8, ...){
  if(length(layers)==1) neur = c(ncol(x),size,ncol(y)) # check for NA
  else neur=layers
  net <- newff(net.type=net.type,
					n.neurons= neur,
               learning.rate.global=1e-2,
               momentum.global=0.5,
                        error.criterium="LMS",
               Stao=NA,
               hidden.layer="tansig", 
#             output.layer="purelin",    # default for amore
               output.layer="tansig",
             method="ADAPTgdwm")
  result <- train(net, x, y, error.criterium="LMS",
                  report=TRUE, show.step=maxit/10, n.shows=10 )
  class(result) = "amorenet"
  result
}

predict.amorenet = function(anet,x,...){
  if(!is.matrix(x) && !is.data.frame(x)) x = t(x)            #turn vec x into matrix.
  sim.MLPnet(anet$net,x,...)}

amorenetNorm = function(x,y,...){
# take means and sd col by col
  xmn = numeric(0); xsd = numeric(0)
  for(i in 1:ncol(x)) xmn = c(xmn,mean(x[,i]))
  for(i in 1:ncol(x)) xsd = c(xsd,sqrt(var(x[,i])))
  xnorm = x
  for(i in 1:ncol(x))
    xnorm[,i] = (x[,i]-xmn[i])/xsd[i]
  ymn = numeric(0); ysd = numeric(0)
  for(i in 1:ncol(y)) ymn = c(ymn,mean(y[,i]))
  for(i in 1:ncol(y)) ysd = c(ysd,sqrt(var(y[,i])))
  ynorm = y
  for(i in 1:ncol(y))
    ynorm[,i] = (y[,i]-ymn[i])/ysd[i]
  anet = amorenet(xnorm,ynorm,...)
  res = list(xmn = xmn,
    xsd = xsd,
    ymn = ymn,
    ysd = ysd,
    anet = anet)
  class(res) = "amorenetNorm"
  res
}
  
predict.amorenetNorm = function(anet,x,...){
  if(!is.matrix(x) && !is.data.frame(x)) x = t(x)            #turn vec x into matrix.
  if(ncol(x)!=length(anet$xmn))
     stop("predict.amorenetNorm:  data not compatible")
  xnorm = x
  for(i in 1:ncol(x)) xnorm[,i] = (x[,i]-anet$xmn[i])/anet$xsd[i]
  yraw = sim.MLPnet(anet$anet$net,xnorm,...)
  ynorm = yraw
  for(i in 1:ncol(yraw)) 
    ynorm[,i] = yraw[,i]*anet$ysd[i] + anet$ymn[i]
  ynorm
}


nnetNorm = function(x,y,...){
# take means and sd col by col
  xmn = numeric(0); xsd = numeric(0)
  for(i in 1:ncol(x)) xmn = c(xmn,mean(x[,i]))
  for(i in 1:ncol(x)) xsd = c(xsd,sqrt(var(x[,i])))
  xnorm = x
  for(i in 1:ncol(x))
    xnorm[,i] = (x[,i]-xmn[i])/xsd[i]
#  ymn = numeric(0); ysd = numeric(0)
#  for(i in 1:ncol(y)) ymn = c(ymn,mean(y[,i]))
#  for(i in 1:ncol(y)) ysd = c(ysd,sqrt(var(y[,i])))
#  normalize y to be in 0-1 for nnet's logistic output
  ynorm = y
  ymn=numeric(ncol(y)); yslop=numeric(ncol(y))
  for(i in 1:ncol(y)){
    ymn[i] = min(y[,i]); ymx = max(y[,i]); yslop[i] = 1/(ymx-ymn)
    ynorm[,i] = yslop[i]*(y[,i]-ymn)
  }
  net = nnet(xnorm,ynorm,...)
  res = list(xmn = xmn,
    xsd = xsd,
    ymn = ymn,
    yslop = yslop,
    net = net)
  class(res) = "nnetNorm"
  res
}
  
predict.nnetNorm = function(net,x,...){
  if(!is.matrix(x) && !is.data.frame(x)) x = t(x)            #turn vec x into matrix.
  if(ncol(x)!=length(net$xmn))
     stop("predict.nnetNorm:  data not compatible")
  xnorm = x
  for(i in 1:ncol(x)) xnorm[,i] = (x[,i]-net$xmn[i])/net$xsd[i]
  yraw = predict(net$net,xnorm,...)
  ynorm = yraw
#  browser()
  for(i in 1:ncol(yraw)) 
    ynorm[,i] = yraw[,i]/net$yslop[i] + net$ymn[i]
  ynorm
}


