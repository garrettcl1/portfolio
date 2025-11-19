

## R-code for Evolutionary Design of Experiments (EDoE)
## Copyright N. Packard, ProtoLife, May 2008


#  modeltst... = old versions, using kriging stuff
#  tstmodel = newest version
#  uses data generator objects as in datGauss.R, cliffGauss.R
#  model object must have predict method.

source('datGauss.R')
source('cliffGauss.R')
source('scanData.R')

testModel = function(model,datagen,points,errmx = NA,new=0,Ngrid=50,doplot=F)
{                                       # Assume model has predict method
    rtn = list()                        # to stash results
    if(!is.matrix(points))
        points = t(points)                  # vec -> matrix
    if(ncol(points) != datagen$d)
        stop("dimension mismatch")
                                        # compute errors, etc
    preds = numeric(nrow(points))
    targs = numeric(nrow(points))
    err = numeric(nrow(points))
    for(i in 1:nrow(points)){
        preds[i] = predict(model,t(points[i,]))
        targs[i] = sampleData(datagen,t(points[i,]))
        tst = abs(targs[i] - preds[i])
        if(length(tst)!=1) browser()
        err[i] = tst
    }
    if(is.na(errmx))  errmx = max(err)
    cat("max err = ",max(err),"\n")
    cat("mean err = ",mean(err),"\n")
    rankcor = cor(rank(targs),rank(preds))
    cat("rank cor = ",rankcor,"\n")
    ## stash return values
    rtn$maxerr = max(err)
    rtn$meanerr = mean(err)
    rtn$rankcor = rankcor
    rtn$err = err
    rtn$targs = targs
                                        # magic numbers to scale plots right for overlay...  
    xl=c(0.033,0.967); yl=c(0.033,0.967)
    if(doplot){
                                        # draw response surface
        scanData(datagen,Ngrid=Ngrid)
                                        # show points with circle proportional to error
        par(new=T)
        for(i in 1:nrow(points)){
            if(i==1){                   # do 1st plot to set plot pars
                plot(points[i,2],points[i,1],xlim=xl,ylim=yl,xlab="",ylab="",cex=err[i]*10/errmx,col=4)
                points(points[i,2],points[i,1],pch=20,col=4)
            }
            else{
                if(i<nrow(points)-new)  # new in red old in blue
                    curcol = 4 else curcol= 2
                points(points[i,2],points[i,1],cex=err[i]*10/errmx,col=curcol)
                points(points[i,2],points[i,1],pch=20,col=curcol)
            }
        }
        
                                        # model response surface...
        if(datagen$d>2)
            vec = rep(.5,datagen$d-2)
        else
            vec = numeric(0)
        par(new=T)
        xx = (1:Ngrid)/Ngrid
        yy = (1:Ngrid)/Ngrid
        res = matrix(ncol=Ngrid,nrow=Ngrid)
        r=numeric(0)
        for(i in 1:Ngrid)
            for(j in 1:Ngrid){
                res[i,j] = predict(model,t(c(xx[j],yy[i],vec)))
            }
        contour(res,col=4,labcex=1.5,xlim=xl,ylim=yl,xlab="",ylab="",add=T)
    }
    invisible(rtn)
}

############################################################
##  utilities

# create all pairs of c(a[i],b[j])
xpair = function(a,b){                  
  res = matrix(nrow=(length(a))*(length(a)),ncol=2)
  k = 1

  for(i in 1:length(a))
    for(j in 1:length(a)){
      res[k,] = c(a[i],b[j])
      k = k+1
    }
  res
}

Dim = function(ob,...) UseMethod("Dim",ob)
Dim.Krig = function(kk){ncol(kk$x)}
Dim.gaussdat = function(gd){  gd$d}
Dim.cliffdat = function(cd){  ncol(cd$ctr)}
Dim.nnet = function(nn){  nn$n[1]}
Dim.nnetNorm = function(nn){nn$net$n[1]}
Dim.amorenet = function(nn){length(nn$net$input)}
Dim.amorenetNorm = function(nn){length(nn$anet$net$input)}


##  sampsq sample values of data generator over unit square

sampsq = function(dat,Ngrid=Ngrid){
  xx = (1:Ngrid)/Ngrid; yy = xx
  allpts = xpair(xx,yy)
  if(dat$d>2)
    for(i in 1:(dat$d-2)) allpts = cbind(allpts,rep(.5,nrow(allpts)))
  res = sampleData(dat,allpts)
  res = matrix(res,ncol=Ngrid,nrow=Ngrid)
  res
}

## predict using model over unit square, return prediction surface
## need Dim (above)
predsq = function(model,Ngrid=50){
  xx = (1:Ngrid)/Ngrid; yy = xx
  allpts = xpair(xx,yy)
  d = Dim(model)
  if(d>2)
    for(i in 1:(d-2)) allpts = cbind(allpts,rep(.5,nrow(allpts)))
  res = predict(model,allpts)
  res = matrix(res,ncol=Ngrid,nrow=Ngrid)
  res
}

## only for Krig models...
predsesq = function(model,Ngrid=50){
  xx = (1:Ngrid)/Ngrid; yy = xx
  allpts = xpair(xx,yy)
  d = Dim(model)
  if(d>2)
    for(i in 1:(d-2)) allpts = cbind(allpts,rep(.5,nrow(allpts)))
  res = predict.se(model,allpts)
  res = matrix(res,ncol=Ngrid,nrow=Ngrid)
  res
}



# for(i in 1:100) targ = rbind(targ,sampleData(dat,trset[i,]))

############################################################
##  old model test functions


modeltst = function(d=2,Ninit=10,Ntry=1000)
{
  gdat = rangauss(d=d)
  tstdat = matrix(runif(d*Ninit),ncol=d)
  for(i in 1:nrow(tstdat)){
    if(i==1)
      plot(tstdat[i,2],tstdat[i,1],xlim=xl,ylim=yl,xlab="",ylab="",pch=15)
    else
      points(tstdat[i,2],tstdat[i,1],pch=15)
  }
  par(new=T)
  scandat(gdat,col=4)
# magic numbers to scale plots right...  
  xl=c(0.033,0.967); yl=c(0.033,0.967)

  obsdat = numeric(0)
  for(i in 1:Ninit)
    obsdat[i] = sampleData(gdat,tstdat[i,])
  k = krige(tstdat,obsdat)
  nn = nnet(tstdat,obsdat,size=4)

  par(new=T)
  xx = (1:100)/100
  yy = (1:100)/100
  res = matrix(ncol=100,nrow=100)
  r=numeric(0)
  for(i in 1:100)
  for(j in 1:100){
	res[i,j] = predict(nn,c(xx[j],yy[i]))
  }
  contour(res,col=3,xlim=xl,ylim=yl,xlab="",ylab="",add=T)

  par(new=T)
  xx = (1:100)/100
  yy = (1:100)/100
  res = matrix(ncol=100,nrow=100)
  r=numeric(0)
  for(i in 1:100)
  for(j in 1:100){
	res[i,j] = quality(k,c(xx[j],yy[i]),predict(nn,c(xx[j],yy[i])))
  }
  contour(res,col=4,xlim=xl,ylim=yl,xlab="",ylab="",add=T)

  qx = matrix(nrow=Ntry,ncol=d)
  qq = numeric(Ntry)
  for(i in 1:Ntry){
    qx[i,]=runif(d)
    ob = predict(nn,qx[i,])
    qq[i] = quality(k,qx[i,],ob)
  }
  qqord = rev(order(qq))
  qqbest = qqord[1:(Ntry/10)]
  par(new=T)
  for(i in 1:length(qqbest)){
    ii = qqord[i]
    if(i==1)
      plot(qx[ii,2],qx[ii,1],xlim=xl,ylim=yl,xlab="",ylab="",pch=17,col=4)
    else
      points(qx[ii,2],qx[ii,1],pch=17,col=4)
  }
}

modeltst2 = function(d=2,Ninit=10,Ntry=1000,Niter=16)
{
  gdat = rangauss(d=d)
  tstdat = matrix(runif(d*Ninit),ncol=d)

  obsdat = numeric(0)
  for(i in 1:Ninit)
    obsdat[i] = sampleData(gdat,tstdat[i,])
# magic numbers to scale plots right...  
  xl=c(0.033,0.967); yl=c(0.033,0.967)

  for(ij in 1:Niter){
    obsdat = (obsdat-min(obsdat))/(max(obsdat)-min(obsdat))
    k = krige(tstdat,obsdat)
    nn = nnet(tstdat,obsdat,size=5)
#    nn = Rbf(tstdat,obsdat,size=5)
browser()

    for(i in 1:nrow(tstdat)){
      if(i==1)
        plot(tstdat[i,2],tstdat[i,1],xlim=xl,ylim=yl,xlab="",ylab="",pch=15)
      else
        points(tstdat[i,2],tstdat[i,1],pch=15)
    }
    par(new=T)
    scandat(gdat,col=4,labcex=2)
    par(new=T)
    xx = (1:100)/100
    yy = (1:100)/100
    res = matrix(ncol=100,nrow=100)
    r=numeric(0)
    for(i in 1:100)
      for(j in 1:100){
	res[i,j] = predict(nn,c(xx[j],yy[i]))
      }
    contour(res,col=3,labcex=2,xlim=xl,ylim=yl,xlab="",ylab="",add=T)

    par(new=T)
    xx = (1:100)/100
    yy = (1:100)/100
    res = matrix(ncol=100,nrow=100)
    r=numeric(0)
    for(i in 1:100)
      for(j in 1:100){
	res[i,j] = quality(k,c(xx[j],yy[i]),predict(nn,c(xx[j],yy[i])))
      }
    contour(res,col=4,labcex=2,xlim=xl,ylim=yl,xlab="",ylab="",add=T)

    err = numeric(nrow(tstdat))
    for(i in 1:nrow(tstdat))
      err[i] = abs(sampleData(gdat,tstdat[i,]) - predict(nn,tstdat[i,]))
    errmx = max(err)
    par(new=T)
    for(i in 1:nrow(tstdat)){
      if(i==1)
        plot(tstdat[i,2],tstdat[i,1],xlim=xl,ylim=yl,xlab="",ylab="",cex=err[i]*10/errmx)
      else
        points(tstdat[i,2],tstdat[i,1],cex=err[i]*10/errmx)
    }

    qx = matrix(nrow=Ntry,ncol=d)
    qq = numeric(Ntry)
    for(i in 1:Ntry){
      qx[i,]=runif(d)
      ob = predict(nn,qx[i,])
      qq[i] = quality(k,qx[i,],ob)
    }
    qqord = rev(order(qq))
    qqbest = qqord[1:(Ntry/10)]
    par(new=T)
    ii = qqord[1]
    plot(qx[ii,2],qx[ii,1],xlim=xl,ylim=yl,xlab="",ylab="",pch=17,col=4)
    tstdat = rbind(tstdat,qx[ii,])
    obsdat = c(obsdat,sampleData(gdat,qx[ii,]))
  }

}


modeltst3 = function(d=2,Ninit=10,Ntry=1000,Niter=16)
{
  gdat = rangauss(d=d)
  tstdat = matrix(runif(d*Ninit),ncol=d)
# magic numbers to scale plots right...  
  xl=c(0.033,0.967); yl=c(0.033,0.967)

  obsdat = numeric(0)
  for(i in 1:Ninit)
    obsdat[i] = sampleData(gdat,tstdat[i,])

  mse = function(x,y){mean((x-y)*(x-y))}
  for(ij in 1:Niter){
    obsdat = (obsdat-min(obsdat))/(max(obsdat)-min(obsdat))
    k = krige(tstdat,obsdat)
#    nn = GGNet(tstdat,obsdat) 
    nn = GGNet(tstdat,obsdat,PredQuality=mse) 
#    nn = nnet(tstdat,obsdat,size=5)
#    nn = Rbf(tstdat,obsdat,size=5)
browser()						# stop to see graphics
obmx = max(obsdat)
    for(i in 1:nrow(tstdat)){
      if(i==1)
        plot(tstdat[i,2],tstdat[i,1],xlim=xl,ylim=yl,xlab="",ylab="",pch=0,cex=1+obsdat[i]*10/obmx )
      else
        points(tstdat[i,2],tstdat[i,1],pch=0,cex=1+obsdat[i]*10/obmx)
    }
    par(new=T)
#    scandat(gdat,col=4,labcex=1.5)
    par(new=T)
    xx = (1:100)/100
    yy = (1:100)/100
    res = matrix(ncol=100,nrow=100)
    r=numeric(0)
    for(i in 1:100)
      for(j in 1:100){
	res[i,j] = predict(nn,c(xx[j],yy[i]))
      }
    contour(res,col=3,labcex=1.5,xlim=xl,ylim=yl,xlab="",ylab="",add=T)

    err = numeric(nrow(tstdat))
    for(i in 1:nrow(tstdat))
      err[i] = abs(sampleData(gdat,tstdat[i,]) - predict(nn,tstdat[i,]))
    errmx = max(err)
    par(new=T)
   for(i in 1:nrow(tstdat)){
      if(i==1)
        plot(tstdat[i,2],tstdat[i,1],xlim=xl,ylim=yl,xlab="",ylab="",cex=err[i]*10/errmx)
      else
        points(tstdat[i,2],tstdat[i,1],cex=err[i]*10/errmx)
    }

    par(new=T)
    xx = (1:100)/100
    yy = (1:100)/100
    res = matrix(ncol=100,nrow=100)
    r=numeric(0)
    for(i in 1:100)
      for(j in 1:100){
	res[i,j] = quality(k,c(xx[j],yy[i]),predict(nn,c(xx[j],yy[i])))
      }
    contour(res,col=4,labcex=1.5,xlim=xl,ylim=yl,xlab="",ylab="",add=T)

    qx = matrix(nrow=Ntry,ncol=d)
    qq = numeric(Ntry)
    for(i in 1:Ntry){
      qx[i,]=runif(d)
      ob = predict(nn,qx[i,])
      qq[i] = quality(k,qx[i,],ob)
    }
    qqord = rev(order(qq))
    qqbest = qqord[1:(Ntry/10)]
    par(new=T)
    ii = qqord[1]
    plot(qx[ii,2],qx[ii,1],xlim=xl,ylim=yl,xlab="",ylab="",pch=17,col=4)
    tstdat = rbind(tstdat,qx[ii,])
    obsdat = c(obsdat,sampleData(gdat,qx[ii,]))
  }

}

