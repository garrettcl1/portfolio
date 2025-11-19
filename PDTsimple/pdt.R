source('amorenet.R')
source('cliffGauss.R')
source('choosePoints.R')
source('testModel.R')

## R-code for Evolutionary Design of Experiments (EDoE)
## Copyright N. Packard, ProtoLife, May 2008

#############################################################################
# pdt loop

dopdt = function(model,                 # model 
    modelpar,                           # model params
    chooser,                            # point chooser
    chooserpar,                         # point chooser params
    dat,                                # data generator object e.g. from rancliff()
    Npop=100,
    Ngen=10,
    pause=T,
    doplot=T,
    saveplot=F){
    res = list()                        # to accumulate results
    
    ## First generation:
    ## random sample of points
    chooserpar$N = Npop
    chooser = do.call(chooser,chooserpar)
    d = dat$d
    train = matrix(runif(d*Npop),ncol=d) # random points to train on (GGG: gen 1)
    targ = sampleData(dat,train)            # sample data generator 
    modelpar$x = train
    modelpar$y = targ
    curmodel = do.call(model,modelpar)
    cat("---------------------------------Finished generation ",1,": ")
    res[[1]] = testModel(curmodel,dat,train,errmx=0.2)
    if(saveplot)
      dev.copy2pdf(file='pdfplot1.pdf')
    ## First gen done...
    ## main loop
    for(g in 2:Ngen){
        if(pause){
            cat("Enter 'c' to continue, 'Q' to exit.\n")
            browser()                       #pause ...
        }
        ## choose new pts from top modeled
        new = choosePoints(chooser,curmodel) # choose new data points using appropriate choosePoints method for chooser
        train = rbind(train,new)
        targ = rbind(targ,sampleData(dat,new))
        ## train the latest model
        modelpar$x = train
        modelpar$y = targ
        curmodel = do.call(model,modelpar)
        ## take a look
        cat("---------------------------------Finished generation ",g,":\n")
        res[[g]] = testModel(curmodel,dat,train,errmx=0.2,new=Npop,doplot=doplot)
        if(saveplot)
          dev.copy2pdf(file=paste('pdfplot',g,'.pdf',sep=""))
    }
    res
}


## Usage:
    
## datgen = rancliff()
## mychooser = chooseRandFit
## mychooserparams = list(fac=100)
## mymodel=amorenetNorm
## mymodelparams = list(layers=c(2,10,10,1),maxit=500)           #train, targ will be added

## result = dopdt(mymodel,mymodelparams,mychooser,mychooserparams,datgen,Npop=10)




## utility to get targets (fitness) of last generation
get.targs = function(pdt){
    Npop = length(pdt[[1]]$targs)
    targs = lapply(pdt,function(x){x$targs[(length(x$targs)-Npop):length(x$targs)]})
    targall=numeric(0)
    for(tt in targs){
        foo = tt[order(tt)]
        targall = c(targall,foo)
    }
    targall
}
## usage:
## after pdtresult = dopdt(...)
## barplot(get.targs(pdtresult))

get.meanerrs = function(pdt){
    foo = lapply(pdt,function(x){x$meanerr})
    foo = unlist(foo)
    foo
}


