
########################################################
## scan a data generator,
## make a plot of the surface
## for a data generator with d>2, sample hyperplane of 1st two coordinates
## all other coords set to 0.5
scanData = function(datgen,Ngrid=50,docontour=T,...)
{
    xx = (1:Ngrid)/Ngrid
    yy = xx
                                        # magic numbers to scale plots right... 
    xl=c(0.033,0.967); yl=c(0.033,0.967)

    vec = rep(.5,datgen$d - 2)             # to fill out dimensions > 2

    res = matrix(ncol=Ngrid,nrow=Ngrid)
    r=numeric(0)
    for(i in 1:Ngrid)
        for(j in 1:Ngrid){
            res[i,j] = sampleData(datgen,c(xx[j],yy[i],vec))
        }
    image(res, col = terrain.colors(100), axes = T, ...)
    if(docontour)
        contour(res,xlim=xl,ylim=yl,add=T,...)
                                        #  filled.contour(res,xlim=xl,ylim=yl,...)
    par(new=T)
    for(i in 1:nrow(datgen$ctr)){
        if(i==1)
            plot(datgen$ctr[i,2],datgen$ctr[i,1],xlim=xl,ylim=yl,xlab="",ylab="",xaxt='n',yaxt='n')
        else
            points(datgen$ctr[i,2],datgen$ctr[i,1])
    }
}

## wireframe plot of data generator
wiredat = function(datgen,               # synthetic data generator
    Ngrid=100,                            # number of grid points each axis
    ...)
{
    xx = (1:Ngrid)/Ngrid
    yy = xx
                                        # magic numbers to scale plots right... 
    xl=c(0.033,0.967); yl=c(0.033,0.967)

    vec = rep(.5,datgen$d - 2)             # to fill out dimensions > 2

    res = matrix(ncol=Ngrid,nrow=Ngrid)
    r=numeric(0)
    for(i in 1:Ngrid)
        for(j in 1:Ngrid){
            res[i,j] = sampleData(datgen,c(xx[j],yy[i],vec))
        }
    persp(xx,yy,res,phi=30,theta=100,axes=F,box=F)
                                        #  persp(xx,yy,res,phi=30,theta=100,axes=T,box=T,xlab="x1",ylab='x2',zlab='Response',...)
}

## example:
## source('cliffGauss.R')
## source('scanData.R')
## datgen = rancliff()
## scanData(datgen)
