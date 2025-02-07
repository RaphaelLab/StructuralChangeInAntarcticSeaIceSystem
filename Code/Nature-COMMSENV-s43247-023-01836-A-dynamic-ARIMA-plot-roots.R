plot.armaroots <- function(x, type, xlab="Real", ylab="Imaginary",
    main=paste("Inverse roots of AR and MA characteristic polynomial"),
    ...)
{
  oldpar <- par(pty='s')
  on.exit(par(oldpar))
  if(type=="AR"){
    plot(c(-1.5,1.5), c(-1.5,1.5), xlab=xlab, ylab=ylab,
       type="n", bty="n", xaxt="n", yaxt="n",
       sub="cross is MA; outside red", main=main, ...)
    axis(1, at=c(-1,0,1), line=0.5, tck=-0.025)
    axis(2, at=c(-1,0,1), label=c("-i","0","i"),
      line=0.5, tck=-0.025)
    circx <- seq(-1,1,l=501)
    circy <- sqrt(1-circx^2)
    lines(c(circx,circx), c(circy,-circy), col='gray')
    lines(c(-2,2), c(0,0), col='gray')
    lines(c(0,0), c(-2,2), col='gray')
  }
  if(length(x) > 0)
  {
    inside <- abs(x) > 1
    points(1/x[inside], pch=ifelse(type=="AR",19,4), col='black')
    if(sum(!inside) > 0)
      points(1/x[!inside], pch=ifelse(type=="AR",19,4), col='red')
  }
}
