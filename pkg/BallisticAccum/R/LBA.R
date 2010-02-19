
LBAsingle = function(rt,resp,Nresp,muB=1,sig2B=1,mu0=-1,sig20=1,a0=1,b0=1,Niters=100,sdMetRT=1,sdMetMuD=.1,sdMetsig2D=.1,sdMetB=.1,sdMett0=2)
  {
    .Call("LBAsingleC",rt,as.integer(resp),as.integer(Nresp),as.integer(length(rt)),muB,sig2B,mu0,sig20,a0,b0,Niters,sdMetRT,sdMetB,sdMetMuD,sdMetsig2D,sdMett0,package="BallisticAccum")
  }

LBAlike = function(rt,muD=-1,sig2D=1,bThresh=2,t0=1)
  {
    .Call("RlogLikelihood",rt, muD, sig2D, bThresh, t0, package="BallisticAccum")
  }

LBAlikeMu = function(mu,rt,j,sig2,bThresh,t0,mu0,sig20)
{
  .Call("RlogCondPostMuD",rt,mu,as.integer(j),sig2,bThresh,t0,mu0,sig20,as.integer(dim(rt)[1]),package="BallisticAccum")
}

LBAsampleTruncRT = function(cut,start,mu,sig2,bThresh,t0,sigMet=.5)
  {
      .Call("RsampleTruncRT",cut,start,mu,sig2,bThresh,t0,sigMet,package="BallisticAccum")
    }


LBAfullLike = function(rt,muD=-1,sig2D=.01,bThresh=2,t0=2)
  {
    .Call("RfullLogLikelihood",rt, muD, sig2D, bThresh, t0, as.integer(dim(rt)[1]), as.integer(dim(rt)[2]), package="BallisticAccum")
  }


