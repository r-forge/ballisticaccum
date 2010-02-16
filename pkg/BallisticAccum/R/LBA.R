
LBAsingle = function(rt,resp,Nresp,muB=1,sig2B=1,mu0=1,sig20=1,a0=1,b0=1,Niters=100,sdMetRT=.1,sdMetMuD=.1,sdMetsig2D=.1,sdMetB=.1,sdMett0=.1)
  {
    .Call("LBAsingleC",rt,as.integer(resp),as.integer(Nresp),as.integer(length(rt)),muB,sig2B,mu0,sig20,a0,b0,Niters,sdMetRT,sdMetB,sdMetMuD,sdMetsig2D,sdMetB,sdMett0,package="BallisticAccum")
  }
