/* LBA of Brown and Heathcote (2008) */
/* C code by Richard D. Morey */


#include <R.h>
#include <Rmath.h>  
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R_ext/Utils.h>
#include <Rversion.h>
#include <Rconfig.h>
#include <R_ext/Constants.h>
#include <R_ext/Lapack.h>
#include <R_ext/Random.h>
#include <R_ext/BLAS.h>



void R_CheckUserInterrupt(void);


SEXP LBAsingleC(SEXP rt, SEXP resp, SEXP Nresp, SEXP Ntrials, SEXP muB, SEXP sig2B, SEXP mu0, SEXP sig20, SEXP a0, SEXP b0, SEXP Niters, SEXP sigMetT, SEXP sigMetB, SEXP sigMetMuD, SEXP sigMetsig2D, SEXP sigMett0);

SEXP alloc3Darray(SEXPTYPE mode, int nrow, int ncol, int nface);
int compare_doubles (const void *a, const void *b);
double expXPhiA(double x,double a);
double logPhiAminusPhiB(double a, double b);
double expXphiAminusphiB(double x, double a, double b, int returnlog);
double samplet0(double *allRT, double *muD, double sig2D, double bThresh, int Ntrials, int Nresp, double minRT, double start, double sigMet, int *acc);
double sampleBthresh(double *allRT, double *muD, double sig2D, double t0, double muB, double sig2B,int Ntrials, int Nresp, double start, double sigMet,int *acc);
double samplesig2D(double *allRT, double *muD, double bThresh, double t0, double a0, double b0, int Ntrials, int Nresp, double start, double sigMet, int *acc);
double sampleMuD(double *allRT, int j, double sig2D, double bThresh, double t0, double mu0, double sig20, int Ntrials, double start, double sigMet, int *acc);
double sampleTruncRT(double rtcut, double muD, double sig2D, double bThresh, double t0, double start, double sigMet, int *acc, double *like);
double logCondPostMuD(double *allRT, double muD, int j, double sig2D, double bThresh, double t0, double mu0, double sig20, int Ntrials);
double logCondPostB(double *allRT, double *muD, double sig2D, double bThresh, double t0, double muB, double sig2B, int Nresp, int Ntrials);
double logCondPostsig2D(double *allRT, double *muD, double sig2D, double bThresh, double t0, double a0, double b0, int Nresp, int Ntrials);
double logCondPostt0(double *allRT, double *muD, double sig2D, double bThresh, double t0, double minRT, int Nresp, int Ntrials);
double logLikelihood(double tij, double muD, double sig2D, double bThresh, double t0);

SEXP RlogCondPostMuD(SEXP allRTR, SEXP muDR, SEXP jR, SEXP sig2DR, SEXP bThreshR, SEXP t0R, SEXP mu0R, SEXP sig20R, SEXP NtrialsR);
SEXP RlogLikelihood(SEXP tijR, SEXP muDR, SEXP sig2DR, SEXP bThreshR, SEXP t0R);
SEXP RsampleTruncRT(SEXP cutR, SEXP startR, SEXP muDR, SEXP sig2DR, SEXP bThreshR, SEXP t0R, SEXP sigMetR);

/**
 * Allocate a 3-dimensional array
 *
 * @param mode The R mode (e.g. INTSXP)
 * @param nrow number of rows
 * @param ncol number of columns
 * @param nface number of faces
 *
 * @return A 3-dimensional array of the indicated dimensions and mode
 */
SEXP alloc3Darray(SEXPTYPE mode, int nrow, int ncol, int nface)
{
    SEXP s, t;
    int n;

    if (nrow < 0 || ncol < 0 || nface < 0)
	error("negative extents to 3D array");
    if ((double)nrow * (double)ncol * (double)nface > INT_MAX)
	error("alloc3Darray: too many elements specified");
    n = nrow * ncol * nface;
    PROTECT(s = allocVector(mode, n));
    PROTECT(t = allocVector(INTSXP, 3));
    INTEGER(t)[0] = nrow;
    INTEGER(t)[1] = ncol;
    INTEGER(t)[2] = nface;
    setAttrib(s, R_DimSymbol, t);
    UNPROTECT(2);
    return s;
}


int compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  
  return (*da > *db) - (*da < *db);
}

double expXPhiA(double x,double a){
  double ret=0;
  switch(a<-5){
  case 1:
    ret=exp(x-pow(a,2)/2-log(-a)-.94/pow(a,2)-.5*log(2*M_PI));
    break;
  case 0:
    Rprintf("\n****Tried to use approximation inappropriately. x=%f, a=%f\n\n",x,a);
    break;
  } 
  return(ret);
}

double logPhiAminusPhiB(double a, double b){
  int i;
  double c0=.2316419;
  double c[5]={.319381530,-.356563782,1.781477937,-1.821255978,1.330274429};
  double da=0,db=0,ta,tb,fa,fb,g,lza,lzb;
  if(a<0||b<0){     
    Rprintf("\n****Tried to use approximation inappropriately. a=%f, b=%f\n\n",a,b);
  }
  ta=1/(1+c0*a);
  tb=1/(1+c0*b);  
  lza=-pow(a,2)/2;
  lzb=-pow(b,2)/2;
  for(i=0;i<5;i++){
    da+=c[i]*pow(ta,i+1);
    db+=c[i]*pow(tb,i+1);
  }
  fa=exp(lza+log(da));
  fb=exp(lzb+log(db));
  g=log(fb-fa);
  return(-.5*log(2*M_PI)+g);
}

double expXphiAminusphiB(double x, double a, double b, int returnlog){
  double ret;
  switch(abs(a)>5&&abs(b)>5){
  case 1:
    if(a>5&&b>5){
      ret=x+logPhiAminusPhiB(a,b);
    }else if(a<-5&&b<-5){
      ret=x+logPhiAminusPhiB(-b,-a);
    }else{
      ret=x+log(pnorm(a,0,1,1,0)-pnorm(b,0,1,1,0));
    }
    break;
  default:
    ret=x+log(pnorm(a,0,1,1,0)-pnorm(b,0,1,1,0));
    break;
  }
  if(!returnlog){
    return(exp(ret));
  }else{
    return(ret);
  }
}



SEXP LBAsingleC(SEXP rt, SEXP resp, SEXP Nresp, SEXP Ntrials, SEXP muB, SEXP sig2B, SEXP mu0, SEXP sig20, SEXP a0, SEXP b0, 
				SEXP Niters, SEXP sigMetT, SEXP sigMetB, SEXP sigMetMuD, SEXP sigMetsig2D, SEXP sigMett0)
{
	
	int NrespC=0, NtrialsC=0, *respC, NitersC=0,i=0,j=0,m=0,accT=0,accB=0,accMuD=0,accsig2D=0,acct0=0;
	double *rtC, muBC=0, sig2BC=0, mu0C=0, sig20C=0, a0C=0, b0C=0,minRT=DBL_MAX;
	double bThresh=0, sig2D=0,t0=0,*allRTp,*chainsp, *accRatesp;
	double sigMetTC=0, sigMetBC=0, sigMett0C=0, sigMetMuDC=0, sigMetsig2DC=0, *likeTotalp, likeSingle=0;
	SEXP allRT, chains, returnList, accRates, likeTotal;


	NrespC = INTEGER_VALUE(Nresp);
	NtrialsC = INTEGER_VALUE(Ntrials);
	NitersC = INTEGER_VALUE(Niters);
	respC = INTEGER_POINTER(resp);
	muBC = REAL(muB)[0];
	sig2BC = REAL(sig2B)[0];
	mu0C = REAL(mu0)[0];
	sig20C = REAL(sig20)[0];
	a0C = REAL(a0)[0];
	b0C = REAL(b0)[0];
	rtC = REAL(rt);
	sigMetTC = REAL(sigMetT)[0];
	sigMetBC = REAL(sigMetB)[0];
	sigMett0C = REAL(sigMett0)[0];
	sigMetsig2DC = REAL(sigMetsig2D)[0];
	sigMetMuDC = REAL(sigMetMuD)[0];
	
	double muD[NrespC];
	
	PROTECT(chains = allocMatrix(REALSXP, NitersC, 3 + NrespC));
	PROTECT(allRT = alloc3Darray(REALSXP,NtrialsC,NrespC,NitersC));
	PROTECT(returnList = allocVector(VECSXP, 4));
	PROTECT(accRates = allocVector(REALSXP,5)); 
	PROTECT(likeTotal = allocVector(REALSXP,NitersC-1));
	
	allRTp = REAL(allRT);
	chainsp = REAL(chains);
	accRatesp = REAL(accRates);
	likeTotalp = REAL(likeTotal);
	
	GetRNGstate();
	
	/* Initialize Starting Values */
	for(j=0;j<NrespC;j++){
		muD[j] = -1;
		chainsp[j*NitersC] = muD[j];
		for(i=0;i<NtrialsC;i++){
			if(j==(respC[i]-1)){
				allRTp[j*NtrialsC + i] = rtC[i];
			}else{
				allRTp[j*NtrialsC + i] = rtC[i]*1.3;
			}
		}
	}
	
	for(i=0;i<NtrialsC;i++){
		if(rtC[i]<minRT) minRT = rtC[i];
	}
	
	bThresh = 2;
	sig2D = 0.1;
	t0 = 0.95*minRT;

	chainsp[(NrespC + 0)*NitersC] = sig2D;
	chainsp[(NrespC + 1)*NitersC] = bThresh;
	chainsp[(NrespC + 2)*NitersC] = t0;
	
	/* Start MCMC chain */
	for(m=1;m<NitersC;m++){
			
			/* Sample RTs */
			likeTotalp[m-1]=0;
			for(i=0;i<NtrialsC;i++){
				for(j=0;j<NrespC;j++){
					if((respC[i]-1)==j){
						allRTp[m*(NrespC*NtrialsC) + j*NtrialsC + i] = rtC[i];
					}else{
					  allRTp[m*(NrespC*NtrialsC) + j*NtrialsC + i] = sampleTruncRT(rtC[i],muD[j],sig2D,bThresh,t0,allRTp[(m-1)*(NrespC*NtrialsC) + j*NtrialsC + i],sigMetTC,&accT,&likeSingle);
					  likeTotalp[m-1] += likeSingle;
					}
				}
			}
			
			/*  Sample  muD[j] */
			for(j=0;j<NrespC;j++){
			  //muD[j] = sampleMuD(allRTp,j,sig2D,bThresh,t0,mu0C,sig20C,NtrialsC,muD[j],sigMetMuDC,&accMuD);
				chainsp[j*NitersC + m] = muD[j];
			}
			
			//sig2D = samplesig2D(allRTp,muD,bThresh,t0,a0C,b0C,NtrialsC,NrespC,sig2D,sigMetsig2DC,&accsig2D);
			//bThresh = sampleBthresh(allRTp,muD,sig2D,t0,muBC,sig2BC,NtrialsC,NrespC,bThresh,sigMetBC,&accB);
			//t0 = samplet0(allRTp,muD,sig2D,bThresh,NtrialsC,NrespC,minRT,t0,sigMett0C,&acct0);
		
			chainsp[(NrespC + 0)*NitersC + m] = sig2D;
			chainsp[(NrespC + 1)*NitersC + m] = bThresh;
			chainsp[(NrespC + 2)*NitersC + m] = t0;
		
	}
	
	accRatesp[0] = (double)(accT)/((NrespC-1)*NitersC*NtrialsC);
	accRatesp[1] = (double)(accMuD)/(NrespC*NitersC);
	accRatesp[2] = (double)(accsig2D)/NitersC;
	accRatesp[3] = (double)(accB)/NitersC;
	accRatesp[4] = (double)(acct0)/NitersC;
			
	
	SET_VECTOR_ELT(returnList, 0, chains);
	SET_VECTOR_ELT(returnList, 1, allRT);
	SET_VECTOR_ELT(returnList, 2, accRates);
	SET_VECTOR_ELT(returnList, 3, likeTotal);
	
	UNPROTECT(5);
	
	PutRNGstate();
	
	return(returnList);
	
}

SEXP RlogLikelihood(SEXP tijR, SEXP muDR, SEXP sig2DR, SEXP bThreshR, SEXP t0R)
{
	double a=0,b=0;
	double tij = REAL(tijR)[0];
	double muD = REAL(muDR)[0];
	double sig2D = REAL(sig2DR)[0];
	double bThresh = REAL(bThreshR)[0];
	double t0 = REAL(t0R)[0];
	SEXP ansR;
	
	PROTECT(ansR = allocVector(REALSXP, 1));
	double *ans = REAL(ansR);
	
	if(bThresh<=0){
		ans[0] = log(pnorm((muD-log(tij-t0)-sig2D)/sqrt(sig2D),0,1,1,0));
	}else{
		a = (log(tij-t0) - log(bThresh - 1) - muD + sig2D)/sqrt(sig2D);
		b = (log(tij-t0) - log(bThresh) - muD + sig2D)/sqrt(sig2D);
	    ans[0]  = expXphiAminusphiB(0.0, a, b, 1);	
	}
	ans[0] += -muD + sig2D/2;

	UNPROTECT(1);
	
	return(ansR);
}

SEXP RsampleTruncRT(SEXP cutR, SEXP startR, SEXP muDR, SEXP sig2DR, SEXP bThreshR, SEXP t0R, SEXP sigMetR)
{
  double a=0,b=0;
  double cut = REAL(cutR)[0];
  double muD = REAL(muDR)[0];
  double sig2D = REAL(sig2DR)[0];
  double bThresh = REAL(bThreshR)[0];
  double t0 = REAL(t0R)[0];
  double start = REAL(startR)[0];
  double sigMet = REAL(sigMetR)[0];
  double temp=0;
  int acc=0;

  SEXP ansR;

  PROTECT(ansR = allocVector(REALSXP, 1));
  double *ans = REAL(ansR);

  GetRNGstate();
  ans[0] = sampleTruncRT(cut, muD, sig2D, bThresh, t0, start, sigMet, &acc, &temp);
  PutRNGstate();

  UNPROTECT(1);

  return(ansR);
}


SEXP RlogCondPostMuD(SEXP allRTR, SEXP muDR, SEXP jR, SEXP sig2DR, SEXP bThreshR, SEXP t0R, SEXP mu0R, SEXP sig20R, SEXP NtrialsR)
{
  double *ans;
  int i=0;
  double *allRT, muD, sig2D, bThresh, t0, mu0, sig20;
  int j,Ntrials;
  SEXP ansR;

  PROTECT(ansR = allocVector(REALSXP,1));
  ans = REAL(ansR);
  allRT = REAL(allRTR);
  muD = REAL(muDR)[0];
  j = INTEGER_VALUE(jR)-1;
  sig2D = REAL(sig2DR)[0];
  bThresh = REAL(bThreshR)[0];
  t0 = REAL(t0R)[0];
  mu0 = REAL(mu0R)[0];
  sig20 = REAL(sig20R)[0];
  Ntrials = INTEGER_VALUE(NtrialsR);

  ans[0] = logCondPostMuD(allRT, muD, j, sig2D, bThresh, t0, mu0, sig20, Ntrials);


  UNPROTECT(1);

  return(ansR);
}


double logLikelihood(double tij, double muD, double sig2D, double bThresh, double t0)
{
	double ans=0,a=0,b=0;
	
	if(bThresh<=0){
		ans = log(pnorm((muD-log(tij-t0)-sig2D)/sqrt(sig2D),0,1,1,0));
	}else{
		a = (log(tij-t0) - log(bThresh - 1) - muD + sig2D)/sqrt(sig2D);
		b = (log(tij-t0) - log(bThresh) - muD + sig2D)/sqrt(sig2D);
	    ans  = expXphiAminusphiB(0.0, a, b, 1);	
	}
	return(ans);
}

double logCondPostt0(double *allRT, double *muD, double sig2D, double bThresh, double t0, double minRT, int Nresp, int Ntrials)
{
	double ans=-DBL_MAX;
	int i=0,j=0;
	
	if(t0<0 || t0>minRT) return(ans);
	ans = 0;
	for(i=0;i<Ntrials;i++){
		for(j=0;j<Nresp;j++){
			ans += logLikelihood(allRT[j*Ntrials + i], muD[j], sig2D, bThresh, t0);
		}
	}
	return(ans);
}

double logCondPostsig2D(double *allRT, double *muD, double sig2D, double bThresh, double t0, double a0, double b0, int Nresp, int Ntrials)
{
	double ans=-DBL_MAX;
	int i=0,j=0;
	
	if(sig2D<0) return(ans);
	ans = 0;
	for(i=0;i<Ntrials;i++){
		for(j=0;j<Nresp;j++){
			ans += logLikelihood(allRT[j*Ntrials + i], muD[j], sig2D, bThresh, t0);
		}
	}
	
	ans += Ntrials*Nresp*sig2D/2 - (a0+1)*log(sig2D) - b0/sig2D;
	return(ans);
}

double logCondPostB(double *allRT, double *muD, double sig2D, double bThresh, double t0, double muB, double sig2B, int Nresp, int Ntrials)
{
	double ans=-DBL_MAX;
	int i=0,j=0;
	
	if(bThresh<0) return(ans);
	
	ans = 0;
	for(i=0;i<Ntrials;i++){
		for(j=0;j<Nresp;j++){
			ans += logLikelihood(allRT[j*Ntrials + i], muD[j], sig2D, bThresh, t0);
		}
	}
	
	ans += -pow(log(bThresh)-muB,2)/(2*sig2B);
	return(ans);
}


double logCondPostMuD(double *allRT, double muD, int j, double sig2D, double bThresh, double t0, double mu0, double sig20, int Ntrials)
{
	double ans=0;
	int i=0;
	
	for(i=0;i<Ntrials;i++){
			ans += logLikelihood(allRT[j*Ntrials + i], muD, sig2D, bThresh, t0);
			//Rprintf("i:%d j:%d t:%f\n",i,j,allRT[j*Ntrials + i]);
	}
	
	ans += -Ntrials*muD - 0.5/sig20*pow(muD - mu0,2);
	return(ans);
}

double sampleTruncRT(double rtcut, double muD, double sig2D, double bThresh, double t0, double start, double sigMet, int *acc, double *like)
{
	//double cand = rexp(1)*sigMet+rtcut;
	double cand = rnorm(0,sigMet) + start;
	double u=0;
	double likeCand = 0, likeStart = 0;
	if(cand<rtcut) return(start);
	u = runif(0,1);
	
	likeCand = logLikelihood(cand, muD, sig2D, bThresh, t0);
	likeStart = logLikelihood(start, muD, sig2D, bThresh, t0);
	if(log(u) < (likeCand - likeStart))
		{
			*like = likeCand - muD + sig2D/2;
			*acc = *acc + 1;
			return(cand);
		}else{
			*like = likeStart - muD + sig2D/2;
			return(start);
		}
}

double sampleMuD(double *allRT, int j, double sig2D, double bThresh, double t0, double mu0, double sig20, int Ntrials, double start, double sigMet, int *acc)
{
	double cand = start + rnorm(0,sigMet);
	double u=0;
	double likeCand = 0, likeStart = 0;

	u = runif(0,1);
	
	likeCand = logCondPostMuD(allRT, cand, j, sig2D, bThresh, t0, mu0, sig20, Ntrials);
	likeStart = logCondPostMuD(allRT, start, j, sig2D, bThresh, t0, mu0, sig20, Ntrials);
	if(log(u) < (likeCand - likeStart))
		{
			*acc = *acc + 1;
			return(cand);
		}else{
			return(start);
		}
}

double samplesig2D(double *allRT, double *muD, double bThresh, double t0, double a0, double b0, int Ntrials, int Nresp, double start, double sigMet, int *acc)
{
	double cand = start + rnorm(0,sigMet);
	double u=0;
	double likeCand = 0, likeStart = 0;
	if(cand<0) return(start);
	u = runif(0,1);
	
	likeCand = logCondPostsig2D(allRT, muD, cand, bThresh, t0, a0, b0, Nresp, Ntrials);
	likeStart = logCondPostsig2D(allRT, muD, start, bThresh, t0, a0, b0, Nresp, Ntrials);
	if(log(u) < (likeCand - likeStart))
		{
			*acc = *acc + 1;
			return(cand);
		}else{
			return(start);
		}	
}

double sampleBthresh(double *allRT, double *muD, double sig2D, double t0, double muB, double sig2B,int Ntrials, int Nresp, double start, double sigMet,int *acc)
{
	double cand = start + rnorm(0,sigMet);
	double u=0;
	double likeCand = 0, likeStart = 0;
	if(cand<0) return(start);
	u = runif(0,1);
	
	likeCand = logCondPostB(allRT, muD, sig2D, cand, t0, muB, sig2B, Nresp, Ntrials);
	likeStart = logCondPostB(allRT, muD, sig2D, start, t0, muB, sig2B, Nresp, Ntrials);
	if(log(u) < (likeCand - likeStart))
		{
			*acc = *acc + 1;
			return(cand);
		}else{
			return(start);
		}	
}

double samplet0(double *allRT, double *muD, double sig2D, double bThresh, int Ntrials, int Nresp, double minRT, double start, double sigMet, int *acc)
{
	double cand = start + rnorm(0,sigMet);
	double u=0;
	double likeCand = 0, likeStart = 0;
	
	if(cand<0 || cand>minRT) return(start);
	u = runif(0,1);
	
	
	likeCand = logCondPostt0(allRT, muD, sig2D, bThresh, cand, minRT, Nresp, Ntrials);
	likeStart = logCondPostt0(allRT, muD, sig2D, bThresh, start, minRT, Nresp, Ntrials);
	if(log(u) < (likeCand - likeStart))
		{
			*acc = *acc + 1;
			return(cand);
		}else{
			return(start);
		}	
}


