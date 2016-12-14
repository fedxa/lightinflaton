#include "chidecays.h"
#include <math.h>
#include <complex.h>
#include <malloc.h>
#include <gsl/gsl_spline.h>


/* Bad style -- global variables */
/* double muF2; */
/* double tau; */

const double GF = 1.16637e-5;
const double me = 0.510998910e-3;
const double mmu = 105.6583668e-3;
const double mtau = 1776.84e-3;
const double mpi = 139.57018e-3;
const double mpi0 = 134.9766e-3;
const double mW = 80.398;
const double mZ = 91.1876;
const double md = 5.04e-3;
const double mu = 2.55e-3;
const double ms = 105e-3;
const double mc = 1.27;
const double mb = 4.20;
const double mt = 171.3;

const double mK = 493.677e-3;
const double mK0 = 497.614e-3;

const double mDmeson = 1865e-3;
const double mBmeson = 5279e-3;

const double alpha = 1. / 137.035999679;
const double alphaSMZ = 0.1176;


/* Quark MSbar masses */
double cfunc (double mu) {
	double x = alphaS(mu)/M_PI;
	if (mu < mc)
		return ( pow(9.0/2*x, 4.0/9) *
			 (1 + x * (0.895 + x * (1.371 + x * 1.952))) );
	else if (mu < mb)
		return ( pow(25.0/6*x, 12.0/25) *
			 (1 + x * (1.014 + x * (1.389 + x * 1.091))) );
	else if (mu < mt)
		return ( pow( (23.0/6 * x), 12.0/23) *
			 (1 + x * ( 1.175+ x * ( 1.501 + x*0.1725))) );
	else
		return ( pow( (7.0/2 * x), 4.0/7) *
			 (1 + (x * (1.398 + (x * (1.793 + (x * -0.6834)))))) );
}
double msbar(double mu) {
	return (104e-3 * (cfunc(mu) /  cfunc (mtau)));
}
double mcbar (double mu) {
	return (mc * ( cfunc(mu) / cfunc(mc)));
}
double mbbar (double mu) {
	return (mb * (cfunc(mu) / cfunc(mb)));
}


/* Lepton decays */
double gammall(double mchi, double beta, double ml) {
	if (mchi > (2*ml))
		return ((beta * ml * ml *
			 pow( (1 - ((4 * ml * ml) / mchi / mchi)), 3.0/2)) /
			(4*M_PI*mchi) );
	else
		return 0;
}

double gammaee (double mchi, double beta) {
	return gammall( mchi, beta, me);
}
double gammamumu (double mchi, double beta) {
	return gammall( mchi, beta, mmu);
}
double gammatautau (double mchi, double beta) {
	return gammall( mchi, beta, mtau);
}
double gammass (double mchi, double beta) {
	double a = (alphaS(mchi) / M_PI);
	return ((3 * gammall(mchi, beta, msbar(mchi))) *
		(1 + (a * (5.67))));
}
double gammassKcut (double mchi, double beta) {
	if (mchi > ( 2 * mK))
		return gammass( mchi, beta);
	else
		return 0;
}
double gammass0 (double mchi, double beta) {
	return (3 * gammall( mchi, beta, ms));
}
double gammacc (double mchi, double beta) {
	if (mchi > (2 * mDmeson)) {
		double a = (alphaS(mchi)/ M_PI);
		return ((3 * gammall( mchi, beta, mcbar(mchi))) *
			(1 + (a * (5.67))));
	} else
		return 0;
}
double gammacc0 (double mchi, double beta) {
	if (mchi > (2 * mDmeson))
		return (3 * gammall( mchi, beta, mc));
	else
		return 0;
}
double gammabb (double mchi, double beta) {
	if (mchi > (2 * mBmeson)) {
		double a = (alphaS( mchi) / M_PI);
		return ((3 * gammall( mchi, beta, mbbar(mchi))) *
			(1 + ( a * (5.67))));
	} else
		return 0;
}
double gammatt (double mchi, double beta) {
	return (3 * gammall( mchi, beta, mt));
}


/* Meson decays */
double gammaf0f0 (double mchi, double beta, double mf) {
	if (mchi > (2 * mf))
		return (beta * mchi *
			pow( (2.0/9 + ((11.0/9 * mf * mf) / mchi / mchi)), 2) *
			sqrt( (1 - ((4 * mf * mf) / mchi / mchi)) ) /
			16 / M_PI);
	else
		return 0;
}

double gammapipi (double mchi, double beta) {
	return (2 * gammaf0f0( mchi, beta, mpi));
}
double gammapi0pi0 (double mchi, double beta) {
	return gammaf0f0( mchi, beta, mpi0);
}
double gammaKKchiral (double mchi, double beta) {
	return (2 * gammaf0f0( mchi, beta, mK));
}
double gammaK0K0chiral (double mchi, double beta) {
	return (2 * gammaf0f0( mchi, beta, mK0));
}


/* These are just for show -- Do not use, not real formulas!!!! */
double gammaKK (double mchi, double beta) {
	double x = gammaKKchiral( mchi, beta);
	double y = ((gammass( mchi, beta) / 2));
	return ((x * y) / (x + y));
}
double gammaK0K0 (double mchi, double beta) {
	double x = gammaK0K0chiral( mchi, beta);
	double y = ( gammass(mchi, beta)/ 2);
	return ((x * y)/ (x + y));
}


double ratiopipimumu_q[] = { 3.00895E-1, 4.19882E-1, 6.00949E-1, 7.97556E-1, 9.42801E-1, 9.74610E-1, 9.87942E-1, 9.95567E-1, 1.01006E+0, 1.05375E+0, 1.09500E+0, 1.16214E+0, 1.24225E+0 };
double ratiopipimumu_v[] = { 1.58601E+0, 7.60993E+0, 1.21441E+1, 1.97691E+1, 7.78151E+1, 1.85547E+2, 2.41353E+2, 2.22743E+2, 8.00972E+1, 4.20843E+1, 2.42281E+1, 9.45603E+0, 1.65230E+0 };

double ratioKKmumu_q[] = { 9.86220E-1, 9.88912E-1, 9.91682E-1, 9.96949E-1, 1.00738E+0, 1.05910E+0, 1.09527E+0, 1.16755E+0, 1.24247E+0, 1.32257E+0, 1.39754E+0, 1.5E+0, 2E+0 };
double ratioKKmumu_v[] = { 1.81770E+0, 1.65447E+1, 4.21244E+1, 5.52994E+1, 6.76958E+1, 6.68872E+1, 6.14375E+1, 4.35612E+1, 3.11097E+1, 2.17555E+1, 1.70559E+1, 1E+1, 0.41E+1 };


struct ratiosplines {
	gsl_interp_accel *accpi;
	gsl_spline *splinepi;
	gsl_interp_accel *accK;
	gsl_spline *splineK;
};


/* This should be called to init the spline structures before calling the KK and pipi functions */
struct ratiosplines* ratio_spline_init() {
	struct ratiosplines *rs = (struct ratiosplines *)malloc(sizeof(struct ratiosplines));
	rs->accpi = gsl_interp_accel_alloc();
	rs->splinepi = gsl_spline_alloc (gsl_interp_linear, sizeof(ratiopipimumu_q)/sizeof(double));
	gsl_spline_init (rs->splinepi, ratiopipimumu_q, ratiopipimumu_v, sizeof(ratiopipimumu_q)/sizeof(double));
	rs->accK = gsl_interp_accel_alloc();
	rs->splineK = gsl_spline_alloc (gsl_interp_linear, sizeof(ratioKKmumu_q)/sizeof(double));
	gsl_spline_init (rs->splineK, ratioKKmumu_q, ratioKKmumu_v, sizeof(ratioKKmumu_q)/sizeof(double));
	return rs;
}


double gammapipitot(double mchi, double beta, ratiosplines* rs) {
  if (mchi < 3.00895E-1)
	  return (gammapipi(mchi, beta) + gammapi0pi0( mchi, beta));
  else if (mchi < 1.24225E+0)
	  return ((gsl_spline_eval(rs->splinepi, mchi, rs->accpi) * gammamumu( mchi, beta)));
  else
	  return -1;
/*    (t (gammagg mchi beta)))) */
}


/* Do not believe this above mchi>1.2!!!! */
double gammaKKtot (double mchi, double beta, ratiosplines* rs) {
  if (mchi <= (2 * mK))
	  return 0;
  else if ((mchi > 9.86220E-1) && (mchi < 2.0))
	  return (gsl_spline_eval(rs->splineK, mchi, rs->accK) * gammamumu( mchi, beta));
  else
	  return (gammaKK ( mchi, beta) + gammaK0K0( mchi, beta));
}




/* Photon decay  */
double complex xfunc (double y) {
	if (y > 1)
		return (atan(1.0 / sqrt (y - 1)));
	else
		return (1.0/2 * (M_PI + ((0 + 1*I) * log( ( (1 + (sqrt (1 - y)))/
							    ( 1- (sqrt (1 - y))))))));
}
double complex x2func (double y) {
	return cpow( xfunc(y), 2);
}
double complex Fw (double y) {
	return (2 + (3 * y * (1 + ((2 - y) * x2func(y)))));
}
double complex Ff (double y) {
	return (-2 * y * (1 + ((1 -y) * x2func(y))));
}
double complex Fgamma (double mchi) {
#define y(m) ((4 * (m) * (m)) / mchi / mchi)
#define term(f, s, t) ((f) * cpow((s), 2) * Ff( y((t))))
	double complex res = Fw(y(mW));
	res += term(1, 1 ,me);
	res += term(1, 1 ,mmu);
	res += term(1, 1 ,mtau);
	res += term(3, -1./3 ,md);
	res += term(3, 2./3 ,mu);
	res += term(3, -1./3 ,ms);
	res += term(3, 2./3 ,mc);
	res += term(3, -1./3 ,mb);
	res += term(3, 2./3 ,mt);
	return res;
#undef y
#undef term
}

/* Be carefull -- it is ok only for small mh! (check at okun') ? */
double gammagammagamma (double mchi, double beta) {
	double complex F = Fgamma(mchi);
	return ( ( beta * mchi * alpha * alpha *
		   creal((F * conj(F))) ) /
		 (128*M_PI*M_PI*M_PI) );
}

double alphaSlambda (double mu, double l, int nf) {
	double b0 = (11 - (2./3 * nf));
	double b1 = (51 - (19./3 * nf));
	double b2 = (2857 + (-5033./9 * nf) + (325./27 * nf * nf));
	double lm = (2 * log( (mu / l)));
	double llm = log(lm);
	return ( (4 * M_PI * (1 +
			    ((-2 * b1 * llm) / b0 / b0 / lm) +
			    ( (4 * b1 * b1 * (pow( (llm - 1./2), 2) +
					      ((b2 * b0) / 8 / b1 / b1) +
					      -5./4)) /
			      pow(b0, 4) / lm / lm))) /
		 b0 / lm);
}

double alphaSrun (double mu, double mu0, double as0, double nf) {
	return (as0/
		(1 + (as0 * ((33 - (2 * nf))/ 6 /M_PI) * log( (mu /mu0)))));
}
/* (defparameter +alphaSmt+ (alphaSrun +mt+ +MZ+ +alphaSMZ+ 5)) */
const double alphaSmt = 0.10784314311759013;
/* (defparameter +alphaSmb+ (alphaSrun +mb+ +MZ+ +alphaSMZ+ 5)) */
const double alphaSmb = 0.21062098036870402;
/* (defparameter +alphaSmc+ (alphaSrun +mc+ +mb+ +alphaSmb+ 4)) */
const double alphaSmc = 0.3163024362555688;
/* Rough approximation for alphaS at scale mu */
/* IMPORTANT: mu is the scale in GeV, while in alphas_f90 the argument is _sqared_ scale! */
double alphaS (double mu) {
  if (mu < mc)
	  return alphaSrun( mu, mc, alphaSmc, 3);
  else if ( mu < mb)
	  return alphaSrun( mu, mb, alphaSmb, 4);
  else if ( mu < mt)
	  return alphaSrun( mu, mt, alphaSmt, 5);
  else
	  return alphaSrun( mu, mt, alphaSmt, 6);
}

double complex Fg (double mchi) {
#define y(m) ((4 * (m) * (m)) / mchi / mchi)
#define term(m) (Ff(y((m))))
	return term(md)+term(mu)+term(ms)+term(mc)+term(mb)+term(mt);
}
double gammagg (double mchi, double beta) {
	double complex F = Fg(mchi);
	return ( ( beta * mchi * pow( alphaS( mchi), 2) *
		   creal( ( F * conj( F))) ) /
		 64 / M_PI / M_PI / M_PI );
}

double gammatotle (double mchi, double beta, ratiosplines* rs) {
	return ( gammaee(mchi, beta) + gammamumu(mchi, beta) + gammatautau(mchi, beta) +
		 gammapipitot(mchi, beta, rs) +
		 gammaKKtot(mchi, beta, rs) +
		 gammagammagamma(mchi, beta) );
}

double gammatothe (double mchi, double beta) {
	return ( gammaee(mchi, beta) + gammamumu(mchi, beta) + gammatautau(mchi, beta) +
		 gammagg(mchi, beta) + gammassKcut(mchi, beta) +
		 gammacc(mchi, beta) + gammabb(mchi, beta) + gammatt(mchi, beta)
		);
}

double gammatot (double mchi, double beta, ratiosplines* rs) {
  if (mchi < 1.23)
	  return gammatotle( mchi, beta, rs);
  else if ( mchi< 1.76)
	  /*  Very wild hand made fit in the intermediate region of masses */
	  return (1.1 * gammaKKtot(mchi, beta, rs));
  else
	  return gammatothe( mchi, beta);
}
