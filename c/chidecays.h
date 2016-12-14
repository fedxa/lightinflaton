#ifndef __CHIDECAYS_H__
#define __CHIDECAYS_H__

/* safegueard for C99 */
#ifndef M_PI
# define M_PI		3.14159265358979323846	/* pi */
#endif

const double GF;
const double me;
const double mmu;
const double mtau;
const double mpi;
const double mpi0;
const double mW;
const double mZ;
const double md;
const double mu;
const double ms;
const double mc;
const double mb;
const double mt;

const double mK;
const double mK0;

const double mDmeson;
const double mBmeson;

const double alpha;
const double alphaSMZ;


/* Lepton decays */
double gammall(double mchi, double beta, double ml);

double gammaee (double mchi, double beta);
double gammamumu (double mchi, double beta);
double gammatautau (double mchi, double beta);
double gammass (double mchi, double beta);
double gammassKcut (double mchi, double beta);
double gammass0 (double mchi, double beta);
double gammacc (double mchi, double beta);
double gammacc0 (double mchi, double beta);
double gammabb (double mchi, double beta);
double gammatt (double mchi, double beta);


/* Meson decays */
double gammaf0f0 (double mchi, double beta, double mf);

double gammapipi (double mchi, double beta);
double gammapi0pi0 (double mchi, double beta);
double gammaKKchiral (double mchi, double beta);
double gammaK0K0chiral (double mchi, double beta);


/* These are just for show -- Do not use, not real formulas!!!! */
double gammaKK (double mchi, double beta);
double gammaK0K0 (double mchi, double beta);


typedef struct ratiosplines ratiosplines;


/* This should be called to init the spline structures before calling the KK and pipi functions */
ratiosplines* ratio_spline_init();

double gammapipitot(double mchi, double beta, ratiosplines* rs);

/* Do not believe this above mchi>1.2!!!! */
double gammaKKtot (double mchi, double beta, ratiosplines* rs);




/* Photon decay  */
/* Be carefull -- it is ok only for small mh! (check at okun') ? */
double gammagammagamma (double mchi, double beta);

/* IMPORTANT: mu is the scale in GeV, while in alphas_f90 the argument is _sqared_ scale! */
double alphaS (double mu);

double gammagg (double mchi, double beta);



double gammatotle (double mchi, double beta, ratiosplines* rs);

double gammatothe (double mchi, double beta);

double gammatot (double mchi, double beta, ratiosplines* rs);


#define bree(mchi, beta, rs) (gammaee((mchi), (beta))/gammatot((mchi), (beta), (rs)))
#define brmumu(mchi, beta, rs) (gammamumu((mchi), (beta))/gammatot((mchi), (beta), (rs)))
#define brtautau(mchi, beta, rs) (gammatautau( (mchi), (beta))/gammatot((mchi), (beta), (rs)))
#define brpipitot(mchi, beta, rs) (gammapipitot( (mchi), (beta), (rs))/gammatot((mchi), (beta), (rs)))
#define brkktot(mchi, beta, rs) (gammaKKtot( (mchi), (beta), (rs))/gammatot((mchi), (beta), (rs)))
#define brgammagamma(mchi, beta, rs) (gammagammagamma( (mchi), (beta))/gammatot((mchi), (beta), (rs)))

#define brgg(mchi, beta, rs) (gammagg( (mchi), (beta))/gammatot((mchi), (beta), (rs)))
#define brss(mchi, beta, rs) (gammass( (mchi), (beta))/gammatot((mchi), (beta), (rs)))
#define brsskcut(mchi, beta, rs) (gammassKcut( (mchi), (beta))/gammatot((mchi), (beta), (rs)))
#define brcc(mchi, beta, rs) (gammacc( (mchi), (beta))/gammatot((mchi), (beta), (rs)))
#define brbb(mchi, beta, rs) (gammabb( (mchi), (beta))/gammatot((mchi), (beta), (rs)))
#define brtt(mchi, beta, rs) (gammatt( (mchi), (beta))/gammatot((mchi), (beta), (rs)))


#endif /* !__CHIDECAYS_H__ */
