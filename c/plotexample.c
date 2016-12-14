#include "chidecays.h"
#include <math.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

	double beta = 1e-13;
	
	double minm = 0.1;
	double maxm = 10;
	int n = 100;

	ratiosplines* rs = ratio_spline_init();
	
	for (int i=0; i<=n; i++) {
		double m = exp(log(minm)+i*(log(maxm)-log(minm))/n);
		printf("%g", m);
		printf("\t%g", gammatot(m, beta, rs));
		printf("\t%g", gammaee(m, beta));
		printf("\t%g", gammamumu(m, beta));
		printf("\t%g", gammatautau(m, beta));
		printf("\t%g", gammapipitot(m, beta, rs));
		printf("\t%g", gammaKKtot(m, beta, rs));
		printf("\t%g", gammagammagamma(m, beta));
		printf("\t%g", gammagg(m, beta));
		printf("\t%g", gammassKcut(m, beta));
		printf("\t%g", gammacc(m, beta));
		printf("\t%g", gammabb(m, beta));
		printf("\t%g", gammatt(m, beta));
		printf("\n");
	}
	return 0;
}
