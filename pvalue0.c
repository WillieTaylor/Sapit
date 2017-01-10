/*
cc pvalue.c -o pvalue util/wt/util.o util/wt/sort.o -lm

*/
#include "util/wt/incl/util.h"
#include <math.h>


main(argc,argv) int argc; char *argv[];
{
int i, j, k, m, n, maln,naln,maxaln; 
float r, s, t, t0s,t1s,t2s,t3s,tset, smin,smax, range;
float fit, sum, ssum, wsum, mean, stdv, ymax, zmax;
float *a, b, rmsd, x,y,z;
char line[222], file[222];
FILE *dat, *plt;
float	damp = 30.0*30.0; // should be adjusted to fit data
int	margin = 10;
int	better = 0;
	a = (float*)alloca(sizeof(float)*100000);
	sscanf(argv[2],"%f", &damp);
	sscanf(argv[3],"%d", &margin);
	Pr(damp) Pi(margin) NL
	damp *= damp;
	// read overall.dat (full length native)
	strcpy(file,argv[1]);
	strcat(file,"/");
	strcat(file,"overall.dat");
	dat = fopen(file,"r");
	read_line(dat,line);
	sscanf(line,"%d %f %f %f %f %f %f", &naln, &rmsd, &s, &t0s, &t1s, &t2s, &t3s);
	x = (float)naln;
	b = rmsd/(sqrt(x)*(1-exp(-x*x/damp)));
	fclose(dat);
	// read random.dat
	strcpy(file,argv[1]);
	strcat(file,"/");
	strcat(file,"random.dat");
	dat = fopen(file,"r");
	i = 0;
	while (1)
	{ float s0, t0, t1, t2, t3;
	  int	io = read_line(dat,line);
		if (io<=0) break;
		sscanf(line,"%d %f %f %f %f %f %f", &n, &r, &s0, &t0, &t1, &t2, &t3);
		if (n<naln-margin || n>naln+margin) continue;
		x = (float)n;
		// r = a*sqrt(x) : a = r/sqrt(N)
		// r = a*sqrt(x)*(1-exp(-x*x/(30*30))) : a = r/(sqrt(N)*(1-exp(-x*x/(30*30))))
		a[i] = r/(sqrt(x)*(1-exp(-x*x/damp)));
		i++;
	}
	fclose(dat);
	n = i;
	Pi(n)
	sum = ssum = 0.0;
	for (i=0; i<n; i++) {
		sum += a[i];
		ssum += a[i]*a[i];
		if (a[i] < b) better++;
	}
	mean = sum/(float)n;
	stdv = sqrt((ssum-mean*mean*(float)n)/(float)(n-1));
	Pr(mean) Pr(stdv)
	x = (float)naln;
	y = rmsd/(sqrt(x)*(1-exp(-x*x/damp)));
	Pr(y) Pr(mean-y)
	y = mean - y;
	z = y/stdv;
	b = (float)better/(float)n;
	Pi(n) Pi(better) Pr(b) Pr(z) NL
}
