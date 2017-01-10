/*
cc stest.c -o stest util/wt/util.o util/aa/stutest.o util/wt/sort.o -lm

*/
#include "util/wt/incl/util.h"
#include "util/aa/incl/student.h"
#include <math.h>


main(argc,argv) int argc; char *argv[];
{
int i, j, k, m, n, in1,in2, maln,naln,maxaln; 
float r, s, t, t0s,t1s,t2s,t3s,tset, smin,smax, range;
float fit, sum, ssum, mean1, mean2, stdv1, stdv2;
float *a, rmsd, x,y,z, fin;
char line[222], file[222];
FILE *dat, *plt;
float	damp = 30.0;
int	margin = 10;
	a = (float*)alloca(sizeof(float)*100000);
	sscanf(argv[2],"%f", &damp);
	sscanf(argv[3],"%d", &margin);
	Ps(argv[1]) Pr(damp) NL
	damp *= damp;
	// read overall.dat (full length native)
	strcpy(file,argv[1]);
	strcat(file,"/");
	strcat(file,"overall.dat");
	dat = fopen(file,"r");
	maln = 999;
	naln = -99;
	in1 = 0;
	while (1)
	{ int	io = read_line(dat,line);
		if (io<=0) break;
		sscanf(line,"%d %f %f %f %f %f %f", &n, &rmsd, &s, &t0s, &t1s, &t2s, &t3s);
		if (n < maln) maln = n;
		if (n > naln) naln = n;
		x = (float)n;
		// r = a*sqrt(x) : a = r/sqrt(N)
		// r = a*sqrt(x)*(1-exp(-x*x/(30*30))) : a = r/(sqrt(N)*(1-exp(-x*x/(30*30))))
		a[in1] = rmsd/(sqrt(x)*(1-exp(-x*x/damp)));
		in1++;
	}
	fclose(dat);
	if (naln-maln < margin*2) {
		n = (maln+naln)/2;
		maln = n - margin;
		naln = n + margin;
	}
	Pi(in1) Pi(maln) Pi(naln)
	sum = ssum = 0.0;
	for (i=0; i<in1; i++) {
		sum += a[i];
		ssum += a[i]*a[i];
	}
	fin = (float)in1;
	mean1 = sum/fin;
	stdv1 = sqrt((ssum-mean1*mean1*fin)/(fin-1.0));
	Pr(mean1) Pr(stdv1)
	// read random.dat
	strcpy(file,argv[1]);
	strcat(file,"/");
	strcat(file,"randoms.dat");
	dat = fopen(file,"r");
	in2 = 0;
	while (1)
	{ float s0, t0, t1, t2, t3;
	  int	io = read_line(dat,line);
		if (io<=0) break;
		sscanf(line,"%d %f %f %f %f %f %f", &n, &r, &s0, &t0, &t1, &t2, &t3);
		if (n<maln || n>naln) continue;
		x = (float)n;
		a[in2] = r/(sqrt(x)*(1-exp(-x*x/damp)));
		in2++;
	}
	fclose(dat);
	Pi(in2)
	sum = ssum = 0.0;
	for (i=0; i<in2; i++) {
		sum += a[i];
		ssum += a[i]*a[i];
	}
	fin = (float)in2;
	mean2 = sum/fin;
	stdv2 = sqrt((ssum-mean2*mean2*fin)/(fin-1.0));
	Pr(mean2) Pr(stdv2)
	NL
	stutest(mean1,mean2,stdv1*stdv1,stdv2*stdv2,in1,in2);
}
