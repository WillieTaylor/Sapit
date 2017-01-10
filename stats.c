/*
cc stats.c -o stats util/wt/util.o util/aa/student.o -lm -m32

*/
#include "util/wt/incl/util.h"
#include "util/aa/incl/student.h"

main (argc, argv) int argc; char *argv[]; {
int	i, j, k, m, n1,n2;
float	sum, sumsq, mean1,stdev1, mean2,stdev2;
char	line[222];
FILE	*dat, *out;
int	freq[200];
	for (i=0; i<200; i++) freq[i] = 0;
	dat = fopen("decoy.dat","r");
	sum = sumsq = 0.0;
	n1 = 0;
	while (1) { int io; float x;
		io = read_line(dat,line);
		if (io<=0) break;
		sscanf(line,"%f", &x);
		if (x<0.0) continue;
		sum += x; sumsq += x*x;
		n1++;
		if (x<0.0 || x>199.99) continue;
		freq[(int)x]++;
	}
	mean1 = sum/(float)n1;
	stdev1 = sqrt((sumsq-2.0*mean1*sum+(float)n1*mean1*mean1)/(float)(n1-1));
Pi(n1) Pr(mean1) Pr(stdev1) NL
	out = fopen("decoy.plot.dat","w");
	for (i=0; i<200; i++) fprintf(out,"%d %d\n", i,freq[i]);
	fclose(out);
	fclose(dat);
	for (i=0; i<200; i++) freq[i] = 0;
	dat = fopen("native.dat","r");
	sum = sumsq = 0.0;
	n2 = 0;
	while (1) { int io; float x;
		io = read_line(dat,line);
		if (io<=0) break;
		sscanf(line,"%f", &x);
		if (x<0.0) continue;
		sum += x; sumsq += x*x;
		n2++;
		if (x<0.0 || x>199.99) continue;
		freq[(int)x]++;
	}
	mean2 = sum/(float)n2;
	stdev2 = sqrt((sumsq-2.0*mean2*sum+(float)n2*mean2*mean2)/(float)(n2-1));
Pi(n2) Pr(mean2) Pr(stdev2) NL
	out = fopen("native.plot.dat","w");
	for (i=0; i<200; i++) fprintf(out,"%d %d\n", i,freq[i]);
	mean2 = sum/(float)n2;
	stutest(mean1,mean2, stdev1*stdev1,stdev2*stdev2, n1,n2);
}
