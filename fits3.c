/*
cc fits3.c -o fits3 util/wt/util.o util/wt/sort.o -lm -m32

check mean std values with
fits2 3chy5nul A 0
cat 3chy5nul/random.dat | awk '{n++; s+=$2; ss+=$2*$2; if(n>1) {m=s/n; z=sqrt((ss-2*m*s+n*m*m)/(n-1)); print n,m,z}}' | tail -1
*/
#include "util/wt/incl/util.h"
#include <math.h>

#define N 200
#define R 2.0
#define T 0.0001
#define U 0.0002
#define RMST0 0
#define RMSV0 20
#define RMSW0 20
#define BASE 2


main(argc,argv) int argc; char *argv[];
{
int EVD;
int i, j, k, m, n, maln,naln,maxaln; 
float x, y, z, a, b, c;
float r, s, t, t0s,t1s,t2s,t3s,tset, smin,smax, range;
float fit, sum, ssum, wsum, mean, stdv, ymax, zmax, z0;
float aopt, bopt, copt, best;
float *d,*w, e[N], f[N], g[N], h[N];
float set, damp = 2.0;
int   freq[N];
int   invert, imax, nmax;
float lave, mave, nave;
float rlp, ss;
int   *p;
char line[222], file[222], mode;
FILE *dat, *plt;
	d = (float*)alloca(sizeof(float)*100000);
	w = (float*)alloca(sizeof(float)*100000);
	p = (int*)alloca(sizeof(int)*100000);
	invert = 0;
	mode = *argv[2];
	if (islower(mode)) { mode = toupper(mode); invert = 1; }
	EVD = 0;
	if (argc>3 && *argv[3]=='1') EVD = 1;
	Pi(EVD) NL
	strcpy(line,argv[1]);
	strcat(line,"/");
	strcat(line,"pdb2pdb1.dat");
	dat = fopen(line,"r");
	read_line(dat,line);
	sscanf(line,"%d %f %f %f %f %f %f", &naln, &r, &s, &t0s, &t1s, &t2s, &t3s);
	fclose(dat);
	maln = (1*naln)/2;
	Pi(maln) Pi(naln) NL
	smin = 99999.9;
	sum = ssum = smax = wsum = 0.0;
	for (i=0; i<N; i++) f[i] = 0;
	strcpy(file,argv[1]);
	strcat(file,"/");
	strcat(file,"nat1rev2.dat");
	dat = fopen(file,"r");
	n = 0;
	nmax = 0;
	mave = 0.0;
	while (1)
	{ float s0, t0, t1, t2, t3;
	  int	io = read_line(dat,line);
		if (io<=0) break;
		sscanf(line,"%d %f %f %f %f %f %f", &i, &r, &s0, &t0, &t1, &t2, &t3);
		if (i>nmax) nmax = i;
		mave += (float)i;
		n++;
	}
	mave /= (float)n;
Pi(naln) Pi(maln) Pi(nmax) Pi(n) Pr(mave)
	fclose(dat);
	dat = fopen(file,"r");
	n = 0;
	nave = 0.0;
	while (1)
	{ float s0, t0, t1, t2, t3;
	  int	io = read_line(dat,line);
		if (io<=0) break;
		sscanf(line,"%d %f %f %f %f %f %f", &i, &r, &s0, &t0, &t1, &t2, &t3);
		if ((float)i<mave) {
			lave += (float)i;
		} else {
			nave += (float)i;
			n++;
		}
	}
	nave /= (float)n;
	lave /= (float)n;
Pi(n) Pr(lave) Pr(nave) NL
	fclose(dat);
	n = 0;
	dat = fopen(file,"r");
	while (1)
	{ float wa, wb, fi, fm=(float)maln, fn=nave;
	  float s0, t0, t1, t2, t3;
	  int	io = read_line(dat,line);
		if (io<=0) break;
		sscanf(line,"%d %f %f %f %f %f %f", &i, &r, &s0, &t0, &t1, &t2, &t3);
		if (t0<0.000001) t0=0.000001;
		if (t1<0.000001) t1=0.000001;
		if (t2<0.000001) t2=0.000001;
		if (t3<0.000001) t3=0.000001;
		fi = (float)i;
		wa = fn-fi;
		if (fi>fn) w[n]=1.0; else w[n] = exp(-wa*wa*0.01)+0.01;
		if (mode=='A') s = r;
		if (mode=='B') s = s0;
		if (mode=='C') s = t0;
		if (mode=='D') s = t2;
		if (mode=='R') s = 100.0-fi/(BASE+r);
		if (mode=='S') s = r/fi;
		if (mode=='T') { fi=fi-RMST0; s = r/(1.0-exp(-fi*fi*T)); }
		if (mode=='U') s = r/((1.0-exp(-fi*fi*U))*sqrt(fi));
		if (mode=='V') {
/*			fi=fi-RMSV0; s = log(100/t0)/fi; */
			s = (100-t0)/((1.0-exp(-fi*fi*U))*sqrt(fi));
		} /* L+G'98 */
		if (mode=='W') {
/*			fi=fi-RMSW0; s = log(100/t2)/fi; */
			s = (100-t2)/((1.0-exp(-fi*fi*U))*sqrt(fi));
		} /* Gauss */
		if (invert) s = 100.0 - 100.0/(BASE+s*10.0);
		if (mode=='X') s = 10.0-log(s0);
		if (mode=='Y') s = 200.0-sqrt(s0);
		d[n] = s;
		sum += s;
		wsum += d[n]*w[n];
		ssum += s*s;
		if (s>smax) smax=s;
		if (s<smin) smin=s;
		n++;
	}
	fclose(dat);
	if (!n) {
		y = ymax = z = zmax = -9.9;
		Pr(y) Pr(ymax) Pr(z) Pr(zmax) NL
		exit(1);
	}
	range = smax-smin;
	mean = sum/(float)n;
	stdv = sqrt((ssum-mean*mean*(float)n)/(float)(n-1));
Pi(n) Pr(mean) Pr(stdv) NL
Pr(smin) Pr(smax) Pr(range) NL
Pr(sum) Pr(wsum) NL
	sort(0,d,0,p,n,1);
	ssum = 0.0;
	for (i=n-1; i>-1; i--) {
		set = d[p[i]];
		s = d[p[i]]*w[p[i]];
		ssum += s;
		x = 100.0*ssum/wsum;
		if (x>5.0) break;
	}
Pr(set) NL
	strcpy(line,argv[1]);
	strcat(line,"/");
	strcat(line,"pdb1pdb2.dat");
	dat = fopen(line,"r");
	z = zmax = -9.9;
	while (1) { float fi, s0; int io = read_line(dat,line);
		if (io<=0) break;
		sscanf(line,"%d %f %f", &i, &r, &s0);
		if (i<40) continue;
		if (i<(int)lave) continue;
		fi = (float)i;
		if (mode=='A') s = r;
		if (mode=='B') s = s0;
		if (mode=='R') {
			s = 100.0-fi/(BASE+r);
		}
		if (mode=='S') {
			s = r/fi;
		}
		if (mode=='T') {
			if (fi < RMST0) continue;
			fi = fi - RMST0;
			s = r/(1.0-exp(-fi*fi*T));
		}
		if (mode=='U') {
			s = r/((1.0-exp(-fi*fi*U))*sqrt(fi));
		}
		if (mode=='V') {
/*
			fi = fi - RMSV0;
			s = log(100/t0s)/fi;
*/
			r = 100-t0s;
			s = r/((1.0-exp(-fi*fi*U))*sqrt(fi));
		}
		if (mode=='W') {
/*
			fi = fi - RMSW0;
			s = log(100/t2s)/fi;
*/
			r = 100-t2s;
			s = r/((1.0-exp(-fi*fi*U))*sqrt(fi));
		}
		if (invert) s = 100.0 - 100.0/(BASE+s*10.0);
		if (mode=='X') s = 10.0-log(s0);
		if (mode=='Y') s = 200.0-sqrt(s0);
		if (s>set) continue;
		Pi(i) Pr(r) Pr(s) Pr(smin) Pr(set) NL
		z = 1.0 - s/set;
		x = (float)i/(float)nave;
		z *= 1.0-exp(-x*x*damp);
		x=exp(-x*x*damp); Pr(x) 
		if (z>zmax) { imax=i; zmax=z; }
	}
	Pr(z) Pr(zmax) Pi(imax) NL
/*
	if (zmax<0.0) { 
		y = ymax = z = zmax = -9.9;
		Pr(y) Pr(ymax) Pr(z) Pr(zmax) NL
		exit(1);
	}
*/
	fclose(dat);
	Pi(n) Pr(mean) Pr(stdv) NL
	Pr(smin) Pr(smax) Pr(range) NL
	sum = ssum = wsum = 0.0;
	for (i=0; i<N; i++) f[i] = 0;
	for (i=0; i<n; i++) {
		s = 50.0+(d[i]-smin)*100.0/range;
		m = (int)s;
		sum += s*w[i]; ssum += s*s*w[i];
		if (m<N && m>-1) f[m] += w[i];
		wsum += w[i];
	}
	mean = sum/wsum;
	stdv = sqrt((ssum-mean*mean*wsum)/wsum);
	Pi(n) Pr(mean) Pr(stdv) NL
	e[0] = f[0];
	for (i=1; i<N; i++) e[i] = e[i-1]+f[i];
/* weighting? 
	for (i=1; i<N; i++) {
		if (i<N/2) e[i] = e[i-1]+f[i];
		      else e[i] = e[i-1]+f[i]*(1.0-0.01*(float)i);
	}
*/
	best = 9999999.9;
	for (j=0; j<=100; j++) { float bs = 0.01*(float)(j-50);
		b = mean;
		if (bs<0.0) b /= 1.0-bs; else b *= 1.0+bs;
		for (k=0; k<=100; k++) { float cs = 0.04*(float)(k-50); 
			c = stdv;
			if (cs<0.0) c /= 1.0-cs; else c *= 1.0+cs;
			if (EVD) {
				for (i=0; i<N; i++) {
		                      	x = ((float)i-b)/c;
					g[i] = 500.0*exp(x)*exp(-exp(x))/c;
					h[i] = 500.0*(1.0-exp(-exp(x)));
				}
			} else {
                        	for (i=0; i<N; i++) { float cc=2.0*c*c;
				x=(float)i-b; g[i] = 20.0*exp(-x*x/cc);
				}
                        	h[0] = g[0];
                        	for (i=1; i<N; i++) h[i] = g[i]+h[i-1];
			}
			a = wsum/h[N-1];
			for (i=0; i<N; i++) h[i] = h[i]*a;
			a *= 10.0;
			fit = 0.0;
			for (i=0; i<N; i++) { float dif = e[i]-h[i];
				fit += dif*dif*(1.0-e[i]/e[N-1]);;
			}
			fit = sqrt(fit);
			if (fit<best) { best=fit; aopt=a; bopt=b; copt=c; }
		}
	}
	Pr(aopt) Pr(bopt) Pr(copt) Pr(best) NL
	if (EVD) {
		for (i=0; i<N; i++) {
	               	x = ((float)i-bopt)/copt;
			g[i] = 500.0*exp(x)*exp(-exp(x))/copt;
			h[i] = 500.0*(1.0-exp(-exp(x)));
		}
	} else {
        	for (i=0; i<N; i++) { float cc=2.0*copt*copt;
		x=(float)i-bopt; g[i] = 20.0*exp(-x*x/cc);
		}
        	h[0] = g[0];
        	for (i=1; i<N; i++) h[i] = g[i]+h[i-1];
	}
	a = wsum/h[N-1];
	for (i=0; i<N; i++) h[i] = h[i]*a;
	fit = 0.0;
	for (i=0; i<N; i++) { float dif = e[i]-h[i];
		fit += dif*dif*(1.0-e[i]/e[N-1]);;
	}
	fit = sqrt(fit);
	Pr(aopt) Pr(bopt) Pr(copt) Pr(fit) NL
	plt = fopen("dist.out","w");
	for (i=0; i<N; i++) fprintf(plt,"%5d %f %f %f %f\n", i, f[i],0.1*e[i], g[i],0.1*h[i]);
	fclose(plt);
	strcpy(line,argv[1]);
	strcat(line,"/");
	strcat(line,"pdb1pdb2.dat");
	dat = fopen(line,"r");
	y = ymax = -9.9;
	while (1)
	{ float fi, s0; int io = read_line(dat,line);
	  double u;
		if (io<=0) break;
		sscanf(line,"%d %f %f", &i, &r, &s0);
		if (i<40) continue;
		if (i<(int)lave) continue;
		fi = (float)i;
		if (mode=='A') s = r;
		if (mode=='B') s = s0;
		if (mode=='R') {
			s = 100.0-fi/(BASE+r);
		}
		if (mode=='S') {
			s = r/fi;
		}
		if (mode=='T') {
			if (fi < RMST0) continue;
			fi = fi - RMST0;
			s = r/(1.0-exp(-fi*fi*T));
		}
		if (mode=='U') {
			s = r/((1.0-exp(-fi*fi*U))*sqrt(fi));
		}
		if (mode=='V') {
/*
			fi = fi - RMSV0;
			s = log(100/t0s)/fi;
*/
			s = (100-t0s)/((1.0-exp(-fi*fi*U))*sqrt(fi));
		}
		if (mode=='W') {
/*
			fi = fi - RMSW0;
			s = log(100/t2s)/fi;
*/
			s = (100-t2s)/((1.0-exp(-fi*fi*U))*sqrt(fi));
		}
		if (invert) s = 100.0 - 100.0/(BASE+s*10.0);
		if (mode=='X') s = 10.0-log(s0);
		if (mode=='Y') s = 100.0-sqrt(s0);
/*
		if (s>set) continue;
*/
		Pi(i) Pr(r) Pr(s) Pr(smin) Pr(set) 
		ss = s = 50.0+(s-smin)*100.0/range;
		if (EVD) {
		        y=(s-bopt)/copt;
			y = 500.0*(1.0-exp(-exp(y)))*a;
			y = y/h[N-1];
			y = -log10(y);
		} else {
			y = (bopt-s)/copt;
		}
		x = (float)i/(float)nave;
		y *= 1.0-exp(-x*x*damp);
		if (y>ymax) { ymax=y; imax=i; }
	}
	Pr(y) Pr(ymax)
	z = 50.0+(0.0-smin)*100.0/range;
	z0 = z;
	rlp = 100.0*(ss-z0)/(bopt-z0);
	plt = fopen("native.dat","w");
	fprintf(plt,"%f\n", rlp);
	fclose(plt);
	if (EVD) {
		z=(z-bopt)/copt;
		z = 500.0*(1.0-exp(-exp(z)))*a;
		z = z/h[N-1];
		z = -log10(z);
	} else {
		z = (bopt-z)/copt;
	}
	z *= 1.0-exp(-x*x*damp);
	zmax = ymax/z;
	z = y/z;
	Pr(z) Pr(zmax) NL
	for (i=0; i<N; i++) freq[i] = 0;
	strcpy(file,argv[1]);
	strcat(file,"/");
	strcat(file,"nat1rev2.dat");
	dat = fopen(file,"r");
	plt = fopen("decoy.dat","w");
	n = 0;
	while (1)
	{ float wa, wb, fi, fm=(float)maln, fn=nave;
	  float s0, t0, t1, t2, t3;
	  int	io = read_line(dat,line);
		if (io<=0) break;
		sscanf(line,"%d %f %f %f %f %f %f", &i, &r, &s0, &t0, &t1, &t2, &t3);
		if (t0<0.000001) t0=0.000001;
		if (t1<0.000001) t1=0.000001;
		if (t2<0.000001) t2=0.000001;
		if (t3<0.000001) t3=0.000001;
		fi = (float)i;
		wa = fn-fi;
		if (fi>fn) w[n]=1.0; else w[n] = exp(-wa*wa*0.01)+0.01;
		if (mode=='A') s = r;
		if (mode=='B') s = s0;
		if (mode=='C') s = t0;
		if (mode=='D') s = t2;
		if (mode=='R') s = 100.0-fi/(BASE+r);
		if (mode=='S') s = r/fi;
		if (mode=='T') { fi=fi-RMST0; s = r/(1.0-exp(-fi*fi*T)); }
		if (mode=='U') s = r/((1.0-exp(-fi*fi*U))*sqrt(fi));
		if (mode=='V') {
/*			fi=fi-RMSV0; s = log(100/t0)/fi; */
			s = (100-t0)/((1.0-exp(-fi*fi*U))*sqrt(fi));
		} /* L+G'98 */
		if (mode=='W') {
/*			fi=fi-RMSW0; s = log(100/t2)/fi; */
			s = (100-t2)/((1.0-exp(-fi*fi*U))*sqrt(fi));
		} /* Gauss */
		s = 50.0+(s-smin)*100.0/range;
		rlp = 100.0*(s-z0)/(bopt-z0);
		fprintf(plt,"%f\n", rlp);
		freq[(int)s]++;
		if (EVD) {
		        y=(s-bopt)/copt;
			y = 500.0*(1.0-exp(-exp(y)))*a;
			y = y/h[N-1];
			y = -log10(y);
		} else {
			y = (bopt-s)/copt;
		}
		rlp = y/z0;
		if (rlp<0.0) continue;
		n++;
	}
	fclose(plt);
	fclose(dat);
	plt = fopen("freq.dat","w");
	for (i=0; i<N; i++) fprintf(plt,"%4d %4d\n",i,freq[i]);
	fclose(plt);
}
