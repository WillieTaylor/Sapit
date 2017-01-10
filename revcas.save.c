/*
cc revcas.c -o revcas util/wt/util.o util/wt/geom.o util/wt/sort.o util/aa/pdbprot.o util/aa/cones.o util/aa/matrix.o -lm

*/
#include <stdlib.h>
#include <alloca.h>
#include "util/wt/incl/util.h"
#include "util/wt/incl/geom.h"
#include "util/aa/incl/pdbprot.h"
#include "util/aa/incl/matplot.h"
#include "util/aa/incl/matrix.h"

double  drand48();

#define NOUT 10

main(argc,argv)
int argc; char *argv[];
{
int	i, j, m, n, nres, nout, breaks, anyseq = -1;
int	more=2;
Vec	*nat;
char	*res;
float	*acc, *acs;
int	*at, *in;
FILE	*out;
char	mode;
Pdbentry_ *prot;
long    rseed = 1234;
	if (argc>3) {
		sscanf(argv[3],"%d", &rseed);
		Pi(rseed) NL
	}
        srand48(rseed);
	prot = get_pdb(argv[1],1,1);
	cones(prot);
	nres = prot->Chains[0].Aano;
	m = nres+more+2;
	nat = (Vec*)alloca(sizeof(Vec)*m*2);
	res = (char*)alloca(sizeof(char)*m*2);
	acc = (float*)alloca(sizeof(float)*m*2);
	acs = (float*)alloca(sizeof(float)*m+2);
	at = (int*)alloca(sizeof(int)*m+2);
	in = (int*)alloca(sizeof(int)*m+2);
	printf("COMPND   %s\n",prot->Compound);
	printf("SOURCE   %s\n",prot->Source);
	for (i=0; i<nres; i++)
	{ Atom_ *Atoms = prot->Chains[0].Atoms;
               	nat[i+1].x = Atoms[i].X;
               	nat[i+1].y = Atoms[i].Y;
               	nat[i+1].z = Atoms[i].Z;
               	acc[i+1] = Atoms[i].Bfact;
		res[i+1] = Atoms[i].Aa;
	}
	mode = *argv[2];
	if (islower(mode)) {
		mode = toupper(mode);
		flips(nat,acc,nres,more);
	}
	tidyup(nat,nres,20);
	cycle(nat,acc,nres,more);
	nout = nres+1;
	acc[0] = acc[nout] = acc[nout+1] = 0.0;
	res[0] = res[nout] = res[nout+1] = 'B';
	for (i=1; i<=nout; i++) {
		vcopy(nat[i],nat+i+nout+1);
		acc[i+nout+1] = acc[i];
		res[i+nout+1] = res[i];
	}
	out = fopen("circle.out","w");
	for (i=0; i<=nres+1 && mode=='N'; i++) {
		fprintf(out,"ATOM%7d  CA  %3s A%4d    %8.3f%8.3f%8.3f  0.00%6.2f\n",
               		i, aa_code13(res[i]), i, nat[i].x, nat[i].y, nat[i].z, acc[i]);
	}
	for (i=nres+1; i>=0 && mode=='R'; i--) {
		fprintf(out,"ATOM%7d  CA  %3s A%4d    %8.3f%8.3f%8.3f  0.00%6.2f\n",
               		i, aa_code13(res[i]), nres-i+1, nat[i].x, nat[i].y, nat[i].z, acc[i]);
	}
	fprintf(out,"TER\n");
	fclose(out);
	out = fopen("test.out","w");
	for (i=1; i<=nres && mode=='N'; i++) {
		fprintf(out,"ATOM%7d  CA  %3s Z%4d    %8.3f%8.3f%8.3f  0.00%6.2f\n",
               		i, aa_code13(res[i]), i, nat[i].x, nat[i].y, nat[i].z, acc[i]);
	}
	for (i=nres; i>=1 && mode=='R'; i--) {
		fprintf(out,"ATOM%7d  CA  %3s Z%4d    %8.3f%8.3f%8.3f  0.00%6.2f\n",
               		i, aa_code13(res[i]), nres-i+1, nat[i].x, nat[i].y, nat[i].z, acc[i]);
	}
	fprintf(out,"TER\n");
	n = 0;
	for (j=2; j<nres-1; j++) acs[n++] = 2.0*acc[j-1]+acc[j]+acc[j+1]+acc[j+2]*2.0;
	sort(0,acs,0,at,n,1);
	m = 1;
	for (i=0; i<=nres+1; i++) in[i] = 0;
	for (i=n-1; i>=0; i--) { int skip, rati = at[i]+2;
		if (in[rati]<0) continue;
		if (acs[at[i]]>0.5) continue;
		skip = 0;
		for (j=rati-2; j<=rati+3; j++) {
			if (in[j]>0) skip = 1;
		}
		if (skip) continue;
		for (j=rati-2; j<=rati+3; j++) in[j] = -1;
		in[rati] = m++;
	}
	nout = 0;
	for (j=2; j<nres-1; j++)
	{ char	c; int k;
		if (in[j]<1) continue;
		if (in[j]>NOUT) continue;
		c = (char)('A'+nout);
		if (nout>25) c = (char)('a'+nout-25);
		for (i=j+1, k=1; i<=j+nres && mode=='N'; i++, k++) {
			fprintf(out,"ATOM%7d  CA  %3s %c%4d    %8.3f%8.3f%8.3f  0.00%6.2f\n",
                		i, aa_code13(res[i]), c, k, nat[i].x, nat[i].y, nat[i].z, acc[i]);
		}
		for (i=j+nres, k=1; i>=j+1 && mode=='R'; i--, k++) {
			fprintf(out,"ATOM%7d  CA  %3s %c%4d    %8.3f%8.3f%8.3f  0.00%6.2f\n",
                		i, aa_code13(res[i]), c, k, nat[i].x, nat[i].y, nat[i].z, acc[i]);
		}
		fprintf(out,"TER\n");
		nout++;
	}
	Pi(nout) NL
	fclose(out);
}

typedef struct { int a,b; float s; char c; } Pairs;

cycle (nat,acc,n,more) Vec *nat; float *acc; int n, more; {
int	i,j,k,m, iend,jend;
float	acci, accj;
float	cost;
Pairs	*link;
int	*look;
Vec	*old, *cut, guide, cent, mid;
char	*sec;
	Pi(n) NL
	m = n+more+2;
	old = (Vec*)alloca(sizeof(Vec)*m);
	cut = (Vec*)alloca(sizeof(Vec)*m);
	link = (Pairs*)alloca(sizeof(Pairs)*(m*m/4));
	look = (int*)alloca(sizeof(int)*(m*m/4));
	sec = (char*)alloca(sizeof(char)*m);
	k = 0;
	acci = 0.0;
	for (i=1; i<n/4; i++) {
		accj = 0.0;
		for (j=n; j>3*n/4; j--)
		{ float d, s, span;
		  int	lost = i-1+n-j;
			d = vdif(nat[i],nat[j]);
			s = 100.0*(float)lost + d*d + 100.0*(acci+accj);
			accj += acc[j];
			span = sqrt(lost+more)*10.0;
			if (span<d) continue;
			link[k].a = i;
			link[k].b = j;
			link[k].s = s;
			k++;
		}
		acci += acc[i];
	}
	sort(link,0,0,look,k,1);
	iend = link[look[k-1]].a;
	jend = link[look[k-1]].b;
	cost = link[look[k-1]].s;
	Pi(k) Pi(iend) Pi(jend) Pr(cost) NL
	for (i=1; i<=n; i++) sec[i] = '-';
	k = n+more;
	m = k/2+1;
	j = 1;
	for (i=m; i<=k; j++, i++) vcopy(nat[i],old+j);
	for (i=1;  i<m; j++, i++) vcopy(nat[i],old+j);
	for (i=n+1; i<=k; i++) { vinit(old+i); old[i].x = 999.9; }
	for (i=1; i<iend; i++) vinit(nat+i);
	for (i=jend+1; i<=k; i++) vinit(nat+i);
	j = 1;
	for (i=m; i<=k; j++, i++) vcopy(nat[i],cut+j);
	for (i=1;  i<m; j++, i++) vcopy(nat[i],cut+j);
	refine(sec,cut,old,k,200);
	j = 1;
	for (i=m; i<=k; j++, i++) vcopy(cut[i],nat+j);
	for (i=1;  i<m; j++, i++) vcopy(cut[i],nat+j);
	vcopy(nat[k], nat+0);
}

flips (nat,acc,n,more) Vec *nat; float *acc; int n, more; {
int	i,m, a,b,c,d,e,f, spare, need=n/6;
int	best[6]; float ad,eb,cf, max=0.0;
float	dsum, asum, gaps, lost, cost, cut=12.0;
Vec	mid, *new;
	new = (Vec*)alloca(sizeof(Vec)*n);
	for (a=1;  a<n; a++) { for (b=a+1;  b<a+5; b++) {
		for (c=b+10;  c<n; c++) { for (d=c+1;  d<c+5; d++) {
			for (e=d+10; e<n; e++) { for (f=e+1; f<e+5; f++) {
			  int	bc = c-b, de = e-d, fa = a+n-f, low;
				if (bc<need || de<need || fa<need) continue;
			  	ad = vdif(nat[a],nat[d]);
				eb = vdif(nat[e],nat[b]);
				cf = vdif(nat[c],nat[f]);
				if (ad>cut || eb>cut || cf>cut) continue;
				dsum = ad+eb+cf;
				low = (int)(dsum/6.0);
				spare = (b-a)+(d-c)+(f-e)-3;
				if (spare < low ) continue;
				if (spare > 9 ) continue;
				lost = cbrt((float)((b-a)*(d-c)*(f-e)));
				gaps = cbrt((float)(bc*de*fa));
				asum = acc[a]+acc[b]+acc[c]+acc[d]+acc[e]+acc[f];
				asum = 6.0-asum;
				cost = (gaps*asum)/(dsum*lost);
				if (cost > max) {
					best[0]=a; best[1]=b;
					best[2]=c; best[3]=d;
					best[4]=e; best[5]=f;
					max = cost;
				}
			} }
		} }
	} }
	a=best[0];b= best[1];
	c=best[2]; d=best[3];
	e=best[4]; f=best[5];
	spare = (b-a)+(d-c)+(f-e)-3;
        ad = vdif(nat[a],nat[d]);
        eb = vdif(nat[e],nat[b]);
        cf = vdif(nat[c],nat[f]);
	dsum = ad+eb+cf;
Pi(a) Pi(b) Pi(c) Pi(d) Pi(e) Pi(f) NL
Pr(ad) Pr(eb) Pr(cf) NL
Pr(dsum) Pi(spare) NL
	m = 0;
	for (i=0; i<=a; i++) { vcopy(nat[i],new+m); m++; }
	vave(nat[a],nat[d],&mid);
	if (spare==6 || (spare==5 && ad>eb) || (spare==4 && ad>eb && ad>cf)) {
		vave(nat[a],mid,new+m); m++;
		if (spare==9 || (spare==8 && ad>eb) || (spare==7 && ad>eb && ad>cf)) {
			vcopy(mid,new+m); m++;
		}
		vave(nat[d],mid,new+m); m++;
	} else {
		vcopy(mid,new+m); m++;
	}
	for (i=d; i<=e; i++) { vcopy(nat[i],new+m); m++; }
	vave(nat[e],nat[b],&mid);
	if (spare==6 || (spare==5 && eb>ad) || (spare==4 && eb>ad && eb>cf)) {
		vave(nat[e],mid,new+m); m++;
		if (spare==9 || (spare==8 && eb>ad) || (spare==7 && eb>ad && eb>cf)) {
			vcopy(mid,new+m); m++;
		}
		vave(nat[b],mid,new+m); m++;
	} else {
		vcopy(mid,new+m); m++;
	}
	for (i=b; i<=c; i++) { vcopy(nat[i],new+m); m++; }
	vave(nat[c],nat[f],&mid);
	if (spare==6 || (spare==5 && cf>ad) || (spare==4 && cf>ad && cf>eb)) {
		vave(nat[c],mid,new+m); m++;
		if (spare==9 || (spare==8 && cf>ad) || (spare==7 && cf>ad && cf>eb)) {
			vcopy(mid,new+m); m++;
		}
		vave(nat[f],mid,new+m); m++;
	} else {
		vcopy(mid,new+m); m++;
	}
	for(i=f; i<n; i++) vcopy(nat[i], new+i);
	for(i=0; i<n; i++) vcopy(new[i], nat+i);
}

refine (sec,ca,cb,len,m) char *sec; Vec *ca, *cb; int len, m;
{
int	i, j, k, l, n;
Vec	*work, *cab;
int	*real, *jump;
int	Nter, Cter;
Vec	cent, step, atom, cext, mid, aim;
float	extend = 3.0;
int	pacsum, lenn, ntry = 50;
	real = (int*)alloca(sizeof(int)*(len+2));
	jump = (int*)alloca(sizeof(int)*(len+2));
	work = (Vec*)alloca(sizeof(Vec)*(len+2)*2);
	cab = (Vec*)alloca(sizeof(Vec)*(len+2)*2);
	n = 0;
	lenn = len*2;
        vinit(&cent);
        for (j=1; j<=len; j++) { 
		if (vsqr(ca[j]) < 0.001) {
			real[j] = 0; 
		} else {
			real[j] = 1;
			vsum(ca[j], &cent);
			n++;
	 	}
		if (vdif(cb[j-1],cb[j])>4.5) jump[j-1]=jump[j]=1; else jump[j]=0;
	}
	vdiv(&cent,(float)n);
        for (i=1; i<=len; i++) {
                if (!real[i]) { Nter = i-1; break; }
	}
        for (i=len; i; i--) {
                if (!real[i]) { Cter = i+1; break; }
	}
	n = Cter-Nter;
	vave(ca[Cter],ca[Nter],&mid);
	vsub(mid,cent,&aim);
	vnorm(&aim); vmul(&aim,100.0);
	vadd(mid,aim,&cext);
	vsub(ca[Cter],ca[Nter],&step);
	vdiv(&step,(float)n);
	vcopy(ca[Nter],&atom);
        for (k=Nter+1; k<Cter; k++)
	{float	nn = 0.5*(float)n,
		x = fabs((float)(k-Nter)-nn),
		h = 1.0 + 3.0*sqrt(nn-x);
		vadd(atom,step,&atom);
		vsub(cext,atom,&aim);
		vnorm(&aim); vmul(&aim,extend);
		ca[k].x = atom.x + h*aim.x + 3.0*(drand48()-0.5);
		ca[k].y = atom.y + h*aim.y + 3.0*(drand48()-0.5);
		ca[k].z = atom.z + h*aim.z + 3.0*(drand48()-0.5);
	}
        for (i=0; i<m; i++) {
		for (j=Nter+1; j<Cter-1; j++) { int bump;
			bump = 0;
			for (k=1; k<=len; k++) {
				if (k>=j-1 && k<=j+1) continue;
				if (k>Nter && k<Cter) {
					if (vddif(ca[j],ca[k])<30.0) bump++;
				} else {
					if (vddif(ca[j],ca[k])<60.0) bump++;
				}
			}
			if (bump) continue;
			vsub(cent,ca[j],&aim); vmul(&aim,0.01); 
			vadd(ca[j],aim,ca+j);
		}
        	for (j=1; j<=len; j++) vcopy(ca[j],work+j);
        	for (j=2; j<=len; j++) {
			if ((j<Nter || j>Cter) && (jump[j-1] && jump[j])) continue;
			separate(ca,work,3.8,0.5,j-1,j);
		}
        	for (j=1; j<=len; j++) vcopy(work[j],ca+j);
		for (j=Nter+1; j<Cter-2; j++)
                { float hca = phand(ca[j-1],ca[j],ca[j+1],ca[j+2]),
                        hcb = phand(cb[j-1],cb[j],cb[j+1],cb[j+2]),
                        d = hca-hcb; d = d*d;
			if (jump[j-1]||jump[j]||jump[j+1]||jump[j+2]) continue;
                        if (d>100.0) reset(ca,cb,j);
		}
        	for (j=Nter; j<=Cter; j++) {
        		for (k=1; k<=len; k++) vcopy(ca[k],work+k);
			for (k=j-4; k<=j+4; k++) { float d = vdif(cb[j],cb[k]);
				if (d<7.0) separate(ca,work,d,0.5,j,k);
			}
        		for (k=1; k<=len; k++) vcopy(work[k],ca+k);
		}
        	for (j=1; j<len-1; j++) {
        		for (k=1; k<=len; k++) vcopy(ca[k],work+k);
			for (k=j+2; k<=len; k++) {
				if (vddif(ca[j],ca[k])<20.0) separate(ca,work,4.5,0.5,j,k);
			}
        		for (k=1; k<=len; k++) vcopy(work[k],ca+k);
		}
	}
	for (i=0; i<50; i++) {
        	for (j=1; j<len-1; j++) {
        		for (k=1; k<=len; k++) vcopy(ca[k],work+k);
			for (k=j+2; k<=len; k++) {
				if (vddif(ca[j],ca[k])<20.0) separate(ca,work,4.5,0.5,j,k);
			}
        		for (k=1; k<=len; k++) vcopy(work[k],ca+k);
		}
        	for (j=1; j<=len; j++) vcopy(ca[j],work+j);
        	for (j=2; j<=len; j++) {
			if ((j<Nter || j>Cter) && (jump[j-1] && jump[j])) continue;
			separate(ca,work,3.8,0.5,j-1,j);
		}
		separate(ca,work,3.8,0.5,1,len);
        	for (j=1; j<=len; j++) vcopy(work[j],ca+j);
	}
}

tidyup (ca,len,m) Vec *ca; int len, m;
{
int	i, j, k, l, n;
Vec	*work;
	work = (Vec*)alloca(sizeof(Vec)*(len+2)*2);
	for (i=0; i<m; i++) {
        	for (j=1; j<len-1; j++) {
        		for (k=1; k<=len; k++) vcopy(ca[k],work+k);
			for (k=j+2; k<=len; k++) {
				if (vddif(ca[j],ca[k])<20.0) separate(ca,work,4.5,0.5,j,k);
			}
        		for (k=1; k<=len; k++) vcopy(work[k],ca+k);
		}
	}
	for (i=0; i<m; i++) {
        	for (j=1; j<=len; j++) vcopy(ca[j],work+j);
        	for (j=1; j< len; j++) separate(ca,work,3.8,0.5,j,j+1);
        	for (j=1; j<=len; j++) vcopy(work[j],ca+j);
	}
}

reset (new,old,n) Vec *new, *old; int n;
{
Vec	a,b,c,d, p;
Vec	e,f,g,h, q;
Mat	bpc, fqg;
	vcopy(new[n-1],&a); vcopy(new[n-0],&b); vcopy(new[n+1],&c); vcopy(new[n+2],&d);
	vcopy(old[n-1],&e); vcopy(old[n-0],&f); vcopy(old[n+1],&g); vcopy(old[n+2],&h);
	vave(b,c,&p); vave(f,g,&q);
	setframe(a,p,d,&bpc);
	setframe(e,q,h,&fqg);
	vsub(a,p,&a); vsub(d,p,&d);
	vsub(b,p,&b); vsub(c,p,&c);
	VmulM(&bpc,a,&a); VmulM(&bpc,d,&d);
	VmulM(&bpc,b,&b); VmulM(&bpc,c,&c);
	vsub(e,q,&e); vsub(h,q,&h);
	vsub(f,q,&f); vsub(g,q,&g);
	VmulM(&fqg,e,&e); VmulM(&fqg,h,&h);
	VmulM(&fqg,f,&f); VmulM(&fqg,g,&g);
	vave(e,a,&a); vave(h,d,&d);
	vcopy(f,&b); vcopy(g,&c);
	MmulV(&bpc,a,&a); MmulV(&bpc,d,&d);
	MmulV(&bpc,b,&b); MmulV(&bpc,c,&c);
	vadd(a,p,&a); vadd(d,p,&d);
	vadd(b,p,&b); vadd(c,p,&c);
	vcopy(a,new+n-1); vcopy(d,new+n+2);
	vcopy(b,new+n+0); vcopy(c,new+n+1);
}

separate (old,new,dist,wt,m,n)
Vec    *old, *new;
float   dist, wt;
int     m, n;
{
        Vec     disp;
        float   gap, shift;
        if (wt>1.0) wt=1.0;
        gap = vdif(old[m],old[n]);
        if (gap < 0.1) {
                new[m].x += drand48();
                new[m].y += drand48();
                new[m].z += drand48();
                new[n].x -= drand48();
                new[n].y -= drand48();
                new[n].z -= drand48();
                return;
        } else {
                vsub(old[n],old[m],&disp);
        }
        shift = wt*(gap-dist)/(1.0+gap*2.0);
        vmul(&disp,shift);
        vadd(new[m],disp,&(new[m]));
        vsub(new[n],disp,&(new[n]));
}

setframe (a, b, c, frame)
    Vec a, b, c;
    Mat *frame;
{
    int    i;
    Vec    x, y, z ;
        vsub(c,a,&x);
        vave(c,a,&c);
        vsub(c,b,&y);
        vprod(y,x,&z);
        vprod(z,x,&y);
        vnorm(&x);
        vnorm(&y);
        vnorm(&z);
        VtoM(x,y,z,frame);
}
