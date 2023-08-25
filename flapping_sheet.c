#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tag.h"
#include "tension.h"
//#include "view.h"

#define hl 0.02
#define hg 10*hl
#define rho_l 1.
#define rho_g 1.2e-3
#define U_liq 1.
#define Ug (Re*mu_g)/(rho_g*hl)
// #define Ul 4.84 //Dimotakis
#define Ul (sqrt(rho_l)*U_liq+sqrt(rho_g)*Ug)/(sqrt(rho_l)+sqrt(rho_g))
#define Re 300.
#define mu_g 1.63e-5
#define mu_l mu_g*58.82
#define stension 1.84
#define length 0.150
#define delta_g hl/2.

int MAXL = 12;
double vol_cut = 3.0518e-5;
scalar f0[];

u.n[left]  = dirichlet(y<=(-hl/2-hg) ? 0.05*Ug*tanh(-2*(y+(hl/2+hg))/delta_g) : (y>=(-hl/2-hg) && y<=-hl/2) ? Ug*tanh(-2*(y-(-hl/2-hg))/delta_g)*tanh(2*(y-(-hl/2))/delta_g) : (y>=-hl/2 && y<=hl/2) ? Ul : (y>=hl/2 && y<=(hl/2+hg)) ? Ug*tanh(-2*(y-hl/2)/delta_g)*tanh(2*(y-(hl/2+hg))/delta_g) : y>=(hl/2+hg) ? 0.05*Ug*tanh(2*(y-(hl/2+hg))/delta_g) : 0. );
u.t[left] = dirichlet(0);
#if dimension > 2
u.r[left] = dirichlet(0);
#endif
p[left] = neumann(0);
f[left] = f0[];

//since code is second order we can set backflow to zero (sponge method)
//can also coarsen mesh close to boundary
u.n[right] = u.n[] > 0 ? neumann(0) : 0;
//u.n[right] = neumann(0);
p[right] = dirichlet(0);

int main(){
	L0=1;
	origin (0, -L0/2, -L0/2);
	periodic(back);
	rho1 = rho_l;
	rho2 = rho_g;
	mu1 = mu_l;
	mu2 = mu_g;
	f.sigma = stension;
	TOLERANCE = 1e-4;
	init_grid(128);
	size (L0);
	run();
}

event init(t=0){
	if (!restore (file = "restart")) {

		refine (y >= -5*hl*L0 && y <= 5*hl*L0 && level < 10);
		restore (file = "init");
		scalar m1[];
		fraction(m1, hl/2. - y);
		scalar m2[];
		fraction(m2, hl/2. + y);
		//scalar m3[];
		//fraction (m3, sq(x-(length-0.08-R0)) + sq (z) - sq(R0));
		scalar m4[];
		fraction(m4, sq(hl/2.)-sq(x-length)-sq(y));

		foreach()
			f0[] += m1[]*m2[]*(x<=length);
			f0.refine = f0.prolongation = fraction_refine;
			restriction ({f0});

		foreach() {
			f[] += union (f0[],m4[]);
		}
		boundary({f,f0});
	}}

event logfile (i++) {
	if (i == 0)
		fprintf(ferr, "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
	fprintf(ferr, "%g\t %g\t %d\t %d\t %d\t %ld\t %g\t %g\n", t, dt, mgp.i, mgpf.i, mgu.i, grid->tn, perf.t, perf.speed);
}

event gfssnapshot (t=0.; t+= 1e-5; t<=1.) {
	char name[80];
	sprintf (name, "gfs-snapshot-%3.5f.gfs", t);
	scalar omega[];
 	vorticity (u, omega);
 	scalar omgx[];
 	foreach()
   		omgx[] = (u.z[0,1,0] - u.z[0,-1,0]
             		- u.y[0,0,1] + u.y[0,0,-1])/(2.*Delta);
   	boundary((scalar *){omgx});
	output_gfs (file = name, translate = true, list = {f, u, p, omega, omgx});
	//sprintf (name, "gfs-restart");
	//dump (name);
}

event snapshot (t = 0.; t += 0.0001; t<= 3.8) {
	char name[80];
	sprintf (name, "bv-snapshot-%g",t);
	dump (name);
	scalar pid[];
	foreach()
		pid[] = fmod(pid()*(npe() + 37), npe());
	boundary ({pid});
	sprintf (name, "restart");
	dump (name);
}

event adapt (i++) {
	scalar f1[];
	foreach()
		f1[] = f[];
	boundary ({f1});
	//adapt on the velocity field, try changing 0.1 to a larger value
	adapt_wavelet ((scalar *){f1}, (double[]){1e-3}, minlevel = 7, maxlevel = MAXL);
	unrefine (x > L0*0.8 && level > 5);
}

event droplets (t+=2e-6)
{
	scalar m[];
	foreach()
		m[] = f [] > 1e-3;
	int n = tag (m);
	double v[n];
	coord b[n];
	for (int j=0; j < n; j++)
		v[j] = b[j].x = b[j].y = b[j].z = 0.;
	foreach_leaf()
		if (m[] > 0) {
			int j = m[] - 1;
			v[j] += dv()*f[];
			coord p = {x,y,z};
			foreach_dimension()
				b[j].x += dv()*f[]*p.x;
		}
	#if _MPI
		MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#endif
	for (int j=0; j < n; j++)
		fprintf (fout, "%d %g %d %g %g %g\n", i, t, j, v[j], b[j].x/v[j], b[j].y/v[j]);
	fflush (fout);
}

event remove_drops ( i+=10 )
{
	scalar m[];
	foreach()
		m[] = f[] > 0.;
	int n = tag (m);
	double v[n];
	int remove_flag[n];
	for (int j = 0; j < n; j++) {
		v[j] = 0.;
		remove_flag[j] = 0;
	}
	foreach_leaf()
		if (m[] > 0) {
			int j = m[] - 1;
			v[j] += dv()*f[];
		}

	#if _MPI
		MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#endif

	for (int j = 0; j < n; j++)
		if (v[j] < vol_cut)
			remove_flag[j] = 1;

	foreach()
		if (m[] > 0) {
			int j = m[] - 1;
			if ( remove_flag[j] == 1 )
				f[] = 0.;
		}
}
