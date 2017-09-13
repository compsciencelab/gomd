#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define DIM 3
#define XX  0   
#define YY  1  
#define ZZ  2
#define TIMEFACTOR 48.88821
#define BOLTZMAN 0.001987191

typedef float real;
//typedef double real;
typedef struct { real x, y, z; } rvec;


static inline real rvec_norm2(const rvec a) {
  return a.x*a.x+a.y*a.y+a.z*a.z;
}
static inline real rvec_prod(const rvec a,const rvec b) {
  return (a.x*b.x+a.y*b.y+a.z*b.z);
}
static inline void rvec_sub(rvec* a, const rvec b, const rvec c) {
  a->x=b.x-c.x; a->y=b.y-c.y; a->z=b.z-c.z;
}
static inline void rvec_inc(rvec* a, const rvec b) {
  a->x+=b.x; a->y+=b.y; a->z+=b.z;
}
static inline void rvec_dec(rvec* a, const rvec b) {
  a->x-=b.x; a->y-=b.y; a->z-=b.z;
}
static inline void rvec_smul(rvec* v1, const real a, const rvec v2) {
 v1->x = a*v2.x; v1->y = a*v2.y; v1->z = a*v2.z;
}

void read_input(char* inputfile, rvec* pos, rvec* vel, int natoms) {
	FILE* fp = fopen (inputfile,"r");
	if (fp == NULL) {perror("ljmd");exit(1);}
	int i;
	for (i = 0; i< natoms; i++) {
		float tpos[DIM],tvel[DIM];
		int res = fscanf(fp, "%f %f %f %f %f %f",&tpos[XX],&tpos[YY],&tpos[ZZ],&tvel[XX],&tvel[YY],&tvel[ZZ]);
		if (res != 6) {printf("Cannot load input file %s\n",inputfile);exit(1);}
		pos[i].x = (real) tpos[XX];
		pos[i].y = (real) tpos[YY];
		pos[i].z = (real) tpos[ZZ];
		vel[i].x = (real) tvel[XX];
		vel[i].y = (real) tvel[YY];
		vel[i].z = (real) tvel[ZZ];
	}	
}

void write_output(char* outputfile, rvec* pos, rvec* vel, int natoms) {
	FILE* fp = fopen (outputfile,"w");
        if (fp == NULL) {perror("ljmd");exit(1);}
	int i;
	for (i = 0; i< natoms; i++) {//it does not work if I change to double
		fprintf(fp,"%f %f %f %f %f %f\n",pos[i].x, pos[i].y, pos[i].z, vel[i].x, vel[i].y, vel[i].z);
	}

}

void first_VV(rvec* pos, rvec* vel, rvec* force, int natoms, real mass, real dt) {
	real imass = 1.0/mass;
	int i;
	for (i=0; i<natoms; i++) { 
		real coeff = 0.5*dt*dt*imass;
		pos[i].x +=  vel[i].x*dt + coeff*force[i].x;
		pos[i].y +=  vel[i].y*dt + coeff*force[i].y;
		pos[i].z +=  vel[i].z*dt + coeff*force[i].z;
		coeff =  0.5*dt*imass; 
		vel[i].x +=  coeff*force[i].x; 
		vel[i].y +=  coeff*force[i].y; 
		vel[i].z +=  coeff*force[i].z; 
	}	
}

void second_VV(rvec* pos, rvec* vel, rvec* force, int natoms, real mass, real dt, real* Ekin) {
	real imass = 1.0/mass;
	*Ekin = 0.0;
	int i,a;
	for (i=0; i<natoms; i++) {//this is crap and worse than C++ for efficency
		real coeff =  0.5*dt*imass;
		vel[i].x +=  coeff*force[i].x;
		vel[i].y +=  coeff*force[i].y;
		vel[i].z +=  coeff*force[i].z;
		*Ekin += rvec_prod(vel[i],vel[i])*mass;  
	}
	*Ekin *= 0.5;	
}

real ljtable[2];
void nonbonded_compute(const rvec* pos,  rvec* force, int natoms,  real *Epot) {
	*Epot = 0.0;
	int i,j;
	for (i = 0; i<natoms; i++) {
	rvec forcei = {0,0,0};
	  for (j = i+1; j<natoms; j++) {
		rvec rji;
		rvec_sub(&rji, pos[j],  pos[i]);//rji = pos[j] - pos[i]
		real dist2 = rvec_norm2(rji);
		real r = sqrt(dist2);
                real r_1 = 1.0/r;
                real r_2 = r_1*r_1;
                real r_6 = r_2*r_2*r_2;
                real r_12 = r_6*r_6;
		real vdwA = ljtable[0];//mimic reading lj parameters
                real vdwB = ljtable[1];
                real AmBterm = (vdwA * r_6 - vdwB)*r_6;
                *Epot += AmBterm;
                real force_r = (6.0 * (vdwA*r_12 + AmBterm) * r_1 - AmBterm )*r_1;

		rvec forceij;
		rvec_smul(&forceij, -force_r, rji);//forceij = -force_r*r_ji
		rvec_inc(&forcei,forceij);//forcei += forceij;
		rvec_dec(&force[j],forceij);//force[j] -= forceij;//store on force[j]
	  }
	rvec_inc(&force[i],forcei);//store on force[i]
	} 
}


void periodic_boundary(rvec* pos, int natoms, real lx, real ly, real lz, rvec O) {
   int i;
   rvec delta;
   real ilx = 1.0/lx;
   real ily = 1.0/ly;
   real ilz = 1.0/lz;
   for (i = 0; i < natoms; i++) {
     delta.x = pos[i].x - O.x;
     delta.y = pos[i].y - O.y;
     delta.z = pos[i].z - O.z;
     delta.x -= lx * rint( ilx * delta.x );
     delta.y -= ly * rint( ily * delta.y );
     delta.z -= lz * rint( ilz * delta.z );
     pos[i].x = delta.x + O.x;
     pos[i].y = delta.y + O.y;
     pos[i].z = delta.z + O.z;
   } 
}

rvec delta(const Vector& p1,const  Vector& p2) const {
  rvec delta;
  delta.x  = p1.x - p2.x;
  delta.y  = p1.y - p2.y;
  delta.z  = p1.z - p2.z;
  delta.x -= lx * rint( ilx * delta.x );
  delta.y -= ly * rint( ily * delta.y );
  delta.z -= lz * rint( ilz * delta.z );
  return delta;
}

real kinetic_to_temp(real Ekin, int nfree) {
  return  2.0*Ekin/(3.0*nfree*BOLTZMAN);
}

void berendsden_thermostat(rvec* vel, int natoms, real iT, real tT, real rel) {
  int i;
  real lambda = sqrt( 1.0 + rel*( tT/iT - 1.0));
  for (i = 0; i < natoms; i++) vel[j] *= lambda;
}

int main(int argc, char* argv[]) {
	printf("Start\n");
	char inputfile[] = "input.dat";
	char outputfile[] = "output.dat";
	int natoms = atoi(argv[1]);
	int niter = atoi(argv[2]);
	int io = (argc==3)?io=1:atoi(argv[3]);//inputoutput
	real dt = 1.0/TIMEFACTOR;
	real mass = 39.948;
	real sigma = 3.4;
	real eps = 0.238;
	ljtable[0] = 4.0*eps*pow(sigma,12.0);
	ljtable[1] = 4.0*eps*pow(sigma,6.0); 
	rvec *pos = (rvec*) malloc(natoms*sizeof(rvec));
	rvec *vel = (rvec*) malloc(natoms*sizeof(rvec));
	rvec *force = (rvec*) malloc(natoms*sizeof(rvec));
	memset((void *)force, 0, natoms*sizeof(rvec));

	if (io) read_input(inputfile, pos, vel, natoms);

	real Ekin = 0.0, Epot = 0.0;
	nonbonded_compute(pos,force,natoms, &Epot);
	printf("%d %f %f %f\n", 0, Epot, Ekin, Epot+Ekin);
        real T;
	int n;
	for (n = 0; n < niter; n++) {
		first_VV(pos,vel,force,natoms,mass,dt);
		memset((void *)force, 0, natoms*sizeof(rvec));
		nonbonded_compute(pos,force,natoms, &Epot);
		second_VV(pos,vel,force,natoms,mass,dt, &Ekin);
                T = kinetic_to_temp(Ekin, natoms);
	        printf("%d %f %f %f %f\n", n, Epot, Ekin, Epot+Ekin, T);
	}
	if (io) write_output(outputfile, pos, vel, natoms);

	free(pos);
	free(vel);
	free(force);	
	printf("Finish\n");
	return 0;
}


