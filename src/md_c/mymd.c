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

typedef float  real;
typedef struct { real x, y, z; } rvec;
// Lee el fichero de entrada, Input.dat,tambien psodicion, velocidad, numero de atomos
void read_input(char* inputfile, rvec* pos, rvec* vel, int natoms) {
	FILE* fp = fopen (inputfile,"r");//(fichero,"opciones")
	if (fp == NULL) {perror("Error in reading file\n");exit(1);}
	int i;
	for (i = 0; i< natoms; i++) {
		float tpos[DIM],tvel[DIM];
		int res = fscanf(fp, "%f %f %f %f %f %f",&tpos[XX],&tpos[YY],&tpos[ZZ],&tvel[XX],&tvel[YY],&tvel[ZZ]);
		if (res != 6) {printf("Cannot parse input file %s\n",inputfile);exit(1);}
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
	if (fp == NULL) {perror("Cannot write output file");exit(1);}
	int i;
	for (i = 0; i< natoms; i++) {
		fprintf(fp,"%f %f %f %f %f %f\n",pos[i].x, pos[i].y, pos[i].z, vel[i].x, vel[i].y, vel[i].z);
	}
}

void first_VV(rvec* pos, rvec* vel, rvec* force, int natoms, real mass, real dt) {
    int i;
    real imass = 1.0/mass;
    for (i = 0; i < natoms; i++) {
            real coeff = 0.5*dt*dt*imass;
            pos[i].x += dt*vel[i].x+coeff*force[i].x;
            pos[i].y += dt*vel[i].y+coeff*force[i].y;
            pos[i].z += dt*vel[i].z+coeff*force[i].z;
            
            real coeffvel = 0.5*dt*imass;
            vel[i].x += coeffvel*force[i].x;
            vel[i].y += coeffvel*force[i].y;
            vel[i].z += coeffvel*force[i].z;
        }
}

void second_VV(rvec* vel, rvec* force, int natoms, real mass, real dt,float *Ekin) {
    int i;
    real imass = 1.0/mass;
    real tmpE = 0.0;
    //Velocity update, second part//
    for (i = 0; i < natoms; i++) {
            real coeffvel = 0.5*dt*imass;
            vel[i].x += coeffvel*force[i].x;
            vel[i].y += coeffvel*force[i].y;
            vel[i].z += coeffvel*force[i].z;
            tmpE += vel[i].x*vel[i].x + vel[i].y*vel[i].y + vel[i].z*vel[i].z;
        }
    *Ekin = tmpE*0.5*mass;
}

real A,B;

void nonbonded_compute(const rvec* pos, rvec* force, int natoms, real *Epot) {
    *Epot = 0.0;
    real rcut2 = 12.0*12.0;
    int i,j;
    for (i = 0; i<natoms; i++) {
        for (j = i+1; j<natoms; j++) {
            rvec rij;
	    //Calcular la componente del radio en cada eje
            rij.x = pos [i].x - pos[j].x;
            rij.y = pos [i].y - pos[j].y;
            rij.z = pos [i].z - pos[j].z;
            real r2 = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
	    //Para aligerar el calculo se calcula para los atomos en un radio menor que rcut2//
            if (r2<rcut2) {
                real r = sqrt(r2);
                real r_1 = 1.0/r;
                real r_2 = r_1*r_1;
                real r_6 = r_2*r_2*r_2;
                real r_12 = r_6*r_6;
                
                real AmbTerm = (A*r_6 - B)*r_6;
                *Epot += AmbTerm;
                real force_r = (6.0*(A*r_12 + AmbTerm)*r_1 - AmbTerm)*r_1;
                
                rvec forceij;
                forceij.x = force_r*rij.x;
		        force[i].x += forceij.x;
		        force[j].x -= forceij.x;
		        forceij.y = force_r*rij.y;
		        force[i].y += forceij.y;
		        force[j].y -= forceij.y;
		        forceij.z = force_r*rij.z;
                force[i].z += forceij.z;
		        force[j].z -= forceij.z;
            }
        }
    }
    
}


real Kinetic_to_Temp(real Ekin, int natoms){
		//Calculamos la Temperatura para cata Ekin, ya que T depende de Ekin
		real Temp=Ekin*2.0/(3.0*BOLTZMAN*natoms);
		return Temp;
}

void berendsen_thermostat (rvec* vel, int natoms, real dt){
	int i;
	real lambda=sqrt(1.0+0.01*(300.0/dt-1.0));
	for (i = 0; i < natoms; i ++){
		vel[i].x*=lambda;
		vel[i].y*=lambda;
		vel[i].z*=lambda;
	}
}

void add_ext_forces (rvec* pos, int natoms, rvec* force){
	int i;
	real k=-1.0;
	for (i=0; i<natoms;i++){
		force[i].x+=k*pos[i].x;
		force[i].y+=k*pos[i].y;
		force[i].z+=k*pos[i].z;
	}
}


void analyze_gdr(rvec* pos, int natoms){
	real r2_inner_shell=5.0*5.0;
	FILE* fileout = fopen ("gdr.dat","a");
	if (fileout == NULL){
		perror("Cannot write output file");
		exit(1);
	}
	int i,j;
	for (i=0;i<natoms;i++){
		real r2=pos[i].x*pos[i].x+pos[i].y*pos[i].y+pos[i].z*pos[i].z;
		if (r2<r2_inner_shell){
			for (j=0;j<natoms;j++){
				rvec dist;
				dist.x=pos[i].x-pos[j].x;
				dist.y=pos[i].y-pos[j].y;
				dist.z=pos[i].z-pos[j].z;
				real d=sqrt(dist.x*dist.x+dist.y*dist.y+dist.z*dist.z);
				fprintf(fileout,"%f\n",d);
			}
		}
	}
	fclose(fileout);
}
int main(int argc, char* argv[]) {
	printf("Start\n");
	char inputfile[] = "input.dat";
	char outputfile[] = "output.dat";
	int natoms = atoi(argv[1]);
	int niter = atoi(argv[2]);
	rvec *pos = (rvec*) malloc(natoms*sizeof(rvec));
	rvec *vel = (rvec*) malloc(natoms*sizeof(rvec));
	rvec *force = (rvec*) malloc(natoms*sizeof(rvec));
	memset((void *)force, 0, natoms*sizeof(rvec));
	read_input(inputfile, pos, vel, natoms);
	//MD skeleton. You have all the coordinates loaded up and ready to be used.
	int n;
	real sigma = 3.4;
	real eps = 0.238;
	A = 4.0*eps*pow(sigma,12.0);
	B = 4.0*eps*pow(sigma,6.0);
	real dt = 1.0/TIMEFACTOR;
	real mass = 39.948;
	real imass = 1.0/mass;
	real Ekin = 0.0, Epot = 0.0;
	real Temp = 0.0;
	nonbonded_compute(pos,force,natoms, &Epot);
	printf("%d %f %f %f\n", 0, Epot, Ekin, Epot+Ekin);
	for (n = 0; n < niter; n++) {
	        first_VV(pos, vel, force, natoms, mass, dt);
	        memset((void *)force, 0, natoms*sizeof(rvec));
	        nonbonded_compute(pos,force,natoms, &Epot);
		add_ext_forces(pos,natoms,force);
	        //
		second_VV(vel,force,natoms,mass,dt,&Ekin);
		Temp=Kinetic_to_Temp(Ekin,natoms);
		berendsen_thermostat(vel,natoms,dt);
	        //real coeff = 0.5*dt*dt*imass;
		if (n%10==0) {
			printf("%d %f %f %f %f %d\n",n,Ekin,Epot,Ekin+Epot,Temp,0);
			analyze_gdr(pos,natoms);
		}
	}


//End
	write_output (outputfile,pos,vel,natoms);
	free(pos);
	free(vel);
	free(force);
	printf("%d %f %f %f\n", n, Epot, Ekin, Epot+Ekin,0);
	printf("Finish\n");
	return 0;
}