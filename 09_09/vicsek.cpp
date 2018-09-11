//To complile this code, type following command.
// gcc vicsek.c -o vicsek -lm && ./vicsek
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include "random.h"

using namespace std;

const long long int N=400;
int steps=100000, near_count;
int time_max = 100;
double x[N], y[N], vx[N], vy[N], theta[N], t[N], theta_avg=0.;
double L=10,vc=0.3, dt=1, sigma=0.15, r=1,eta=0.,lambda = 0.1;

//double unirand(double a, double b); // fn to return uniform random number in (a,b)
//double dist(double x1, double x2, double y1, double y2, double L); // fn to compute the distance between i th particle and j th particle at periodic boundary condition
//void pdist(double D[], double x[], double y[], int n, double L); // fn to compute the pairwise distance a periodic boundary condition
//void cal_mean_angle(double D[], double mean_ang[], double theta[], int N, double r); // fn to compute mea angle
//void filewrite(double x[], double y[], double vx[], double vy[], int iter, int N); // fn to write down the positions and the velocities

double unirand(double a, double b){
    return a + (b-a) * Uniform();
}

int uniint(int a, int b){
    return a + rand()%(b-a);
}


double dist(double x1, double x2, double y1, double y2){
    double dx, dy;
    dx = fabs(x1 - x2);
    if (dx > L-dx)
        dx = L-dx;
    dy = fabs(y1 - y2);
    if (dy > L-dy)
        dy = L-dy;
    return sqrt(dx*dx + dy*dy);
}

double bdist(int i,int j){
    return(dist(x[i], x[j], y[i], y[j]));
}

void update(){

    int i = uniint(0,N);
    int j = uniint(0,N);

    if (bdist(i,j)<= r){
        theta_avg = atan2(cos(theta[i])+cos(theta[j]),sin(theta[i])+sin(theta[j]));
        theta[i] = theta_avg+ rand_normal(0,sigma);
        theta[j] = theta_avg+ rand_normal(0,sigma);
    }
}

double order(){

    double mean_cos = 0.;
    double mean_sin = 0.;

    for(int i = 0;i<N;i++)
    {
        mean_cos += cos(theta[i]);
        mean_sin += sin(theta[i]);
    }

    return(sqrt(mean_cos*mean_cos+mean_sin*mean_sin)/N);

}

void evaluate()
{
    for(int i=0; i<N; i++){
        theta[i] = unirand(-M_PI, M_PI);
        x[i] = unirand(0.0, L);
        y[i] = unirand(0.0, L);
        vx[i] = vc * cos(theta[i]);
        vy[i] = vc * sin(theta[i]);
    }

    //
    // start vicsek model
    //
    for (int iter=0; iter<steps; iter++){

        // calculate mean angle
       update();


        // update positions and velocities
        for (int i=0; i<N; i++){
            vx[i] = vc * cos(theta[i]); vy[i] = vc * sin(theta[i]);
            x[i] += vx[i] * dt; y[i] += vy[i] * dt;
        }

        // periodic boundary condition
        for (int i=0; i<N; i++){
            x[i] = fmodl(x[i] + L, L);
            y[i] = fmodl(y[i] + L, L);
        }
    }


}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
int main(void){
    FILE *fp;
    fp = fopen("boltzmann_100.txt", "w");

    clock_t begin = clock();
	/*
	 * N: the number of particle
	 * steps: the number of run steps
	 * x[i], y[i]: the x,y-position of i th particle
	 * vx[i], vy[i]: the vx,vy-velosity of i th particle
	 * theta[i]: the angle of i th particle
	 * D[i][j]: pairwise distance between i th particle and j th particle
	 *
	 * L: system size
	 * vc: speeds of the particles
	 * dt: time interval
	 * eta: order
	 * r: interaction radius
	 *
	 */


	// make "data" directory to store positions and velocities data.
	mkdir("data", 0755);


    // seed
    init_genrand((unsigned)time(NULL));
    srand (time(NULL));

    double order_value = 0.;

    for(sigma=0;sigma<=2.;sigma+=0.1)
    {
        order_value = 0.;
        for(int i=0;i<100;i++)
        {
            evaluate();
            order_value += order();
        }

        fprintf(fp, "%lf\t%lf\n",sigma,order_value/100.);

    }

    fclose(fp);

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cout<<"Time: "<<elapsed_secs<<"\n";
    return 0;
}
