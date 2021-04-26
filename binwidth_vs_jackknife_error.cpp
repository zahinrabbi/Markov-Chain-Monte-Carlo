#include<iostream>
#include<cmath>
#include<complex>
#include<fstream>
#include<ctime>

using namespace std;

int niter = 50000;
double step_size = 0.5;
double x, backup_x, dx;
double action_init, action_fin;
double metropolis;
int naccept;


double jackknife(int niter, int width, double sample[]){
    double sum = 0;
    double average;
    for(int i=0; i<niter; i++){
        sum += sample[i];
    }
    average = sum/((double)niter);
    double bin_variance, bin_variance_squared, jacked_sum = 0;
    double jacked_avg[niter];
    double temp;
    double k;
    //for(int w=1; w<=(iter/2); w++){
        k = ((double)niter)/((double)width);
        bin_variance_squared = 0;
        for(int m=0; m<(int)k; m++){
            jacked_sum = sum;
            for(int n=0; n<width; n++){
                jacked_sum -= sample[(m*width)+n];
            }
            jacked_avg[m] = jacked_sum/((double)(niter-width));
            temp = jacked_avg[m] - average;
            bin_variance_squared = temp * temp;
        }
        bin_variance = sqrt(((double)(niter-1)*bin_variance_squared)/((double)niter));
    //}
    return bin_variance;
}

int main(){
    ofstream fout("jackknife_error.txt");
    srand((unsigned)time(NULL));
    x = 0;
    double sum_xx = 0;
    int naccept = 0;
    double sample[niter];
    for(int width=1; width<=100; width++){
        for(int iter=0; iter<niter; iter++){
            backup_x = x;
            action_init = 0.5*x*x;
            dx = (double)rand()/RAND_MAX;
            dx = (dx - 0.5)*step_size*2;
            x = x + dx;
            action_fin = 0.5*x*x;
            metropolis = (double)rand()/RAND_MAX;
            if(exp(action_init-action_fin) > metropolis){
                naccept = naccept + 1;
            }
            else{
                x = backup_x;
            }
            sum_xx += x*x;
            sample[iter] = sum_xx/((double)iter+1);
        }
        fout << width << "        " << sample[niter-1]/((double)width) << "        " << jackknife(niter, width, sample)/((double)width) << endl;
    }
    return 0;
}
