#include<iostream>
#include<cmath>
#include<complex>
#include<fstream>
#include<ctime>

using namespace std;

const int N = 10; //size of the matrices
const int ntau = 10; //final 'time step' in the leap frog method
const int niter = 7500; //sample size
const double dtau = 0.01; //step size in the leap frog method
double mu_squared;

int box_muller(double& p, double& q){
    double r,s;
    r = (double)rand()/RAND_MAX;
    s = (double)rand()/RAND_MAX;
    p = sqrt(-2*log(r))*sin(2*M_PI*s);
    q = sqrt(-2*log(r))*cos(2*M_PI*s);
    return 0;
}

double action(const complex<double> phi[N][N], double mu_squared){
    //takes initialized phi matrix as argument
    complex<double> phi2[N][N];
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<N; k++){
                phi2[i][j] += phi[i][k]*phi[k][j]; //we multiply the matrix with itself to create phi squared matrix
            }
        }
    }
    double action = 0;
    //our action has a mu_squared term multiplied to it
    for(int i=0; i<N; i++){
        action += real(0.5*mu_squared*phi2[i][i]); //calculates the trace of phi squared matrix
    }
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            action += 0.25*real(phi2[i][j]*phi2[j][i]); //calculates the trace of phi quartic matrix and adds trace of phi squared matrix
        }
    }
    return action*((double)N); //returns the total trace calculated
}

double phi2value(const complex<double> phi[N][N], double mu_squared){
    //this function calculates trace of just the phi^2 matrix
    complex<double> phi2[N][N];
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<N; k++){
                phi2[i][j] += phi[i][k]*phi[k][j];
            }
        }
    }
    double phi2value = 0;
    for(int i=0; i<N; i++){
        phi2value += real(0.5*mu_squared*phi2[i][i]); //trace of phi^2 matrix with the mu_squared term
    }
    return phi2value*((double)N);
}

double hamiltonian(const complex<double> phi[N][N],const complex<double> P[N][N]){
    double hmltn;
    hmltn = action(phi, mu_squared);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            hmltn += 0.5*real(P[i][j]*P[j][i]); //calculates the trace of the pseudo hamiltonian H = phi^2 + phi^4 + randomly genrated momentum^2
        }
    }
    return hmltn;
}

int force(complex<double> (&delh)[N][N], const complex<double> phi[N][N]){
    //referenced argument, the derivative of the hamiltonian to be calculated in this function
    complex<double> phi2[N][N], phi3[N][N];
    // the derivative is phi + phi^3 as dH/dx = dS/dx
    for(int i=0; i<N; i++){ //initializing the phi^2 and phi^3 matrices
        for(int j=0; j<N; j++){
            phi2[i][j] = complex<double>(0.0, 0.0);
            phi3[i][j] = complex<double>(0.0, 0.0);
        }
    }
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<N; k++){
                phi2[i][j] += phi[i][k]*phi[k][j];
            }
        }
    }
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<N; k++){
                phi3[i][j] += phi2[i][k]*phi[k][j];
            }
        }
    }
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            delh[i][j] = (phi[i][j] + phi3[i][j])*((double)N); //we add up the terms for the derivative matrix
        }
    }
    return 0;
}

int Molecular_Dynamics(complex<double> (&phi)[N][N], double& hmltn_i, double& hmltn_f, int ntau){
    complex<double> P[N][N];
    complex<double> delh[N][N];
    double r1,r2;
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            box_muller(r1, r2); //calling the box_muller function defined above to randomly generate the auxiliary momentum
            P[i][j] = r1/sqrt(2) + r2/sqrt(2)*complex<double>(0, 1);
            P[j][i] = r1/sqrt(2) - r2/sqrt(2)*complex<double>(0, 1);
        }
    }
    for(int i=0; i<N; i++){
        box_muller(r1, r2);
        P[i][i] = complex<double>(r1, 0.0);
    }
    hmltn_i = hamiltonian(phi, P); //initial hamiltonian before leap frog method
    //now we update our phi and momentum matrices using leap frog method which in turn updates our hamiltonian
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            phi[i][j] += P[i][j]*0.5*dtau; //first step of the leap frog method
        }
    }
    for(int step=1; step<ntau; step++){
        force(delh, phi); //function called to use the derivative of the hamiltonian below
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
                P[i][j] = P[i][j] - delh[i][j]*dtau; //intermediate steps of the leap frog method
                phi[i][j] += P[i][j]*dtau;
            }
        }
    }
    force(delh, phi);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            P[i][j] = P[i][j] - delh[i][j]*dtau; //final step of the leap frog method
            phi[i][j] += P[i][j]*0.5*dtau;
        }
    }
    hmltn_f = hamiltonian(phi,P); //final hamiltonian
    return 0;
}

double jackknife(int iter, double sample[]){
    //the jackknife function which takes the sample and it's size as arguments
    double sum = 0;
    double average;
    for(int i=0; i<iter; i++){
        sum += sample[i]; //sum all the values in the sample
    }
    average = sum/((double)iter);
    double bin_variance, bin_variance_squared, jacked_sum = 0;
    double jacked_avg[iter];
    double temp;
    double k;
    for(int w=1; w<=(iter/2); w++){
        //we split the sample into bins and vary the number of the size of the bin from 1 to half of the sample size
        k = ((double)iter)/((double)w); //number of bins created depending on the bin width
        bin_variance_squared = 0;
        for(int m=0; m<(int)k; m++){
            //we run a loop for the distinct bins created
            jacked_sum = sum;
            for(int n=0; n<w; n++){
                jacked_sum -= sample[(m*w)+n]; //we remove the bin from the sum
            }
            jacked_avg[m] = jacked_sum/((double)(iter-w)); //average of the sample with the bin removed
            temp = jacked_avg[m] - average;
            bin_variance_squared = temp * temp;
        }
        bin_variance = sqrt(((double)(iter-1)*bin_variance_squared)/((double)iter)); //formula for the jackknife error
    }
    return bin_variance;
}

int main(){
    ofstream fout("phi_model_mu_dependence.txt");
    complex<double> phi[N][N];
    complex<double> backup_phi[N][N];
    double hmltn_i, hmltn_f, metropolis;
    double sample[niter];
    srand((unsigned)time(NULL));
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            phi[i][j] = complex<double>(0.0, 0.0); //initializing the phi matrix
        }
    }
    int naccept = 0;
    double sum = 0;
    for(mu_squared=-10.0; mu_squared<=10.0; mu_squared=mu_squared+.25){
        for(int iter=0; iter<niter; iter++){
            //we start a loop with increasing sample size
            for(int i=0; i<N; i++){
                for(int j=0; j<N; j++){
                    backup_phi[i][j] = phi[i][j]; //initializing the matrix which will hold the same value of the phi matrix if it fails the metropolis test
                }
            }
            Molecular_Dynamics(phi, hmltn_i, hmltn_f, ntau);
            metropolis = (double)rand()/RAND_MAX;
            if(exp(hmltn_i-hmltn_f) > metropolis){
                naccept=naccept+1; //metropolis test, the new matrix is accepted
            }
            else{
                for(int i=0; i<N; i++){
                    for(int j=0; j<N; j++){
                        phi[i][j] = backup_phi[i][j]; //the new matrix is rejected
                    }
                }
            }
            sum += phi2value(phi, mu_squared); //trace of just the phi^2 term
            sample[iter] = sum; //the trace is stored in the sample array to be jackknifed. the size of the array increases with each iteration
            if(iter == niter-1){
                //we take the sample with maximum data points and calculate the error with jackknife function
                //we plot mu_squared vs the trace of the phi^2 term
                fout << mu_squared << "        " << sum/((double)(iter+1)*N*N*ntau) << "        " <<
                jackknife(iter+1, sample)/((double)(iter+1)*N*N) << endl;
            }
        }
    }
    return 0;
}

