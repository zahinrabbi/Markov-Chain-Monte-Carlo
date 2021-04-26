#include<iostream>
#include<cmath>
#include<complex>
#include<fstream>
#include<ctime>

using namespace std;

const int N = 20;
const int niter = 10000;
int ntau;
const double dtau = 0.05;

int box_muller(double& p, double& q){
    double r,s;
    r = (double)rand()/RAND_MAX;
    s = (double)rand()/RAND_MAX;
    p = sqrt(-2*log(r))*sin(2*M_PI*s);
    q = sqrt(-2*log(r))*cos(2*M_PI*s);
    return 0;
}

double action(const complex<double> phi[N][N]){
    complex<double> phi2[N][N];
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<N; k++){
                phi2[i][j] += phi[i][k]*phi[k][j];
            }
        }
    }
    double action = 0;
    for(int i=0; i<N; i++){
        action += real(0.5*phi2[i][i]);
    }
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            action += 0.25*real(phi2[i][j]*phi2[j][i]);
        }
    }
    return action*((double)N);
}

double hamiltonian(const complex<double> phi[N][N],const complex<double> P[N][N]){
    double hmltn;
    complex<double> store;
    hmltn = action(phi);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            store = P[i][j]*P[j][i];
            hmltn += 0.5*real(store);
        }
    }
    return hmltn;
}

int force(complex<double> (&delh)[N][N], const complex<double> phi[N][N]){
    complex<double> phi2[N][N], phi3[N][N];
    for(int i=0; i<N; i++){
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
            delh[i][j] = (phi[i][j] + phi3[i][j])*((double)N);
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
            box_muller(r1, r2);
            P[i][j] = r1/sqrt(2) + r2/sqrt(2)*complex<double>(0, 1);
            P[j][i] = r1/sqrt(2) - r2/sqrt(2)*complex<double>(0, 1);
        }
    }
    for(int i=0; i<N; i++){
        box_muller(r1, r2);
        P[i][i] = complex<double>(r1, 0.0);
    }
    hmltn_i = hamiltonian(phi, P);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            phi[i][j] += P[i][j]*0.5*dtau;
        }
    }
    for(int step=1; step<ntau; step++){
        force(delh, phi);
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
                P[i][j] = P[i][j] - delh[i][j]*dtau;
                phi[i][j] += P[i][j]*dtau;
            }
        }
    }
    force(delh, phi);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            P[i][j] = P[i][j] - delh[i][j]*dtau;
            phi[i][j] += P[i][j]*0.5*dtau;
        }
    }
    hmltn_f = hamiltonian(phi,P);
    return 0;
}

double jackknife(int iter, double sample[]){
    double sum = 0;
    double average;
    for(int i=0; i<iter; i++){
        sum += sample[i];
    }
    average = sum/((double)iter);
    double bin_variance, bin_variance_squared, jacked_sum = 0;
    double jacked_avg[iter];
    double temp;
    double k;
    for(int w=1; w<=(iter/2); w++){
        k = ((double)iter)/((double)w);
        bin_variance_squared = 0;
        for(int m=0; m<(int)k; m++){
            jacked_sum = sum;
            for(int n=0; n<w; n++){
                jacked_sum -= sample[(m*w)+n];
            }
            jacked_avg[m] = jacked_sum/((double)(iter-w));
            temp = jacked_avg[m] - average;
            bin_variance_squared = temp * temp;
        }
        bin_variance = sqrt(((double)(iter-1)*bin_variance_squared)/((double)iter));
    }
    return bin_variance;
}

int main(){
    ofstream fout("ntau_vs_error.txt");
    complex<double> phi[N][N];
    complex<double> backup_phi[N][N];
    double hmltn_i, hmltn_f, metropolis;
    double sample[niter];
    srand((unsigned)time(NULL));
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            phi[i][j] = complex<double>(0.0, 0.0);
        }
    }
    int naccept = 0;
    double sum = 0;
    int ntau;
    for(int i=1; i<=20; i++){
        for(int iter=0; iter<niter; iter++){
            for(int i=0; i<N; i++){
                for(int j=0; j<N; j++){
                    backup_phi[i][j] = phi[i][j];
                }
            }
            Molecular_Dynamics(phi, hmltn_i, hmltn_f, i);
            metropolis = (double)rand()/RAND_MAX;
            if(exp(hmltn_i-hmltn_f) > metropolis){
                naccept=naccept+1;
            }
            else{
                for(int i=0; i<N; i++){
                    for(int j=0; j<N; j++){
                        phi[i][j] = backup_phi[i][j];
                    }
                }
            }
            sum += action(phi);
            sample[iter] = sum;
            //cout << iter*ntau << "        " << sum/((double)(iter+1)*N*N) << "        " << ((double)naccept)/((double)(iter+1)) << "        "
            //<< jackknife(iter+1, sample)/((double)(iter+1)*N*N) << endl;
            if(iter == niter-1){
                fout << i << "        " << sum/((double)(iter+1)*N*N*i) << "        " << jackknife(iter+1, sample)/((double)(iter+1)*N*N)
                << endl;
            }
        }
    }
    fout.close();
    return 0;
}
