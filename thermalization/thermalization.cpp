#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>

using namespace std;

int main(){
    ofstream fout("thermalization.txt");
    int iter, niter = 5000;
    double step_size = 0.5;
    double x, backup_x, dx;
    double action_init, action_fin;
    double metropolis;
    double sum = 0;

    srand((unsigned)time(NULL));

    x = 100.0;
    int naccept = 0;

    for(iter=1; iter<=niter; iter++){
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
        sum += x;
        fout << iter << "        " << sum/iter <<  endl;
    }

    cout << naccept << endl;

    return 0;
}
