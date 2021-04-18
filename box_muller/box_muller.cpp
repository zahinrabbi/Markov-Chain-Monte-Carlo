#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<ctime>

using namespace std;

ofstream fout("box_muller_1.txt");
ofstream file("box_muller_2.txt");

void histogram(double x[], double y[], int n){
    double temp;
    for(int i=0; i<n-1; i++){
        for(int j=i+1; j<n; j++){
            if(x[j] < x[i]){
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }
    }
    for(int i=0; i<n-1; i++){
        for(int j=i+1; j<n; j++){
            if(y[j] < y[i]){
                temp = y[i];
                y[i] = y[j];
                y[j] = temp;
            }
        }
    }
    for(int k=0; k<n; k++){
        fout << fixed << setprecision(2) << x[k] << endl;
    }

    for(int k=0; k<n; k++){
        file << fixed << setprecision(2) << y[k] << endl;
    }

    cout << endl;
}

int main(){

    double x[1000];
    double y[1000];
    double p, q;

    srand((unsigned)time(NULL));

    for(int i=0; i<1000; i++){
        p = (double)rand()/RAND_MAX;
        q = (double)rand()/RAND_MAX;
        x[i] = sqrt(-2*log(p))*sin(2*M_PI*q);
        y[i] = sqrt(-2*log(p))*cos(2*M_PI*q);
    }

    cout << endl;
    histogram(x, y, 1000);

    return 0;
}
