#include<iostream>
#include<fstream>
#include<iomanip>
#include<ctime>

using namespace std;

ofstream fout("histogram.txt");

void histogram(double data[], int n){
    double temp;
    for(int i=0; i<n-1; i++){
        for(int j=i+1; j<n; j++){
            if(data[j] < data[i]){
                temp = data[i];
                data[i] = data[j];
                data[j] = temp;
            }
        }
    }
    for(int k=0; k<n; k++){
        fout << fixed << setprecision(2) << data[k] << endl;
    }
    cout << endl;
}

int main(){

    double data[1000];

    srand((unsigned)time(NULL));

    for(int i=0; i<1000; i++){
        data[i] = (double)rand()/RAND_MAX;
    }

    /*for(int j=0; j<1000; j++){
        cout << data[j] << "    ";
    }*/

    cout << endl;
    histogram(data, 1000);

    return 0;
}
