#include<iostream>
#include<complex>
#include<ctime>

using namespace std;

int main(){
    int n;
    double a, b;

    cout << "Choose the dimension of the square matrix: " << endl;
    cin >> n;
    complex<double> matrix[n][n];

    srand((unsigned)time(NULL));

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            a = (double)rand()/RAND_MAX;
            b = (double)rand()/RAND_MAX;
            matrix[i][j] = complex<double>(a, b);
        }
    }

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout << matrix[i][j] << "     ";
        }
        cout << endl;
    }

    return 0;
}
