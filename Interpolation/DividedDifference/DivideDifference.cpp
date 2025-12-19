#include<bits/stdc++.h>
using namespace std;

double DivideDiff(const vector<double>& x, const vector<double>& y, double value, int n){
    vector<vector<double>> f(n, vector<double>(n));

    for(int i = 0; i < n; i++){
        f[i][0] = y[i];
    }

    for(int j = 1; j < n; j++){
        for(int i = 0; i < n - j; i++){
            f[i][j] = (f[i + 1][j - 1] - f[i][j - 1]) / (x[i + j] - x[i]);
        }
    }

    double result = f[0][0];

    for(int i = 1; i < n; i++){
        double term = f[0][i];
        for(int j = 0; j < i; j++){
            term *= (value - x[j]);
        }
        result += term;
    }

    return result;
}

int main(){
    ifstream in("input.txt");
    ofstream out("output.txt");

    int n;
    in >> n;

    vector<double> x(n), y(n);
    for(int i = 0; i < n; i++){
        in >> x[i];
    }
    for(int i = 0; i < n; i++){
        in >> y[i];
    }

    double value;
    in >> value;

    out << fixed << setprecision(6);
    out <<"The result is: "<< DivideDiff(x, y, value, n) << endl;

    return 0;
}
