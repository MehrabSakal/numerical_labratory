#include<bits/stdc++.h>
using namespace std;

double u_cal(double u, int n){
    double temp = u;
    for(int i = 1; i < n; i++){
        temp *= (u - i);
    }
    return temp;
}

double fact(int n){
    double temp = 1;
    for(int i = 1; i <= n; i++){
        temp *= i;
    }
    return temp;
}

int main(){
    ifstream in("input.txt");
    ofstream out("output.txt");

    out << fixed << setprecision(3);

    int n;
    in >> n;

    vector<double> x(n);
    vector<vector<double>> y(n, vector<double>(n));

    for(int i = 0; i < n; i++){
        in >> x[i];
    }
    for(int i = 0; i < n; i++){
        in >> y[i][0];
    }

    for(int i = 1; i < n; i++){
        for(int j = 0; j < n - i; j++){
            y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
        }
    }

    for(int i = 0; i < n; i++){
        out << x[i] << " ";
        for(int j = 0; j < n; j++){
            out << y[i][j] << " ";
        }
        out << endl;
    }

    double value;
    in >> value;

    double h = x[1] - x[0];
    double u = (value - x[0]) / h;

    double sum = y[0][0];

    for(int i = 1; i < n; i++){
        sum += (u_cal(u, i) * y[0][i]) / fact(i);
    }

    out <<"\nThe result is: "<<sum<<endl;

    return 0;
}
