#include <bits/stdc++.h>
using namespace std;

double func(double x){
    return x*x*x + 2*x*x + x + 1;
}

double dfunc(double x){
    return 3*x*x + 4*x + 1;
}

double ddfunc(double x){
    return 6*x + 4;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if(!fin){
        cout<<"Error: input.txt not found!"<<endl;
        return 0;
    }
    if(!fout){
        cout<<"Error: Can't open output.txt!"<<endl;
        return 0;
    }

    fout<<fixed<<setprecision(3);
    double upperLimit, lowerLimit, p;
    int n;
    fin>>upperLimit>>lowerLimit>>n>>p;

    double h = (upperLimit - lowerLimit) / n;

    double x[n], y[n][n];
    for(int i=0; i<n; i++){
        x[i] = lowerLimit + i*h;
    }
    for(int i=0; i<n; i++){
        y[i][0] = func(x[i]);
    }

    for(int j=1; j<n; j++){
        for(int i=0; i<n-j; i++){
            y[i][j] = y[i+1][j-1] - y[i][j-1];
        }
    }

    fout<<"\nThe difference table is : \n";
    for(int i=0; i<n; i++){
        fout<<x[i]<<"\t";
        for(int j=0; j<n-i; j++){
            fout<<y[i][j]<<"\t";
        }
        fout<<endl;
    }

    double d = x[1] - x[0];
    double u = (p - x[0]) / d;

    double y1 = ( y[0][1] + ((2*u - 1)/2)*y[0][2] + ((3*u*u - 6*u+ 2)/6)*y[0][3] + ((4*u*u*u - 18*u*u - 14*u - 6)/24)*y[0][4] ) / d;
    fout<<"\nThe value of f'(p) is = "<<y1<<endl;

    double y2 = ( y[0][2] + (u-1)*y[0][3] + ((12*u*u - 36*u - 14)/24)*y[0][4] ) / (d*d);
    fout<<"\nThe value of f''(p) is = "<<y2<<endl;

    double err1 = ( fabs(dfunc(p) - y1) / fabs(dfunc(p)) ) * 100;
    double err2 = ( fabs(ddfunc(p) - y2) / fabs(ddfunc(p)) ) * 100;

    fout<<"\nThe relative error of f'(p) = "<<err1<<"%\n";
    fout<<"\nThe relative error of f''(p) = "<<err2<<"%\n";

    fin.close();
    fout.close();
    return 0;
}