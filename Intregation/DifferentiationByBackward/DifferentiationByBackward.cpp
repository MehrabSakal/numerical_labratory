#include <bits/stdc++.h>
using namespace std;

double func(double x){
    return ( pow(sin(x),5) + 4*pow(sin(x),4) + 1 );
}

double dfunc(double x){
    return ( 5*pow(sin(x),4)*cos(x) + 16*pow(sin(x),3)*cos(x) );
}

double ddfunc(double x){
    return ( -5*pow(sin(x),5) + 20*pow(sin(x),3)*pow(cos(x),2) - 16*pow(sin(x),4) + 48*pow(sin(x),2)*pow(cos(x),2) );
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
        for(int i=n-1; i>=j; i--){
            y[i][j] = y[i][j-1] - y[i-1][j-1];
        }
    }

    fout<<"\nThe difference table is : \n";
    for(int i=0; i<n; i++){
        fout<<x[i]<<"\t";
        for(int j=0; j<=i; j++){
            fout<<y[i][j]<<"\t";
        }
        fout<<endl;
    }

    double d = x[1] - x[0];
    double v = (p - x[n-1]) / d;

    double y1 = ( y[n-1][1] + ((2*v + 1)/2)*y[n-1][2] + ((3*v*v + 6*v+ 2)/6)*y[n-1][3] + ((4*v*v*v + 18*v*v + 22*v + 6)/24)*y[n-1][4] ) / d;
    fout<<"\nThe value of f'(p) is = "<<y1<<endl;

    double y2 = ( y[n-1][2] + (v+1)*y[n-1][3] + ((12*v*v + 36*v + 22)/24)*y[n-1][4] ) / (d*d);
    fout<<"\nThe value of f''(p) is = "<<y2<<endl;

    double err1 = ( abs(dfunc(p) - y1) / abs(dfunc(p)) ) * 100;
    double err2 = ( abs(ddfunc(p) - y2) / abs(ddfunc(p)) ) * 100;
    fout<<"\nThe relative error of f'(p) = "<<err1<<"%\n";
    fout<<"\nThe relative error of f''(p) = "<<err2<<"%\n";

    fin.close();
    fout.close();
    return 0;
}