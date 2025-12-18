#include<bits/stdc++.h>
using namespace std;

double f(double x){
    return (x*x*x*x)-(5*x*x)+4;
}

double f_prime(double x){
    return (4*x*x*x)-10*x;
}

vector<double> guesses(int n){
    vector<double> a;
    for(int i = -n; i < n - 1; i++){
        double x = i;
        double y = i + 1;
        if(f(x) == 0){
            a.push_back(x);
        }
        if(f(x) * f(y) < 0){
            a.push_back(x);
        }
    }
    return a;
}

void newton_raphson(double a, ofstream &out){
    double c = a;
    int itr = 0;
    out << fixed << setprecision(4);
    while(itr != 1000){
        if(f_prime(a) == 0){
            out << "Derivative is zero, no further iteration possible\n";
            return;
        }
        c = a - (f(a) / f_prime(a));
        if(fabs(c - a) < 0.0001){
            break;
        }
        out << "c:" << a << " c(next):" << c << "\n";
        a = c;
        itr++;
    }
    out << "Root:" << c << " Iteration:" << itr + 1 << "\n";
    out << "----------------------------\n";
}

int main(){
    ifstream in("input.txt");
    ofstream out("output.txt");

    int n;
    in >> n;

    vector<double> choose = guesses(n);
    for(double x : choose){
        newton_raphson(x, out);
    }
    return 0;
}
