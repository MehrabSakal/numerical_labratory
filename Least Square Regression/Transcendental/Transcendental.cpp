#include <bits/stdc++.h>
#include<fstream>

using namespace std;

int main() {
      ifstream in("input.txt");
    ofstream out("output.txt");
    int n;
    cout << "Enter number of data points: ";
    cin >> n;

    vector<double> x(n), y(n);

    cout << "Enter x values:\n";
    for (int i = 0; i < n; i++) cin >> x[i];

    cout << "Enter y values:\n";
    for (int i = 0; i < n; i++) cin >> y[i];


    vector<double> Y(n);
    for (int i = 0; i < n; i++) {
        if (y[i] <= 0) {
            cout << "Error: y values must be positive for exponential regression.\n";
            return 0;
        }
        Y[i] = log(y[i]);
    }


    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;

    for (int i = 0; i < n; i++) {
        sumX  += x[i];
        sumY  += Y[i];
        sumXY += x[i] * Y[i];
        sumX2 += x[i] * x[i];
    }


    double b = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    double C = (sumY - b * sumX) / n;


    double a = exp(C);

    cout << "\nExponential Regression Model:\n";
    cout << "y = " << a << " * e^(" << b << "x)\n";


    double x_new;
    cout << "\nEnter value of x to estimate y: ";
    cin >> x_new;

    double y_est = a * exp(b * x_new);

    cout << "Estimated y = " << y_est << endl;

    return 0;
}
