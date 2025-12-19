#include <bits/stdc++.h>
#include<fstream>
using namespace std;


vector<double> gaussElimination(vector<vector<double>> A, vector<double> B) {
    int n = A.size();

    for (int i = 0; i < n; i++) {

        for (int k = i + 1; k < n; k++) {
            if (fabs(A[i][i]) < fabs(A[k][i])) {
                swap(A[i], A[k]);
                swap(B[i], B[k]);
            }
        }


        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = 0; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            B[k] -= factor * B[i];
        }
    }


    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = B[i];
        for (int j =0; j < n; j++) {
                if(i!=j)x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
    return x;
}


int main() {
      ifstream in("input.txt");
    ofstream out("output.txt");
    int n, degree;
    cout << "Enter number of data points: ";
    cin >> n;

    vector<double> x(n), y(n);

    cout << "Enter x values:\n";
    for (int i = 0; i < n; i++) cin >> x[i];

    cout << "Enter y values:\n";
    for (int i = 0; i < n; i++) cin >> y[i];

    cout << "Enter degree of polynomial: ";
    cin >> degree;

    int m = degree + 1;

    vector<vector<double>> A(m, vector<double>(m, 0));
    vector<double> B(m, 0);
    vector<double> powX(2 * degree + 1, 0);


    for (int i = 0; i < n; i++) {
        double v = 1;
        for (int j = 0; j <= 2 * degree; j++) {
            powX[j] += v;
            v *= x[i];
        }
    }


    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            A[i][j] = powX[i + j];
        }
    }


    for (int i = 0; i < n; i++) {
        double v = 1;
        for (int j = 0; j < m; j++) {
            B[j] += v * y[i];
            v *= x[i];
        }
    }


    vector<double> coeff = gaussElimination(A, B);

    cout << "\nPolynomial Coefficients:\n";
    for (int i = 0; i < m; i++) {
        cout << "a" << i << " = " << coeff[i] << endl;
    }


    cout << "\nPolynomial Equation:\n";
    cout << "y = ";

    for (int i = 0; i < m; i++) {

        if (i == 0) {
            cout << coeff[i];
        } else {
            if (coeff[i] >= 0)
                cout << " + " << coeff[i];
            else
                cout << " - " << fabs(coeff[i]);
        }


        if (i >= 1) {
            cout << "x";
            if (i > 1) cout << "^" << i;
        }
    }
    cout << endl;


    double x_new;
    cout << "\nEnter value of x to estimate y: ";
    cin >> x_new;

    double y_est = 0, p = 1;
    for (int i = 0; i < m; i++) {
        y_est += coeff[i] * p;
        p *= x_new;
    }

    cout << "Estimated y = " << y_est << endl;

    return 0;
}
