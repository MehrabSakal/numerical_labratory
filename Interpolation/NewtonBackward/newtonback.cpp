#include <bits/stdc++.h>
#include<fstream>
using namespace std;

long long fact(int n) {
    long long f = 1;
    for (int i = 1; i <= n; i++)
        f *= i;
    return f;
}

double newtonBackward(vector<double>& x, vector<double>& y, double value) {
    int n = x.size();
    double h = x[1] - x[0];

    vector<vector<double>> diff(n, vector<double>(n, 0));


    for (int i = 0; i < n; i++)
        diff[i][0] = y[i];


    for (int j = 1; j < n; j++) {
        for (int i = n - 1; i >= j; i--) {
            diff[i][j] = diff[i][j - 1] - diff[i - 1][j - 1];
        }
    }

    cout << "\nBackward Difference Table:\n";
    for (int i = 0; i < n; i++) {
        cout << x[i] << "\t";
        for (int j = 0; j <= i; j++)
            cout << diff[i][j] << "\t";
        cout << endl;
    }

    double p = (value - x[n - 1]) / h;
    double result = diff[n - 1][0];

    double p_term = 1;
    for (int k = 1; k < n; k++) {
        p_term *= (p + (k - 1));
        result += (p_term * diff[n - 1][k]) / fact(k);
    }

    return result;
}

int main() {
      ifstream in("input.txt");
    ofstream out("output.txt");
    int n;
    cout << "Enter number of data points: ";
    cin >> n;

    vector<double> x(n), y(n);

    cout << "Enter x values (equally spaced):\n";
    for (int i = 0; i < n; i++)
        cin >> x[i];

    cout << "Enter corresponding y values:\n";
    for (int i = 0; i < n; i++)
        cin >> y[i];

    double x1, x2;
    cout << "Enter first value (x1): ";
    cin >> x1;

    cout << "Enter second value (x2): ";
    cin >> x2;

    double F_x1 = newtonBackward(x, y, x1);
    double F_x2 = newtonBackward(x, y, x2);

    cout << "\n============================================";
    cout << "\n       RESULTS USING NEWTON BACKWARD";
    cout << "\n============================================";
    cout << "\nF(" << x1 << ") = " << F_x1;
    cout << "\nF(" << x2 << ") = " << F_x2;
    cout << "\nStudents between " << x1 << " and " << x2 << " = "
         << (F_x2 - F_x1);

    return 0;
}
