#include <bits/stdc++.h>
#include<fstream>
using namespace std;


double f(const vector<double>& a, double x) {
    double sum = 0;
    int n = a.size() - 1;
    for (int i = 0; i <= n; i++)
        sum += a[i] * pow(x, n - i);
    return sum;
}


double f1(const vector<double>& a, double x) {
    double sum = 0;
    int n = a.size() - 1;
    for (int i = 0; i < n; i++)
        sum += a[i] * (n - i) * pow(x, n - i - 1);
    return sum;
}


double f2(const vector<double>& a, double x) {
    double sum = 0;
    int n = a.size() - 1;
    for (int i = 0; i < n - 1; i++)
        sum += a[i] * (n - i) * (n - i - 1) * pow(x, n - i - 2);
    return sum;
}


void printPolynomial(const vector<double>& a) {
    int n = a.size() - 1;
    cout << "f(x) = ";
    for (int i = 0; i <= n; i++) {
        if (a[i] == 0) continue;

        cout << a[i];
        if (n - i > 0) cout << "x^" << (n - i);
        if (i != n) cout << " + ";
    }
    cout << endl;
}


int main() {
    ifstream in("input.txt");
    ofstream out("output.txt");
    int n;
    cout << "Enter degree: ";
    cin >> n;

    vector<double> a(n + 1);
    cout << "Enter coefficients (a_n to a_0): ";
    for (int i = 0; i <= n; i++)
        cin >> a[i];

    double L, U;
    cout << "Enter lower limit: ";
    cin >> L;
    cout << "Enter upper limit: ";
    cin >> U;

    int intervals;
    cout << "Enter number of intervals (multiple of 3): ";
    cin >> intervals;

    if (intervals % 3 != 0) {
        cout << "Error: Number of intervals must be a multiple of 3 for Simpson's 3/8 Rule.\n";
        return 0;
    }

    double p;
    cout << "Enter value of p: ";
    cin >> p;

    cout << "\n===========================\n";
    printPolynomial(a);


    vector<double> x(intervals + 1), y(intervals + 1);
    double h = (U - L) / intervals;

    for (int i = 0; i <= intervals; i++) {
        x[i] = L + i * h;
        y[i] = f(a, x[i]);
    }


    cout << "\nForward Difference Table:\n";
    vector<vector<double>> diff(intervals + 1, vector<double>(intervals + 1, 0));

    for (int i = 0; i <= intervals; i++)
        diff[i][0] = y[i];

    for (int j = 1; j <= intervals; j++)
        for (int i = 0; i <= intervals - j; i++)
            diff[i][j] = diff[i + 1][j - 1] - diff[i][j - 1];

    for (int i = 0; i <= intervals; i++) {
        for (int j = 0; j <= intervals - i; j++)
            cout << setw(12) << diff[i][j] << " ";
        cout << endl;
    }

    double simpson = y[0] + y[intervals];

    for (int i = 1; i < intervals; i++) {
        if (i % 3 == 0)
            simpson += 2 * y[i];
        else
            simpson += 3 * y[i];
    }

    simpson *= (3.0 * h / 8.0);

    cout << "\nIntegral (Simpson 3/8 Rule) = " << simpson << endl;
    cout << "f'(p) = " << f1(a, p) << endl;
    cout << "f''(p) = " << f2(a, p) << endl;

    return 0;
}
