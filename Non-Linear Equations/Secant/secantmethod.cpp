#include <bits/stdc++.h>
using namespace std;

double f(double x, vector<double>& a, int degree) {
    double sum = 0;
    for (int i = 0; i <= degree; i++)
        sum += a[i] * pow(x, degree - i);
    return sum;
}

bool secant(double x0, double x1,
            vector<double>& a, int degree,
            double& root, int& iterations) {

    double tol = 1e-6;
    int maxIter = 100;
    iterations = 0;

    double f0 = f(x0, a, degree);
    double f1 = f(x1, a, degree);

    for (int i = 0; i < maxIter; i++) {
        if (fabs(f1 - f0) < 1e-12)
            return false;

        double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        double f2 = f(x2, a, degree);

        iterations++;

        if (fabs(f2) < tol) {
            root = x2;
            return true;
        }

        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f2;
    }
    return false;
}

int main() {
    int degree;
    cout << "Enter degree of polynomial: ";
    cin >> degree;

    vector<double> coeff(degree + 1);
    cout << "Enter coefficients (highest degree to constant):\n";
    for (int i = 0; i <= degree; i++)
        cin >> coeff[i];

    double xmin, xmax, step;
    cout << "Enter search range [xmin xmax]: ";
    cin >> xmin >> xmax;

    cout << "Enter step size: ";
    cin >> step;

    vector<double> roots;
    double tol = 1e-6;

    cout << "\nDetected roots:\n";

    for (double x = xmin; x < xmax; x += step) {
        double x1 = x;
        double x2 = x + step;

        double f1 = f(x1, coeff, degree);
        double f2 = f(x2, coeff, degree);


        if (fabs(f1) < tol) {
            roots.push_back(x1);
            cout << "Root = " << x1 << " (exact)\n";
        }


        else if (f1 * f2 < 0) {
            double root;
            int iter;
            if (secant(x1, x2, coeff, degree, root, iter)) {
                bool duplicate = false;
                for (double r : roots)
                    if (fabs(r - root) < 1e-4)
                        duplicate = true;

                if (!duplicate) {
                    roots.push_back(root);
                    cout << "Root = " << root
                         << "   Iterations = " << iter << endl;
                }
            }
        }
    }

    if (roots.empty())
        cout << "No roots found in given range.\n";

    return 0;
}
