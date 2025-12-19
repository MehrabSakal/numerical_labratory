#include <bits/stdc++.h>

using namespace std;

double f(double x, double y) {
    return x * y + y;
}

int main() {
    ifstream in("input.txt");
    ofstream out("output.txt");

    if (!in.is_open()) {
        return 1;
    }

    double x0, y0, xn, h;

    if (in >> x0 >> y0 >> xn >> h) {
        int n = round((xn - x0) / h);
        double x = x0;
        double y = y0;

        out << fixed << setprecision(6);
        out << x << " " << y << endl;

        for (int i = 0; i < n; i++) {
            double k1 = h * f(x, y);
            double k2 = h * f(x + h / 2.0, y + k1 / 2.0);
            double k3 = h * f(x + h / 2.0, y + k2 / 2.0);
            double k4 = h * f(x + h, y + k3);

            y = y + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
            x = x + h;

            out << x << " " << y << endl;
        }
    }

    in.close();
    out.close();

    return 0;
}