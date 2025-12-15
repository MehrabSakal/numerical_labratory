#include<bits/stdc++.h>
using namespace std;

double E = 1e-4;

double f(vector<double> cof, double x, int n)
{
    int deg = n;
    double sum = 0;
    for(int i = 0; i <= n; i++)
    {
        sum += cof[i] * pow(x, deg);
        deg--;
    }
    return sum;
}

void Bisection(int n, vector<double> cof, double lower, double upper, ofstream &out)
{
    double mid = lower;
    int iteration = 0;

    if(fabs(f(cof, mid, n)) < 1e-12)
    {
        out << "Root: " << mid << "\nIteration: " << iteration << "\n";
        return;
    }

    while(fabs(upper - lower) > E)
    {
        iteration++;
        mid = lower - f(cof, lower, n) *
              ((upper - lower) / (f(cof, upper, n) - f(cof, lower, n)));

        if(fabs(f(cof, mid, n)) < 1e-12)
        {
            out << "Root: " << mid << "\nIteration: " << iteration << "\n";
            return;
        }

        if(f(cof, lower, n) * f(cof, mid, n) < 0)
            upper = mid;
        else
            lower = mid;
    }

    out << "Approximate root: " << mid
        << "\nIteration: " << iteration << "\n";
}

int main()
{
    ifstream in("input.txt");
    ofstream out("output.txt");

    if(!in.is_open() || !out.is_open())
    {
        cerr << "File error\n";
        return 1;
    }

    int n;
    in >> n;

    vector<double> cof(n + 1);
    for(int i = 0; i <= n; i++)
        in >> cof[i];

    double xmax = abs(sqrt(pow((cof[1] / cof[0]), 2) - 2 * (cof[2] / cof[0])));
    double start = -xmax;
    double end = xmax;
    double step = 0.5;

    while(start <= end)
    {
        double y1 = start;
        double y2 = start + step;

        if(f(cof, y1, n) * f(cof, y2, n) < 0)
        {
            Bisection(n, cof, y1, y2, out);
            out << "\n";
        }
        start += step;
    }

    in.close();
    out.close();
}
