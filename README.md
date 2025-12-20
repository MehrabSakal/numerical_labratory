# Table of Contents

- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
  - [Matrix Inversion](#matrix-inversion)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)

- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)
  - [Newton Rapson Method](#newton-rapson-method)
    - [Theory](#newton-rapson-theory)
    - [Code](#newton-rapson-code)
    - [Input](#newton-rapson-input)
    - [Output](#newton-rapson-output)
  - [Secant Method](#secant-method)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)
  
- [Least Square Regression](#least-square-regression)
  - [Linear Regression Method](#linear-regression-method)
    - [Theory](#linear-regression-theory)
    - [Code](#linear-regression-code)
    - [Input](#linear-regression-input)
    - [Output](#linear-regression-output)
  - [Transcendental Regression Method](#transcendental-regression-method)
    - [Theory](#transcendental-regression-theory)
    - [Code](#transcendental-regression-code)
    - [Input](#transcendental-regression-input)
    - [Output](#transcendental-regression-output)
  - [Polynomial Regression Method](#polynomial-regression-method)
    - [Theory](#polynomial-regression-theory)
    - [Code](#polynomial-regression-code)
    - [Input](#polynomial-regression-input)
    - [Output](#polynomial-regression-output)

- [Interpolation Methods](#interpolation-methods)
  - [Newton Forward Interpolation](#newton-forward-interpolation)
    - [Theory](#newton-forward-theory)
    - [Code](#newton-forward-code)
    - [Input](#newton-forward-input)
    - [Output](#newton-forward-output)
  - [Newton Backward Interpolation](#newton-backward-interpolation)
    - [Theory](#newton-backward-theory)
    - [Code](#newton-backward-code)
    - [Input](#newton-backward-input)
    - [Output](#newton-backward-output)
  - [Newton Divided Difference Interpolation](#newton-divided-difference-interpolation)
    - [Theory](#newton-divided-difference-theory)
    - [Code](#newton-divided-difference-code)
    - [Input](#newton-divided-difference-input)
    - [Output](#newton-divided-difference-output)
- [Integration and Differentiation](#integration-and-differentiation)
  - [Simpson 1/3 Rule](#simpson-13-rule)
    - [Theory](#simpson-13-rule-theory)
    - [Code](#simpson-13-rule-code)
    - [Input](#simpson-13-rule-input)
    - [Output](#simpson-13-rule-output)
  - [Simpson 3/8 Rule](#simpson-38-rule)
    - [Theory](#simpson-38-rule-theory)
    - [Code](#simpson-38-rule-code)
    - [Input](#simpson-38-rule-input)
    - [Output](#simpson-38-rule-output)
  - [Differentiation by forward](#differentiation-by-forward)
    - [Theory](#differentiation-by-forward-theory)
    - [Code](#differentiation-by-forward-code)
    - [Input](#differentiation-by-forward-input)
    - [Output](#differentiation-by-forward-output)
  - [Differentiation by backward](#differentiation-by-backward)
    - [Theory](#differentiation-by-backward-theory)
    - [Code](#differentiation-by-backward-code)
    - [Input](#differentiation-by-backward-input)
    - [Output](#differentiation-by-backward-output)   
- [Ordinary Differential Equation](#ordinary-differential-equation)
  - [Runge Kutta Method](#runge-kutta-method)
    - [Theory](#runge-kutta-theory)
    - [Code](#runge-kutta-code)
    - [Input](#runge-kutta-input)
    - [Output](#runge-kutta-output)
---

### Solution of Linear Equations

### Gauss Elimination Method

#### Gauss Elimination Theory
[Add your theory content here]

#### Gauss Elimination Code
```cpp
#include<bits/stdc++.h>
using namespace std;

vector<double> gauss_elimination(vector<vector<double>>& mat, int n, ofstream& out){
    const double eps = 1e-9;

    for(int k = 0; k < n - 1; k++){
        int i_max = k;
        for(int i = k + 1; i < n; i++){
            if(fabs(mat[i][k]) > fabs(mat[i_max][k])){
                i_max = i;
            }
        }

        if(fabs(mat[i_max][k]) < eps){
            return {};
        }

        swap(mat[i_max], mat[k]);

        for(int i = k + 1; i < n; i++){
            double factor = mat[i][k] / mat[k][k];
            for(int j = k; j <= n; j++){
                mat[i][j] -= mat[k][j] * factor;
            }
        }
    }

    out << fixed << setprecision(3);
    for(int i = 0; i < n; i++){
        for(int j = 0; j <= n; j++){
            out << mat[i][j] << " ";
        }
        out <<endl;
    }
    out <<endl;

    vector<double> x(n);
    for(int i = n - 1; i >= 0; i--){
        if(fabs(mat[i][i]) < eps){
            return {};
        }
        x[i] = mat[i][n];
        for(int j = i + 1; j < n; j++){
            x[i] -= mat[i][j] * x[j];
        }
        x[i] /= mat[i][i];
    }

    return x;
}

int main(){
    ifstream in("input.txt");
    ofstream out("output.txt");

    int n;
    in >> n;

    vector<vector<double>> mat(n, vector<double>(n + 1));
    for(int i = 0; i < n; i++){
        for(int j = 0; j <= n; j++){
            in >> mat[i][j];
        }
    }

    vector<double> solution = gauss_elimination(mat, n, out);

    if(solution.empty()){
        out << "No unique solution";
        return 0;
    }

    out << fixed << setprecision(3);
    for(double v : solution){
        out <<"Root: "<< v <<endl;
    }

    return 0;
}

```

#### Gauss Elimination Input
```
5 
2 1 -1 3 2 9 
1 3 2 -1 1 8 
3 2 4 1 -2 20 
2 1 3 2 1 17 
1 -1 2 3 4 15
```

#### Gauss Elimination Output
```
3.000 2.000 4.000 1.000 -2.000 20.000 
0.000 2.333 0.667 -1.333 1.667 1.333 
0.000 0.000 -3.571 2.143 3.571 -4.143 
0.000 0.000 0.000 2.400 7.000 7.960 
0.000 0.000 0.000 0.000 -1.083 -1.283 

Root: 5.154
Root: -1.000
Root: 2.262
Root: -0.138
Root: 1.185

```

---

### Gauss Jordan Elimination Method

#### Gauss Jordan Theory
[Add your theory content here]

#### Gauss Jordan Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int gauss_jordan(vector<vector<double>>& mat, int n, vector<double>& x) {
    const double eps = 1e-12;
    vector<int> where(n, -1);

    for (int col = 0, row = 0; col < n && row < n; ++col) {
        int sel = row;
        for (int i = row; i < n; ++i) {
            if (fabs(mat[i][col]) > fabs(mat[sel][col])) sel = i;
        }
        if (fabs(mat[sel][col]) < eps) continue;
        swap(mat[sel], mat[row]);
        where[col] = row;

        for (int i = 0; i < n; ++i) {
            if (i != row) {
                double factor = mat[i][col] / mat[row][col];
                for (int j = col; j <= n; ++j)
                    mat[i][j] -= mat[row][j] * factor;
            }
        }
        ++row;
    }

    x.assign(n, 0);
    for (int i = 0; i < n; ++i) {
        if (where[i] != -1)
            x[i] = mat[where[i]][n] / mat[where[i]][i];
    }

    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < n; ++j)
            sum += mat[i][j] * x[j];
        if (fabs(sum - mat[i][n]) > eps) return 0; 
    }

    for (int i = 0; i < n; ++i)
        if (where[i] == -1) return 2;

    return 1;
}

int main() {
    ifstream in("input.txt");
    ofstream out("output.txt");

    int n;
    in >> n;

    vector<vector<double>> mat(n, vector<double>(n + 1));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= n; ++j)
            in >> mat[i][j];

    vector<double> solution;
    int res = gauss_jordan(mat, n, solution);

    if (res == 0) out << "No solution";
    else if (res == 2) out << "Infinitely many solutions";
    else {
        out << fixed << setprecision(3);
        for (double v : solution)
            out <<"Root: "<<v<< endl;
    }

    return 0;
}

```

#### Gauss Jordan Input
```
5 
2 1 -1 3 2 9 
1 3 2 -1 1 8 
3 2 4 1 -2 20 
2 1 3 2 1 17 
1 -1 2 3 4 15
```

#### Gauss Jordan Output
```
Root: 5.154
Root: -1.000
Root: 2.262
Root: -0.138
Root: 1.185

```

---

### LU Decomposition Method

#### LU Decomposition Theory
[Add your theory content here]

#### LU Decomposition Code
```cpp
#include<bits/stdc++.h>
using namespace std;
void LU()
{
    int n;
    cin >> n;
    vector<vector<double>> a(n, vector<double>(n));
    vector<vector<double>> l(n, vector<double>(n, 0));
    vector<vector<double>> u(n, vector<double>(n, 0));
    vector<double> b(n), y(n), x(n);

    int index = -1;

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            cin >> a[i][j];
        cin >> b[i];
    }

    for(int i = 0; i < n; i++)
    {
        for(int k = i; k < n; k++)
        {
            double sum = 0;
            for(int j = 0; j < i; j++)
                sum += l[i][j] * u[j][k];

            u[i][k] = a[i][k] - sum;
        }

        if(fabs(u[i][i]) < 1e-9)
            index = i;

        for(int k = i; k < n; k++)
        {
            if(k == i)
                l[i][i] = 1;
            else
            {
                double sum = 0;
                for(int j = 0; j < i; j++)
                    sum += l[k][j] * u[j][i];

                l[k][i] = (a[k][i] - sum) / u[i][i];
            }
        }
    }

    for(int i = 0; i < n; i++)
    {
        double sum = 0;
        for(int j = 0; j < i; j++)
            sum += l[i][j] * y[j];
        y[i] = b[i] - sum;
    }

    if(index != -1)
    {
        if(fabs(y[index]) < 1e-9)
            cout << "Infinite Solution" << endl;
        else
            cout << "No Solution" << endl;
        return;
    }

    cout << "Unique Solution exists:" << endl;

    for(int i = n - 1; i >= 0; i--)
    {
        double sum = 0;
        for(int j = i + 1; j < n; j++)
            sum += u[i][j] * x[j];
        x[i] = (y[i] - sum) / u[i][i];
    }

    for(int i=0;i<x.size();i++)
    {
        cout<<"x"<<i+1<<":"<<x[i]<<endl;
    }


}

int main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    while(1)
    {
        string choice;//if you want to solve equation the type yes unless no
        cin>>choice;
        if(choice=="yes")
         LU();
        else
         break;
    }

    return 0;
}

```
#### LU Decomposition Input
```
yes
5
2 1 -1 3 2 9
1 3 2 -1 1 8
3 2 4 1 -2 20
2 1 3 2 1 17
1 -1 2 3 4 15
yes
2
1 1 2
2 2 4
yes
2
1 1 2
2 2 5
no
```

#### LU Decomposition Output
```
Unique Solution exists:
x1:5.15385
x2:-1
x3:2.26154
x4:-0.138462
x5:1.18462
Infinite Solution
No Solution
```

---

### Matrix Inversion

#### Matrix Inversion Theory
[Add your theory content here]

#### Matrix Inversion Code
```cpp
#include <bits/stdc++.h>
#include<fstream>

using namespace std;


void getCofactor(const vector<vector<double>>& A, vector<vector<double>>& temp,
                 int p, int q, int n)
{
    int i = 0, j = 0;

    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {


            if (row != p && col != q) {
                temp[i][j++] = A[row][col];

                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}


double determinant(const vector<vector<double>>& A, int n)
{
    if (n == 1)
        return A[0][0];

    double det = 0;
    int sign = 1;
    vector<vector<double>> temp(n, vector<double>(n));

    for (int f = 0; f < n; f++) {
        getCofactor(A, temp, 0, f, n);
        det += sign * A[0][f] * determinant(temp, n - 1);
        sign = -sign;
    }
    return det;
}


void adjoint(const vector<vector<double>>& A, vector<vector<double>>& adj)
{
    int n = A.size();

    if (n == 1) {
        adj[0][0] = 1;
        return;
    }

    int sign = 1;
    vector<vector<double>> temp(n, vector<double>(n));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {


            getCofactor(A, temp, i, j, n);

            sign = ((i + j) % 2 == 0 ? 1 : -1);

            adj[j][i] = sign * determinant(temp, n - 1);
        }
    }
}


bool inverse(const vector<vector<double>>& A, vector<vector<double>>& inv)
{
    int n = A.size();
    double det = determinant(A, n);

    if (det == 0) {
        cout << "Matrix is singular — inverse does not exist.\n";
        return false;
    }

    vector<vector<double>> adj(n, vector<double>(n));
    adjoint(A, adj);


    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = adj[i][j] / det;

    return true;
}

int main()
{
    ifstream in("input.txt");
    ofstream out("output.txt");
    int n;
    cout << "Enter size of matrix n (n x n): ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<vector<double>> inv(n, vector<double>(n));

    cout << "\nEnter matrix elements:\n";
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            cin >> A[i][j];

    if (inverse(A, inv)) {
        cout << "\nInverse Matrix (by Adjoint Method):\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                cout << fixed << setprecision(6) << inv[i][j] << "  ";
            cout << endl;
        }
    }

    return 0;
}

```

#### Matrix Inversion Input
```
2
4 7
2 6
```

#### Matrix Inversion Output
```
Inverse Matrix (by Adjoint Method):
0.600000  -0.700000
-0.200000  0.400000
```

---

### Solution of Non-Linear Equations

### Bisection Method

#### Bisection Theory
[Add your theory content here]

#### Bisection Code
```cpp
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
        mid = (lower+upper)/2;

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

```

#### Bisection Input
```
4
1 0 -5 0 4
```

#### Bisection Output
```
Approximate root: -1.99999
Iteration: 13

Approximate root: -0.999985
Iteration: 13

Approximate root: 1.00001
Iteration: 13

Approximate root: 2.00001
Iteration: 13
```

---

### False Position Method

#### False Position Theory
[Add your theory content here]

#### False Position Code
```cpp
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

void FalsePosition(int n, vector<double> cof, double lower, double upper, ofstream &out)
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
            FalsePosition(n, cof, y1, y2, out);
            out << "\n";
        }
        start += step;
    }

    in.close();
    out.close();
}

```

#### False Position Input
```
4
1 0 -5 0 4
```

#### False Position Output
```
Root: -2
Iteration: 20

Approximate root: -1
Iteration: 3

Root: 1
Iteration: 7

Root: 2
Iteration: 30
```

---
### Newton Rapson Method

### Newton Rapson Theory

### Newton Rapson Code
```cpp
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

```
### Newton Rapson Input
```
5
```
### Newton Rapson Output
```
Root:-2.0000 Iteration:1
----------------------------
Root:-1.0000 Iteration:1
----------------------------
Root:1.0000 Iteration:1
----------------------------
Root:2.0000 Iteration:1
----------------------------

```
---

### Secant Method

### Secant Theory

### Secant Code
```cpp

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
```
### Secant Input
```
3
1 0 -1 -1
1 2
0.2


3
1 -6 11 -6
0 4
0.45
```
### Secant Output
```
Detected roots:
Root = 1.32472   Iterations = 4


Detected roots:
Root = 1   Iterations = 6
Root = 2   Iterations = 3
Root = 3   Iterations = 
```
---
### Least Square Regression

### Linear Regression Method

### Linear Regression Theory

### Linear Regression Code
```cpp
#include <bits/stdc++.h>
#include <fstream>
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

    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;

    for (int i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }


    double m = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    double c = (sumY - m * sumX) / n;

    cout << "\nLinear Regression Line:\n";
    cout << "y = " << m << "x + " << c << endl;


    double x_new;
    cout << "\nEnter value of x to estimate y: ";
    cin >> x_new;

    double y_est = m * x_new + c;

    cout << "Estimated y = " << y_est << endl;

    return 0;
}
```
### Linear Regression Input
```
5
1 2 3 4 5
3 5 7 9 11
6
```
### Linear Regression Output
```
Linear Regression Line:
y = 2x + 1
Estimated y = 13
```
---

### Transcendental Regression Method

### Transcendental Regression Theory

### Transcendental Regression Code
```cpp
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
```
### Transcendental Regression Input
```
5
0 1 2 3 4
2.0 2.7 3.6 4.9 6.7
2.5
```
### Transcendental Regression Output
```
Exponential Regression Model:
y = 1.99163 * e^(0.30139x)
Estimated y = 4.23097
```
---

### Polynomial Regression Method

### Polynomial Regression Theory

### Polynomial Regression Code
```cpp
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
```
### Polynomial Regression Input
```
5
0 1 2 3 4
1 4 9 16 25
2
2.5
```
### Polynomial Regression Output
```
Polynomial Coefficients:
a0 = 1
a1 = 2
a2 = 1

Polynomial Equation:
y = 1 + 2x + 1x^2
Estimated y = 12.25
```
---


### Interpolation Methods

### Newton Forward Interpolation

### Newton Forward Theory

### Newton Forward Code
```cpp
code here
```
### Newton Forward Input
```
input here
```
### Newton Forward Output
```
output here
```
---
### Newton Backward Interpolation

### Newton Backward Theory

### Newton Backward Code
```cpp
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
```
### Newton Backward Input
```
5
35 45 55 65 75
31 42 51 35 31
40
45
```
### Newton Backward Output
```
Backward Difference Table:
35      31
45      42      11
55      51      9       -2
65      35      -16     -25     -23
75      31      -4      12      37      60

Backward Difference Table:
35      31
45      42      11
55      51      9       -2
65      35      -16     -25     -23
75      31      -4      12      37      60

============================================
       RESULTS USING NEWTON BACKWARD
============================================
F(40) = 32.9688
F(45) = 42
Students between 40 and 45 = 9.03125
```
---

### Newton Divided Difference Interpolation

### Newton Divided Difference Theory

### Newton Divided Difference Code
```cpp
code here
```
### Newton Divided Difference Input
```
input here
```
### Newton Divided Difference Output
```
output here
```
---

### Integration and Differentiation

### Simpson 1/3 Rule

### Simpson 1/3 Rule Theory

### Simpson 1/3 Rule Code
```cpp
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
    for (int i = 0; i <= n; i++) cin >> a[i];

    double L, U;
    cout << "Enter lower limit: ";
    cin >> L;
    cout << "Enter upper limit: ";
    cin >> U;

    int intervals;
    cout << "Enter number of intervals (even): ";
    cin >> intervals;

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
        if (i % 2 == 0)
            simpson += 2 * y[i];
        else
            simpson += 4 * y[i];
    }
    simpson *= (h / 3.0);

    cout << "\nIntegral (Simpson 1/3 Rule) = " << simpson << endl;
    cout << "f'(p) = " << f1(a, p) << endl;
    cout << "f''(p) = " << f2(a, p) << endl;

    return 0;
}
```
### Simpson 1/3 Rule Input
```
2
1 2 4
0
4
4
1
```
### Simpson 1/3 Rule Output
```
===========================
f(x) = 1x^2 + 2x^1 + 4

Forward Difference Table:
           4            3            2            0            0
           7            5            2            0
          12            7            2
          19            9
          28

Integral (Simpson 1/3 Rule) = 53.3333
f'(p) = 4
f''(p) = 2
```
---

### Simpson 3/8 Rule

### Simpson 3/8 Rule Theory

### Simpson 3/8 Rule Code
```cpp

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
```
### Simpson 3/8 Rule Input
```
3
1 0 1 1
0
3
3
2
```
### Simpson 3/8 Rule Output
```
===========================
f(x) = 1x^3 + 1x^1 + 1

Forward Difference Table:
           1            2            6            6
           3            8           12
          11           20
          31

Integral (Simpson 3/8 Rule) = 27.75
f'(p) = 13
f''(p) = 12
```
---

### Differentiation by forward

### Differentiation by forward Theory

### Differentiation by forward Code
```cpp
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
```
### Differentiation by forward Input
```
5 0 6 0.5
```
### Differentiation by forward Output
```

The difference table is : 
0.000	1.000	2.801	6.250	3.472	0.000	-0.000	
0.833	3.801	9.051	9.722	3.472	-0.000	
1.667	12.852	18.773	13.194	3.472	
2.500	31.625	31.968	16.667	
3.333	63.593	48.634	
4.167	112.227	

The value of f'(p) is = 3.750

The value of f''(p) is = 7.000

The relative error of f'(p) = 0.000%

The relative error of f''(p) = 0.000%

```
---
### Differentiation by backward

### Differentiation by backward Theory

### Differentiation by backward Code
```cpp
code here
```
### Differentiation by backward Input
```
input here
```
### Differentiation by backward Output
```
output here
```
---

### Ordinary Differential Equation

### Runge Kutta Method

### Runge Kutta Theory

### Runge Kutta Code
```cpp
code here
```
### Runge Kutta Input
```
input here
```
### Runge Kutta Output
```
output here
```
---
