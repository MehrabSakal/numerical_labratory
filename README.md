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
  - [Differentiation](#differentiation)
    - [Theory](#differentiation-theory)
    - [Code](#differentiation-code)
    - [Input](#differentiation-input)
    - [Output](#differentiation-output)
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
```python
# Add your code here
```

#### Gauss Elimination Input
```
[Add your input format here]
```

#### Gauss Elimination Output
```
[Add your output format here]
```

---

### Gauss Jordan Elimination Method

#### Gauss Jordan Theory
[Add your theory content here]

#### Gauss Jordan Code
```python
# Add your code here
```

#### Gauss Jordan Input
```
[Add your input format here]
```

#### Gauss Jordan Output
```
[Add your output format here]
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
code here
```
### Newton Rapson Input
```
input here
```
### Newton Rapson Output
```
output here
```
---

### Secant Method

### Secant Theory

### Secant Code
```cpp
code here
```
### Secant Input
```
input here
```
### Secant Output
```
output here
```
---
### Least Square Regression

### Linear Regression Method

### Linear Regression Theory

### Linear Regression Code
```cpp
code here
```
### Linear Regression Input
```
input here
```
### Linear Regression Output
```
output here
```
---

### Transcendental Regression Method

### Transcendental Regression Theory

### Transcendental Regression Code
```cpp
code here
```
### Transcendental Regression Input
```
input here
```
### Transcendental Regression Output
```
output here
```
---

### Polynomial Regression Method

### Polynomial Regression Theory

### Polynomial Regression Code
```cpp
code here
```
### Polynomial Regression Input
```
input here
```
### Polynomial Regression Output
```
output here
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
code here
```
### Simpson 1/3 Rule Input
```
input here
```
### Simpson 1/3 Rule Output
```
output here
```
---

### Simpson 3/8 Rule

### Simpson 3/8 Rule Theory

### Simpson 3/8 Rule Code
```cpp
code here
```
### Simpson 3/8 Rule Input
```
input here
```
### Simpson 3/8 Rule Output
```
output here
```
---

### Differentiation

### Differentiation Theory

### Differentiation Code
```cpp
code here
```
### Differentiation Input
```
input here
```
### Differentiation Output
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
