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
```python
# Add your code here
```

#### Matrix Inversion Input
```
[Add your input format here]
```

#### Matrix Inversion Output
```
[Add your output format here]
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
