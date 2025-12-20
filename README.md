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

### 1. Introduction
Gauss Elimination is a direct algorithm used to solve systems of linear equations. It transforms a system into an "upper triangular" matrix (Row Echelon Form), making it easy to solve the variables one by one from the bottom up.

It consists of two main stages:
1.  **Forward Elimination:** Converting the matrix into an upper triangular form.
2.  **Back Substitution:** Solving for the unknowns starting from the last equation.



### 2. Mathematical Principle
Given a system of equations represented in matrix form $Ax = b$, we construct an **Augmented Matrix** $[A|b]$.

The goal is to use Elementary Row Operations to transform the matrix $A$ into an Upper Triangular Matrix $U$, where all elements below the main diagonal are zero:

$$
\begin{bmatrix}
a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33}
\end{bmatrix}
\rightarrow
\begin{bmatrix}
a_{11} & a_{12} & a_{13} \\
0 & u_{22} & u_{23} \\
0 & 0 & u_{33}
\end{bmatrix}
$$

### 3. The Algorithm Steps

**Stage 1: Forward Elimination**
For each column $j$ (from 1 to $n-1$):
1.  Find a **Pivot Row** (usually the current row $i$).
2.  For every row $k$ below the pivot row:
    * Calculate the multiplier: $m = \frac{A_{kj}}{A_{jj}}$
    * Subtract the pivot row scaled by $m$ from row $k$:
      $Row_k = Row_k - m \times Row_j$
    * This makes the element $A_{kj}$ zero.

**Stage 2: Back Substitution**
Once the matrix is upper triangular:
1.  Solve for the last variable $x_n$ directly:
    $x_n = \frac{b_n}{A_{nn}}$
2.  Move upwards to row $i = n-1, n-2, \dots, 1$:
    * Substitute the known values of $x$ found so far.
    * Solve for $x_i$:
      $$x_i = \frac{b_i - \sum_{j=i+1}^{n} A_{ij} x_j}{A_{ii}}$$

### 4. Complexity and Pitfalls
* **Time Complexity:** $O(n^3)$. It is computationally expensive for very large systems (e.g., thousands of equations).
* **Division by Zero:** If a pivot element ($A_{ii}$) is zero, the calculation fails.
* **Round-off Error:** In floating-point arithmetic, small errors can accumulate, leading to inaccurate results.

### 5. Advantages vs. Disadvantages

**Advantages**
* **General Purpose:** Works for any system of linear equations (provided a unique solution exists).
* **Exact:** Theoretically produces the exact answer (unlike iterative methods like Gauss-Seidel).

**Disadvantages**
* **Sensitivity:** Without "Partial Pivoting" (swapping rows to get the largest pivot), it is very sensitive to small errors.
* **Destructive:** It alters the original matrix $A$. If you need $A$ later, you must copy it first.

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

### 1. Introduction
The Gauss-Jordan Elimination method is a variation of Gauss Elimination. While Gauss Elimination stops when the matrix is in "Upper Triangular" form (Row Echelon Form), Gauss-Jordan continues the process to convert the matrix into a **Diagonal Matrix** (Reduced Row Echelon Form).

This means that after the operations are complete, the values of the unknowns can be read directly from the rightmost column without needing any "Back Substitution."



### 2. Mathematical Principle
Given an augmented matrix $[A|b]$, the goal is to transform the matrix $A$ into the **Identity Matrix** $I$.

$$
\begin{bmatrix}
a_{11} & a_{12} & a_{13} & | & b_1 \\
a_{21} & a_{22} & a_{23} & | & b_2 \\
a_{31} & a_{32} & a_{33} & | & b_3
\end{bmatrix}
\rightarrow
\begin{bmatrix}
1 & 0 & 0 & | & x_1 \\
0 & 1 & 0 & | & x_2 \\
0 & 0 & 1 & | & x_3
\end{bmatrix}
$$

Once the left side becomes the Identity Matrix, the right side automatically becomes the solution vector $x$.

### 3. The Algorithm Steps
1.  **Augment:** Form the augmented matrix $[A|b]$.
2.  **Iterate:** For each column $j$ (from 1 to $n$):
    * **Normalize Pivot:** Divide the entire row $j$ by the pivot element $A_{jj}$ so that the pivot becomes 1.
        $$Row_j = \frac{Row_j}{A_{jj}}$$
    * **Eliminate Other Rows:** For **every other** row $i$ (where $i \neq j$):
        * Subtract a multiple of the pivot row from row $i$ to make the element $A_{ij}$ zero.
        * $Row_i = Row_i - (A_{ij} \times Row_j)$

3.  **Result:** The final column contains the solution values $x_1, x_2, \dots, x_n$.

### 4. Complexity and Analysis
* **Time Complexity:** $O(n^3)$.
* **Comparison:** It requires approximately **50% more operations** than standard Gauss Elimination because you must eliminate entries *above* the diagonal as well as below it.
* **Matrix Inversion:** This method is the standard way to find the **Inverse of a Matrix**. If you augment $A$ with the Identity Matrix $[A|I]$ and apply Gauss-Jordan, the result is $[I|A^{-1}]$.

### 5. Advantages vs. Disadvantages

**Advantages**
* **Direct Solution:** No back substitution step is required.
* **Matrix Inversion:** It is excellent for calculating the inverse of a matrix.

**Disadvantages**
* **Inefficient for Systems:** For simply solving equations ($Ax=b$), it is slower than Gauss Elimination.
* **More Round-off Error:** Since it performs more arithmetic operations, there is a higher chance of accumulated floating-point errors.

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

## LU Decomposition Method

### LU Decomposition Theory

### 1. Introduction
LU Decomposition is a method used to solve systems of linear equations by factoring a matrix into two simpler triangular matrices. It transforms a square matrix $A$ into the product of a **Lower triangular matrix ($L$)** and an **Upper triangular matrix ($U$)**.

This method is particularly useful when you need to solve the same system of equations for multiple different result vectors, as the decomposition only needs to be done once.



### 2. Mathematical Principle
Given a square matrix $A$, we aim to find two matrices $L$ and $U$ such that:

$$A = L \cdot U$$

Where:
* **$L$ (Lower Triangular):** A matrix where all elements *above* the main diagonal are zero. (Usually, the diagonal elements are set to 1).
* **$U$ (Upper Triangular):** A matrix where all elements *below* the main diagonal are zero.

**The System of Equations:**
To solve $Ax = b$:
1.  Substitute $A$ with $LU$:  $LUx = b$
2.  Let $Ux = y$. Then the equation becomes: $Ly = b$

This breaks the hard problem into two easy steps:
1.  Solve $Ly = b$ for $y$ (Forward Substitution).
2.  Solve $Ux = y$ for $x$ (Backward Substitution).

### 3. The Algorithm Steps
There are different ways to decompose the matrix (Doolittle, Crout, Cholesky). The **Doolittle Algorithm** is the most common for general use:

1.  **Decomposition (Find L and U):**
    For each row $i$ and column $j$:
    * **Find U:** Calculate the upper triangle elements.
        $$U_{ij} = A_{ij} - \sum_{k=0}^{i-1} L_{ik} U_{kj}$$
    * **Find L:** Calculate the lower triangle elements.
        $$L_{ji} = \frac{1}{U_{ii}} (A_{ji} - \sum_{k=0}^{i-1} L_{jk} U_{ki})$$

2.  **Forward Substitution:**
    Solve $Ly = b$ to find $y$. Since $L$ is lower triangular, you start from the top ($y_1$) and work down.

3.  **Backward Substitution:**
    Solve $Ux = y$ to find $x$. Since $U$ is upper triangular, you start from the bottom ($x_n$) and work up.

### 4. Complexity and Conditions
* **Computational Cost:** $O(n^3)$ for the decomposition, but only $O(n^2)$ for the substitution steps.
* **Pivot Elements:** The method can fail if a pivot element (diagonal element $U_{ii}$) becomes zero. To prevent this, **Partial Pivoting** (swapping rows) is often added to the algorithm.

### 5. Advantages vs. Disadvantages

**Advantages**
* **Efficiency for Multiple Inputs:** If you are solving $Ax = b$ for many different vectors $b$, LU is much faster than Gaussian Elimination because you decompose $A$ only once.
* **Inversion:** It is a very efficient way to calculate the inverse of a matrix.

**Disadvantages**
* **Zero Pivots:** Basic LU decomposition fails if a diagonal element is zero (requires pivoting logic).
* **Memory:** Requires storing two new matrices ($L$ and $U$), though in-place algorithms exist to save space.

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

## Matrix Inversion

### Matrix Inversion Theory

Matrix Inversion is a numerical method used to solve a system of linear equations of the form:

**A·X = B**

where:
* **A** is the coefficient matrix
* **X** is the vector of unknowns
* **B** is the constant vector

If the inverse of matrix **A** exists, the solution is:

**X = A⁻¹·B**

## Condition for Inverse

The inverse of a square matrix exists only if its determinant is non-zero:

**det(A) ≠ 0**

If **det(A) = 0**, the matrix is singular and the system has either no solution or infinitely many solutions.

## Determinant

The determinant of a matrix is calculated using cofactor expansion:

in plain text:

det(A) = Σ(j=1 to n) (-1)^(1+j) · a₁ⱼ · det(M₁ⱼ)


## Cofactor

The cofactor of an element **aᵢⱼ** is defined as:

**Cᵢⱼ = (-1)^(i+j) · det(Mᵢⱼ)**

where **Mᵢⱼ** is the minor matrix obtained by deleting row i and column j.

## Adjoint Matrix

The adjoint of a matrix is the transpose of the cofactor matrix:

**adj(A) = [Cᵢⱼ]ᵀ**

## Inverse Matrix

The inverse of matrix **A** is given by:

**A⁻¹ = (1/det(A)) · adj(A)**

## Solving the System Using Inverse Matrix

Given:

**A·X = B**

Multiplying both sides by **A⁻¹**:

**X = A⁻¹·B**

This method provides a unique solution only when **det(A) ≠ 0**.

## Augmented Matrix Representation

The system is provided as an augmented matrix:

[A|B] = [a₁₁  a₁₂  ⋯  a₁ₙ | b₁]
        [a₂₁  a₂₂  ⋯  a₂ₙ | b₂]
        [ ⋮    ⋮   ⋱   ⋮  | ⋮ ]
        [aₙ₁  aₙ₂  ⋯  aₙₙ | bₙ]


* First **n** columns represent matrix **A**
* Last column represents vector **B**

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

## Solution of Non-Linear Equations

## Bisection Method

## Bisection Theory

### 1. Introduction
The Bisection Method is one of the fundamental numerical techniques used to find the root of a non-linear equation f(x) = 0. It is classified as a "bracketing method" because it requires two initial guesses that bracket (enclose) the root. It is widely used due to its simplicity and guaranteed convergence for continuous functions.

### 2. Mathematical Principle
The method is based on the **Intermediate Value Theorem**. 

The theorem states:
If a function f(x) is continuous on an interval [a, b], and the signs of f(a) and f(b) are opposite (meaning f(a) * f(b) < 0), then there is at least one root 'c' strictly between 'a' and 'b' such that f(c) = 0.

Essentially, if the function goes from positive to negative, it must cross the x-axis (zero) somewhere in between.

### 3. The Algorithm Steps
Given a function f(x) and an interval [a, b] such that f(a) and f(b) have opposite signs:

1.  **Calculate Midpoint:** Compute the middle point 'c' of the interval: 
    c = (a + b) / 2

2.  **Evaluate Function:** Calculate the value of f(c).

3.  **Check for Root:**
    * If f(c) is 0 (or very close to 0 within a specified tolerance), then 'c' is the root. Stop here.

4.  **Update Interval:**
    * If f(a) and f(c) have **opposite signs** (f(a) * f(c) < 0), the root lies in the left half. 
        -> Set b = c (The new interval is [a, c])
    * If f(a) and f(c) have the **same sign**, the root lies in the right half. 
        -> Set a = c (The new interval is [c, b])

5.  **Repeat:** Repeat steps 1-4 until the interval size (b - a) is smaller than your desired error tolerance.

### 4. Convergence Analysis
The Bisection Method is known for being reliable but relatively slow compared to other methods.

* **Convergence Type:** Linear Convergence.
* **Convergence Rate:** The interval size is halved in every iteration. This means the method gains approximately one decimal digit of accuracy for every 3.3 iterations.
* **Error Estimate:** The absolute error after 'n' iterations is guaranteed to be less than:
    Error <= (Original Interval Size) / 2^n

### 5. Stopping Criteria
In a computer program, we cannot iterate forever. The loop typically stops when one of these conditions is met:

1.  **Absolute Error:** The width of the interval |b - a| is less than a small number epsilon (e.g., 0.00001).
2.  **Function Value:** The value of |f(c)| is extremely close to zero.
3.  **Iteration Limit:** A safety counter (e.g., 100 iterations) is reached to prevent infinite loops if something goes wrong.

### 6. Advantages vs. Disadvantages

**Advantages**
* **Guaranteed Convergence:** As long as the function is continuous and the interval is valid, it will ALWAYS find a root.
* **Simple:** It requires no knowledge of derivatives (unlike Newton-Raphson).
* **Error Bound:** You can easily calculate exactly how many iterations are needed to reach a specific accuracy before you even start.

**Disadvantages**
* **Slow:** It converges slowly. If you need high precision, it takes many iterations.
* **No Complex Roots:** It cannot find imaginary or complex roots (e.g., x^2 + 1 = 0).
* **Multiple Roots:** If the interval contains multiple roots, it will only find one of them, and it's unpredictable which one it will find.

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

## False Position Method

#### False Position Theory

### 1. Introduction
The False Position Method (or *Regula Falsi*) is a root-finding algorithm that combines features of the Bisection Method and the Secant Method. Like Bisection, it is a "bracketing method" (requires two initial guesses with opposite signs), but instead of blindly taking the midpoint, it uses a linear interpolation to estimate the root position. This often results in faster convergence.

### 2. Mathematical Principle
The method relies on the concept that if a function is continuous, the root is likely closer to the point where the function value is smaller (closer to zero).

Instead of cutting the interval in half, we draw a straight line (a secant line) connecting the points `(a, f(a))` and `(b, f(b))`. The point where this line crosses the x-axis is our new estimate `c`.

**The Formula:**
The new estimate `c` is calculated using:

$$c = \frac{a \cdot f(b) - b \cdot f(a)}{f(b) - f(a)}$$



### 3. The Algorithm Steps
Given a function f(x) and an interval [a, b] such that f(a) and f(b) have opposite signs:

1.  **Calculate Estimate:** Compute the weighted average point 'c' using the formula above.

2.  **Evaluate Function:** Calculate the value of f(c).

3.  **Check for Root:**
    * If f(c) is 0 (or within tolerance), stop. 'c' is the root.

4.  **Update Interval:**
    * If f(a) and f(c) have **opposite signs**, the root is in the left side.
        -> Set b = c
    * If f(a) and f(c) have the **same sign**, the root is in the right side.
        -> Set a = c

5.  **Repeat:** Repeat steps 1-4 until the error is sufficiently small.

### 4. Convergence Analysis
* **Convergence Type:** Linear (often faster than Bisection but slower than Newton-Raphson).
* **Behavior:** It typically converges faster than Bisection because it uses the magnitude of the function values to "aim" for the root. However, for functions with significant curvature (very convex or concave), one end of the interval can get "stuck," causing convergence to slow down significantly.

### 5. Stopping Criteria
Common conditions to stop the loop:
1.  **Function Tolerance:** |f(c)| < epsilon
2.  **Step Tolerance:** |c_new - c_old| < epsilon (Note: We use this instead of |b-a| because the interval width might not go to zero in False Position).
3.  **Max Iterations:** To prevent infinite loops.

### 6. Advantages vs. Disadvantages

**Advantages**
* **Guaranteed Convergence:** Like Bisection, it will always find a root if the initial interval brackets one.
* **Faster than Bisection:** usually converges quicker because it considers the values of f(x), not just the signs.

**Disadvantages**
* **One-Sided Approach:** If the function is curved, one of the interval endpoints (a or b) might stay fixed for many iterations, causing the method to slow down.
* **Complex Formula:** The formula is slightly more expensive to calculate than the simple average used in Bisection.

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

The Newton-Raphson Method is a numerical technique used to determine the root of a non-linear equation of the form f(x) = 0. Unlike the Secant or False Position methods, the Newton-Raphson Method requires the derivative of the function to compute successive approximations. It uses a single initial guess that is close to the root.

### Initial Guess

Let the initial guess be:

x₀

### Iterative Formula

The next approximation of the root is computed using the formula:

xₙ₊₁ = xₙ - f(xₙ) / f'(xₙ)

Here, xₙ₊₁ is obtained by finding the x-intercept of the tangent line to the curve y = f(x) at the point (xₙ, f(xₙ)).

### Iteration Process

The process is repeated iteratively, updating the approximation as follows:

xₙ ← xₙ₊₁

until the difference between successive approximations becomes smaller than a prescribed tolerance:

|xₙ₊₁ - xₙ| < ε

or the absolute value of the function at the current approximation satisfies:

|f(xₙ₊₁)| < ε

### Algorithm Steps

1. Start with an initial guess x₀
2. Compute f(xₙ) and f'(xₙ)
3. Calculate the next approximation:

   xₙ₊₁ = xₙ - f(xₙ) / f'(xₙ)

4. Check convergence criteria:
   - If |xₙ₊₁ - xₙ| < ε or |f(xₙ₊₁)| < ε, stop
   - Otherwise, set xₙ = xₙ₊₁ and repeat from step 2

### Convergence Characteristics

The Newton-Raphson Method typically converges faster than the Bisection or Secant Methods, especially when the initial guess is sufficiently close to the root.

**However, it may fail to converge if:**
* The derivative is zero or very small at some iteration: f'(xₙ) ≈ 0
* The initial guess is far from the actual root
* The function has multiple roots or local extrema near the initial guess

### Advantages

* **Fast convergence** - Quadratic convergence rate when conditions are favorable
* **Efficient** - Requires fewer iterations compared to other methods

### Disadvantages

* **Requires derivative** - Must compute or know f'(x)
* **Sensitive to initial guess** - Poor initial guess may lead to divergence
* **May fail** - When f'(x) = 0 or is very small
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

The Secant Method is a numerical technique used to determine the root of a non-linear equation of the form f(x) = 0. Unlike the Bisection or False Position methods, the Secant Method does not require the function to change sign over an interval, and it uses two initial approximations that are close to the root.

### Initial Guesses

Let the initial guesses be:

x₀ and x₁

### Iterative Formula

The next approximation of the root is computed using the formula:

xₙ₊₁ = xₙ - f(xₙ) · (xₙ - xₙ₋₁) / (f(xₙ) - f(xₙ₋₁))

Here, xₙ₊₁ is the intersection of the secant line passing through the points (xₙ₋₁, f(xₙ₋₁)) and (xₙ, f(xₙ)) with the x-axis.

### Iteration Process

The process is repeated iteratively, updating the two previous approximations as follows:

xₙ₋₁ ← xₙ
xₙ ← xₙ₊₁

until the difference between successive approximations becomes smaller than a prescribed tolerance:

|xₙ₊₁ - xₙ| < ε

or the absolute value of the function at the current approximation satisfies:

|f(xₙ₊₁)| < ε

### Algorithm Steps

1. Start with two initial guesses x₀ and x₁
2. Compute f(xₙ₋₁) and f(xₙ)
3. Calculate the next approximation:

   xₙ₊₁ = xₙ - f(xₙ) · (xₙ - xₙ₋₁) / (f(xₙ) - f(xₙ₋₁))

4. Check convergence criteria:
   - If |xₙ₊₁ - xₙ| < ε or |f(xₙ₊₁)| < ε, stop
   - Otherwise, update: xₙ₋₁ = xₙ, xₙ = xₙ₊₁ and repeat from step 2

### Convergence Characteristics

The Secant Method typically converges faster than the Bisection Method, although it may fail to converge if the initial guesses are not sufficiently close to the root. It is particularly useful when derivative information is not available, unlike the Newton-Raphson Method.

**May fail to converge if:**
* The initial guesses are not sufficiently close to the root
* f(xₙ) - f(xₙ₋₁) ≈ 0 (division by near-zero value)
* The function has multiple roots or discontinuities near the initial guesses

### Input Characteristics

1. The first line contains two real numbers for the initial guesses:

   L  R

   representing x₀ and x₁.

2. The second line contains the allowed error (tolerance):

   ε

### Output Characteristics

* The approximate root of the function is displayed.
* For each iteration, the corresponding values and the approximate root are shown.

### Advantages

* **No derivative required** - Unlike Newton-Raphson Method
* **Faster than Bisection** - Super-linear convergence rate
* **No sign change required** - Unlike Bisection or False Position methods

### Disadvantages

* **Requires two initial guesses** - Both should be close to the root
* **May diverge** - If initial guesses are poor or function is problematic
* **Division by zero risk** - When f(xₙ) = f(xₙ₋₁)
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

### 1. Introduction
Newton's Forward Interpolation is a technique used to estimate the value of a function for an intermediate point ($x$) based on a set of given data points. This method is specifically designed for datasets where the interval between points is **equal** (equidistant).

It is most accurate when the point you are trying to calculate is located near the **beginning** of the dataset ($x$ is close to $x_0$).

### 2. Mathematical Principle
The method uses "Forward Differences". A forward difference is calculated by subtracting the current value from the next value:
$$\Delta y_0 = y_1 - y_0$$



**The General Formula:**
To find $y$ at a point $x$:

$$y(x) = y_0 + u \Delta y_0 + \frac{u(u-1)}{2!} \Delta^2 y_0 + \frac{u(u-1)(u-2)}{3!} \Delta^3 y_0 + \dots$$

Where:
* $y_0$: The first value in the $y$ column.
* $\Delta y_0, \Delta^2 y_0$: The top values of the difference columns.
* $u$: A standardized variable representing the distance from $x_0$.

**Calculating 'u':**
$$u = \frac{x - x_0}{h}$$

*(Where $h$ is the constant difference between x-values, i.e., $h = x_1 - x_0$)*

### 3. The Algorithm Steps
1.  **Verify Intervals:** Ensure that the input $x$ values are equally spaced (constant step size $h$).
2.  **Construct Difference Table:**
    * First Order Differences: $\Delta y_i = y_{i+1} - y_i$
    * Second Order Differences: $\Delta^2 y_i = \Delta y_{i+1} - \Delta y_i$
    * Continue until you reach a single value or negligible error.
3.  **Calculate u:** Determine how far $x$ is from the starting point $x_0$ using $u = (x - x_0) / h$.
4.  **Apply Formula:** Sum the terms using the values from the **top row** of the difference table.

### 4. Conditions
* **Equal Spacing:** The method fails if the difference between $x$ values is not constant.
* **Position:** Best used when the value to be interpolated is near the start of the table ($x < x_n/2$). If the value is near the end, Newton's *Backward* Interpolation should be used instead.

### 5. Advantages vs. Disadvantages

**Advantages**
* **Simpler Calculation:** For polynomial approximation, it is often simpler to compute manually than Lagrange's method.
* **Convergent:** For smooth functions, the error decreases as you add higher-order difference terms.

**Disadvantages**
* **Rigid Input:** Strictly requires equally spaced data points.
* **Location Dependent:** Accuracy degrades if the target point is far from the beginning of the list (use Backward or Central interpolation for other areas).

### Newton Forward Code
```cpp
#include<bits/stdc++.h>
using namespace std;

double u_cal(double u, int n){
    double temp = u;
    for(int i = 1; i < n; i++){
        temp *= (u - i);
    }
    return temp;
}

double fact(int n){
    double temp = 1;
    for(int i = 1; i <= n; i++){
        temp *= i;
    }
    return temp;
}

int main(){
    ifstream in("input.txt");
    ofstream out("output.txt");

    out << fixed << setprecision(3);

    int n;
    in >> n;

    vector<double> x(n);
    vector<vector<double>> y(n, vector<double>(n));

    for(int i = 0; i < n; i++){
        in >> x[i];
    }
    for(int i = 0; i < n; i++){
        in >> y[i][0];
    }

    for(int i = 1; i < n; i++){
        for(int j = 0; j < n - i; j++){
            y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
        }
    }

    for(int i = 0; i < n; i++){
        out << x[i] << " ";
        for(int j = 0; j < n; j++){
            out << y[i][j] << " ";
        }
        out << endl;
    }

    double value;
    in >> value;

    double h = x[1] - x[0];
    double u = (value - x[0]) / h;

    double sum = y[0][0];

    for(int i = 1; i < n; i++){
        sum += (u_cal(u, i) * y[0][i]) / fact(i);
    }

    out <<"\nThe result is: "<<sum<<endl;

    return 0;
}
```
### Newton Forward Input
```
4
1 2 3 4
1 8 27 64
2.5
```
### Newton Forward Output
```
1.000 1.000 7.000 12.000 6.000 
2.000 8.000 19.000 18.000 0.000 
3.000 27.000 37.000 0.000 0.000 
4.000 64.000 0.000 0.000 0.000 

The result is: 15.625
```
---
### Newton Backward Interpolation

### Newton Backward Theory

Newton's Backward Interpolation is used to estimate the value of an unknown variable x which is greater than the middle value of the given data. This method is applicable when the difference between any two consecutive values of x is constant.

### Given Data Points

Let the given data points be x₀, x₁, …, xₙ with corresponding values y₀, y₁, …, yₙ.

### Conditions

The data must satisfy the condition:

xᵢ - xᵢ₋₁ = h (constant), 1 ≤ i ≤ n

and the interpolation point should satisfy:

x > (x₀ + xₙ) / 2

### Backward Difference

The backward difference is defined as:

∇yᵢ = yᵢ - yᵢ₋₁

Higher order backward differences are:

∇²yᵢ = ∇(∇yᵢ)
∇³yᵢ = ∇(∇²yᵢ)

and so on.

### Interpolation Parameter

Let

u = (x - xₙ) / h

### Newton's Backward Interpolation Formula

Then the Newton's Backward Interpolation formula is given by:

y(x) = yₙ + u∇yₙ + [u(u+1) / 2!]∇²yₙ + [u(u+1)(u+2) / 3!]∇³yₙ + ⋯

### Algorithm Steps

1. Verify that the data points have equal spacing: xᵢ - xᵢ₋₁ = h
2. Verify that x > (x₀ + xₙ) / 2
3. Construct the backward difference table:
   - Calculate first differences: ∇yᵢ = yᵢ - yᵢ₋₁
   - Calculate second differences: ∇²yᵢ = ∇yᵢ - ∇yᵢ₋₁
   - Continue for higher order differences
4. Calculate u = (x - xₙ) / h
5. Apply the interpolation formula using the last row of the difference table

### When to Use

* Use Newton's Backward Interpolation when interpolating near the **end** of the data set
* The interpolation point x should be closer to xₙ than to x₀
* Particularly useful for extrapolation beyond the last data point

### Advantages

* **Efficient for end-point interpolation** - Best suited when x is near xₙ
* **Simple computation** - Uses backward differences from the last row
* **Equal spacing** - Works with equally spaced data points

### Disadvantages

* **Requires equal spacing** - Data points must have constant interval h
* **Limited applicability** - Only effective for x > (x₀ + xₙ) / 2
* **Accuracy decreases** - As x moves further from the data range
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

Newton's Divided Difference Interpolation is used to estimate the value of an unknown variable x when the data points are not equally spaced. This method constructs an interpolating polynomial that passes through all given data points.

### Given Data Points

Let the given data points be x₀, x₁, x₂, …, xₙ with corresponding values y₀, y₁, y₂, …, yₙ, where the values of x are not necessarily equally spaced.

### First Divided Difference

The first divided difference is defined as:

f[xᵢ, xᵢ₊₁] = (yᵢ₊₁ - yᵢ) / (xᵢ₊₁ - xᵢ)

### Second Divided Difference

The second divided difference is defined as:

f[xᵢ, xᵢ₊₁, xᵢ₊₂] = (f[xᵢ₊₁, xᵢ₊₂] - f[xᵢ, xᵢ₊₁]) / (xᵢ₊₂ - xᵢ)

### Higher Order Divided Differences

Higher order divided differences are defined recursively as:

f[xᵢ, xᵢ₊₁, …, xᵢ₊ₖ] = (f[xᵢ₊₁, …, xᵢ₊ₖ] - f[xᵢ, …, xᵢ₊ₖ₋₁]) / (xᵢ₊ₖ - xᵢ)

### Newton's Divided Difference Interpolation Formula

The Newton's Divided Difference Interpolation formula is given by:

y(x) = y₀ + (x - x₀)f[x₀, x₁] + (x - x₀)(x - x₁)f[x₀, x₁, x₂]
     + (x - x₀)(x - x₁)(x - x₂)f[x₀, x₁, x₂, x₃] + ⋯

### Algorithm Steps

1. Construct the divided difference table:
   - First column: y₀, y₁, …, yₙ
   - Second column: First divided differences f[xᵢ, xᵢ₊₁]
   - Third column: Second divided differences f[xᵢ, xᵢ₊₁, xᵢ₊₂]
   - Continue until only one value remains
2. Use the first row of each column in the interpolation formula
3. Evaluate y(x) at the desired point x

### When to Use

This method is suitable for interpolation when the data points are **unequally spaced** and provides a flexible way to construct the interpolating polynomial.

### Advantages

* **Works with unequal spacing** - No requirement for constant interval between data points
* **Flexible** - Can be easily extended to add more data points
* **Efficient computation** - Reuses previously computed differences

### Disadvantages

* **Computational complexity** - Requires construction of complete difference table
* **Numerical errors** - Can accumulate in higher order differences
* **Sensitivity** - Results may be sensitive to the ordering of data points
### Newton Divided Difference Code
```cpp
#include<bits/stdc++.h>
using namespace std;

double DivideDiff(const vector<double>& x, const vector<double>& y, double value, int n){
    vector<vector<double>> f(n, vector<double>(n));

    for(int i = 0; i < n; i++){
        f[i][0] = y[i];
    }

    for(int j = 1; j < n; j++){
        for(int i = 0; i < n - j; i++){
            f[i][j] = (f[i + 1][j - 1] - f[i][j - 1]) / (x[i + j] - x[i]);
        }
    }

    double result = f[0][0];

    for(int i = 1; i < n; i++){
        double term = f[0][i];
        for(int j = 0; j < i; j++){
            term *= (value - x[j]);
        }
        result += term;
    }

    return result;
}

int main(){
    ifstream in("input.txt");
    ofstream out("output.txt");

    int n;
    in >> n;

    vector<double> x(n), y(n);
    for(int i = 0; i < n; i++){
        in >> x[i];
    }
    for(int i = 0; i < n; i++){
        in >> y[i];
    }

    double value;
    in >> value;

    out << fixed << setprecision(6);
    out <<"The result is: "<< DivideDiff(x, y, value, n) << endl;

    return 0;
}

```
### Newton Divided Difference Input
```
4
1 2 3 4
1 8 27 64
2.5

```
### Newton Divided Difference Output
```
The result is: 15.625000

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
#include <bits/stdc++.h>
using namespace std;

double func(double x){
    return ( pow(sin(x),5) + 4*pow(sin(x),4) + 1 );
}

double dfunc(double x){
    return ( 5*pow(sin(x),4)*cos(x) + 16*pow(sin(x),3)*cos(x) );
}

double ddfunc(double x){
    return ( -5*pow(sin(x),5) + 20*pow(sin(x),3)*pow(cos(x),2) - 16*pow(sin(x),4) + 48*pow(sin(x),2)*pow(cos(x),2) );
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
        for(int i=n-1; i>=j; i--){
            y[i][j] = y[i][j-1] - y[i-1][j-1];
        }
    }

    fout<<"\nThe difference table is : \n";
    for(int i=0; i<n; i++){
        fout<<x[i]<<"\t";
        for(int j=0; j<=i; j++){
            fout<<y[i][j]<<"\t";
        }
        fout<<endl;
    }

    double d = x[1] - x[0];
    double v = (p - x[n-1]) / d;

    double y1 = ( y[n-1][1] + ((2*v + 1)/2)*y[n-1][2] + ((3*v*v + 6*v+ 2)/6)*y[n-1][3] + ((4*v*v*v + 18*v*v + 22*v + 6)/24)*y[n-1][4] ) / d;
    fout<<"\nThe value of f'(p) is = "<<y1<<endl;

    double y2 = ( y[n-1][2] + (v+1)*y[n-1][3] + ((12*v*v + 36*v + 22)/24)*y[n-1][4] ) / (d*d);
    fout<<"\nThe value of f''(p) is = "<<y2<<endl;

    double err1 = ( abs(dfunc(p) - y1) / abs(dfunc(p)) ) * 100;
    double err2 = ( abs(ddfunc(p) - y2) / abs(ddfunc(p)) ) * 100;
    fout<<"\nThe relative error of f'(p) = "<<err1<<"%\n";
    fout<<"\nThe relative error of f''(p) = "<<err2<<"%\n";

    fin.close();
    fout.close();
    return 0;
}
```
### Differentiation by backward Input
```
5 0 15 4.5
```
### Differentiation by backward Output
```

The difference table is : 
0.000	1.000	
0.333	1.050	0.050	
0.667	1.675	0.626	0.576	
1.000	3.427	1.752	1.126	0.550	
1.333	5.437	2.010	0.257	-0.869	-1.419	
1.667	5.904	0.467	-1.542	-1.800	-0.931	0.489	
2.000	4.356	-1.548	-2.015	-0.473	1.326	2.257	1.768	
2.333	2.291	-2.065	-0.517	1.499	1.972	0.646	-1.611	-3.380	
2.667	1.195	-1.096	0.969	1.486	-0.013	-1.985	-2.631	-1.019	2.360	
3.000	1.002	-0.193	0.903	-0.066	-1.551	-1.538	0.447	3.078	4.097	1.737	
3.333	1.005	0.003	0.197	-0.706	-0.641	0.910	2.449	2.002	-1.076	-5.173	-6.910	
3.667	1.221	0.216	0.213	0.016	0.722	1.363	0.453	-1.996	-3.997	-2.922	2.251	9.161	
4.000	2.064	0.843	0.627	0.415	0.399	-0.324	-1.687	-2.139	-0.144	3.854	6.775	4.524	-4.637	
4.333	3.288	1.224	0.381	-0.246	-0.661	-1.060	-0.736	0.951	3.090	3.234	-0.620	-7.396	-11.920	-7.283	
4.667	3.989	0.701	-0.523	-0.903	-0.657	0.004	1.064	1.800	0.850	-2.240	-5.474	-4.854	2.542	14.462	21.745	

The value of f'(p) is = 2.298

The value of f''(p) is = -10.493

The relative error of f'(p) = 5.023%

The relative error of f''(p) = 17.384%

```
---

### Ordinary Differential Equation

### Runge Kutta Method

### Runge Kutta Theory
## The Runge-Kutta Method (RK4)

### 1. Introduction
The Runge-Kutta methods are a family of iterative methods used to approximate solutions to Ordinary Differential Equations (ODEs). The most widely used version is the **Fourth-Order Runge-Kutta (RK4)** method. It is far more accurate than the simple Euler's method because it takes four different slope measurements for every single step forward.

### 2. Mathematical Principle
Given a differential equation of the form:
$$\frac{dy}{dx} = f(x, y)$$

We want to find the value of $y$ at a future point. RK4 works by calculating a weighted average of four different slopes ($k_1, k_2, k_3, k_4$) within one step interval $h$.



**The Four Slopes:**
1.  **$k_1$ (Start Slope):** The slope at the beginning of the interval.
2.  **$k_2$ (Midpoint Slope 1):** The slope at the midpoint, using $k_1$ to estimate $y$.
3.  **$k_3$ (Midpoint Slope 2):** Another slope at the midpoint, but using $k_2$ to improve the estimate.
4.  **$k_4$ (End Slope):** The slope at the end of the interval.

### 3. The Formulas
To move from a point $(x_n, y_n)$ to the next point $(x_{n+1}, y_{n+1})$ with a step size $h$:

$$k_1 = h \cdot f(x_n, y_n)$$
$$k_2 = h \cdot f(x_n + \frac{h}{2}, y_n + \frac{k_1}{2})$$
$$k_3 = h \cdot f(x_n + \frac{h}{2}, y_n + \frac{k_2}{2})$$
$$k_4 = h \cdot f(x_n + h, y_n + k_3)$$

**The Final Update Formula:**
$$y_{n+1} = y_n + \frac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

*(Note: The middle slopes $k_2$ and $k_3$ are given more weight because they represent the "average" behavior in the center of the step.)*

### 4. The Algorithm Steps
1.  **Define:** Start with an initial condition $(x_0, y_0)$, a step size $h$, and the target $x$ value.
2.  **Loop:** While $x < target$:
    * Calculate $k_1, k_2, k_3, k_4$ using the formulas above.
    * Calculate the new $y$ value: $y_{new} = y_{old} + \frac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4)$.
    * Update $x$: $x_{new} = x_{old} + h$.
3.  **Repeat:** Continue until the target $x$ is reached.

### 5. Accuracy and Convergence
* **Order of Accuracy:** Fourth-order ($O(h^4)$). This means if you halve the step size $h$, the error decreases by a factor of 16 ($2^4$).
* **Stability:** It is much more stable and accurate than Euler's method or the Midpoint method for the same step size.

### 6. Advantages vs. Disadvantages

**Advantages**
* **High Accuracy:** It provides excellent accuracy without requiring an extremely small step size.
* **Self-Starting:** Unlike multi-step methods (like Adams-Bashforth), RK4 only needs the current point to calculate the next one.

**Disadvantages**
* **Computational Cost:** It requires evaluating the derivative function $f(x,y)$ four times per step (Euler's method only calculates it once).
* **Complexity:** The formula is more complex to implement than basic first-order methods.

### Runge Kutta Code
```cpp
#include<bits/stdc++.h>
using namespace std;

double dydx(double x, double y)
{
    return x + y*y;
}

double rungekutta(double x0, double y0, double xp, double h)
{
    double x = x0;
    double y = y0;

    while (x < xp)
    {
        double k1 = h * dydx(x, y);
        double k2 = h * dydx(x + 0.5*h, y + 0.5*k1);
        double k3 = h * dydx(x + 0.5*h, y + 0.5*k2);
        double k4 = h * dydx(x + h, y + k3);

        y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
        x += h;
    }
    return y;
}

int main()
{
    freopen("input.txt","r",stdin);
    freopen("output.txt","w",stdout);

    double x0, y0, xp, h;
    cin >> x0 >> y0 >> xp >> h;

    cout << "Result: " << rungekutta(x0, y0, xp, h);
}

```
### Runge Kutta Input
```
0 1 0.2 0.1
```
### Runge Kutta Output
```
Result: 1.27356
```
---
