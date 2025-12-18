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
