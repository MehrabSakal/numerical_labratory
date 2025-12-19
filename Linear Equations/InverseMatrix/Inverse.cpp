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
