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
