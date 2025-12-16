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
