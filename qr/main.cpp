#include <iostream>
#include <vector>
#include <cmath>
#include <vector>

double norm (const std::vector<double>& x) { //norma vectora
    double s = 0;
    for (double y:x) {
        s+=y*y;
    }
    return std::sqrt(s);
}

std::vector<std::vector<double>> transp (const std::vector<std::vector<double>>& A) {
    int m =A.size(); //rows
    int n = A[0].size(); //colls
    std::vector<std::vector<double>> T(n, std::vector<double>(m));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            T[j][i]=A[i][j];
        }
    }
    return T;
}

//matrix on vector
std::vector<double> multiply (const std::vector<std::vector<double>>&A, const std::vector<double>& x) {
    int m = A.size();
    int n = A[0].size();
    std :: vector<double> res(m,0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            res[i] += A[i][j] * x[j];
        }
    }
    return res;
}

std::vector<std::vector<double>> multiply (const std::vector<std::vector<double>>&A, const std::vector<std::vector<double>>&B) {
    int m = A.size();
    int n = B[0].size();
    int p = B.size();
    std::vector<std::vector<double>> res(m, std::vector<double>(n,0));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < p; k++) {
                res[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return res;
}

void housh(std::vector<std::vector<double>> A, std::vector<std:: vector<double>>&Q,std::vector<std:: vector<double>>&R) {
    int m = A.size();
    int n = A[0].size();
    Q = std::vector<std::vector<double>>(m, std::vector<double>(m, 0));
    for (int i =0; i<m;i++) {
        Q[i][i]=1;
    }
    for (int k=0 ; k<n &&k<(m-1); k++) {
        std::vector<double> x(m-k); //белаем вектор для части столбцов
        for (int i = k; i < m; i++) {
            x[i-k]=A[i][k];
        }
        double n_x = norm(x);
        if (x[0]>0) {
            n_x=-n_x;
        }
        std::vector<double> y =x;
        y[0] -= n_x;
        double n_y=norm(y);
        for (double &v : y) v /= n_y;

        std::vector<std::vector<double>> P(m, std::vector<double>(m,0));
        for (int i = 0; i < m; i++) {
            P[i][i]=1;
        }
        for (int i = k; i<m; i++) {
            for (int j = k; j <m; j++) {
                P[i][j]-=2*y[i-k]*y[j-k];
            }
        }
        A = multiply(P,A);
        Q= multiply(Q,P);
    }
    R=A;
}

//Rx=b
std::vector<double> back(const std::vector<std::vector<double>>& R,const std::vector<double>& b) {
    int n=R.size();
    std::vector<double> res(n,0);
    for (int i = n-1; i >=0; i--) {
        double sum = b[i];
        for (int j = i + 1; j < n; j++)
            sum -= R[i][j] * res[j];

        res[i] = sum / R[i][i];
    }
    return res;
}

std::vector<double> solveQR(
    std::vector<std::vector<double>> A,
    std::vector<double> f)
{
    std::vector<std::vector<double>> Q, R;
    housh(A, Q, R);
    std::vector<std::vector<double>> Qt = transp(Q);
    std::vector<double> y = multiply(Qt, f);
    return back(R, y);
}

int main()
{
    std::vector<std::vector<double>> A =
    {
        {2, 1},
        {1, 3}
    };
    std::vector<double> f = {1, 2};
    std::vector<double> x = solveQR(A, f);
    for (double v : x)
        std::cout << v << " ";
    return 0;
}
// TIP See CLion help at <a
// href="https://www.jetbrains.com/help/clion/">jetbrains.com/help/clion/</a>.
//  Also, you can try interactive lessons for CLion by selecting
//  'Help | Learn IDE Features' from the main menu.