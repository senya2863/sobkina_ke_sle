#include "tridiagonal.h"

std::vector<double> solve(
    const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    const std::vector<double>& d)
{
    int n=b.size();

    std::vector<double> c1(n);
    std::vector<double> d1(n);
    std::vector<double> x(n);

    //right step
    c1[0]=c[0]/b[0];
    d1[0]=d[0]/b[0];

    for(int i=1;i<n;i++){
        double m = b[i]-a[i]*c1[i-1];
        c1[i] = i<n-1?c[i]/m:0;
        d1[i]=(d[i]-a[i]*d1[i-1])/m;
    }
    //reverse step
        x[n-1]=d1[n-1];
        for (int i=n-2;i>=0;i--){
            x[i]=d1[i]-c1[i]*x[i+1];
        }
        return x;
    }