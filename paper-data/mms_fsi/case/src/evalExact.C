// Reads "fluid x y t" or "solid X Y t" or "interface X t" lines from stdin,
// prints the exact manufactured values. Build:
//   clang++ -O1 -std=c++17 -I ../derive/generated src/evalExact.C -o evalExact
#include <cstdio>
#include <iostream>
#include <string>
#include "mmsFsiExact.H"
int main()
{
    std::string kind;
    while (std::cin >> kind)
    {
        if (kind == "fluid")
        {
            double x, y, t, U[2], p[1], f[2], G[4], gp[2];
            std::cin >> x >> y >> t;
            mmsFsi::fluid(x, y, t, U, p, f, G, gp);
            std::printf("%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n", U[0], U[1], p[0], f[0], f[1], gp[0], gp[1], G[0], G[1], G[2], G[3]);
        }
        else if (kind == "solid")
        {
            double X, Y, t, D[2], U[2], b[2], s[6];
            std::cin >> X >> Y >> t;
            mmsFsi::solid(X, Y, t, D, U, b, s);
            std::printf("%.15g %.15g %.15g %.15g %.15g %.15g\n", D[0], D[1], U[0], U[1], b[0], b[1]);
        }
        else if (kind == "interface")
        {
            double X, t, x[2], tr[2], n[2];
            std::cin >> X >> t;
            mmsFsi::interface(X, t, x, tr, n);
            std::printf("%.15g %.15g %.15g %.15g %.15g %.15g\n", x[0], x[1], tr[0], tr[1], n[0], n[1]);
        }
    }
    return 0;
}
