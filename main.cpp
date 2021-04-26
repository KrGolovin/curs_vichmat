#include <iostream>
#include <cmath>

extern "C" {
#include "cmath.h"
#include "zeroin.h"
#include "nelmin.h"
#include "quanc.h"
}


double T = 1000;

double y(double x)
{
  return sqrt(x) - 2 * cos(M_PI * x / 2);
}

double f(int n, double *zp)
{
  double z = *zp;
  return std::exp(z) * (2 * z * z - 4) +
         std::pow(2 * z * z - 1, 2) +
         std::exp(2 * z) -
         3 * std::pow(z, 4);

}

double integralF(double x)
{
  return 1 / (std::pow(x, 5) * (std::exp(1.432 / T / x) - 1));
}

int main()
{
  int flag;
  double valueY = zeroin(0, 0.00000001, y, 0.0000001, &flag);
  double lowerBound = 5.548752e-5 * valueY;
  
  std::cout << lowerBound << std::endl;
  double xMin = 0.1;
  double ynewlo;
  double step = 0.001;
  int icount = 0;
  int numres = 0;
  int ifault = 0;
  double reltol = 1.0e-10;
  double abstol = 0;
  nelmin(f, 1, &xMin, &ynewlo, 0.0001, &step, 10000, &icount, 900, &numres, &ifault, reltol, abstol);
  double upperBound = 13.02892e-5 * xMin;
  
  std::cout << upperBound << std::endl;
  double a, b, epsabs, epsrel, errest, posn;
  int nfe;
  epsrel = 1.0e-10;
  epsabs = 0.0;
  double resultR[9] = {};
  for (int i = 0; i < 9; ++i, T += 1000) {
    quanc8(integralF, lowerBound, upperBound, epsabs, epsrel, resultR + i, &errest,
           &nfe, &posn, &flag);
    resultR[i] *= (64.77 / std::pow(T, 4));
    std::cout << "T = " << T <<" | " << resultR[i] << std::endl;
  }
  return 0;
}
