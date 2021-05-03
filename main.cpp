#include <iostream>
#include <cmath>

extern "C" {
#include "cmath.h"
#include "zeroin.h"
#include "fmin.h"
#include "quanc.h"
}

double T = 1000;

double y(double x) {
  return sqrt(x) - 2 * cos(M_PI * x / 2);
}

double f(double z) {
  return std::exp(z) * (2 * z * z - 4) + std::pow(2 * z * z - 1, 2) + std::exp(2 * z) - 3 * std::pow(z, 4);

}

double integralF(double x) {
  return 1 / (std::pow(x, 5) * (std::exp(1.432 / T / x) - 1));
}

int main() {
  double epsabs, epsrel, errest, posn;
  epsrel = 1.0e-10;
  epsabs = 0.0;
  for (int i = 0; i < 2; ++i) {
    int flag;
    double valueY = zeroin(0, 0.001, y, epsrel, &flag);
    double lowerBound = 5.548752e-5 * valueY;

    std::cout << "lambda 1 =  " << lowerBound << std::endl;

    double xMin = FMin(f, 0.1, 1, epsrel);
    double upperBound = 13.02892e-5 * xMin;

    std::cout << "lambda 2 = " << upperBound << std::endl;
    int nfe;
    epsabs = 0.0;
    double resultR[9] = {};
    double absresult[9] = {3.08707e6, 1.9908e11, 1.03247e13, 8.32523e13, 3.07104e14, 7.54059e14, 1.45732e15, 2.4176e15,
                           3.61697e15};
    std::cout << "            Quanc8 |       AbsValue   |    AbsErr " << std::endl;
    for (int i = 0; i < 9; ++i, T += 1000) {
      quanc8(integralF, lowerBound, upperBound, epsabs, epsrel, resultR + i, &errest, &nfe, &posn, &flag);
      resultR[i] *= (64.77 / std::pow(T, 4));
      absresult[i] *= (64.77 / std::pow(T, 4));
      std::cout << "T = " << T << " | " << resultR[i] << " | " << absresult[i] << " | "
                << abs(resultR[i] - absresult[i]) << std::endl;
    }
    epsrel = 1.0e-4;
    T = 1000;
  }
  return 0;
}
