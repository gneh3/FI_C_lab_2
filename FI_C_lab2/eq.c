#include <math.h>
#include "eq.h"

double F27(double x)
{
    double t = log(x);
    return sin(t) - cos(t) + 2 * t;
}

double dF27(double x)
{
    double t = log(x);
    return (cos(t) + sin(t) + 2) / x;
}

double d2F27(double x)
{
    double t = log(x);
    return -2 * (sin(t) + 1) / (x * x);
}

double f27(double x)
{
    double t = log(x);
    return exp((cos(t) - sin(t)) / 2);
}

double df27(double x)
{
    double t = log(x);
    return f27(x) * (-sin(t) - cos(t)) / (2 * x);
}

double F28(double x)
{
    return x - 2 + sin(1 / x);
}

double dF28(double x)
{
    return 1 - cos(1 / x) / (x * x);
}

double d2F28(double x)
{
    return -sin(1 / x) / pow(x, 4) + 2 * cos(1 / x) / pow(x, 3);
}

double f28(double x)
{
    return 2 - sin(1 / x);
}

double df28(double x)
{
    return cos(1 / x) / (x * x);
}

double F1(double x)
{
    return exp(x) + log(x) - 10 * x;
}

double dF1(double x)
{
    return exp(x) + 1 / x - 10;
}

double d2F1(double x)
{
    return exp(x) - 1 / (x * x);
}

double f1(double x)
{
    return log(10 * x - log(x));
}

double df1(double x)
{
    return (10 - 1 / x) / (10 * x - log(x));
}

const expression_data tasks[3] = {
    { "sin(ln x) - cos(ln x) + 2 ln x = 0", 1.0, 3.0, 1.3749,
      F27, dF27, d2F27, f27, df27 },
    { "x - 2 + sin(1/x) = 0",               1.2, 2.0, 1.3077,
      F28, dF28, d2F28, f28, df28 },
    { "e^x + ln x - 10x = 0",               3.0, 4.0, 3.5265,
      F1,  dF1,  d2F1,  f1,  df1 },
};