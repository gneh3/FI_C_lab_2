#include <math.h>
#include "eq.h"

double F27(double x)
{
    double t = log(x);
    return sin(t) - cos(t) + 2 * t;
}

double f27(double x)
{
    double t = log(x);
    return exp((cos(t) - sin(t)) / 2);
}

double F28(double x)
{
    return x - 2 + sin(1 / x);
}

double f28(double x)
{
    return 2 - sin(1 / x);
}

double F1(double x)
{
    return exp(x) + log(x) - 10 * x;
}

double f1(double x)
{
    return log(10 * x - log(x));
}

const expression_data tasks[3] = {
    { "sin(ln x) - cos(ln x) + 2 ln x = 0", 1.0, 3.0, 1.3749, F27, f27 },
    { "x - 2 + sin(1/x) = 0",               1.2, 2.0, 1.3077, F28, f28 },
    { "e^x + ln x - 10x = 0",               3.0, 4.0, 3.5265, F1,  f1  },
};
