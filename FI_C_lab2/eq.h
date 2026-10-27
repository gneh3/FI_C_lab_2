#pragma once

#include <math.h>

typedef struct {
    const char* name;
    double a, b;
    double ideal_result;
    double (*F)(double);   // F(x) = 0
    double (*f)(double);   // x = f(x), для метода итераций
} expression_data;

double F27(double x);
double f27(double x);

double F28(double x);
double f28(double x);

double F1(double x);
double f1(double x);

extern const expression_data tasks[3];
