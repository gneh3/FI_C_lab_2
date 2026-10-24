#pragma once

#include <math.h>

typedef struct {
    const char* name;          
    double a, b;               
    double table_root;         
    double (*F)(double);       // уравнение 
    double (*dF)(double);      // F'
    double (*d2F)(double);     // F'' 
    double (*f)(double);       // x = f(x)
    double (*df)(double);      // f'
} expression_data;

double F27(double x); 

double dF27(double x);

double d2F27(double x);

double f27(double x);

double df27(double x);

double F28(double x);

double dF28(double x);

double d2F28(double x);

double f28(double x);

double df28(double x);

double F1(double x);

double dF1(double x);

double d2F1(double x);

double f1(double x);

double df1(double x);

extern const expression_data tasks[3];
