#include <stdio.h>
#include <math.h>
#include "eq.h"

#define EPS 1e-6
#define H 1e-5            // шаг для численной производной
#define ITERS_LIMIT 100  

// F'(x) ~ (F(x+h) - F(x-h)) / 2h
double d(double (*F)(double), double x)
{
	return (F(x + H) - F(x - H)) / (2 * H);
}

// F''(x)~(F(x+h)-2F(x)+F(x-h))/h^2
double d_2(double (*F)(double), double x)
{
	return (F(x + H) - 2 * F(x) + F(x - H)) / (H * H);
}

double dichotomy_method(double (*F)(double), double a, double b, double eps,
	int* iterations, double* accuracy)
{
	*iterations = 0;
	*accuracy = 0;

	double fa = F(a);
	if (fa * F(b) >= 0)
	{
		// на концах один знак
		*iterations = -1;
		return a;
	}

	while (b - a > eps)
	{
		double middle = (a + b) / 2;
		double fmiddle = F(middle);
		if (fmiddle * fa > 0)
		{
			a = middle;     // корень в правой половине
			fa = fmiddle;
		}
		else
			b = middle;     // в левой
		(*iterations)++;
	}

	*accuracy = (b - a) / 2;
	return (a + b) / 2;
}

double iterations_method(double (*f)(double), double a, double b, double eps, int iters_limit,
	int* iterations, double* accuracy)
{
	double x = (a + b) / 2;
	double x_prev;
	*iterations = 0;

	do
	{
		// |f'(x)|<1
		if (fabs(d(f, x)) >= 1)
		{
			*iterations = -1;
			return x;
		}

		x_prev = x;
		x = f(x);
		(*iterations)++;

		if (*iterations == iters_limit || isnan(x) || isinf(x))
		{
			*iterations = 0;
			*accuracy = fabs(x - x_prev);
			return -1;
		}
	} while (fabs(x - x_prev) >= eps);

	*accuracy = fabs(x - x_prev);
	return x;
}

double newton_method(double (*F)(double), double a, double b, double eps, int iters_limit,
	int* iterations, double* accuracy)
{
	double x = (a + b) / 2;
	double x_prev;
	*iterations = 0;

	do
	{
		double fx = F(x);
		double d1 = d(F, x);
		double d2 = d_2(F, x);

		// |F * F''|< (F')^2
		if (fabs(fx * d2) >= d1 * d1)
		{
			*iterations = -1;
			return x;
		}

		x_prev = x;
		x = x - fx / d1;
		(*iterations)++;

		if (*iterations == iters_limit || isnan(x) || isinf(x))
		{
			*iterations = 0;
			*accuracy = fabs(x - x_prev);
			return -1;
		}
	} while (fabs(x - x_prev) >= eps);

	*accuracy = fabs(x - x_prev);
	return x;
}

void print_result(const char* method, double result, double ideal_result, double accuracy, int iterations)
{
	if (iterations == -1)
		printf("%-11s condition failed. x = %.4f\n", method, result);
	else
		printf("%-11s x = %.7f table: %.4f accuracy: %.2e iterations: %d\n",
			method, result, ideal_result, accuracy, iterations);
}

int main(void)
{
	for (int i = 0; i < 3; i++)
	{
		expression_data t = tasks[i];

		int it;
		double acc;
		double res;

		printf("%s,  [%.1f, %.1f]\n", t.name, t.a, t.b);

		res = dichotomy_method(t.F, t.a, t.b, EPS, &it, &acc);
		print_result("dichotomy", res, t.ideal_result, acc, it);

		res = iterations_method(t.f, t.a, t.b, EPS, ITERS_LIMIT, &it, &acc);
		print_result("iterations", res, t.ideal_result, acc, it);

		res = newton_method(t.F, t.a, t.b, EPS, ITERS_LIMIT, &it, &acc);
		print_result("newton", res, t.ideal_result, acc, it);

		printf("\n");
	}

	return 0;
}
