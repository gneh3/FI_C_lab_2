#include <stdio.h>
#include <math.h>

#include "eq.h"

#define EPS 1e-6          // требуемая точность
#define ITERS_LIMIT 1000  // защита от зависания, если метод расходится

double max_on_segment(double (*g)(double), double a, double b)
{
	double m = 0;
	for (int i = 0; i <= 1000; i++)
	{
		double v = fabs(g(a + (b - a) * i / 1000));
		if (v > m) m = v;
	}
	return m;
}

double newton_condition_max(double (*F)(double), double (*dF)(double), double (*d2F)(double),
	double a, double b)
{
	double m = 0;
	for (int i = 0; i <= 1000; i++)
	{
		double x = a + (b - a) * i / 1000;
		double v = fabs(F(x) * d2F(x)) / (dF(x) * dF(x));
		if (v > m) m = v;
	}
	return m;
}

// все методы при неудаче возвращают -1 и *iterations = 0.

double dichotomy_method(double (*F)(double), double a, double b, double eps,
	int* iterations, double* accuracy)
{
	*iterations = 0;
	*accuracy = 0;

	while (b - a > eps)
	{
		double fa = F(a);
		double fb = F(b);
		if (fa * fb > 0)
		{
			// метод неприменим
			return -1;
		}

		double middle = (a + b) / 2;
		double fmiddle = F(middle);
		if (fmiddle * fa > 0)
			a = middle;     // корень в правой половине
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
	double x_prev = (a + b) / 2;
	double x = f(x_prev);
	*iterations = 1;

	while (fabs(x - x_prev) >= eps)
	{
		x_prev = x;
		x = f(x);
		(*iterations)++;
		if (*iterations == iters_limit || isnan(x) || isinf(x))
		{
			// не сошлось
			*iterations = 0;
			*accuracy = fabs(x - x_prev);
			return -1;
		}
	}

	*accuracy = fabs(x - x_prev);
	return x;
}

double newton_method(double (*F)(double), double (*dF)(double), double a, double b, double eps,
	int iters_limit, int* iterations, double* accuracy)
{
	double x_prev = (a + b) / 2;
	double x = x_prev - F(x_prev) / dF(x_prev);
	*iterations = 1;

	while (fabs(x - x_prev) >= eps)
	{
		x_prev = x;
		x = x_prev - F(x_prev) / dF(x_prev);
		(*iterations)++;
		if (*iterations == iters_limit || isnan(x) || isinf(x))
		{
			// не сошлось
			*iterations = 0;
			*accuracy = fabs(x - x_prev);
			return -1;
		}
	}

	*accuracy = fabs(x - x_prev);
	return x;
}

// одна строка таблицы: метод, корень, корень из таблицы, точность, итерации
void print_result(const char* method, double root, double table_root, double accuracy, int iterations)
{
	if (iterations == 0)
		printf("%-11s broken\n", method);
	else
		printf("%-11s x = %.7f table: %.4f accuracy: %.2e iterations: %d\n",
			method, root, table_root, accuracy, iterations);
}

int main(void)
{
	for (int i = 0; i < 3; i++)
	{
		expression_data t = tasks[i];

		int it;
		double acc;
		double root;

		printf("%s,  [%.1f, %.1f]\n", t.name, t.a, t.b);

		double q_iter = max_on_segment(t.df, t.a, t.b);
		double q_newton = newton_condition_max(t.F, t.dF, t.d2F, t.a, t.b);
		printf("check: F(a)*F(b) = %.3f (%s)\n", t.F(t.a) * t.F(t.b),
			t.F(t.a) * t.F(t.b) < 0 ? "< 0, OK" : "NOT < 0");
		printf("check: max|f'(x)| = %.3f (%s)\n", q_iter,
			q_iter < 1 ? "< 1" : "> 1");
		printf("check: max|F*F''|/(F')^2 = %.3f (%s)\n", q_newton,
			q_newton < 1 ? "< 1" : "> 1, newton condition failed on [a,b]");

		root = dichotomy_method(t.F, t.a, t.b, EPS, &it, &acc);
		print_result("dichotomy", root, t.table_root, acc, it);

		root = iterations_method(t.f, t.a, t.b, EPS, ITERS_LIMIT, &it, &acc);
		print_result("iterations", root, t.table_root, acc, it);

		root = newton_method(t.F, t.dF, t.a, t.b, EPS, ITERS_LIMIT, &it, &acc);
		print_result("newton", root, t.table_root, acc, it);

		printf("\n");
	}

	return 0;
}
