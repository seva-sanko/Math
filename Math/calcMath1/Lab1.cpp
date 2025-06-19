#include <stdio.h>
#include <iostream>
#include <cmath>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

//Исходная функция.
double function(double x) {
	return cos(x) / (1 + x);
}

//Расчет полинома Лагранжа.
double lagrange(int degree, double x[], double y[], double xk) {
	double result = 0.0;

	for (int i = 0; i < degree; i++) {
		double yi = y[i];
		for (int j = 0; j < degree; j++) {
			if (i != j) {
				yi *= (xk - x[j]) / (x[i] - x[j]);
			}
		}
		result += yi;
	}
	return result;
}

//Исходный интеграл.
double integrand(double x, void* parameters) {
	double m = *(double *)parameters;
	return pow(abs(x * x + 2 * x - 3), m);
}

int main() {
	const int degree = 9;
	double x[degree], y[degree];

	for (int k = 0; k < degree; k++) {
		x[k] = 0.2 * k;
		y[k] = function(x[k]);
	}

	const int points = 8;
	double xi[points];
	double yf[points];
	double yl[points];
	double ys[points];

	//Получение значений сплайн-функции через библиотеку gsl.
	gsl_interp_accel* acc = gsl_interp_accel_alloc();
	gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, degree);

	gsl_spline_init(spline, x, y, degree);

	for (int i = 0; i < points; i++) {
		xi[i] = 0.1 + 0.2 * i;
		yf[i] = function(xi[i]);
		yl[i] = lagrange(degree, x, y, xi[i]);
		ys[i] = gsl_spline_eval(spline, xi[i], acc);
	}

	//Вывод результатов.
	double differenceL = 0.0;
	double differenceS = 0.0;
	std::cout << "\n Values of three functions at specific points: " << std::endl;
	for (int i = 0; i < points; i++) {
		double value = yf[i];
		double valueL = yl[i];
		double valueS = ys[i];
		std::cout << " x = " << xi[i] << "; F(x) = " << value << "; L(x) = " << valueL << "; S(x) = " << valueS << ";" << std::endl;	
		differenceL += abs(value - valueL);
		differenceS += abs(value - valueS);
	}

	std::cout << "\n Differences: " << std::endl;
	std::cout << " Lagrange: " << differenceL << std::endl;
	std::cout << " Spline: " << differenceS << std::endl;

	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);

	std::cout << "\n===================================================================================" << std::endl;

	double result, err;
	double lower = 0.0;
	double upper = 2.1;
	double ms[] = { -1.0, -0.5 };

	//Расчет интеграла с помощью QUANC8(cquad с лимитом в 30 делений) через библиотеку gsl.
	for (double m : ms) {
		gsl_function Func;
		Func.function = &integrand;
		Func.params = &m;
		gsl_integration_cquad_workspace* cqws = gsl_integration_cquad_workspace_alloc(30);
		size_t nevals;

		gsl_integration_cquad(&Func, lower, upper, 0, 1e-10, cqws, &result, &err, &nevals);

		std::cout << "\n x integral for |x^2 + 2x - 3|^(" << m;
		std::cout<< ") in the range from " << lower << " to " << upper << " equals: " << result << std::endl;

		gsl_integration_cquad_workspace_free(cqws);
	}
	std::cout << "\n===================================================================================" << std::endl;
	return 0;
}