#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "RKF45.hpp"

const double relerr = 1.0e-13;
const double abserr = 1.0e-13;
const double t_start = 1.0;
const double t_end = 2.0;

double exact_solution(double t) {
	return exp(2 * t);
}

void equations(const double t, double* x, double* dx) {
	dx[0] = x[1];
	dx[1] = (t + 1.0) * x[1] / t + 2.0 * (t - 1.0) * x[0] / t;
}

double localRKF[4][2];

void rkf(double h, int count) {
	const int n = 2;
	const int AMT = 500;
	double x[2] = {exp(2), 2 * exp(2)};
	double dx[2];

	double* t = new double (t_start);
	double* tout = new double (t_start);
	double* re = new double (relerr);
	double* ae = new double (abserr);
	int flag = 1;

	double data[AMT][2];
	int i = 0;

	for (*tout = t_start; *tout <= t_end + h; *tout += h) {
		r8_rkf45(equations, n, x, dx, t, *tout, re, *ae, flag);
		data[i][0] = *tout;
		data[i++][1] = x[0];
	}
	delete t, tout, re, ae, flag;

	std::cout << "\n  * RKF45 method with h = " << std::setprecision(4) << h << " *\n" << std::endl;
	std::cout << "  t\tRKF45 method\t\t\tExact value\t\t\tDifference" << std::endl;

	double num = 0;
	double global = 0;
	for (int k = 0; k < 101; k++) {
		double hprint = data[k][0];
		double result = data[k][1];
		double exact = exact_solution(hprint);

		if (hprint - num * 0.1 - 1 > -0.000001) {
			std::cout << std::fixed << std::setprecision(1) << hprint << "\t";
			std::cout << std::fixed << std::setprecision(20) << result << "\t\t";
			std::cout << std::fixed << std::setprecision(20) << exact << "\t\t";
			std::cout << std::fixed << std::setprecision(20) << abs(exact - result) << std::endl;
			global += abs(exact - result);
			num += 1.0;
		}
	}
	localRKF[count][0] = abs(exact_solution(data[1][0]) - data[1][1]);
	localRKF[count][1] = global;
	std::cout << "==============================================================================================\n";
}

double localTrap[4][2];

void trap(double h, int count) {
	double t = t_start;
	const int AMT = 500;
	double x = exp(2);
	double data[AMT][2];

	std::cout << "\n  * Trapezoidal method with h = " << std::setprecision(4) << h << " *\n" << std::endl;
	std::cout << "  t\tTrapezoidal method\t\tExact value\t\t\tDifference" << std::endl;

	int i = 0;
	while (t <= t_end + h) {
		data[i][0] = t;
		data[i][1] = x;

		double x_next = x + h * (x * (t + h + 1) / (t + h) + 2 * (t - 1) * x / t);
		x += (h / 2) * (x * (t + h + 1) / t + x_next * (t + h + 1) / (t + h) + 2 * (t - 1) * (x + x_next) / (2 * t + h));
		t += h;
		i++;
	}

	double num = 0;
	double global = 0;
	for (int k = 0; k < 101; k++) {
		double hprint = data[k][0];
		double result = data[k][1];
		double exact = exact_solution(hprint);

		if (hprint - num * 0.1 - 1 > -0.000001) {
			std::cout << std::fixed << std::setprecision(1) << hprint << "\t";
			std::cout << std::fixed << std::setprecision(20) << result << "\t\t";
			std::cout << std::fixed << std::setprecision(20) << exact << "\t\t";
			std::cout << std::fixed << std::setprecision(20) << abs(exact - result) << std::endl;
			global += abs(exact - result);
			num += 1.0;
		}
	}
	localTrap[count][0] = abs(exact_solution(data[1][0]) - data[1][1]);
	localTrap[count][1] = global;

	std::cout << "==============================================================================================\n\n";
}

int main() {
	double hint[4] = { 0.1, 0.05, 0.025, 0.0125 };
	int count = 0;

	for (auto h : hint) {
		rkf(h, count);
		count++;
	}
	count = 0;
	for (auto h : hint) {
		trap(h, count);
		count++;
	}


	for (int i = 0; i < 4; i++) {
		std::cout << "For h = " << std::fixed << std::setprecision(4) << hint[i] << ":\n";
		std::cout << "\tTrapezoidal method: \n";
		std::cout << "\t\tLocal error on the first step equals:\t" 
			<< std::fixed << std::setprecision(20) << localTrap[i][0] << "\n";
		std::cout << "\t\tGlobal error equals:\t\t\t" << std::fixed << std::setprecision(20) << localTrap[i][1] << "\n";
		std::cout << "\tRKF45 method: \n";
		std::cout << "\t\tLocal error on the first step equals:\t"
			<< std::fixed << std::setprecision(20) << localRKF[i][0] << "\n";
		std::cout << "\t\tGlobal error equals:\t\t\t" << std::fixed << std::setprecision(20) << localRKF[i][1] << "\n\n";
	}
}