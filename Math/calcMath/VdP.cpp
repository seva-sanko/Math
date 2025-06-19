#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

double w = 0.9517245; 
double tau = 0.0576371;
const double T = 10.0;
const double mu = 0.097;

int func(double t, const double y[], double f[], void* params) {
    f[0] = y[1]; //V'
    f[1] = -w * w * y[0] + 2 * mu * (1 - y[2]) * y[1] - 2 * mu * tau * y[2]; //V''
    f[2] = y[0] * y[0] - y[2] / tau; //(dz/dt)
    return GSL_SUCCESS;
}

//Интегрирование
double integral(double a, double b) {
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function F;
    F.function = [](double x, void* params) { return 1.0 / (sqrt(x) * (1 + pow(x, 1.0 / 3))); };
    gsl_integration_qags(&F, a, b, 0, 1e-7, 1000, workspace, &result, &error);
    gsl_integration_workspace_free(workspace);
    return result;
}

//Метод Гаусса для решения системы уравнений
void solveSystem(double& A, double& B, double& C) {
    double matrix[3][4] = { {16, 24, 18, 50}, {-24, 46, -42, -90}, {18, -42, 49, 85} };
    int n = 3;

    for (int i = 0; i < n; i++) {
        double div = matrix[i][i];
        for (int j = i; j <= n; j++)
            matrix[i][j] /= div;
        for (int k = 0; k < n; k++) {
            if (k == i) continue;
            double mul = matrix[k][i];
            for (int j = i; j <= n; j++)
                matrix[k][j] -= mul * matrix[i][j];
        }
    }
    A = matrix[0][3];
    B = matrix[1][3];
    C = matrix[2][3];
}

double f(double x) {
    return tan(M_PI * x / 4) - x - 3;
}

double f_prime(double x) {
    double sec = 1.0 / pow(cos(M_PI * x / 4), 2);
    return (M_PI / 4) * sec - 1;
}

//Метод Ньютона для нахождения наимаеньшего положительного корня
double newton(double x0, double epsilon) {
    double x = x0;
    while (true) {
        double x_next = x - f(x) / f_prime(x);
        if (abs(x_next - x) < epsilon) {
            return x_next;
        }
        x = x_next;
    }
}

int main() {
    double epsilon = 1e-10;
    double x0 = 1.5; //Начальное приближение для нахождения x*
    double x_star = newton(x0, epsilon);
    double integ = integral(1.0, 5.0);
    double A, B, C;

    tau *= x_star;
    w *= integ;

    w = 1.03;


    std::cout << "Constant values:\n\n";
    std::cout << std::fixed << std::setprecision(10) << "tau =  " << tau << " ( x* = " << x_star << " )" << std::endl;

    solveSystem(A, B, C);
    std::cout << " A  =  " << A << "\n B  =  " << B << "\n C  =  " << C << "\n T  = " << T << "\n mu =  " << mu << "\n w  =  " << w <<
        " ( integral value = " << integ << " )" << "\n";
    std::cout << "======================================================\n\n";

    const double h = 0.01;
    const int numPoints = T / h + 1;
    double t = 0.0;
    double y[3] = { A, B, C };
    double V, V_prime;

    //Решение СУ через RKF45
    gsl_odeiv2_system system = { func, nullptr, 3, nullptr };
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);

    std::cout << "(V; t):\t\t\t\t(V; V'):\n";
    for (int i = 0; i < numPoints; i++) {
        double ti = i * h;
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);
        if (status != GSL_SUCCESS) {
            std::cerr << "Integration failed!" << std::endl;
            break;
        }
        V = y[0];
        V_prime = y[1];
        //std::cout << "(" << V << "; " << ti << ")\t"; //Вывод точек для графика (V, t)
        std::cout << "" << V_prime << "\n";// << V_prime << ")" << std::endl; //Вывод точек для графика (V, V')
    }
    gsl_odeiv2_driver_free(driver);
    return 0;
}