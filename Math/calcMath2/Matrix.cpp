#include <iostream>
#include <vector>
#include <gsl/gsl_linalg.h>


//Создание матрицы А с определенным элементом alpha, создание вектора zi (i = 1, i <= N)
void calculate(double alpha, int N, gsl_matrix *A, gsl_vector *z) {
	//Пробелы для вывода на экран
	std::string str = " ";
	double alp = alpha;
	int alpInt = 0;
	while (true) {
		alp *= 10;
		alpInt = alp;
		if (alpInt % 10 == 0) {
			break;
		}
		str += " ";
	}
	double element;
	//Вывод матрицы А
	std::cout << "Matrix with alpha = " << alpha << ":\n";
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			element = (j >= i) ? j - i + 1.0 : alpha;
			gsl_matrix_set(A, i, j, element);
			std::cout << element << " ";
			if (alpha != 0 && element != alpha) {
				std::cout << str;
			}
		}
		std::cout << "\n";
	}
	//Вывод вектора z
	std::cout << "\nVector z:\n";
	for (int i = 1; i <= N; i++) {
		double sum1 = 0.0, sum2 = 0.0, gk = 0.0;
		for (int k = 1; k <= i - 1; k++) {
			gk = pow(2, k - 2);
			sum1 += gk;
		}
		for (int k = i; k <= N; k++) {
			gk = pow(2, k - 2);
			sum2 += (k - i + 1) * gk;
		}
		double res = alpha * sum1 + sum2;
		gsl_vector_set(z, i - 1, res);
		std::cout << res << "\n";
	}
}

//Норма Фробениуса
double frobeniusNorm(int N, gsl_matrix* A) {
	double norm = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			norm += pow((gsl_matrix_get(A, i, j)), 2);
		}
	}
	return sqrt(norm);
}

//Рассчитывание числа обусловленности через нормы Фробениуса
double matrixCond(int N, gsl_matrix* A) {
	double cond;
	int signum;
	gsl_permutation* perm = gsl_permutation_alloc(N);
	gsl_linalg_LU_decomp(A, perm, &signum);
	gsl_matrix* invertedA = gsl_matrix_alloc(N, N);
	gsl_linalg_LU_invert(A, perm, invertedA);
	cond = frobeniusNorm(N, A) * frobeniusNorm(N, invertedA);

	gsl_permutation_free(perm);
	gsl_matrix_free(invertedA);
	return cond;
}

int main() {
	int N = 6;
	double alpha[] = { 0, 0.25, 0.49, 0.499 };
	for (double a : alpha) {
		//Выделение памяти под матрицу А и вектор z
		gsl_matrix *A = gsl_matrix_alloc(N, N);
		gsl_vector* z = gsl_vector_alloc(N);
		calculate(a, N, A, z);

		gsl_permutation* perm = gsl_permutation_alloc(N);
		int signum;
		//Рассчитывание вектора х через программы DECOMP, SOLVE
		gsl_vector* x = gsl_vector_alloc(N);
		gsl_linalg_LU_decomp(A, perm, &signum);
		gsl_linalg_LU_solve(A, perm, z, x);

		std::cout <<"\nSolution x:\n";
		gsl_vector_fprintf(stdout, x, "%g");

		double cond = matrixCond(N, A);
		std::cout << "\nCondition (Frobenius) number cond(A) = " << cond << ".\n";
		
		//Поиск максимума по k из gk
		gsl_vector* g = gsl_vector_alloc(N);
		for (int k = 1; k <= N; k++) {
			gsl_vector_set(g, k - 1, pow(2, k - 2));
		}
		//Вычитание вектора g от вектора x с записью результата в вектор x
		gsl_vector_sub(x, g);
		//Расчет погрешности
		double err = gsl_blas_dnrm2(x) / gsl_blas_dnrm2(g);
		std::cout << "\nRelative error ||x - g|| / ||g||: " << err << "\n\n";
		std::cout << "==========================================================\n";

		gsl_matrix_free(A);
		gsl_vector_free(z);
		gsl_vector_free(x);
		gsl_vector_free(g);
		gsl_permutation_free(perm);
	}
	return 0;
}