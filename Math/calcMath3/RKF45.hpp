void r4_fehl(void f(float t, float y[], float yp[]), int neqn,
	float y[], float t, float h, float yp[], float f1[], float f2[], float f3[],
	float f4[], float f5[], float s[]);
int r4_rkf45(void f(float t, float y[], float yp[]), int neqn,
	float y[], float yp[], float* t, float tout, float* relerr, float abserr,
	int flag);
float r4_sign(float x);

void r8_fehl(void f(double t, double y[], double yp[]), int neqn,
	double y[], double t, double h, double yp[], double f1[], double f2[], double f3[],
	double f4[], double f5[], double s[]);
int r8_rkf45(void f(double t, double y[], double yp[]), int neqn,
	double y[], double yp[], double* t, double tout, double* relerr, double abserr,
	int flag);
double r8_sign(double x);

void timestamp();
