#include <cmath>
#include <armadillo>
#include <iostream>
#include <iomanip>
//using namespace arma;

double** sparse(int m, int n) {
	double ** A = new double *[n];
	for (int i = 0; i < n; i++) {
		A[i] = new double [n] ();
	}
	return A;
}

double* zeros(int n) {
	double * v = new double [n] ();
	return v;
}

//void* print(double** & A) {
//	std::cout << std::fixed << std::setprecision(8) << std::setw(8);
//	for (int i = 0; i < n; i++) {
//		for (int j = 0; j < n; j++) {
//			std::cout << A[i][j] << " ";
//		}
//		std::cout << std::endl;
//	}
//	std::cout << std::endl;
//}

int main(){
	double pi = 3.14159265358979323846264338328; double epsilon = 0.00000000000000001;
	int N = 5; double h = 1.0 / N; double hh = h * h;
	double** A = sparse(N - 1, N - 1); double* lambda = zeros(N - 1);

    // approximate double derivative
	for (int i = 0; i < N - 1; i++) {A[i][i] = 2.0/hh;}
	for (int i = 0; i < N - 2; i++) {A[i+1][i] = -1.0/hh; A[i][i+1] = -1.0/hh;}

	int k; int l; 
	double tau;
	double t; double t1; double t2; //tan(theta)
	double c; //cos(theta)
	double c2; //cos^2(theta)
	double s; //sin(theta)
	double s2; //sin^2(theta)
	double cs; //sin(theta)*cos(theta)
	double two_Akl;
	double root;
	double maxval;
	double Aij2; double Aik; double Ail; double Akk; double All; double Akl; double Alk; //matrix values used in code

	int count = 0;
	bool condition = true;
	while (condition) {
		maxval = 0.0;
		for (int i = 0; i < N - 1; i++) {
			for (int j = 0; j < N - 1 && i!=j; j++) {
				Aij2 = A[i][j] * A[i][j];
				if (Aij2 > maxval) {
					maxval = Aij2;
					k = i;
					l = j;
				}
			}
		}
		condition = (maxval > epsilon);
		count += 1;

		// calculate values for tau, tan(theta), cos(theta), and sin(theta)
		if (A[k][l] != 0.0) {
			tau = (A[l][l] - A[k][k]) / (2.0 * A[k][l]);
			root = sqrt(1 + tau * tau);
			if (tau > 0) { 
				t = 1.0 / (tau + root); 
			}
			else { 
				t = -1.0 / (-tau + root);
			}

			c = 1.0 / sqrt(1 + t * t); s = t * c;
			c2 = c * c; s2 = s * s; cs = c * s; //trig functions update
		}
		else {
			c = 1.0; s = 0.0; c2 = 1.0; s2 = 0.0; cs = 0.0; //trig functions update
		}

		// update A
		for (int i= 0; i < N - 1 && i != k && i != l; i++) {
			Ail = A[i][l]; Aik = A[i][k];
			A[i][k] = Aik * c - Ail * s; 
			A[k][i] = A[i][k];
			A[i][l] = Ail * c + Aik * s; 
			A[l][i] = A[i][l];
		}
		Akk = A[k][k]; All = A[l][l]; Akl = A[k][l]; Alk = A[l][k]; two_Akl = 2 * Akl;
		A[k][k] = Akk * c2 - two_Akl * cs + All * s2;
		A[l][l] = Akk * s2 + two_Akl * cs + All * c2;
		A[k][l] = 0; //(Akk - All) * cs + c2 * Akl - s2 * Alk;
		A[l][k] = 0; //(Akk - All) * cs + c2 * Alk - s2 * Akl;
	}

	// get eigenvalues
	for (int i = 0; i < N - 1; i++) { lambda[i] = A[i][i]; }

	// analytical eigenvalues
	double* lambda_analytical = zeros(N - 1);
	for (int j = 0; j < N - 1; j++) {
		lambda_analytical[j] = (2.0  - 2.0 * cos((j + 1)*pi / N))/hh;
		std::cout << lambda[j] << ' ' << lambda_analytical[j] << std::endl;
		//std::cout << abs(lambda_analytical[j]- lambda[j])/abs(lambda_analytical[j]) << std::endl;
	}

	std::cout << "success!";
	std::cout << std::cin.get();

};








