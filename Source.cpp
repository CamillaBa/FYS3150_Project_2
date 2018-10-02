#include <armadillo>
#include <iostream>
#include <fstream>
#include <ctime>
#include <ratio>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>
#include <cmath>

double ** sparse( int n) {
	/* This function takes as argument an integer n and returns a new
	2-dimensional array of doubles of shape n times n
	*/

	double ** A = new double *[n];
	for (int i = 0; i < n; i++) {
		A[i] = new double [n] ();
	}
	return A;
}

double* zeros(int n) {
	/* This function takes as argument an integer n and returns a new 1-dimensional
	array of doubles of size n.
	*/

	double * v = new double[n]();
	return v;
}

void clear_memory(double ** & A, int n) {
	/* This function deletes 2 dimensional double arrays.
	*/

	for (int i = 0; i < n; i++) {
		delete[] A[i];
	}
	delete[] A;
}

void get_max(double** &A, int n, int &k, int &l, double &maxval) {
	/* This function takes as argument a double** A and of shape n times n, 
	and sets (k,l) so that Akl is the entry in A with largest absolute value. This value
	is then stored as the double "maxval".
	*/

	double Aij2;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i != j) {
				Aij2 = A[i][j] * A[i][j];
				if (Aij2 > maxval) {
					maxval = Aij2;
					k = i;
					l = j;
				}
			}
		}
	}
}

void print(double** & A, int n) { // figure out a way to avoid passing N as argument
	/* A simple function that prints 2d double arrays to the terminal.
	*/

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void jacobi_rot_algorithm(double ** & A, double ** & U, double * & lambda, int N) { // figure out a way to avoid passing N as argument
	/* This function implements Jacobi's rotiational algorithm to find the eigenvalues and eigenvectors
	of the matrix (2d-array) A. The eigenvalues are then stored in the vector (1d-array) lambda, and
	the eigenvectors are stored as the colum vectors of the (2d-array) U.
	*/

	int k; int l; double epsilon = 0.000000001;
	double tau;
	double t; // tan(theta)
	double c; // cos(theta)
	double c2; // cos^2(theta)
	double s; // sin(theta)
	double s2; // sin^2(theta)
	double cs; // sin(theta)*cos(theta)
	double root;
	double Aik; double Ail; double Akk; double All; double Akl; double Alk; double two_cs_Akl; //matrix values used in code
	double Uik; double Uil;
	double maxval;

	bool condition = true;

	
	int count = 0;
	while (condition) {
		maxval = 0.0;
		get_max(A, (N - 1), k, l, maxval);
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
				t = 1.0 / (tau - root);
			}
			c = 1.0 / sqrt(1 + t * t); s = t * c;
			c2 = c * c; s2 = s * s; cs = c * s; // trig functions update
		}
		else {
			c = 1.0; s = 0.0; c2 = 1.0; s2 = 0.0; cs = 0.0; // trig functions update
		}

		// update A
		Akk = A[k][k]; All = A[l][l]; Akl = A[k][l]; Alk = A[l][k]; two_cs_Akl = 2 * cs * Akl;
		A[k][k] = Akk * c2 - two_cs_Akl  + All * s2;
		A[l][l] = Akk * s2 + two_cs_Akl + All * c2;
		A[k][l] = 0; //(Akk - All) * cs + c2 * Akl - s2 * Alk;
		A[l][k] = 0; //(Akk - All) * cs + c2 * Alk - s2 * Akl;
		for (int i = 0; i < N - 1; i++) {
			if (i != k && i != l) {
				Aik = A[i][k]; Ail = A[i][l];
				A[i][k] = Aik * c - Ail * s;
				A[k][i] = A[i][k];
				A[i][l] = Ail * c + Aik * s;
				A[l][i] = A[i][l];
			}

			// solve for eigenvectors
			Uik = U[i][k]; Uil = U[i][l];
			U[i][k] = c * Uik - s * Uil;
			U[i][l] = c * Uil + s * Uik;
		}
	}

	// enable coment below to print number of interations needed to converge
	//std::cout <<"matrix size: "<< N-1 <<" steps: " << count << "\n\n";

	// obtain eigenvalues
	for (int i = 0; i < N - 1; i++) { lambda[i] = A[i][i]; }
}

void unit_test1(void) {
	/* unit test 1 (search for largest element). Testing the function "get_max".
	*/

	double maxval = 0.0;
	double ** R = sparse(5);
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			R[i][j] = rand() % 10;
		}
	}
	int k, l;

	std::cout << "Unit test: seach for entry with largest absolute value in matrix" << std::endl << std::endl;
	print(R, 5);
	get_max(R, 5, k, l, maxval); // returns integers k,l such that |Rkl| is maximal
	std::cout << "row: " << k << " colum: " << l << std::endl << "Largest absolute value: " << R[k][l] << std::endl << std::endl;
	clear_memory(R, 5);
}

void unit_test2(void) {
	/* unit test 2 (analytical eigenvalues) Checking whether the function "jacobi_rot_algorithm" finds
	the correct eigenvalues and eigenvectors.
	*/

	double n = 4; int n_minus_1;
	double** A = sparse(n); double ** U = sparse(n);  double* lambda = zeros(n);

	// initializing tridiagonal toeplitz matrix
	for (int i = 0; i <n; i++) { A[i][i] = 2.0; }
	for (int i = 0; i <n-1; i++) { A[i + 1][i] = -1.0; A[i][i + 1] = -1.0; }

	// initialize U as the identity matrix
	for (int i = 0; i < n; i++) { U[i][i] = 1; }

	// print matrix we solve the eigenvalue/ eigenvector problem for:
	std::cout << "Unit test: finding eigenvalues and eigenvectors for tridiagonal toeplitz matrix " << std::endl << std::endl;
	print(A, n);

	// jacobi's rotation algorithm (finding eigenvalues and eigenvectors)
	jacobi_rot_algorithm(A, U, lambda, n+1);

	// printing eigenvalues
	std::cout << "yields eigenvalues: " << std::endl << std::endl;
	for (int j = 0; j < n; j++) {
		std::cout << lambda[j] << std::endl;
	}
	std::cout << std::endl;
	std::cout << "and eigenvectors: " << std::endl << std::endl;
	//eigenvectors
	for (int j = 0; j < n; j++) {
		std::cout << '(';
		for (int i = 0; i < n-1; i++) {
			std::cout << U[i][j] << ',';
		}
		n_minus_1 = n - 1;
		std::cout << U[n_minus_1][j] << ')' << std::endl;
	}
	std::cout << std::endl;
	std::cout <<"Analytical solution: " << "http://www.wolframalpha.com/input/?i=%7B%7B2,-1,0,0%7D,%7B-1,2,-1,0%7D,%7B0,-1,2,-1%7D,%7B0,0,-1,2%7D%7D" << std::endl;
	std::cout << std::endl;

	// clear memory
	clear_memory(A, n); clear_memory(U, n);  delete[] lambda;
}

void comparison(int N) {
	// setting up time variables:
	using namespace std::chrono;
	high_resolution_clock::time_point t1, t2; duration<double> time_span;

	// setting up variables for arrays
	double** A = sparse(N - 1); double ** U = sparse(N - 1); double * lambda = zeros(N - 1);

	// setting up armadillo matrices and vectors
	arma::sp_mat  B(N - 1, N - 1); 	arma::vec eigval; 	arma::mat eigvec; arma::mat C(B);
	B.diag() += 2; B.diag(1) += -1; B.diag(-1) += -1;

	// determining number of experiments
	int number_experiments=10;

	// open file for writing out times:
	std::ofstream myfile;
	std::string filename = "jacobi_vs_armadillo_N_equals_"; 
	filename += std::to_string(N);
	filename += ".txt";
	myfile.open(filename);

	std::cout << "Time used by our implementation of Jacobi's algorithm vs armadillo for eigenvalue/eigenvector problem. \n \n";
	std::cout << "Jacobi: "<< '	' << "armadillo: \n";
	for (int j = 0; j < number_experiments; j++) {
		// initialize A as tridiagonal toeplitz matrix
		for (int i = 0; i < N - 1; i++) { A[i][i] = 2.0; }
		for (int i = 0; i < N - 2; i++) { A[i + 1][i] = -1.0; A[i][i + 1] = -1.0; }

		// initialize U as the identity matrix
		for (int i = 0; i < N - 1; i++) { U[i][i] = 1; }

		t1 = high_resolution_clock::now();
		jacobi_rot_algorithm(A, U, lambda, N);
		t2 = high_resolution_clock::now();
		time_span = duration_cast<duration<double>>(t2 - t1);

		// printing message
		std::cout << time_span.count() << '	' ;
		myfile << time_span.count() << '	';

		t1 = high_resolution_clock::now();
		eig_sym(eigval, eigvec, C);
		t2 = high_resolution_clock::now();
		time_span = duration_cast<duration<double>>(t2 - t1);

		// printing message
		std::cout << time_span.count() << std::endl;
		myfile << time_span.count() << std::endl;
	}
	std::cout << std::endl;

	// clearing memory and closing file
	clear_memory(A, N - 1); clear_memory(U, N - 1);  delete[] lambda;
	myfile.close();
}

void count_iterations(int n) {
	// setting up variables for arrays
	double** A = sparse(n); double ** U = sparse(n); double * lambda = zeros(n);

	// initializing tridiagonal toeplitz matrix
	for (int i = 0; i < n; i++) { A[i][i] = 2.0; }
	for (int i = 0; i < n - 1; i++) { A[i + 1][i] = -1.0; A[i][i + 1] = -1.0; }

	// initialize U as the identity matrix
	for (int i = 0; i < n; i++) { U[i][i] = 1; }

	// jacobi's algorithm
	jacobi_rot_algorithm(A, U, lambda, n);
}


int main(){
	/* Uncomment unit_test1() and unit_test2() for unit tests.
	Uncomment comparison(50) for comparison of armadillo solver for eigenvectors/eigenvalues vs our implementation of Jacobi's method.
	The result is printed to "jacobi_vs_armadillo_N_equals_50.txt".

	Uncomment the count_iterations(i) loop and the part in the function "jacobi_rot_algorith" that prints to terminal,
	to see the steps needed to converge for matrices of different size.
	*/

	//unit_test1();
	//unit_test2();
	//comparison(50);
	//for (int i = 11; i <= 101; i += 10) {
	//	count_iterations(i);
	//}

	double pi = 3.14159265359;
	int N; double h; double hh; double** A; double** U; double* lambda; double* lambda_analytical;

	//========================================================================================================================

	// comparison of tridiagonal toeplitz matrix eigenvalues with analytical eigenvalues (buckling beam)
	N = 10; h = 1.0 / N; hh = h * h;
	A = sparse(N - 1); U = sparse(N - 1); lambda = zeros(N - 1); lambda_analytical = zeros(N - 1);

	for (int i = 0; i < N - 1; i++) { A[i][i] = 2.0/hh; }
	for (int i = 0; i < N - 2; i++) { A[i + 1][i] = -1.0 / hh; A[i][i + 1] = -1.0 / hh; }

	// Jacobi's rotation algorithm
	jacobi_rot_algorithm(A, U, lambda, N);
	for (int j = 0; j < N - 1; j++) {
		lambda_analytical[j] = (2.0-2.0*cos(pi*(j+1)/N) )/ hh;
	}

	// print relative error vs analytical solution
	std::sort(lambda, lambda + (N - 1));
	std::sort(lambda_analytical, lambda_analytical + (N - 1));
	std::cout << "Relative error of eigenvalues of buckling beam (N = " << N << ")\n\n";
	std::cout.precision(17);
	std::cout << "analytical: " << "	numerical: " << "	relative error:" << std::endl;


	for (int j = 1; j <= N - 1; j++) {
		std::cout << lambda_analytical[j] << '	'<< lambda[j] << '	' << abs(lambda[j] - lambda_analytical[j])/abs(lambda_analytical[j]) <<std::endl;
	}
	std::cout << std::endl;

	// clear memory
	clear_memory(A, N - 1); clear_memory(U, N - 1);  delete[] lambda; delete[] lambda_analytical;

	//========================================================================================================================

	// schrödingers equation (harmonic oscillator)
	N = 150; h = 10.0 / N; hh = h * h;
	A = sparse(N - 1); U = sparse(N - 1); lambda = zeros(N - 1);

	// approximate double derivative + potential term
	for (int i = 0; i < N - 1; i++) { A[i][i] = 2.0 / hh+i*i*hh; }
	for (int i = 0; i < N - 2; i++) { A[i + 1][i] = -1.0 / hh; A[i][i + 1] = -1.0 / hh; }

	// Jacobi's rotation algorithm
	jacobi_rot_algorithm(A, U, lambda, N);

	// printing 10 smallest eigenvalues to terminal
	std::sort(lambda, lambda + (N - 1));
	std::cout << "Solving Schrödingers for harmonic oscillator potential with (N=" << N << " and infinity = " << h*N << ")" << std::endl;
	std::cout << "The 10 smalles eigenvalues (energies): \n\n";
	for (int j = 0; j < 10; j++) {
		std::cout << lambda[j] << std::endl;
	}
	std::cout << std::endl;

	// clear memory
	clear_memory(A, N - 1); clear_memory(U, N - 1);  delete[] lambda;

	//========================================================================================================================

	// schrödingers equation (two particles)
	N = 150; h = 10.0 / N; hh = h * h; double omega_r = 5.0; double ih;
	A = sparse(N - 1); U = sparse(N - 1); lambda = zeros(N - 1);
	
	// approximate double derivative + potential term
	for (int i = 0; i < N - 1; i++) {
		ih = i * h;
		A[i][i] = 2.0 / hh + omega_r * omega_r* ih*ih + 1/ih;
	}
	for (int i = 0; i < N - 2; i++) { A[i + 1][i] = -1.0 / hh; A[i][i + 1] = -1.0 / hh; }

	// Jacobi's rotation algorithm
	jacobi_rot_algorithm(A, U, lambda, N);

	// printing 10 smallest eigenvalues to terminal
	std::sort(lambda, lambda + (N - 1));
	std::cout << "Solving Schrödingers for 2-particles with (N=" << N << " and infinity = " << h * N << " and omega_r="<< omega_r <<" )" << std::endl;
	std::cout << "The 10 smalles eigenvalues (energies): \n\n";
	for (int j = 0; j < 10; j++) {
		std::cout << lambda[j] << std::endl;
	}

	// clear memory
	clear_memory(A, N - 1); clear_memory(U, N - 1);  delete[] lambda;

	//========================================================================================================================

	// print success message if program runs
	std::cout << std::endl << "success!" << std::endl;

	std::cout << std::cin.get();
};