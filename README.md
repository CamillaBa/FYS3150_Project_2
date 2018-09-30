# FYS3150_Project_2

This program consists of 9 functions in addition to the main function. We first have
three functions that deal with arrays.

double ** sparse( int n)
double* zeros(int n)
void clear_memory(double ** & A, int n)

The function "sparse" returns an array of dimension n times n. 
The function "zeros" returns an array of dimension n.
The function "clear_memory" deletes an n times n array from memory.

Then there are 3 functions that take a matrix and do something with it.
We have:

void get_max(double** &A, int n, int &k, int &l, double &maxval) 
void print(double** & A, int n)
void jacobi_rot_algorithm(double ** & A, double ** & U, double * & lambda, int N)

The function "get_max" sets the integer maxval equal to the pivot element (entry with greatest absolute value) of A.
Moreover, it sets k and l equal to the row number and column number corresponding to the pivot element.
  The function "print" simply prints the array A of dimension n times n to the terminal.
The function "jacobi_rot_algorithm" finds eigenvectors and associated eigenvalues for a given matrix A.
It then stores the eigenvectors in the matrix U as column vectors, and the eigenvalues are stored in the array lambda.

There are two unit tests:

void unit_test1(void)
void unit_test2(void)

"unit_test1" checks the validity of the function "get_max" and prints the result to the terminal.,
"unit_test2" checks the validity of "jacobi_rot_algorithm".

Lastly, there is a function 

void comparison(int N)

This guy compares our implementation of Jacobi's algorithm to the eigenvalue/eigenvector solver of armadillo
in terms of speed. The results are then stored as files with filenames of the type "jacobi_vs_armadillo_N_equals_100".


The main program consists of three parts, separated by lines of equalsigns (============).
The first part compares the tridiagonal toeplitz matrix eigenvalues found using "jacobi_rot_algorithm" to the analytical ones.
The second part finds the 10 lowest eigenvalues of a dimensionless schrödingers equation with harmonic oscillator potential,
and the last part finds the 10 lowest eigenvalues of a dimensionless schrödingers with two particles with variable oscillation.
The results of each part is printed to the terminal.






